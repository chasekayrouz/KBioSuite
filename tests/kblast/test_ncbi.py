import builtins
import io
import json
import os
import threading
from concurrent.futures import ThreadPoolExecutor
from datetime import UTC, datetime, timedelta
from pathlib import Path
from types import SimpleNamespace

import pytest
import requests
from Bio import SeqIO
from conftest import PROTEIN, blast_xml, genomic_region

from kblast.config import Config
from kblast.models import Run, Settings
from kblast.ncbi import NcbiClient
from kblast.network import RateLimiter, RemoteError, SubmissionUnknown, Transport
from kblast.sequence import parse_query


class VirtualClock:
    def __init__(self):
        self.now = 100.0
        self.delays = []

    def wait(self, delay):
        self.delays.append(delay)
        self.now += delay
        return False

    def is_set(self):
        return False


class Response:
    def __init__(self, status=200, text="", retry_after=""):
        self.status_code = status
        self.ok = status < 400
        self.text = text
        self.headers = {"Retry-After": retry_after} if retry_after else {}


def transport_fixture(tmp_path, monkeypatch, responses):
    clock = VirtualClock()
    monkeypatch.setattr("kblast.lifecycle.time.monotonic", lambda: clock.now)
    transport = Transport(
        Config(api_key="test_credential"),
        clock,
        tmp_path / "requests.jsonl",
        RateLimiter(0.11),
        RateLimiter(10),
    )
    calls = []

    def post(url, data, timeout):
        calls.append((clock.now, url, data, timeout))
        value = responses.pop(0)
        if isinstance(value, Exception):
            raise value
        return value

    monkeypatch.setattr(transport.session, "post", post)
    return transport, clock, calls


@pytest.mark.parametrize("failure", [requests.Timeout, requests.exceptions.ChunkedEncodingError])
def test_transient_retry_and_credentials_never_enter_log(tmp_path, monkeypatch, failure):
    transport, clock, calls = transport_fixture(
        tmp_path,
        monkeypatch,
        [
            Response(429, "error includes test_credential", "5"),
            failure("url?api_key=test_credential"),
            Response(200, ">PROTA.1\nMAAA\n"),
        ],
    )
    assert transport.request(
        "eutils", "efetch.fcgi", {"db": "protein", "id": "PROTA.1"}
    ).startswith(">")
    assert len(calls) == 3
    assert 5 in clock.delays and 4 in clock.delays
    assert all(call[2]["api_key"] == "test_credential" for call in calls)
    log = (tmp_path / "requests.jsonl").read_text()
    assert "test_credential" not in log and "api_key" not in log
    assert all("?" not in json.loads(line)["endpoint"] for line in log.splitlines())


def test_retry_is_bounded(tmp_path, monkeypatch):
    transport, _, calls = transport_fixture(tmp_path, monkeypatch, [Response(503)] * 3)
    with pytest.raises(RemoteError, match="3 attempt"):
        transport.request("eutils", "efetch.fcgi", {"id": "X"})
    assert len(calls) == 3


def test_very_long_retry_after_stops_without_early_retry(tmp_path, monkeypatch):
    transport, _, calls = transport_fixture(
        tmp_path, monkeypatch, [Response(429, retry_after="900")]
    )
    with pytest.raises(RemoteError):
        transport.request("eutils", "efetch.fcgi", {"id": "X"})
    assert len(calls) == 1


def test_blast_submission_is_sent_once_and_query_not_logged(tmp_path, monkeypatch):
    transport, _, calls = transport_fixture(
        tmp_path, monkeypatch, [requests.Timeout("test_credential")]
    )
    with pytest.raises(SubmissionUnknown):
        transport.request("blast", "", {"CMD": "Put", "QUERY": PROTEIN}, submit=True)
    assert len(calls) == 1
    assert "api_key" not in calls[0][2]
    log = (tmp_path / "requests.jsonl").read_text()
    assert PROTEIN not in log and "test_credential" not in log
    assert "query_sha256" in log


def test_blast_get_retries_keep_one_minute_spacing(tmp_path, monkeypatch):
    transport, _, calls = transport_fixture(
        tmp_path,
        monkeypatch,
        [Response(503), Response(200, "Status=WAITING"), Response(200, "Status=READY")],
    )
    transport.request("blast", "", {"CMD": "Get", "RID": "EXAMPLE"}, poll=True)
    transport.request("blast", "", {"CMD": "Get", "RID": "EXAMPLE"}, poll=True)
    assert all(b[0] - a[0] >= 60 for a, b in zip(calls, calls[1:], strict=False))


def test_cancelled_transport_sends_no_request(tmp_path, monkeypatch):
    from kblast.network import Cancelled

    transport, _, calls = transport_fixture(tmp_path, monkeypatch, [])
    event = threading.Event()
    event.set()
    transport.cancel = event
    with pytest.raises(Cancelled):
        transport.request("eutils", "efetch.fcgi", {"id": "X"})
    assert not calls


def test_ready_result_retrieval_keeps_ten_second_gate_without_an_extra_minute(
    tmp_path, monkeypatch
):
    transport, _, calls = transport_fixture(
        tmp_path, monkeypatch, [Response(200, "Status=READY"), Response(200, blast_xml())]
    )
    run = Run(
        "a" * 32,
        parse_query(PROTEIN),
        Settings(2, 100),
        rid="EXAMPLE",
        submitted_at=(datetime.now(UTC) - timedelta(minutes=2)).isoformat(),
    )
    assert "BlastXML2" in NcbiClient(transport, tmp_path).wait_for_blast(run, lambda: None)
    assert calls[1][0] - calls[0][0] == 10
    assert all(c[2]["CMD"] == "Get" and c[2]["RID"] == "EXAMPLE" for c in calls)


def test_encoded_post_evidence_omits_query_and_credentials(tmp_path, monkeypatch):
    response = Response(200, ">PROTA.1\nMAAA\n")
    response.request = requests.Request(
        "POST",
        "https://example.test/",
        data={"id": "PROTA.1", "api_key": "test_credential", "email": "test@example.test"},
    ).prepare()
    transport, _, _ = transport_fixture(tmp_path, monkeypatch, [response])
    transport.request("eutils", "efetch.fcgi", {"id": "PROTA.1"})
    log = (tmp_path / "requests.jsonl").read_text()
    assert "test_credential" not in log and "test@example.test" not in log
    received = json.loads(log.splitlines()[-1])
    assert received["api_key_used"] is True
    assert received["sent_parameters"] == {"id": "PROTA.1"}


class StubTransport:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls = []
        self.cancel = VirtualClock()
        self.config = SimpleNamespace(blast_timeout_seconds=3600)

    def request(self, service, endpoint, parameters, **kwargs):
        self.calls.append((service, endpoint, parameters, kwargs))
        return self.responses.pop(0)


def test_record_length_resolution_is_cached_and_versioned(tmp_path):
    summary = {"result": {"uids": ["1"], "1": {"slen": 1000, "accessionversion": "NC_000001.1"}}}
    http = StubTransport([json.dumps(summary)])
    client = NcbiClient(http, tmp_path)
    info = client.nucleotide_info("NC_000001")
    assert info["accession"] == "NC_000001.1" and info["length"] == 1000
    assert client.nucleotide_info("NC_000001") == info
    assert len(http.calls) == 1


def test_unindexed_accession_falls_back_to_efetch(tmp_path):
    raw, _ = genomic_region("NC_000001.1", 1, 1000)
    http = StubTransport(['{"error":"unindexed accession"}', raw])
    info = NcbiClient(http, tmp_path).nucleotide_info("NC_000001.1")
    assert info["length"] == 1000 and "EFetch" in info["method"]
    assert http.calls[-1][2]["rettype"] == "gb"


def test_region_fetch_is_native_validated_and_cached(tmp_path):
    raw, _ = genomic_region("NC_000002.1", 800, 1000)
    http = StubTransport([raw])
    client = NcbiClient(http, tmp_path)
    text, record, stamp = client.fetch_context("NC_000002.1", 800, 1000)
    assert text == raw
    assert len(record) == 201
    assert http.calls[0][2]["strand"] == 1
    assert http.calls[0][2]["rettype"] == "gbwithparts"
    assert client.fetch_context("NC_000002.1", 800, 1000)[2] == stamp
    assert len(http.calls) == 1
    assert list(SeqIO.parse(io.StringIO(text), "genbank"))[0].features[1].location.strand == -1


def test_wrong_region_and_protein_accession_are_rejected(tmp_path):
    raw, _ = genomic_region("NC_000002.1", 800, 1000)
    with pytest.raises(ValueError, match="accession"):
        NcbiClient(StubTransport([raw]), tmp_path).fetch_context("NC_000001.1", 800, 1000)
    with pytest.raises(ValueError, match="incomplete"):
        NcbiClient(StubTransport([raw]), tmp_path).fetch_context("NC_000002.1", 800, 999)
    with pytest.raises(ValueError, match="accession"):
        NcbiClient(StubTransport([">OTHER.1\nMAAA\n"]), tmp_path).fetch_protein("PROTA.1", 1)


def test_polling_uses_same_rid_and_official_xml2_format(tmp_path):
    http = StubTransport(["Status=WAITING", blast_xml()])
    client = NcbiClient(http, tmp_path)
    run = Run(
        "a" * 32,
        parse_query(PROTEIN),
        Settings(2, 100),
        rid="EXAMPLE",
        submitted_at=(datetime.now(UTC) - timedelta(minutes=2)).isoformat(),
        rtoe=9000,
    )
    assert "BlastXML2" in client.wait_for_blast(run, lambda: None)
    assert [call[2]["RID"] for call in http.calls] == ["EXAMPLE", "EXAMPLE"]
    assert all(call[2]["FORMAT_TYPE"] == "XML2_S" for call in http.calls)
    assert any(59.9 <= delay <= 60 for delay in http.cancel.delays)


@pytest.mark.parametrize("status", ["FAILED", "UNKNOWN"])
def test_terminal_and_invalid_blast_statuses_stop_without_resubmission(tmp_path, status):
    http = StubTransport([f"Status={status}"])
    run = Run(
        "a" * 32,
        parse_query(PROTEIN),
        Settings(2, 100),
        rid="EXAMPLE",
        submitted_at=(datetime.now(UTC) - timedelta(minutes=2)).isoformat(),
    )
    with pytest.raises(RemoteError):
        NcbiClient(http, tmp_path).wait_for_blast(run, lambda: None)
    assert len(http.calls) == 1 and http.calls[0][2]["CMD"] == "Get"


def test_first_poll_ignores_large_rtoe_and_ready_without_hits_fetches_xml(tmp_path):
    empty_xml = blast_xml().replace(blast_xml().split("<hits>")[1].split("</hits>")[0], "")
    http = StubTransport(["Status=READY\nThereAreHits=no", empty_xml])
    run = Run(
        "a" * 32,
        parse_query(PROTEIN),
        Settings(2, 100),
        rid="EXAMPLE",
        submitted_at=datetime.now(UTC).isoformat(),
        rtoe=11920,
    )
    assert "BlastXML2" in NcbiClient(http, tmp_path).wait_for_blast(run, lambda: None)
    assert sum(http.cancel.delays) <= 10
    assert http.calls[0][3]["poll"] is True and not http.calls[1][3].get("poll", False)


def test_ready_without_xml_is_an_explicit_retrieval_error(tmp_path):
    http = StubTransport(["Status=READY", "Status=READY"])
    run = Run(
        "a" * 32,
        parse_query(PROTEIN),
        Settings(2, 100),
        rid="EXAMPLE",
        submitted_at=(datetime.now(UTC) - timedelta(minutes=2)).isoformat(),
    )
    with pytest.raises(RemoteError, match="READY.*XML2"):
        NcbiClient(http, tmp_path).wait_for_blast(run, lambda: None)


def opened_files(monkeypatch):
    """Real paths of every file opened through open()/io.open (pathlib and python-dotenv)."""
    opened, real_open = [], io.open

    def spy(file, *args, **kwargs):
        if isinstance(file, (str, bytes, os.PathLike)):
            opened.append(Path(os.path.realpath(os.fsdecode(file))))
        return real_open(file, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", spy)
    monkeypatch.setattr(io, "open", spy)
    return opened


def test_key_discovery_reads_only_dotenv_from_the_working_directory(tmp_path, monkeypatch):
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)
    monkeypatch.setattr("kblast.config.Path.home", lambda: tmp_path / "home")
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    monkeypatch.delenv("NCBI_API_KEY_FILE", raising=False)
    monkeypatch.delenv("NCBI_EMAIL", raising=False)
    # Key-like files in the working directory (the retired checkout-relative developer key and
    # copies of the persistent locations) are never opened; only its .env is.
    for decoy in (
        "private-inputs/NCBI-API_key.txt",
        "NCBI-API_key.txt",
        "api_key",
        ".ncbi/api_key",
        ".config/kblast/ncbi.env",
    ):
        (work / decoy).parent.mkdir(parents=True, exist_ok=True)
        (work / decoy).write_text(
            "NCBI_API_KEY_FILE=decoy\n" if decoy.endswith(".env") else "decoy\n"
        )
    (work / ".env").write_text("NCBI_EMAIL=checkout@example.org\n")
    opened = opened_files(monkeypatch)
    config = Config.discover()
    assert config.api_key == "" and config.email == "checkout@example.org"
    real_work = Path(os.path.realpath(work))
    assert {path for path in opened if path.is_relative_to(real_work)} == {real_work / ".env"}


def test_key_discovery_uses_authorized_local_file_without_repr(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("kblast.config.Path.home", lambda: tmp_path / "home")
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    monkeypatch.delenv("NCBI_API_KEY_FILE", raising=False)
    assert Config.discover().api_key == ""
    explicit = tmp_path / "explicit key"
    explicit.write_text("explicit_test_key\n")
    monkeypatch.setenv("NCBI_API_KEY_FILE", str(explicit))
    config = Config.discover()
    assert config.api_key == "explicit_test_key"
    assert "explicit_test_key" not in repr(config)
    monkeypatch.delenv("NCBI_API_KEY_FILE")
    persistent = tmp_path / "home" / ".ncbi" / "api_key"
    persistent.parent.mkdir(parents=True)
    persistent.write_text("persistent_test_key\n")
    assert Config.discover().api_key == "persistent_test_key"
    unrelated = tmp_path / "unrelated working directory"
    unrelated.mkdir()
    monkeypatch.chdir(unrelated)
    assert Config.discover().api_key == "persistent_test_key"
    assert "persistent_test_key" not in repr(Config.discover())
    monkeypatch.setenv("NCBI_API_KEY", "environment_test_key")
    assert Config.discover().api_key == "environment_test_key"


def test_scientific_looking_error_with_credential_is_not_hashed(tmp_path, monkeypatch):
    transport, _, _ = transport_fixture(
        tmp_path,
        monkeypatch,
        [Response(200, "<IPGReportSet><ERROR>test_credential</ERROR></IPGReportSet>")],
    )
    transport.request("eutils", "efetch.fcgi", {"db": "protein", "rettype": "ipg"})
    events = [json.loads(line) for line in (tmp_path / "requests.jsonl").read_text().splitlines()]
    assert "response_sha256" not in events[-1]
    assert events[-1]["request_seconds"] >= 0 and events[-1]["pace_seconds"] >= 0
    assert "test_credential" not in (tmp_path / "requests.jsonl").read_text()


def test_concurrent_transport_owns_and_closes_one_session_per_thread(tmp_path):
    transport = Transport(
        Config(), threading.Event(), tmp_path / "requests.jsonl", RateLimiter(0.36), RateLimiter(10)
    )
    barrier = threading.Barrier(3)

    def session():
        value = transport.session
        barrier.wait(timeout=5)
        return value

    with ThreadPoolExecutor(max_workers=3) as pool:
        sessions = list(pool.map(lambda _: session(), range(3)))
    assert len({id(s) for s in sessions}) == 3
    assert all(s.headers["User-Agent"].startswith("Kblast") for s in sessions)
    transport.close()
