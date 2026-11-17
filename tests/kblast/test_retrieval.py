"""Deterministic protocol regression tests; no test makes a live request."""

import gzip
import io
import json
import threading
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio import SeqIO
from conftest import HEADER, PROTEIN, blast_xml, genomic_region, ipg_text

from kblast.export import export_bundle
from kblast.models import Hit, Run, Settings
from kblast.ncbi import MappingUnresolved, NcbiClient, genbank_records
from kblast.network import Cancelled, RemoteError
from kblast.parsers import ResponseParseError, ipg_xml_reports, parse_blast_xml, parse_ipg
from kblast.pipeline import Pipeline
from kblast.sequence import parse_query
from kblast.storage import RunStore

FIXTURES = Path(__file__).parent / "fixtures"
HISTORY = """<eLinkResult><LinkSet><DbFrom>protein</DbFrom><IdList><Id>1</Id></IdList>
<LinkSetDbHistory><DbTo>nuccore</DbTo><LinkName>protein_nuccore</LinkName>
<QueryKey>1</QueryKey></LinkSetDbHistory><WebEnv>fixture_history</WebEnv></LinkSet></eLinkResult>"""
EMPTY_LINKS = """<eLinkResult><LinkSet><DbFrom>protein</DbFrom>
<IdList><Id>1</Id></IdList></LinkSet></eLinkResult>"""


class ScriptedHttp:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls, self.logs = [], []
        self.cancel = threading.Event()
        self.config = SimpleNamespace(blast_timeout_seconds=3600)

    def request(self, service, endpoint, parameters, **kwargs):
        self.calls.append((endpoint, parameters))
        response = self.responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response

    def safe_message(self, text):
        return text.replace("fixture_secret", "[redacted]")

    def log(self, entry):
        self.logs.append(entry)

    def close(self):
        pass


def public_fixture(name):
    path = FIXTURES / ("wp_320672066." + name)
    return gzip.decompress(path.read_bytes()).decode() if name.endswith(".gz") else path.read_text()


def test_wp_320672066_captured_ipg_gap_and_exact_linked_context(tmp_path):
    tsv = public_fixture("ipg.tsv")
    xml = public_fixture("ipg.xml")
    full = public_fixture("linked.gb.gz")
    context = public_fixture("context.gb.gz")
    source = SeqIO.read(io.StringIO(full), "genbank")
    target = next(
        f for f in source.features if f.qualifiers.get("protein_id") == ["WP_320672066.1"]
    )
    sequence = target.qualifiers["translation"][0]
    hit = Hit(
        1, "WP_320672066.1", "ref|WP_320672066.1|", "Public regression", "", "", 291, 0, 0, 0, 0, []
    )
    assert parse_ipg(tsv, hit, 20_000) == []  # Reproduces the actual upstream IPG gap.
    assert ipg_xml_reports(xml)[0][2] == tsv
    http = ScriptedHttp(
        [
            f">WP_320672066.1\n{sequence}\n",
            xml,
            HISTORY,
            "<eSearchResult><Count>1</Count></eSearchResult>",
            full,
            context,
        ]
    )
    store = RunStore(tmp_path)
    run = Run("a" * 32, parse_query(PROTEIN), Settings(1, 20_000), hits=[hit])
    store.save(run)

    class MappingPipeline(Pipeline):
        def blast(self):
            pass  # Mapping-only test; query data never leaves the temporary fixture store.

    MappingPipeline(store, run, NcbiClient(http, store.directory(run.id))).execute()
    assert run.status == "success" and hit.mapping_status == "success"
    assert len(run.occurrences) == 1 and not http.responses
    o = run.occurrences[0]
    assert (o.nucleotide_accession, o.cds_start, o.cds_end, o.strand) == (
        "NZ_JAWZSQ010000020.1",
        24826,
        25701,
        "-",
    )
    assert (o.context_start, o.context_end, o.record_length) == (4826, 45701, 107044)
    assert o.mapping_method == "ELink + GenBank CDS" and o.ipg_row == 0
    assert o.target_cds_found is True
    assert (store.directory(run.id) / o.context_file).read_text() == context
    count_request = next(p for endpoint, p in http.calls if endpoint == "esearch.fcgi")
    assert count_request["usehistory"] == "y" and count_request["term"] == "#1"
    assert set(run.timings) >= {
        "blast_submission_wait",
        "blast_result_retrieval_parse",
        "protein_mapping",
        "nucleotide_contexts",
        "local_processing",
    }
    assert store.load(run.id).to_dict() == run.to_dict()
    import zipfile

    with zipfile.ZipFile(io.BytesIO(export_bundle(store.directory(run.id), run))) as archive:
        assert archive.read("provenance/hit_0001_000000.linked.gb").decode() == full


def test_link_history_pages_every_record_and_retains_repeated_cds(tmp_path):
    hit = parse_blast_xml(blast_xml(), 1, 20)[0][0]
    records = []
    for i in range(11):
        _, record = genomic_region(f"NC_{i + 1:06d}.1", 1, 1000)
        record.features.append(deepcopy(record.features[1]))
        records.append(record.format("genbank"))
    http = ScriptedHttp(
        [
            HISTORY,
            "<eSearchResult><Count>11</Count></eSearchResult>",
            "".join(records[:10]),
            records[10],
        ]
    )
    rows = NcbiClient(http, tmp_path).linked_occurrences(hit, 100, {hit.accession})
    assert len(rows) == 22 and len({r.id for r in rows}) == 22
    assert rows[0].cds_start == rows[1].cds_start
    fetches = [p for endpoint, p in http.calls if endpoint == "efetch.fcgi"]
    assert [(p["retstart"], p["retmax"]) for p in fetches] == [(0, 10), (10, 1)]
    assert all(p["query_key"] == "1" for p in fetches)


@pytest.mark.parametrize(
    "raw",
    [
        "",
        "<ERROR>Temporary failure</ERROR>",
        "<eLinkResult><LinkSet><DbFrom>protein</DbFrom></LinkSet></eLinkResult>",
    ],
)
def test_unverified_link_response_cannot_be_no_mapping(tmp_path, raw):
    hit = parse_blast_xml(blast_xml(), 1, 20)[0][0]
    with pytest.raises(ResponseParseError):
        NcbiClient(ScriptedHttp([raw]), tmp_path).linked_occurrences(hit, 100, {hit.accession})


def test_verified_empty_links_and_missing_target_are_different_states(tmp_path):
    hit = parse_blast_xml(blast_xml(), 1, 20)[0][0]
    assert (
        NcbiClient(ScriptedHttp([EMPTY_LINKS]), tmp_path).linked_occurrences(
            hit, 100, {hit.accession}
        )
        == []
    )
    raw, _ = genomic_region("NC_000001.1", 1, 1000)
    http = ScriptedHttp([HISTORY, "<eSearchResult><Count>1</Count></eSearchResult>", raw])
    with pytest.raises(MappingUnresolved, match="no matching"):
        NcbiClient(http, tmp_path).linked_occurrences(hit, 100, {"OTHER.1"})


@pytest.mark.parametrize(
    "failure,status",
    [
        (RemoteError("HTTP 503"), "retrieval_failure"),
        (ResponseParseError("Malformed XML"), "parse_failure"),
        (MappingUnresolved("Missing target CDS"), "unresolved_mapping"),
    ],
)
def test_fallback_failure_never_becomes_no_mapping(run_fixture, failure, status):
    _, run, client, pipeline = run_fixture
    run.settings.top_n = 1
    client.fetch_ipg = lambda *_: HEADER

    def fail(*args, **kwargs):
        raise failure

    client.linked_occurrences = fail
    pipeline.execute()
    assert run.hits[0].mapping_status == status and run.status == "partial_success"
    assert not run.occurrences


def test_legacy_no_mapping_is_rechecked_without_changing_the_original():
    run = Run("a" * 32, parse_query(PROTEIN), Settings(1, 100, "nr_cluster_seq"))
    run.hits = parse_blast_xml(blast_xml(), 1, 20)[0]
    run.hits[0].mapping_status, run.status = "no_mapping", "success"
    original = run.to_dict()
    restored = Run.from_dict(original)
    assert restored.hits[0].mapping_status == "pending" and restored.status == "interrupted"
    assert restored.settings.database == "nr_cluster_seq"
    assert original["hits"][0]["mapping_status"] == "no_mapping"


def test_batch_ipg_preserves_report_boundaries_and_row_multiplicity():
    xml = public_fixture("ipg.xml")
    report = xml.split("<IPGReport ", 1)[1].split("</IPGReport>")[0]
    report = "<IPGReport " + report + "</IPGReport>"
    first = report.replace('product_acc="WP_320672066.1"', 'product_acc="OTHER.1"')
    batched = xml.split("<IPGReport ", 1)[0] + first + report + "</IPGReportSet>"
    reports = ipg_xml_reports(batched)
    assert [r[0] for r in reports] == ["OTHER.1", "WP_320672066.1"]
    assert reports[0][2] == reports[1][2]
    with pytest.raises(ResponseParseError, match="statistics"):
        ipg_xml_reports(xml.replace('prot_count="1"', 'prot_count="2"'))


@pytest.mark.parametrize(
    "text",
    [
        HEADER + "1\tRefSeq\tNC_000001.1\n",
        ipg_text().splitlines()[0] + "\n" + "\t".join(["x"] * 12),
    ],
)
def test_truncated_and_overlong_ipg_rows_are_parse_failures(text):
    hit = parse_blast_xml(blast_xml(), 1, 20)[0][0]
    with pytest.raises(ResponseParseError):
        parse_ipg(text, hit, 100)


def test_batch_proteins_and_ipg_are_attributed_by_accession_not_response_order(tmp_path):
    hits = parse_blast_xml(blast_xml(), 2, 20)[0]
    xml = public_fixture("ipg.xml")
    report = xml[xml.index("<IPGReport ") : xml.index("</IPGReport>") + len("</IPGReport>")]
    batched = (
        xml[: xml.index("<IPGReport ")]
        + "".join(report.replace("WP_320672066.1", h.accession) for h in reversed(hits))
        + "</IPGReportSet>"
    )
    fasta = "".join(f">{h.accession}\n{PROTEIN}\n" for h in reversed(hits))
    http = ScriptedHttp([fasta, batched])
    client = NcbiClient(http, tmp_path)
    client.prefetch_hits(hits)
    assert len(http.calls) == 2
    for h in hits:
        assert client.fetch_protein(h.accession, h.rank) == (PROTEIN, h.accession)
        assert h.accession in client.fetch_ipg(h.accession, h.rank)
    assert len(http.calls) == 2


def test_batched_summary_keeps_exact_versions_and_defers_unindexed_lengths(tmp_path):
    raw, _ = genomic_region("NC_000002.1", 1, 1000)
    summary = {"result": {"uids": ["1"], "1": {"slen": 1000, "accessionversion": "NC_000001.1"}}}
    http = ScriptedHttp([json.dumps(summary)])
    client = NcbiClient(http, tmp_path)
    client.prefetch_nucleotide_info(["NC_000001.1", "NC_000002.1", "NC_000001.1"])
    assert len(http.calls) == 1
    assert client.nucleotide_info("NC_000001.1")["length"] == 1000
    assert client.nucleotide_info("NC_000002.1", allow_unindexed=True) == {
        "accession": "NC_000002.1",
        "length": None,
    }
    assert len(http.calls) == 1


def test_full_context_batch_preserves_bytes_and_native_features(tmp_path):
    raws = [genomic_region(a, 1, 1000)[0] for a in ["NC_000001.1", "NC_000002.1"]]
    http = ScriptedHttp(["".join(raws)])
    keys = [("NC_000001.1", 1, 1000), ("NC_000002.1", 1, 1000)]
    results = dict(NcbiClient(http, tmp_path).context_results(keys, set(keys)))
    assert len(http.calls) == 1
    assert http.calls[0][1]["id"] == "NC_000001.1,NC_000002.1"
    assert "seq_start" not in http.calls[0][1]
    assert [results[k][0] for k in keys] == raws
    assert results[keys[1]][1].features[1].location.strand == -1


def test_partial_full_batch_uses_single_fallback_only_for_missing_record(tmp_path):
    raw, _ = genomic_region("NC_000001.1", 1, 1000)
    http = ScriptedHttp([raw, RemoteError("Second record unavailable")])
    keys = [("NC_000001.1", 1, 1000), ("NC_000002.1", 1, 1000)]
    results = dict(NcbiClient(http, tmp_path).context_results(keys, set(keys)))
    assert results[keys[0]][0] == raw and isinstance(results[keys[1]], RemoteError)
    assert len(http.calls) == 2 and http.calls[1][1]["id"] == "NC_000002.1"


def test_context_pool_is_bounded_and_identical_downloads_are_reused(tmp_path, monkeypatch):
    client = NcbiClient(ScriptedHttp([]), tmp_path)
    keys = [(f"NC_{i:06d}.1", 1, 100) for i in range(1, 4)]
    barrier = threading.Barrier(3)
    lock = threading.Lock()
    calls, active, peak = [], 0, 0

    def fetch(*key, **kwargs):
        nonlocal active, peak
        with lock:
            calls.append(key)
            active += 1
            peak = max(peak, active)
        barrier.wait(timeout=5)
        with lock:
            active -= 1
        return ("fixture", None, "timestamp")

    monkeypatch.setattr(client, "fetch_context", fetch)
    results = list(client.context_results([*keys, keys[0]], set()))
    assert len(results) == 3 and set(calls) == set(keys) and peak == 3


def test_context_region_offset_and_truncated_records_are_rejected(tmp_path):
    raw, _ = genomic_region("NC_000001.1", 1, 169)
    raw = raw.replace("ACCESSION   NC_000001", "ACCESSION   NC_000001 REGION: 2..170")
    with pytest.raises(ValueError, match="REGION"):
        NcbiClient(ScriptedHttp([raw]), tmp_path).fetch_context("NC_000001.1", 1, 169)
    with pytest.raises(ResponseParseError):
        genbank_records(raw.rsplit("//", 1)[0])


def test_cancelled_cached_context_sends_no_request(tmp_path):
    http = ScriptedHttp([])
    http.cancel.set()
    with pytest.raises(Cancelled):
        NcbiClient(http, tmp_path).fetch_context("NC_000001.1", 1, 169)


@pytest.mark.parametrize("start,end", [("", "69"), ("0", "69"), ("10", "invalid")])
def test_invalid_nucleotide_linked_location_is_not_biological_no_mapping(start, end):
    hit = parse_blast_xml(blast_xml(), 1, 20)[0][0]
    text = HEADER + f"1\tRefSeq\tNC_000001.1\t{start}\t{end}\t+\tPROTA.1\tname\torg\t\t\n"
    with pytest.raises(ResponseParseError):
        parse_ipg(text, hit, 100)


def test_unindexed_slice_does_not_download_full_record_or_invent_length(tmp_path):
    raw, _ = genomic_region("NC_000001.1", 1, 169)
    raw = raw.replace("ACCESSION   NC_000001", "ACCESSION   NC_000001 REGION: 1..169")
    http = ScriptedHttp(['{"error":"not indexed"}', raw])
    client = NcbiClient(http, tmp_path)
    client.prefetch_nucleotide_info(["NC_000001.1"])
    result = dict(
        client.context_results([("NC_000001.1", 1, 169)], set(), {("NC_000001.1", 1, 169)})
    )
    assert result[("NC_000001.1", 1, 169)][0] == raw
    assert client.nucleotide_info("NC_000001.1", allow_unindexed=True)["length"] is None
    assert len(http.calls) == 2
    assert http.calls[-1][1]["seq_stop"] == 169


@pytest.mark.parametrize("start,header", [(800, " REGION: 800..1000"), (1, "")])
def test_unindexed_contig_edge_is_proved_by_efetch_response(tmp_path, start, header):
    raw, _ = genomic_region("NC_000002.1", start, 1000)
    raw = raw.replace("ACCESSION   NC_000002", "ACCESSION   NC_000002" + header)
    http = ScriptedHttp([raw])
    client = NcbiClient(http, tmp_path)
    text, record, _ = client.fetch_context("NC_000002.1", start, 1500, allow_clipping=True)
    assert text == raw and len(record) == 1000 - start + 1
    assert client.nucleotide_info("NC_000002.1")["length"] == 1000
    assert len(http.calls) == 1


def test_unindexed_slice_without_native_offset_is_rejected(tmp_path):
    raw, _ = genomic_region("NC_000002.1", 800, 1000)
    with pytest.raises(ValueError, match="native REGION"):
        NcbiClient(ScriptedHttp([raw]), tmp_path).fetch_context(
            "NC_000002.1", 800, 1500, allow_clipping=True
        )


def test_bad_full_batch_record_does_not_lose_another_success(tmp_path):
    first, _ = genomic_region("NC_000001.1", 1, 1000)
    short, _ = genomic_region("NC_000002.1", 800, 1000)
    second, _ = genomic_region("NC_000002.1", 1, 1000)
    http = ScriptedHttp([first + short, second])
    keys = [("NC_000001.1", 1, 1000), ("NC_000002.1", 1, 1000)]
    results = dict(NcbiClient(http, tmp_path).context_results(keys, set(keys)))
    assert results[keys[0]][0] == first and results[keys[1]][0] == second
    assert len(http.calls) == 2


def test_all_link_types_are_combined_without_a_representative_record(tmp_path):
    hit = parse_blast_xml(blast_xml(), 1, 20)[0][0]
    extra = (
        "<LinkSetDbHistory><DbTo>nuccore</DbTo><LinkName>protein_nuccore_wp</LinkName>"
        "<QueryKey>2</QueryKey></LinkSetDbHistory>"
    )
    history = HISTORY.replace("<WebEnv>", extra + "<WebEnv>")
    first, _ = genomic_region("NC_000001.1", 1, 1000)
    second, _ = genomic_region("NC_000002.1", 1, 1000)
    search = (
        "<eSearchResult><Count>2</Count><QueryKey>3</QueryKey>"
        "<WebEnv>fixture_history</WebEnv></eSearchResult>"
    )
    http = ScriptedHttp([history, search, first + second])
    rows = NcbiClient(http, tmp_path).linked_occurrences(hit, 100, {hit.accession})
    assert len(rows) == 2
    assert "linkname" not in http.calls[0][1]
    assert http.calls[1][1]["term"] == "#1 OR #2" and http.calls[2][1]["query_key"] == "3"


def test_partially_incomplete_ipg_supplements_missing_protein_without_removing_rows(run_fixture):
    _, run, client, pipeline = run_fixture
    run.settings.top_n = 1
    original = ipg_text().splitlines()[1]
    missing = "1\tRefSeq\tN/A\t\t\t\tWP_MISSING.1\tname\torg\t\t\n"
    client.fetch_ipg = lambda *_: HEADER + original + "\n" + original + "\n" + missing
    hit = parse_blast_xml(blast_xml(), 1, 20)[0][0]
    extra = parse_ipg(HEADER + ipg_text().splitlines()[3] + "\n", hit, 100)[0]
    extra.id, extra.ipg_row, extra.mapping_method = "h0001_l000001", 0, "ELink + GenBank CDS"
    requested = []
    known = parse_ipg(HEADER + original + "\n", hit, 100)[0]
    repeated = deepcopy(extra)
    repeated.id = "h0001_l000002"

    def supplement(hit, flank, members, **kwargs):
        requested.extend(kwargs["accessions"])
        return [known, extra, repeated]

    client.linked_occurrences = supplement
    pipeline.execute()
    assert requested == ["WP_MISSING.1"]
    assert len(run.occurrences) == 4 and run.status == "success"
    assert sum(o.mapping_method == "IPG" for o in run.occurrences) == 2
    assert sum(o.mapping_method == "ELink + GenBank CDS" for o in run.occurrences) == 2


def test_link_failure_preserves_already_resolved_ipg_occurrences(run_fixture):
    _, run, client, pipeline = run_fixture
    run.settings.top_n = 1
    text = (
        HEADER
        + ipg_text().splitlines()[1]
        + "\n1\tRefSeq\tN/A\t\t\t\tWP_MISSING.1\tname\torg\t\t\n"
    )
    client.fetch_ipg = lambda *_: text

    def fail(*args, **kwargs):
        raise RemoteError("Transient link error")

    client.linked_occurrences = fail
    pipeline.execute()
    assert run.hits[0].mapping_status == "retrieval_failure" and run.status == "partial_success"
    assert len(run.occurrences) == 1 and run.occurrences[0].status == "success"
