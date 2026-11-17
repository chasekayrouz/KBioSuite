import re
import threading
from pathlib import Path

import pytest

from kannotate.gui import AnnotationJob, create_app
from kannotate.search.common import SearchError, run_search
from kannotate.search.rpsblast import database_path


def test_launch_and_status_do_not_set_up_databases(monkeypatch):
    def unexpected(*args, **kwargs):
        raise AssertionError("setup must require an explicit action")

    monkeypatch.setattr("kannotate.gui.setup", unexpected)
    monkeypatch.setattr("kannotate.gui.status", lambda: {"executables": {}, "databases": {}})
    app = create_app(AnnotationJob())
    client = app.test_client()
    index = client.get("/")
    assert index.status_code == 200
    assert b"Kannotate" in index.data
    assert client.get("/api/status").status_code == 403
    token = re.search(r'name="desktop-token" content="([^"]+)"', index.data.decode())[1]
    headers = {"X-Desktop-Token": token}
    assert client.get("/api/status", headers=headers).status_code == 200
    assert (
        client.get("/api/status", headers={**headers, "Origin": "https://example.com"}).status_code
        == 403
    )


def test_only_one_setup_or_annotation_runs_at_a_time():
    entered, release = threading.Event(), threading.Event()
    job = AnnotationJob()

    def operation():
        entered.set()
        release.wait(3)

    job.start(operation)
    assert entered.wait(1)
    with pytest.raises(ValueError, match="current operation"):
        job.start(lambda: None)
    release.set()


def test_space_free_cdd_view_keeps_bytes_and_cleans_links(tmp_path):
    source = tmp_path / "database with spaces"
    source.mkdir()
    original = source / "Cdd.pal"
    original.write_bytes(b"public synthetic database asset")
    with database_path(str(source / "Cdd")) as prefix:
        view = Path(prefix + ".pal")
        assert not any(c.isspace() for c in prefix)
        assert view.is_symlink() and view.read_bytes() == original.read_bytes()
    assert not view.exists()
    assert original.read_bytes() == b"public synthetic database asset"


def test_rps_worker_failure_with_zero_exit_is_not_empty_success(tmp_path):
    program = tmp_path / "rpsblast"
    program.write_text(
        "#!/bin/sh\necho 'Error: [rpsblast] Chunk # 0 exit with error (3)' >&2\nexit 0\n"
    )
    program.chmod(0o755)
    with pytest.raises(SearchError, match="exit with error"):
        run_search([str(program)])


@pytest.mark.parametrize("mode", ["genbank", "csv", "xlsx", "both"])
def test_gui_batch_with_spaces_multiple_formats_and_invalid_input(tmp_path, monkeypatch, mode):
    import time

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    from kannotate import AnnotationEngine

    class NoHits:
        def search(self, *args):
            return []

    engine = AnnotationEngine(NoHits(), NoHits())
    monkeypatch.setattr("kannotate.gui.AnnotationEngine.from_config", lambda: engine)
    inputs = tmp_path / "inputs with spaces"
    inputs.mkdir()
    first, second, bad, genbank = [
        inputs / name
        for name in [
            "first control.faa",
            "second control.fasta",
            "invalid.faa",
            "cluster control.gbk",
        ]
    ]
    first.write_text(">a\nMAVKDLGATNRELVKATGQV\n")
    second.write_text(">b\nMAVKDLGATNRELVKATGQV\n")
    bad.write_text("not a FASTA")
    record = SeqRecord(Seq("ATG" * 20), id="synthetic", description="Offline GUI contract")
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(
            SimpleLocation(0, 60, strand=1),
            type="CDS",
            qualifiers={"locus_tag": ["control"], "translation": ["MAVKDLGATNRELVKATGQV"]},
        )
    ]
    SeqIO.write(record, genbank, "genbank")
    before = genbank.read_bytes()
    destination = tmp_path / "outputs with spaces"
    job = AnnotationJob()
    client = create_app(job).test_client()
    token = re.search(r'name="desktop-token" content="([^"]+)"', client.get("/").data.decode())[1]
    headers = {"X-Desktop-Token": token}
    response = client.post(
        "/api/annotate",
        headers=headers,
        json={
            "inputs": list(map(str, [first, bad, second, genbank])),
            "output": str(destination),
            "threads": 2,
            "genbank_output": mode,
        },
    )
    assert response.status_code == 202
    deadline = time.monotonic() + 5
    while job.snapshot()["busy"] and time.monotonic() < deadline:
        time.sleep(0.01)
    state = job.snapshot()
    assert not state["busy"] and not state["error"]
    assert [r["error"] is None for r in state["results"]] == [True, False, True, True]
    assert (
        sorted(p.suffix for p in destination.iterdir())
        == {
            "genbank": [".csv", ".csv", ".gbk"],
            "csv": [".csv", ".csv", ".csv"],
            "xlsx": [".csv", ".csv", ".xlsx"],
            "both": [".csv", ".csv", ".csv", ".gbk"],
        }[mode]
    )
    assert len(state["results"][-1]["outputs"]) == (2 if mode == "both" else 1)
    assert genbank.read_bytes() == before
    if mode in ("genbank", "both"):
        assert (
            SeqIO.read(destination / "cluster control.kannotate.gbk", "genbank")
            .features[0]
            .location
            == record.features[0].location
        )
    response = client.post(
        "/api/annotate",
        headers=headers,
        json={"inputs": [str(inputs / "missing.faa")], "output": str(destination), "threads": 2},
    )
    assert response.status_code == 400 and "missing" in response.json["error"]


@pytest.mark.parametrize("mode", ["invalid", None, [], {}])
def test_gui_rejects_invalid_output_format(tmp_path, mode):
    job = AnnotationJob()
    client = create_app(job).test_client()
    token = re.search(r'name="desktop-token" content="([^"]+)"', client.get("/").data.decode())[1]
    response = client.post(
        "/api/annotate",
        headers={"X-Desktop-Token": token},
        json={
            "inputs": [str(tmp_path / "source.gbk")],
            "output": str(tmp_path),
            "genbank_output": mode,
        },
    )
    assert response.status_code == 400 and "GenBank output" in response.json["error"]
    assert not job.snapshot()["busy"]


def test_quick_search_cache_defaults_to_installation_and_can_be_isolated(tmp_path, monkeypatch):
    from kannotate.quick_search import parse_queries

    seen = []
    finished = threading.Event()

    def search(proteins, **kwargs):
        seen.append(kwargs["cache_dir"])
        finished.set()
        return {"queries": [], "warnings": [], "provenance": {}, "timings": {}}

    monkeypatch.setattr("kannotate.gui.search", search)
    job = AnnotationJob()
    for cache_dir in (None, tmp_path / "acceptance-cache"):
        finished.clear()
        job.quick_cache_dir = cache_dir
        job.quick_search(parse_queries("MACDEFG"))
        assert finished.wait(5)
        for _ in range(100):
            if not job.snapshot()["busy"]:
                break
            threading.Event().wait(0.02)
    # None selects search()'s persistent <home>/cache/search; acceptance runs pass their own.
    assert seen == [None, tmp_path / "acceptance-cache"]
