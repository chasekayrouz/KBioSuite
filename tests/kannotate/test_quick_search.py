import json
import re
import threading
from dataclasses import replace
from io import StringIO
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from kannotate import AnnotationEngine, Hit
from kannotate.gui import AnnotationJob, create_app
from kannotate.quick_search import import_queries, parse_queries, result_model
from kannotate.search.cache import CachedAdapter, file_identity, sequence_hash
from kannotate.search.common import SearchError


def test_single_multi_raw_input_and_duplicate_names():
    assert parse_queries(">enzyme\nma cdefg\n")[0].sequence == "MACDEFG"
    records = parse_queries(";comment\n>same\nMACDEFG\n>same\nBJZOUX*\n")
    assert [p.header for p in records] == ["same", "same"]
    assert records[0].id != records[1].id
    assert parse_queries("ma cdefg\nHIK")[0].header == "query_1"
    assert parse_queries("ma cdefg\nHIK")[0].sequence == "MACDEFGHIK"


@pytest.mark.parametrize("text", ["", ">\nMACDEFG", ">empty\n", "MACDEFG\n>second\nMACDEFG"])
def test_malformed_input(text):
    with pytest.raises(ValueError):
        parse_queries(text)


@pytest.mark.parametrize("bad", ["-", ".", "1", "!", "_", "é", "ß", ">"])
def test_invalid_protein_input_is_rejected_without_truncation(bad):
    with pytest.raises(ValueError, match="invalid protein character"):
        parse_queries(f">protein\nMACDEFG{bad}HIK\n")


def test_twenty_query_limit_and_fixture_integrity():
    with pytest.raises(ValueError, match="20"):
        parse_queries(">p\nMACDEFG\n" * 21)
    directory = Path(__file__).parent / "fixtures"
    proteins = parse_queries((directory / "quick-search.faa").read_text())
    manifest = json.loads((directory / "quick-search.json").read_text())
    assert len(proteins) == 20
    assert [sequence_hash(p.sequence) for p in proteins] == [
        q["sha256"] for q in manifest["queries"]
    ]


def genbank_input(translations):
    record = SeqRecord(Seq("ATGAAATTT" * len(translations)), id="test", name="test")
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(
            SimpleLocation(i * 9, i * 9 + 9),
            type="CDS",
            qualifiers={
                "locus_tag": ["same"],
                **({"translation": [sequence]} if sequence is not None else {}),
            },
        )
        for i, sequence in enumerate(translations)
    ]
    handle = StringIO()
    SeqIO.write(record, handle, "genbank")
    return handle.getvalue()


def test_genbank_import_uses_cds_translations_and_canonical_derivation(tmp_path):
    text = genbank_input(["ACDEFG*QQ", None, ""])
    path = tmp_path / "original.gbk"
    path.write_text(text)
    proteins, warnings = import_queries("\ufeff" + text, path.name)
    assert [p.sequence for p in proteins] == ["ACDEFG*QQ", "MKF"]
    assert [p.header for p in proteins] == ["test same", "test same"]
    assert proteins[0].id != proteins[1].id
    assert warnings == ["Skipped 1 CDS(s) without a protein translation."]
    assert path.read_text() == text


def test_genbank_import_supports_multiple_records_and_content_over_extension():
    text = genbank_input(["MACDEFG"]) + genbank_input(["HIKLMNP"])
    proteins, warnings = import_queries(text, "unknown.txt")
    assert [p.sequence for p in proteins] == ["MACDEFG", "HIKLMNP"]
    assert warnings == []
    proteins, warnings = import_queries(">protein\nMACDEFG\n", "misnamed.gbk")
    assert proteins == parse_queries(">protein\nMACDEFG\n") and not warnings


def test_genbank_import_rejects_over_limit_empty_and_invalid_translations():
    assert len(import_queries(genbank_input(["MACDEFG"] * 20), "input.GBFF")[0]) == 20
    with pytest.raises(ValueError, match="21 CDS proteins.*up to 20.*Annotate files"):
        import_queries(genbank_input(["MACDEFG"] * 21), "input.gbk")
    with pytest.raises(ValueError, match="no CDS protein translations"):
        import_queries(genbank_input([""]), "input.gbk")
    with pytest.raises(ValueError, match="invalid protein character"):
        import_queries(genbank_input(["ACDEFG!QQ"]), "input.gbk")
    with pytest.raises(ValueError, match="No GenBank records"):
        import_queries("not a GenBank", "input.gbk")


def test_import_api_is_local_token_protected_and_returns_editable_fasta():
    job = AnnotationJob()
    client = create_app(job).test_client()
    token = re.search(r'name="desktop-token" content="([^"]+)"', client.get("/").data.decode())[1]
    headers = {"X-Desktop-Token": token}
    data = {"text": genbank_input(["MACDEFG", ""]), "filename": "input.gbk"}
    assert client.post("/api/quick/import", json=data).status_code == 403
    assert client.post(
        "/api/quick/import", headers={**headers, "Origin": "https://example.com"}, json=data
    ).status_code == 403
    result = client.post("/api/quick/import", headers=headers, json=data)
    assert result.status_code == 200
    assert result.json == {
        "fasta": ">test same\nMACDEFG\n",
        "warnings": ["Skipped 1 CDS(s) without a protein translation."],
    }
    assert not job.snapshot()["busy"] and job.quick is None
    for invalid in [{}, {"text": 12}, {"text": "MACDEFG", "filename": []}, ["text"]]:
        assert client.post("/api/quick/import", headers=headers, json=invalid).status_code == 400


class Adapter:
    def __init__(self, engine="hmm", fail=False):
        self.engine, self.calls, self.fail = engine, [], fail

    def search(self, queries, workdir, threads):
        from kannotate.io.fasta import read_fasta

        proteins = read_fasta(queries)
        self.calls.append([p.sequence for p in proteins])
        if self.fail:
            raise SearchError("database unavailable")
        (workdir / "raw.tsv").write_text("synthetic fixture engine output")
        return [
            Hit(
                self.engine,
                p.header,
                "model",
                "domain",
                "Full description",
                1e-20,
                2,
                len(p.sequence) - 1,
                len(p.sequence),
                bitscore=75,
                spans=((2, 3), (5, len(p.sequence) - 1)),
            )
            for p in proteins
            if p.sequence != "QQQQQQQ"
        ]


def cached_engine(tmp_path, identity=None):
    hmm, rps = Adapter(), Adapter("rps")
    return (
        AnnotationEngine(
            CachedAdapter(hmm, tmp_path / "hmm", identity or {"version": 1}),
            CachedAdapter(rps, tmp_path / "rps", identity or {"version": 1}),
        ),
        hmm,
        rps,
    )


def test_batch_cache_replay_hashes_provenance_zero_hits_and_duplicate_queries(tmp_path):
    engine, hmm, rps = cached_engine(tmp_path)
    proteins = parse_queries(
        ">same\nMACDEFG\n>same\nACDEFGH\n>duplicate sequence\nMACDEFG\n>nohit\nQQQQQQQ"
    )
    first = engine.annotate(proteins, threads=4)
    assert not first.warnings
    assert [len(c) for c in hmm.calls] == [3]
    assert [len(c) for c in rps.calls] == [3]
    renamed = [replace(p, header="new header", id="same user id") for p in reversed(proteins)]
    repeat = engine.annotate(renamed, threads=1)
    assert len(hmm.calls) == len(rps.calls) == 1
    assert [a.protein.header for a in repeat.annotations] == ["new header"] * 4
    assert [a.hmm_hits for a in repeat.annotations[:-1]][0] == ()  # zero-hit record is cached
    assert all(h.query_id == "same user id" for a in repeat.annotations for h in a.hmm_hits)
    assert [a.score for a in repeat.annotations] == [a.score for a in reversed(first.annotations)]
    for entries in repeat.search_provenance.values():
        assert all(p["cached"] and Path(p["raw_output"]).is_dir() for p in entries)
        assert {p["query_sha256"] for p in entries} == {sequence_hash(p.sequence) for p in proteins}
    model = result_model(first)
    assert len({q["id"] for q in model["queries"]}) == 4
    assert len({h["id"] for q in model["queries"] for h in q["hits"]}) == 6
    assert model["queries"][0]["hits"][0]["spans"] == ((2, 3), (5, 6))
    assert model["queries"][0]["hits"][0]["bitscore"] == 75
    assert model["queries"][0]["hits"][0]["evalue"] == 1e-20
    assert {h["source"] for h in model["queries"][0]["hits"]} == {"CDD", "Pfam"}


def test_cache_batches_only_misses_and_invalidates_changed_identity(tmp_path):
    engine, hmm, _ = cached_engine(tmp_path)
    engine.annotate(parse_queries(">a\nMACDEFG"), threads=2)
    engine.annotate(parse_queries(">a\nMACDEFG\n>b\nACDEFGH"), threads=2)
    assert hmm.calls == [["MACDEFG"], ["ACDEFGH"]]
    changed, new_hmm, _ = cached_engine(tmp_path, {"version": 2})
    changed.annotate(parse_queries(">a\nMACDEFG"), threads=2)
    assert new_hmm.calls == [["MACDEFG"]]
    asset = tmp_path / "asset"
    asset.write_text("original")
    first = file_identity(asset)
    asset.write_text("changed")
    assert file_identity(asset) != first


def test_corrupt_cache_researches_and_failed_engines_never_cache(tmp_path):
    engine, hmm, _ = cached_engine(tmp_path)
    proteins = parse_queries(">a\nMACDEFG")
    engine.annotate(proteins, threads=2)
    cache = next(engine.hmm_adapter.directory.glob("*.json"))
    data = json.loads(cache.read_text())
    data["hits"][0]["start"] = 1
    cache.write_text(json.dumps(data))
    engine.annotate(proteins, threads=2)
    assert len(hmm.calls) == 2
    failure = CachedAdapter(Adapter(fail=True), tmp_path / "failure", {"version": 1})
    engine.hmm_adapter = failure
    result = engine.annotate(proteins, threads=2)
    assert "hmmscan: 1 proteins could not be searched" in result.warnings
    assert not list((tmp_path / "failure").rglob("*.json"))
    assert result.annotations[0].rps_hits


def test_unwritable_cache_retains_scientific_result(tmp_path):
    blocked = tmp_path / "not-directory"
    blocked.write_text("block cache creation")
    engine = AnnotationEngine(CachedAdapter(Adapter(), blocked, {}), Adapter("rps"))
    result = engine.annotate(parse_queries(">a\nMACDEFG"), threads=2)
    assert result.annotations[0].hmm_hits
    assert any("cache" in w for w in result.warnings)


def test_quick_api_background_job_security_and_file_results_preserved(monkeypatch):
    entered, release = threading.Event(), threading.Event()
    engine = AnnotationEngine(Adapter(), Adapter("rps"))

    def search(proteins, **kwargs):
        entered.set()
        assert release.wait(3)
        return result_model(engine.annotate(proteins, threads=2))

    monkeypatch.setattr("kannotate.gui.search", search)
    job = AnnotationJob()
    job.update(results=[{"output": "earlier-file.xlsx"}])
    client = create_app(job).test_client()
    token = re.search(r'name="desktop-token" content="([^"]+)"', client.get("/").data.decode())[1]
    headers = {"X-Desktop-Token": token}
    assert client.post("/api/quick/search", json={"text": "MACDEFG"}).status_code == 403
    assert (
        client.post("/api/quick/validate", headers=headers, json={"text": "MAC-DEF"}).status_code
        == 400
    )
    assert client.post("/api/quick/validate", headers=headers, json={"text": "MACDEFG"}).json == {
        "count": 1,
        "lengths": [7],
    }
    assert (
        client.post("/api/quick/search", headers=headers, json={"text": "MACDEFG"}).status_code
        == 202
    )
    try:
        assert entered.wait(1) and job.snapshot()["busy"]
        assert (
            client.post("/api/quick/search", headers=headers, json={"text": "MACDEFG"}).status_code
            == 400
        )
        assert job.snapshot()["results"] == [{"output": "earlier-file.xlsx"}]
        assert client.post("/api/quick/clear", headers=headers, json={}).status_code == 400
    finally:
        release.set()
    # Synchronize with the job's lock rather than asserting wall-time performance.
    import time

    deadline = time.monotonic() + 3
    while job.snapshot()["busy"] and time.monotonic() < deadline:
        time.sleep(0.01)
    response = client.get("/api/quick/results", headers=headers).json
    assert response["revision"] == 1 and len(response["result"]["queries"]) == 1
    assert response["result"]["queries"][0]["name"] == "query_1"
    assert client.post("/api/quick/clear", headers=headers, json={}).json["revision"] == 2
    assert client.get("/api/quick/results", headers=headers).json["result"] is None
    assert job.snapshot()["results"] == [{"output": "earlier-file.xlsx"}]


def test_alignment_evidence_survives_cache_and_result_model(tmp_path):
    from kannotate.models import Alignment

    class AlignedAdapter(Adapter):
        def search(self, *args):
            return [
                replace(h, alignments=(Alignment(2, 6, 1, 5, "ACD-EF", "AC-DEY", 1e-9, 35),))
                for h in super().search(*args)
            ]

    adapter = AlignedAdapter("rps")
    engine = AnnotationEngine(Adapter(), CachedAdapter(adapter, tmp_path, {"schema": 2}))
    proteins = parse_queries(">q\nMACDEFG")
    first = result_model(engine.annotate(proteins, threads=2))
    replay = result_model(engine.annotate(proteins, threads=2))
    assert first["queries"] == replay["queries"]
    assert len(adapter.calls) == 1
    hit = next(h for h in replay["queries"][0]["hits"] if h["source"] == "CDD")
    assert hit["alignments"][0]["query"] == "ACD-EF"
    assert hit["alignments"][0]["model_start"] == 1
    assert hit["alignments"][0]["evalue"] == 1e-9


def test_progress_poll_does_not_probe_resources_and_requires_token(monkeypatch):
    def probe():
        raise AssertionError("Progress polling must not resolve or invoke executables")

    monkeypatch.setattr("kannotate.gui.status", probe)
    job = AnnotationJob()
    client = create_app(job).test_client()
    token = re.search(r'name="desktop-token" content="([^"]+)"', client.get("/").data.decode())[1]
    assert client.get("/api/job").status_code == 403
    for i in range(3):
        job.update(message=f"Batch {i}")
        response = client.get("/api/job", headers={"X-Desktop-Token": token})
        assert response.status_code == 200 and response.json["message"] == f"Batch {i}"
