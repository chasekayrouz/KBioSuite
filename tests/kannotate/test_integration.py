"""Opt-in checks against the locally installed, comprehensive databases."""

import os
from pathlib import Path

import pytest

from kannotate import AnnotationEngine, Protein
from kannotate.io.fasta import read_fasta

pytestmark = pytest.mark.skipif(
    os.environ.get("KANNOTATE_INTEGRATION") != "1", reason="set KANNOTATE_INTEGRATION=1 after setup"
)


def test_real_tools_and_deterministic_annotation_across_thread_counts():
    source = Path(__file__).parents[2] / "src/kannotate/resources/smoke.faa"
    proteins = read_fasta(source) + [Protein("empty", "empty", ""), Protein("short", "short", "M")]
    engine = AnnotationEngine.from_config()
    first = engine.annotate(proteins, threads=1)
    parallel = engine.annotate(proteins, threads=4)
    assert not first.warnings and not parallel.warnings
    assert first.annotations == parallel.annotations
    result = first.annotations[0]
    assert "elongation factor" in result.product.lower()
    assert len(result.hmm_hits) >= 2 and len(result.rps_hits) >= 2
    for hits in (result.hmm_hits, result.rps_hits):
        assert all(h.query_id == result.protein.id for h in hits)
        assert len({h.profile_id for h in hits}) == len(hits)
        assert all(1 <= h.start <= h.end <= len(result.protein.sequence) for h in hits)
    assert all(a.product == "hypothetical protein" for a in first.annotations[1:])


def test_real_quick_search_matches_annotation_pipeline_and_cache(tmp_path):
    from kannotate.quick_search import parse_queries, result_model, search

    fixture = Path(__file__).parent / "fixtures/quick-search.faa"
    proteins = parse_queries(fixture.read_text())[:5]
    baseline = result_model(AnnotationEngine.from_config().annotate(proteins, threads=4))
    first = search(proteins, cache_dir=tmp_path / "cache", threads=4)
    repeat = search(proteins, cache_dir=tmp_path / "cache", threads=1)
    assert first["queries"] == repeat["queries"]
    # The additional display evidence must not change file-workflow annotation.
    import copy

    without_alignments = copy.deepcopy(first["queries"])
    for query in without_alignments:
        for hit in query["hits"]:
            hit["alignments"] = []
    assert baseline["queries"] == without_alignments
    for query in first["queries"]:
        for hit in query["hits"]:
            assert bool(hit["alignments"]) == (hit["source"] == "CDD")
            for alignment in hit["alignments"]:
                assert (
                    alignment["query"].replace("-", "")
                    == query["sequence"][alignment["query_start"] - 1 : alignment["query_end"]]
                )
                assert len(alignment["model"].replace("-", "")) == (
                    alignment["model_end"] - alignment["model_start"] + 1
                )
    assert not first["warnings"] and not repeat["warnings"]
    assert not first["queries"][4]["hits"]
    assert {h["name"] for h in first["queries"][0]["hits"]} == {"MbtH"}
    assert {h["name"] for h in first["queries"][2]["hits"] if h["source"] == "Pfam"} >= {
        "AMP-binding",
        "PP-binding",
        "Condensation",
    }
    assert all(q["cached"] for p in repeat["provenance"].values() for q in p["queries"])


def test_real_parallel_workers_match_full_database_search_and_cache(tmp_path):
    from kannotate.search.hmmer import HmmerAdapter

    proteins = read_fasta(Path(__file__).parent / "fixtures/quick-search.faa")
    baseline = AnnotationEngine.from_config(include_alignments=True)

    class Unpartitioned(HmmerAdapter):
        search = HmmerAdapter._search

    original = baseline.hmm_adapter
    baseline.hmm_adapter = Unpartitioned(original.executable, original.database)
    expected = baseline.annotate(proteins, threads=13)
    engine = AnnotationEngine.from_config(cache_dir=tmp_path / "cache", include_alignments=True)
    actual = engine.annotate(proteins, threads=13)
    assert not expected.warnings and not actual.warnings
    assert expected.annotations == actual.annotations
    partial = AnnotationEngine.from_config(cache_dir=tmp_path / "partial", include_alignments=True)
    partial.annotate(proteins[:5], threads=13)
    mixed = partial.annotate(list(reversed(proteins)), threads=13)
    assert not mixed.warnings and mixed.annotations == list(reversed(expected.annotations))
    for entries in mixed.search_provenance.values():
        assert any(entry["cached"] for entry in entries)
        assert any(not entry["cached"] for entry in entries)
    # Reordered cache hits and duplicates retain alignments and user IDs.
    extra = Protein("extra", "public duplicate", proteins[0].sequence)
    replay = engine.annotate(list(reversed(proteins)) + [extra], threads=13)
    from dataclasses import replace

    assert replay.annotations[:-1] == list(reversed(actual.annotations))
    assert replay.annotations[-1].product == actual.annotations[0].product
    assert replay.annotations[-1].rps_hits == tuple(
        replace(h, query_id="extra") for h in actual.annotations[0].rps_hits
    )
    assert all(
        entry["cached"] for entries in replay.search_provenance.values() for entry in entries
    )


def test_real_resident_session_cache_stop_and_repeated_close(tmp_path, monkeypatch):
    import threading
    import time
    from concurrent.futures import ThreadPoolExecutor

    import psutil

    from kannotate.cancellation import RunCancelled, RunControl
    from kannotate.quick_search import QuickSearchSession, search
    from kannotate.search.resident import ResidentHmmer

    proteins = read_fasta(Path(__file__).parent / "fixtures/quick-search.faa")
    session = QuickSearchSession()
    cache = tmp_path / "cache"
    engine = session.get_engine(None, cache)
    if not isinstance(engine.hmm_adapter.adapter, ResidentHmmer):
        session.close()
        pytest.skip("This runtime has no qualified resident backend")
    try:
        first = search(proteins, cache_dir=cache, threads=13, session=session)
        assert not first["warnings"]
        adapter = session.engine.hmm_adapter.adapter
        pid = adapter.process.pid
        warm = [Protein(p.id, p.header, p.sequence + "A") for p in proteins]
        result = search(warm, cache_dir=cache, threads=13, session=session)
        assert not result["warnings"] and adapter.process.pid == pid
        assert adapter.last_call["load_seconds"] == 0
        before = {str(p): p.read_bytes() for p in cache.rglob("*") if p.is_file()}
        control = RunControl()
        pending = [Protein(p.id, p.header, p.sequence * 4 + "C") for p in proteins]
        entered = threading.Event()
        original_request = adapter._request

        def request(value):
            entered.set()
            return original_request(value)

        monkeypatch.setattr(adapter, "_request", request)
        with ThreadPoolExecutor(max_workers=1) as pool:
            future = pool.submit(
                search, pending, cache_dir=cache, threads=13, control=control, session=session
            )
            assert entered.wait(3)
            control.request_stop()
            with pytest.raises(RunCancelled):
                future.result(timeout=5)
        assert not psutil.pid_exists(pid)
        assert {str(p): p.read_bytes() for p in cache.rglob("*") if p.is_file()} == before
        repeat = search(proteins, cache_dir=cache, threads=13, session=session)
        assert repeat["queries"] == first["queries"]
        assert session.engine.hmm_adapter.adapter.process is None
    finally:
        session.close()
    for index in range(3):
        session = QuickSearchSession()
        try:
            search(proteins, cache_dir=tmp_path / f"cycle-{index}", threads=3, session=session)
            pid = session.engine.hmm_adapter.adapter.process.pid
        finally:
            session.close()
        assert not psutil.pid_exists(pid)

    session = QuickSearchSession()
    control = RunControl()
    with ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(
            search,
            proteins,
            cache_dir=tmp_path / "cold-stop",
            threads=13,
            control=control,
            session=session,
        )
        deadline = time.monotonic() + 3
        child = None
        while time.monotonic() < deadline:
            if session.engine is not None:
                child = session.engine.hmm_adapter.adapter.process
                if child is not None:
                    break
            time.sleep(0.01)
        assert child is not None
        control.request_stop()
        with pytest.raises(RunCancelled):
            future.result(timeout=5)
    assert child.poll() is not None and session.engine is None
    assert not list((tmp_path / "cold-stop").rglob("*.json"))


def test_real_checkpoint_resume_preserves_all_source_rows_and_output_bytes(tmp_path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    from kannotate.batch import annotate_batch
    from kannotate.cancellation import RunCancelled, RunControl

    proteins = read_fasta(Path(__file__).parent / "fixtures/quick-search.faa")[:5]
    records = []
    for index, count in enumerate([260, 260, 0, 513, 2]):
        record = SeqRecord(
            Seq("ATG" * 12), id=f"public_{index}", annotations={"molecule_type": "DNA"}
        )
        record.features = [
            SeqFeature(
                SimpleLocation(0, 21, strand=-1 if i % 2 else 1),
                type="CDS",
                qualifiers={
                    "locus_tag": [f"000{i}"],
                    "product": ['Original, "β" product'],
                    **({"translation": [proteins[i].sequence]} if i < 5 else {"pseudo": [""]}),
                },
            )
            for i in range(count)
        ]
        records.append(record)
    source = tmp_path / "original.public.controls.gbff"
    SeqIO.write(records, source, "genbank")
    original = source.read_bytes()
    engine = AnnotationEngine.from_config()
    kwargs = {"threads": 13, "min_score": 6, "genbank_output": "both", "retain_annotations": False}
    baseline = annotate_batch([source], tmp_path / "baseline", engine, **kwargs)[0]
    assert not baseline.error and not baseline.warnings
    state = tmp_path / "checkpoints"
    control = RunControl()

    def stop(path, count):
        if count:
            control.request_stop()

    with pytest.raises(RunCancelled):
        annotate_batch(
            [source],
            tmp_path / "out",
            engine,
            checkpoint_dir=state,
            control=control,
            on_progress=stop,
            **kwargs,
        )
    before = {str(p): p.read_bytes() for p in state.glob("file-*/group-*/receipt.json")}
    assert len(before) == 1
    resumed = annotate_batch(
        [source], tmp_path / "out", engine, checkpoint_dir=state, resume=True, **kwargs
    )[0]
    assert not resumed.error and not resumed.warnings and resumed.run.protein_count == 1035
    assert [p.read_bytes() for p in resumed.outputs] == [p.read_bytes() for p in baseline.outputs]
    assert all(Path(p).read_bytes() == value for p, value in before.items())
    assert source.read_bytes() == original
