import json
import threading
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from pathlib import Path

import pytest

from kannotate import AnnotationEngine, Hit, Protein
from kannotate.cancellation import RunCancelled, RunControl, check_cancelled
from kannotate.io.fasta import read_fasta
from kannotate.search.cache import CachedAdapter, sequence_hash
from kannotate.search.common import SearchError
from kannotate.search.hmmer import HmmerAdapter
from kannotate.search.rpsblast import load_descriptions
from kannotate.search.scheduling import SearchHits


class FakeHmmer(HmmerAdapter):
    def __init__(self):
        super().__init__("fixture", "fixture")
        self.calls = []

    def _search(self, queries, workdir, threads):
        proteins = read_fasta(queries)
        self.calls.append((tuple(p.header for p in proteins), threads))
        if any(p.sequence == "QQQQQQ" for p in proteins):
            raise SearchError("unusable sequence")
        (workdir / "hmmer.domtbl").write_text("complete synthetic raw evidence")
        return [
            Hit("hmm", p.header, "P1", "family", "description", 1e-20, 1, 6, len(p.sequence))
            for p in proteins
        ]


class NoHits:
    def search(self, *args):
        return []


def proteins(count=24):
    return [Protein(str(i), "duplicate name", "MACDEFG" + "A" * i) for i in range(count)]


def test_shards_preserve_order_and_only_retry_failed_work():
    hmm = FakeHmmer()
    values = proteins()
    values[3] = replace(values[3], sequence="QQQQQQ")
    run = AnnotationEngine(hmm, NoHits()).annotate(values, threads=13)
    assert [a.protein.id for a in run.annotations] == [p.id for p in values]
    assert run.warnings == ["hmmscan: 1 proteins could not be searched"]
    assert not run.annotations[3].hmm_hits
    assert all(a.hmm_hits for i, a in enumerate(run.annotations) if i != 3)
    assert all(threads == 2 for _, threads in hmm.calls)
    initial = [ids for ids, _ in hmm.calls if len(ids) == 6]
    assert len(initial) == 4
    assert all(sum(ids == observed for observed, _ in hmm.calls) == 1 for ids in initial)


def test_cached_nested_artifacts_are_retained_and_deletion_invalidates_entry(tmp_path):
    adapter = FakeHmmer()
    cached = CachedAdapter(adapter, tmp_path / "cache", {"fixture": 1})
    engine = AnnotationEngine(cached, NoHits())
    first = engine.annotate(proteins(), threads=13)
    entry = json.loads(
        (cached.directory / f"{sequence_hash(proteins()[0].sequence)}.json").read_text()
    )
    raw = Path(entry["raw_output"])
    tables = list(raw.rglob("hmmer.domtbl"))
    assert len(tables) == 4
    count = len(adapter.calls)
    repeat = engine.annotate(list(reversed(proteins())), threads=13)
    assert repeat.annotations == list(reversed(first.annotations))
    assert len(adapter.calls) == count
    tables[0].unlink()
    repaired = engine.annotate(proteins(), threads=13)
    assert repaired.annotations == first.annotations
    assert len(adapter.calls) > count


def test_failed_cached_query_is_never_replayed_as_zero_hit(tmp_path):
    class Partial:
        calls = 0

        def search(self, queries, directory, threads):
            self.calls += 1
            items = read_fasta(queries)
            (directory / "raw.tsv").write_text("partial batch")
            return SearchHits(failed_queries=[items[0].header])

    adapter = Partial()
    cached = CachedAdapter(adapter, tmp_path, {"fixture": 1})
    engine = AnnotationEngine(cached, NoHits())
    duplicate = [Protein("same", "a", "MACDEFG"), Protein("same", "b", "MACDEFG")]
    for _ in range(2):
        result = engine.annotate(duplicate, threads=2)
        assert result.warnings == ["hmmscan: 2 proteins could not be searched"]
    assert adapter.calls == 2
    assert not (cached.directory / f"{sequence_hash('MACDEFG')}.json").exists()


def test_all_worker_contexts_observe_stop_before_cleanup():
    lock = threading.Lock()
    entered = threading.Event()
    active = 0

    class Blocked(FakeHmmer):
        def _search(self, *args):
            nonlocal active
            with lock:
                active += 1
                if active == 4:
                    entered.set()
            try:
                while True:
                    check_cancelled()
                    threading.Event().wait(0.01)
            finally:
                with lock:
                    active -= 1

    control = RunControl()
    with ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(
            AnnotationEngine(Blocked(), NoHits()).annotate, proteins(), threads=13, control=control
        )
        try:
            assert entered.wait(3)
            control.request_stop()
            with pytest.raises(RunCancelled):
                future.result(timeout=3)
        finally:
            control.request_stop()
    assert active == 0


def test_description_metadata_cache_tracks_file_identity_and_is_immutable(tmp_path):
    path = tmp_path / "descriptions.tsv"
    path.write_text("1\tP1\tfamily\tfirst description\n")
    first = load_descriptions(path)
    assert load_descriptions(path) is first
    with pytest.raises(TypeError):
        first["1"] = ("changed",)
    path.write_text("1\tP1\tfamily\tsecond description\n")
    assert load_descriptions(path)["1"][2] == "second description"
