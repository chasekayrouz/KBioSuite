"""Canonical annotation API: no FASTA, GenBank or workbook dependencies."""

import math
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import nullcontext
from dataclasses import replace
from itertools import batched
from pathlib import Path

from .cancellation import check_cancelled, run_scope
from .config import default_threads, load_config
from .dependencies import resolve_executable
from .models import Protein, RunResult
from .normalization import is_searchable, normalize_sequence
from .scoring import DEFAULT_MIN_SCORE
from .search.common import SearchError, write_queries
from .search.hmmer import HmmerAdapter
from .search.rpsblast import RpsBlastAdapter
from .search.scheduling import search_budgets
from .selection import select_annotation


class AnnotationEngine:
    def __init__(self, hmm_adapter, rps_adapter, *, database_versions=None):
        self.hmm_adapter = hmm_adapter
        self.rps_adapter = rps_adapter
        self.database_versions = database_versions or {}

    @classmethod
    def from_config(
        cls,
        root: Path | None = None,
        *,
        cache_dir: Path | None = None,
        include_alignments=False,
        pfam_backend="cli",
    ):
        if pfam_backend not in {"cli", "auto"}:
            raise ValueError("Pfam backend must be cli or auto")
        config = load_config(root)
        binaries = {
            name: resolve_executable(name, config["executables"].get(name))
            for name in ("hmmscan", "rpsblast")
        }
        missing = [n for n in ("hmmscan", "rpsblast") if not binaries[n]]
        if missing:
            raise RuntimeError("Missing executables: " + ", ".join(missing))
        databases = config["databases"]
        engine = cls(
            HmmerAdapter(binaries["hmmscan"], databases["pfam"]["path"]),
            RpsBlastAdapter(
                binaries["rpsblast"],
                databases["cdd"]["path"],
                databases["cdd"]["descriptions"],
                include_alignments=include_alignments,
            ),
            database_versions={k: v["version"] for k, v in databases.items()},
        )
        engine.config = config
        if pfam_backend == "auto":
            from .search.identity import reference_identity
            from .search.resident import ResidentHmmer, qualified_identity

            reference = reference_identity(databases["pfam"])
            identity = qualified_identity(engine.hmm_adapter, reference)
            if identity is not None:
                engine.hmm_adapter = ResidentHmmer(engine.hmm_adapter, identity, reference)
        if cache_dir is not None:
            from .search.cache import CachedAdapter

            engine.hmm_adapter = CachedAdapter.from_config(
                engine.hmm_adapter, "pfam", config, cache_dir
            )
            engine.rps_adapter = CachedAdapter.from_config(
                engine.rps_adapter, "cdd", config, cache_dir
            )
        return engine

    def close(self):
        for adapter in (self.hmm_adapter, self.rps_adapter):
            if hasattr(adapter, "close"):
                adapter.close()

    def annotate(
        self,
        proteins,
        *,
        threads=None,
        min_score=DEFAULT_MIN_SCORE,
        control=None,
        evidence_dir=None,
    ) -> RunResult:
        """Return one annotation per protein, preserving order and duplicate user identifiers.

        Search failures are isolated to each engine. On a failed batch, recursively split
        queries to rescue unaffected sequences. No input sequences leave this machine.
        """
        with run_scope(control) as active:
            started = time.perf_counter()
            result = None
            for index, window in enumerate(batched(proteins, 512, strict=False)):
                part = self._annotate(
                    window,
                    threads=threads,
                    min_score=min_score,
                    control=active,
                    evidence_dir=Path(evidence_dir) / str(index) if evidence_dir else None,
                )
                if result is None:
                    result = part
                else:
                    result.annotations.extend(part.annotations)
                    result.protein_count += part.protein_count
                    result.warnings.extend(part.warnings)
                    for name, elapsed in part.timings.items():
                        result.timings[name] = result.timings.get(name, 0.0) + elapsed
                    for name, entries in part.search_provenance.items():
                        result.search_provenance.setdefault(name, []).extend(entries)
            if result is None:
                result = self._annotate([], threads=threads, min_score=min_score, control=active)
            result.timings["total"] = time.perf_counter() - started
            return result

    def _annotate(self, proteins, *, threads, min_score, control, evidence_dir=None):
        started = time.perf_counter()
        for adapter in (self.hmm_adapter, self.rps_adapter):
            if hasattr(adapter, "begin_run"):
                adapter.begin_run()
        threads = default_threads() if threads is None else threads
        if threads < 1 or not math.isfinite(min_score) or min_score < 0:
            raise ValueError("threads must be positive; min_score must be finite and nonnegative")
        originals = list(proteins)
        normalized = []
        for i, p in enumerate(originals):
            control.check()
            normalized.append(
                Protein(f"q{i:08d}", p.header, normalize_sequence(p.sequence).sequence)
            )
        searchable = [p for p in normalized if is_searchable(p.sequence)]
        warnings = []
        timings = {"hmmscan": 0.0, "rpsblast": 0.0}

        def execute(adapter, allocated):
            with control.activate():
                return search(adapter, allocated)

        def search(adapter, allocated):
            failures = []
            calls = 0
            start = time.perf_counter()
            if evidence_dir is not None:
                directory = evidence_dir / ("hmm" if adapter is self.hmm_adapter else "rps")
                directory.mkdir(parents=True, exist_ok=False)
            workspace = (
                nullcontext(directory)
                if evidence_dir is not None
                else tempfile.TemporaryDirectory(prefix="kannotate-")
            )
            with workspace as temp:
                base = Path(temp)

                def attempt(items):
                    control.check()
                    nonlocal calls
                    calls += 1
                    directory = base / str(calls)
                    directory.mkdir()
                    query = directory / "queries.faa"
                    write_queries(items, query)
                    try:
                        hits = adapter.search(query, directory, allocated)
                        control.check()
                        failures.extend(getattr(hits, "failed_queries", ()))
                        return hits
                    except (SearchError, ValueError, OSError) as exc:
                        control.check()
                        # Stop systemic errors instead of retrying every protein.
                        systemic = any(
                            s in str(exc).lower()
                            for s in ("database", "no such file", "permission denied", "timed out")
                        )
                        if len(items) == 1 or systemic:
                            failures.extend(p.id for p in items)
                            return []
                        half = len(items) // 2
                        return attempt(items[:half]) + attempt(items[half:])

                hits = attempt(searchable) if searchable else []
            return hits, time.perf_counter() - start, len(failures)

        if searchable:
            if threads >= 2:
                # Larger batches benefit from multiple two-thread HMMER processes.
                # CachedAdapter partitions only its actual unique cache misses.
                count = len({p.sequence for p in searchable})
                hmm_threads, rps_threads = search_budgets(threads, count)
                with ThreadPoolExecutor(max_workers=2) as pool:
                    futures = {
                        pool.submit(execute, self.hmm_adapter, hmm_threads): 0,
                        pool.submit(execute, self.rps_adapter, rps_threads): 1,
                    }
                    results = [None, None]
                    try:
                        for future in as_completed(futures):
                            results[futures[future]] = future.result()
                    except BaseException:
                        # Also release the other worker on Ctrl-C/unexpected adapter failure.
                        control.request_stop()
                        raise
            else:
                results = [
                    execute(self.hmm_adapter, 1),
                    execute(self.rps_adapter, 1),
                ]
        else:
            results = [([], 0.0, 0), ([], 0.0, 0)]
        by_engine = []
        for name, (hits, elapsed, failures) in zip(("hmmscan", "rpsblast"), results):
            timings[name] = elapsed
            if failures:
                warnings.append(f"{name}: {failures} proteins could not be searched")
            mapping = {}
            for hit in hits:
                check_cancelled()
                mapping.setdefault(hit.query_id, []).append(hit)
            by_engine.append(mapping)
        annotations = []
        for original, p in zip(originals, normalized):
            control.check()
            selected = select_annotation(
                p, by_engine[0].get(p.id, []), by_engine[1].get(p.id, []), min_score
            )
            annotations.append(
                replace(
                    selected,
                    protein=Protein(original.id, original.header, p.sequence),
                    hmm_hits=tuple(replace(h, query_id=original.id) for h in selected.hmm_hits),
                    rps_hits=tuple(replace(h, query_id=original.id) for h in selected.rps_hits),
                )
            )
        timings["total"] = time.perf_counter() - started
        provenance = {}
        for name, adapter in (("hmmscan", self.hmm_adapter), ("rpsblast", self.rps_adapter)):
            warnings.extend(getattr(adapter, "warnings", []))
            provenance[name] = list(getattr(adapter, "provenance", {}).values())
        return RunResult(
            annotations, timings, warnings, threads, self.database_versions.copy(), provenance
        )
