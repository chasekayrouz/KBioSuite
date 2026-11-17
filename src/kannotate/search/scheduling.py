"""Bounded whole-query scheduling; every query still sees the complete database."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import replace

from ..cancellation import current_control, run_scope
from ..io.fasta import read_fasta
from .common import SearchError, write_queries


def search_budgets(threads, count):
    hmm = min(8, 2 * max(1, threads // 3)) if threads >= 6 and count >= 16 else max(1, threads // 2)
    return hmm, threads - hmm


class SearchHits(list):
    def __init__(self, hits=(), failed_queries=()):
        super().__init__(hits)
        self.failed_queries = tuple(failed_queries)


def systemic_error(error):
    return any(
        word in str(error).lower()
        for word in ("database", "no such file", "permission denied", "timed out")
    )


def hmmer_search(search, queries, workdir, threads):
    proteins = [replace(p, id=p.header) for p in read_fasta(queries)]
    workers = (
        min(4, threads // 2, len(proteins) // 4) if threads >= 4 and len(proteins) >= 16 else 1
    )
    per_worker = 2 if workers > 1 else threads
    chunks = [[] for _ in range(workers)]
    lengths = [0] * workers
    for protein in sorted(proteins, key=lambda p: len(p.sequence), reverse=True):
        index = min(range(workers), key=lambda i: lengths[i])
        chunks[index].append(protein)
        lengths[index] += len(protein.sequence)
    with run_scope(current_control()) as control:

        def execute(index, items):
            with control.activate():
                calls = 0

                def attempt(batch):
                    nonlocal calls
                    control.check()
                    calls += 1
                    directory = workdir / f"worker-{index}" / f"attempt-{calls}"
                    directory.mkdir(parents=True)
                    query = directory / "queries.faa"
                    write_queries(batch, query)
                    try:
                        hits = search(query, directory, per_worker)
                        control.check()
                        return SearchHits(hits)
                    except (SearchError, ValueError, OSError) as error:
                        control.check()
                        if len(batch) == 1 or systemic_error(error):
                            return SearchHits(failed_queries=[p.id for p in batch])
                        half = len(batch) // 2
                        left, right = attempt(batch[:half]), attempt(batch[half:])
                        return SearchHits(left + right, left.failed_queries + right.failed_queries)

                return attempt(items) if items else SearchHits()

        results = [None] * workers
        if workers == 1:
            return execute(0, chunks[0])
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(execute, i, chunk): i for i, chunk in enumerate(chunks)}
            try:
                for future in as_completed(futures):
                    results[futures[future]] = future.result()
            except BaseException:
                control.request_stop()
                raise
        return SearchHits(
            [hit for result in results for hit in result],
            [query for result in results for query in result.failed_queries],
        )
