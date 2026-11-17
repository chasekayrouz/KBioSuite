"""Disposable protein search models built on the canonical annotation engine."""

import threading
import time
from dataclasses import asdict
from pathlib import Path

from .cancellation import check_cancelled, run_scope
from .config import default_root
from .engine import AnnotationEngine
from .io.fasta import parse_fasta
from .io.genbank import parse_genbank
from .scoring import DEFAULT_MIN_SCORE
from .search.cache import sequence_hash

MAX_QUERIES = 20


def parse_queries(text):
    proteins = parse_fasta(text, strict=True, allow_raw=True)
    if len(proteins) > MAX_QUERIES:
        raise ValueError(f"Quick Search supports up to {MAX_QUERIES} proteins at a time.")
    return proteins


def import_queries(text, filename=""):
    """Convert a local FASTA/GenBank file to validated, editable protein queries."""
    if not isinstance(text, str) or not isinstance(filename, str):
        raise ValueError("Choose a protein FASTA or GenBank file.")
    text = text.lstrip("\ufeff \t\r\n")
    genbank = not text.startswith(">") and (
        text.startswith("LOCUS")
        or Path(filename).suffix.lower() in {".gb", ".gbk", ".gbff", ".genbank"}
    )
    if not genbank:
        return parse_queries(text), []
    data = parse_genbank(text)
    records = []
    skipped = 0
    for protein, (ri, fi) in zip(data.proteins, data.cds_indices):
        raw = "".join(data.records[ri].features[fi].qualifiers.get("translation", []))
        if not raw.strip():
            skipped += 1
            continue
        # Validate the source translation, not its permissively normalized search copy.
        records.append(f">{protein.header}\n{raw}\n")
    if not records:
        raise ValueError("This GenBank has no CDS protein translations for Quick Search.")
    if len(records) > MAX_QUERIES:
        raise ValueError(
            f"This GenBank contains {len(records)} CDS proteins. "
            f"Quick Search supports up to {MAX_QUERIES}; use Annotate files for larger inputs."
        )
    messages = list(data.warnings)
    if skipped:
        messages.append(f"Skipped {skipped} CDS(s) without a protein translation.")
    return parse_queries("".join(records)), messages


def result_model(run):
    queries = []
    for index, annotation in enumerate(run.annotations):
        check_cancelled()
        protein = annotation.protein
        hits = sorted(
            annotation.hmm_hits + annotation.rps_hits,
            key=lambda h: (-h.score, h.evalue, h.engine, h.profile_id),
        )
        queries.append(
            {
                "id": str(index),
                "name": protein.header,
                "length": len(protein.sequence),
                "sequence": protein.sequence,
                "query_sha256": sequence_hash(protein.sequence),
                "hits": [
                    {
                        "id": f"q{index}-{h.engine}-{i}",
                        "source": "Pfam" if h.engine == "hmm" else "CDD",
                        "name": h.name,
                        "accession": h.profile_id,
                        "description": h.description,
                        "start": h.start,
                        "end": h.end,
                        "evalue": h.evalue,
                        "bitscore": h.bitscore,
                        "score": h.score,
                        "significant": h.score >= DEFAULT_MIN_SCORE,
                        "spans": h.spans or ((h.start, h.end),),
                        "alignments": [asdict(a) for a in h.alignments],
                    }
                    for i, h in enumerate(hits)
                ],
            }
        )
    provenance = {
        name: {
            "identity": entries[0]["identity"] if entries else None,
            "queries": [{k: v for k, v in entry.items() if k != "identity"} for entry in entries],
        }
        for name, entries in run.search_provenance.items()
    }
    return {
        "queries": queries,
        "warnings": run.warnings,
        "timings": run.timings,
        "database_versions": run.database_versions,
        "threads": run.threads,
        "provenance": provenance,
    }


class QuickSearchSession:
    """One serialized app session owns at most one lazy Pfam worker."""

    def __init__(self, *, pfam_backend="auto"):
        self.pfam_backend = pfam_backend
        self.engine = None
        self.identity = None
        self.lock = threading.RLock()

    def close(self):
        with self.lock:
            if self.engine is not None:
                self.engine.close()
                self.engine = None
                self.identity = None

    def get_engine(self, root, cache_dir):
        candidate = AnnotationEngine.from_config(
            root, cache_dir=cache_dir, include_alignments=True, pfam_backend=self.pfam_backend
        )
        identity = [
            [str(a.directory), a.identity] for a in (candidate.hmm_adapter, candidate.rps_adapter)
        ]
        if identity != self.identity:
            self.close()
            self.engine, self.identity = candidate, identity
        else:
            candidate.close()
        return self.engine


def search(proteins, *, root=None, cache_dir=None, threads=None, control=None, session=None):
    """Independent Python API; no GUI, output file or remote service is required."""
    owned = session is None
    session = session or QuickSearchSession()
    try:
        with session.lock, run_scope(control):
            started = time.perf_counter()
            root = Path(root or default_root())
            engine = session.get_engine(root, Path(cache_dir or root / "cache" / "search"))
            initialization = time.perf_counter() - started
            result = result_model(engine.annotate(proteins, threads=threads))
            result["timings"].update(
                initialization=initialization, wall=time.perf_counter() - started
            )
            return result
    except BaseException:
        session.close()
        raise
    finally:
        if owned:
            session.close()
