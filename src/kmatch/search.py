"""Local Foldseek search retaining alignments and target→query transformations."""

import csv
import hashlib
import math
import subprocess
import time
from pathlib import Path

from .database import database_lease, database_status
from .models import DependencyStatus, Progress, SearchResult, write_json
from .process import run_process
from .target_traces import retain_target_traces

FIELDS = (
    "query",
    "target",
    "evalue",
    "prob",
    "bits",
    "fident",
    "alnlen",
    "qlen",
    "tlen",
    "qcov",
    "tcov",
    "qtmscore",
    "ttmscore",
    "alntmscore",
    "lddt",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "qaln",
    "taln",
    "u",
    "t",
)


def hit_rankings(hits, limit=50):
    """Independent display sets; indices always address the unchanged result list."""
    indices = range(len(hits))
    return {
        "qtm": sorted(
            indices, key=lambda i: (-hits[i].query_tm, hits[i].evalue, -hits[i].bits, i)
        )[:limit],
        "significance": sorted(
            indices, key=lambda i: (-hits[i].probability, hits[i].evalue, -hits[i].bits, i)
        )[:limit],
    }


def retained_hit_indices(rankings):
    """Deduplicate overlapping views while preserving saved hit-index order."""
    return sorted({index for ranking in rankings.values() for index in ranking})


def finite(value, field):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"Foldseek {field} is not finite.")
    return number


def parse_foldseek(path: Path):
    results = []
    with path.open(newline="") as stream:
        for line_number, row in enumerate(csv.reader(stream, delimiter="\t"), 1):
            if not row:
                continue
            if len(row) != len(FIELDS):
                raise ValueError(
                    f"Foldseek row {line_number}: expected {len(FIELDS)} fields, "
                    f"received {len(row)}."
                )
            raw = dict(zip(FIELDS, row, strict=True))
            try:
                u = [finite(v, "u") for v in raw["u"].split(",")]
                t = [finite(v, "t") for v in raw["t"].split(",")]
                if len(u) != 9 or len(t) != 3:
                    raise ValueError(
                        "Each hit must retain a 3×3 rotation and 3-vector translation."
                    )
                floats = {
                    key: finite(raw[key], key)
                    for key in (
                        "evalue",
                        "prob",
                        "bits",
                        "fident",
                        "qcov",
                        "tcov",
                        "qtmscore",
                        "ttmscore",
                        "alntmscore",
                    )
                }
                ints = {
                    key: int(raw[key])
                    for key in ("alnlen", "qlen", "tlen", "qstart", "qend", "tstart", "tend")
                }
                if floats["evalue"] < 0 or any(
                    not 0 <= floats[key] <= 1.0001
                    for key in (
                        "prob",
                        "fident",
                        "qcov",
                        "tcov",
                        "qtmscore",
                        "ttmscore",
                        "alntmscore",
                    )
                ):
                    raise ValueError("Invalid Foldseek metric range.")
                if ints["alnlen"] <= 0 or ints["qlen"] <= 0 or ints["tlen"] <= 0:
                    raise ValueError("Invalid Foldseek alignment or structure length.")
                for prefix in ("q", "t"):
                    if (
                        not 1
                        <= ints[prefix + "start"]
                        <= ints[prefix + "end"]
                        <= ints[prefix + "len"]
                    ):
                        raise ValueError("Invalid Foldseek alignment coordinates.")
                if len(raw["qaln"]) != len(raw["taln"]) or len(raw["qaln"]) != ints["alnlen"]:
                    raise ValueError("Foldseek aligned sequences do not match alignment length.")
                lddt = (
                    None
                    if raw["lddt"].lower() in ("", "nan", "na")
                    else finite(raw["lddt"], "lddt")
                )
                result = SearchResult(
                    raw["query"],
                    raw["target"],
                    floats["evalue"],
                    floats["prob"],
                    floats["bits"],
                    floats["fident"],
                    ints["alnlen"],
                    ints["qlen"],
                    ints["tlen"],
                    floats["qcov"],
                    floats["tcov"],
                    floats["qtmscore"],
                    floats["ttmscore"],
                    floats["alntmscore"],
                    lddt,
                    ints["qstart"],
                    ints["qend"],
                    ints["tstart"],
                    ints["tend"],
                    raw["qaln"],
                    raw["taln"],
                    [u[0:3], u[3:6], u[6:9]],
                    t,
                    raw,
                )
                results.append(result)
            except (ValueError, KeyError) as exc:
                raise ValueError(f"Foldseek row {line_number}: {exc}") from exc
    # qTM emphasizes the extent to which a hit accounts for the whole query.
    return sorted(results, key=lambda r: (-r.query_tm, r.evalue, -r.bits))


class FoldseekBackend:
    def __init__(self, config):
        self.config = config

    def executable_status(self):
        if not self.config.foldseek.is_file():
            return DependencyStatus(
                "Foldseek", False, "Foldseek is missing. Install it with Setup."
            )
        try:
            version = subprocess.check_output(
                [str(self.config.foldseek), "version"],
                text=True,
                stderr=subprocess.DEVNULL,
                timeout=5,
            ).strip()
        except (OSError, subprocess.SubprocessError):
            return DependencyStatus(
                "Foldseek", False, "Foldseek cannot execute. Reinstall it with Setup."
            )
        return DependencyStatus("Foldseek", True, "Local Foldseek executable is ready.", version)

    def status(self):
        executable = self.executable_status()
        if not executable.ready:
            return executable
        database = database_status(self.config.pdb_prefix)
        return DependencyStatus(
            "Foldseek / PDB100", database.ready, database.message, executable.version
        )

    def search(self, prediction, output, cancel, progress):
        # Hold through trace retention so an updater cannot retire files still needed here.
        with database_lease(self.config, cancel):
            return self._search(prediction, output, cancel, progress)

    def _search(self, prediction, output, cancel, progress):
        started = time.monotonic()
        state = self.status()
        if not state.ready:
            raise RuntimeError(state.message)
        # Pin this run to one generation even if another process updates `current`.
        database = self.config.pdb_prefix.resolve()
        database_version = database_status(database).version
        query_hash = hashlib.sha256(prediction.structure_path.read_bytes()).hexdigest()
        output.mkdir(parents=True, exist_ok=True)
        tsv, work = output / "foldseek.tsv", output / "work"
        command = [
            self.config.foldseek,
            "easy-search",
            prediction.structure_path,
            database,
            tsv,
            work,
            "--format-output",
            ",".join(FIELDS),
            "--alignment-type",
            "2",
            "--exact-tmscore",
            "1",
            "--threads",
            str(self.config.threads),
            "-s",
            "9.5",
            "--max-seqs",
            "1000",
            "-e",
            "10",
            "--remove-tmp-files",
            "1",
        ]
        progress(Progress("search", "Searching the local PDB100 structure database…"))
        runtime = run_process(
            command,
            output / "foldseek.log",
            cancel,
            progress,
            "search",
            timeout=900,
            scratch_root=self.config.cache / "process",
        )
        results = parse_foldseek(tsv)
        if hashlib.sha256(prediction.structure_path.read_bytes()).hexdigest() != query_hash:
            raise RuntimeError("The query coordinates changed during the Foldseek search.")
        rankings = hit_rankings(results)
        candidates = [results[i] for i in retained_hit_indices(rankings)]
        traces = retain_target_traces(self.config, database, candidates, output, cancel, progress)
        write_json(
            output / "search.json",
            {
                "schema_version": 3,
                "foldseek_version": state.version,
                "database": str(database),
                "database_version": database_version,
                "query_structure_sha256": query_hash,
                "runtime_seconds": time.monotonic() - started,
                "foldseek_runtime_seconds": runtime,
                "target_traces": traces,
                "retained_target_count": len(traces),
                "command": list(map(str, command)),
                "rotation_convention": (
                    "raw Foldseek u/t, row-major; target coordinates → query coordinates"
                ),
                "ranking": "query TM descending, E-value ascending, bit score descending",
                "hit_rankings": rankings,
                "fields": FIELDS,
                "hits": [r.__dict__ for r in results],
            },
        )
        # Keep full saved-result indices stable for clients and alignment export.
        return results
