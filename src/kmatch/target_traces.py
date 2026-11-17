"""Small, run-owned target coordinates survive retirement of the source database."""

import hashlib
import re
import shutil
import tempfile
from pathlib import Path

import gemmi

from .models import Progress
from .process import run_process


def trace_digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def retain_target_traces(config, database, hits, output, cancel, progress):
    """Snapshot displayed hits atomically; never retain Foldseek's symlinked subset DB."""
    from .alignment import ca_chain

    destination = output / "target-traces"
    if destination.exists():
        raise RuntimeError("Target traces already exist; use a new search output directory.")
    records = {}
    if not hits:
        return records
    targets = {hit.target for hit in hits}
    if any(not re.fullmatch(r"[A-Za-z0-9_.-]+", t) or t in {".", ".."} for t in targets):
        raise ValueError("Unsafe target identifier in Foldseek results.")
    progress(Progress("search", "Saving hit coordinates for future alignment exports…"))
    with tempfile.TemporaryDirectory(prefix=".traces-", dir=output) as temporary:
        root = Path(temporary)
        selection, subset, traces = root / "selected.txt", root / "subset", root / "coordinates"
        selection.write_text("\n".join(sorted(targets)) + "\n")
        for command in (
            [config.foldseek, "createsubdb", selection, database, subset, "--id-mode", "1"],
            [
                config.foldseek,
                "convert2pdb",
                subset,
                traces,
                "--pdb-output-format",
                "1",
                "--pdb-output-mode",
                "1",
                "--threads",
                "1",
            ],
        ):
            run_process(
                command + ["-v", "1"],
                output / "foldseek.log",
                cancel,
                progress,
                "target retention",
                timeout=120,
                scratch_root=config.cache / "process",
            )
        keys = set(
            line.split()[0] for line in subset.with_suffix(".index").read_text().splitlines()
        )
        identifiers = {}
        with Path(str(database) + ".lookup").open() as stream:
            for line in stream:
                key, identifier, *_ = line.split()
                if key in keys and identifier in targets:
                    if identifier in identifiers:
                        raise RuntimeError(
                            "A target identifier resolved to multiple database entries."
                        )
                    identifiers[identifier] = key
        if set(identifiers) != targets or len(keys) != len(targets):
            raise RuntimeError("Could not retain exactly one coordinate trace per returned hit.")
        for hit in hits:
            cancel.check()
            path = traces / (hit.target + ".cif")
            structure = gemmi.read_structure(str(path))
            _, sequence, _ = ca_chain(structure, "Retained target")
            observed = hit.target_alignment.replace("-", "").upper()
            if (
                len(sequence) != hit.target_length
                or sequence[hit.target_start - 1 : hit.target_end] != observed
                or any(a.name != "CA" for c in structure[0] for r in c for a in r)
            ):
                raise RuntimeError(
                    "Retained target coordinates disagree with the Foldseek alignment."
                )
            records[hit.target] = {
                "file": "target-traces/" + path.name,
                "sha256": trace_digest(path),
                "database_key": identifiers[hit.target],
            }
        cancel.check()
        # Copy only the verified individual structures, never subset links into the old database.
        retained = root / "retained"
        retained.mkdir()
        for record in records.values():
            name = Path(record["file"]).name
            shutil.copy2(traces / name, retained / name)
        retained.replace(destination)
    return records


def load_target_trace(search_directory, metadata, hit):
    record = metadata.get("target_traces", {}).get(hit.target)
    if record is None:
        return None
    directory = Path(search_directory).resolve()
    path = directory / record["file"]
    if not path.resolve().is_relative_to(directory / "target-traces"):
        raise ValueError("The saved target trace path is outside this run.")
    if not path.is_file() or trace_digest(path) != record["sha256"]:
        raise RuntimeError(
            "This run's retained target trace is missing or changed. Restore the complete run "
            "folder or run a new search; Kmatch will not substitute a different PDB conformation."
        )
    return gemmi.read_structure(str(path)), record["database_key"]
