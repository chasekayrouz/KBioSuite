"""Optional local sequence cache around the canonical profile-search adapters."""

import hashlib
import json
import shutil
import tempfile
from dataclasses import asdict, replace
from pathlib import Path

from ..cancellation import check_cancelled, current_control
from ..io.fasta import read_fasta
from ..models import Alignment, Hit
from .common import write_queries
from .hmmer import HMM_ARGUMENTS
from .identity import file_identity, fingerprint, reference_identity, tool_version
from .rpsblast import RPS_ARGUMENTS
from .scheduling import SearchHits

CACHE_SCHEMA = 4  # Backend implementation/options join the raw-evidence identity.


def sequence_hash(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _raw_manifest(directory):
    result = {}
    for path in sorted(directory.rglob("*")):
        check_cancelled()
        if path.is_symlink():
            raise ValueError("Raw search evidence must not contain symbolic links")
        if path.is_file():
            result[str(path.relative_to(directory))] = [
                path.stat().st_size,
                hashlib.sha256(path.read_bytes()).hexdigest(),
            ]
    return result


def _cache_directory(directory):
    missing = []
    parent = directory
    while not parent.exists():
        missing.append(parent)
        parent = parent.parent
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    control = current_control()
    if control:
        for path in reversed(missing):

            def remove_empty(path=path):
                # Another successful operation may now own data in a shared parent.
                if path.is_dir() and not any(path.iterdir()):
                    path.rmdir()

            control.on_rollback(remove_empty)


def _write_entry(path, entry):
    """Roll back only this replacement; preserve pre-run bytes and concurrent writers."""
    previous = path.read_bytes() if path.exists() else None
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(entry, handle)
        identity = temporary.stat()
        control = current_control()
        if control:

            def restore():
                if not path.exists():
                    return
                stat = path.stat()
                if (stat.st_dev, stat.st_ino) != (identity.st_dev, identity.st_ino):
                    return
                if previous is None:
                    path.unlink()
                else:
                    backup = None
                    try:
                        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
                            backup = Path(handle.name)
                            handle.write(previous)
                        backup.replace(path)
                    finally:
                        if backup:
                            backup.unlink(missing_ok=True)

            control.on_rollback(restore)
        check_cancelled()
        temporary.replace(path)
    finally:
        if temporary:
            temporary.unlink(missing_ok=True)


class CachedAdapter:
    """Batch only missing unique sequences; replay hits under current private IDs.

    Cache failures never turn a successful search into a scientific failure. Failed
    engine calls are never cached. Raw tables/logs stay local beside the cache.
    """

    def __init__(self, adapter, directory, identity):
        self.adapter = adapter
        self.identity = identity
        self.directory = Path(directory) / fingerprint(identity)
        self.provenance = {}
        self.warnings = []

    def begin_run(self):
        self.provenance = {}
        self.warnings = []
        if hasattr(self.adapter, "begin_run"):
            self.adapter.begin_run()

    def close(self):
        if hasattr(self.adapter, "close"):
            self.adapter.close()

    @classmethod
    def from_config(cls, adapter, source, config, directory):
        identity = {
            "schema": CACHE_SCHEMA,
            "source": source,
            **reference_identity(config["databases"][source]),
            "executable": file_identity(adapter.executable),
            "executable_version": tool_version(tuple(file_identity(adapter.executable))),
            "backend": getattr(adapter, "backend_identity", {"backend": "canonical-cli"}),
            "runtime_options": getattr(adapter, "runtime_options", {}),
            "parameters": list(HMM_ARGUMENTS if source == "pfam" else RPS_ARGUMENTS),
            "output_format": getattr(adapter, "output_format", "hmmer-domtbl"),
        }
        return cls(adapter, directory, identity)

    def search(self, queries, workdir, threads):
        # write_queries emits private, nonempty headers; this is the existing parser.
        proteins = [replace(p, id=p.header) for p in read_fasta(queries)]
        groups = {}
        for protein in proteins:
            groups.setdefault(sequence_hash(protein.sequence), []).append(protein)
        entries, missing = {}, []
        checked_raw = {}
        for digest, group in groups.items():
            check_cancelled()
            try:
                entry = json.loads((self.directory / f"{digest}.json").read_text())
                if entry["query_sha256"] != digest or entry["identity"] != self.identity:
                    raise ValueError("Cache identity mismatch")
                if entry["hits_sha256"] != fingerprint(entry["hits"]):
                    raise ValueError("Cache evidence checksum mismatch")
                raw_path = Path(entry["raw_output"])
                if not raw_path.is_dir() or not raw_path.resolve().is_relative_to(
                    self.directory.resolve()
                ):
                    raise ValueError("Cached raw evidence is missing")
                if str(raw_path) not in checked_raw:
                    checked_raw[str(raw_path)] = _raw_manifest(raw_path)
                if not entry["raw_files"] or checked_raw[str(raw_path)] != entry["raw_files"]:
                    raise ValueError("Cached raw evidence is incomplete or changed")
                hits = [
                    Hit(
                        **{
                            **h,
                            "spans": tuple(map(tuple, h["spans"])),
                            "alignments": tuple(Alignment(**a) for a in h["alignments"]),
                        }
                    )
                    for h in entry["hits"]
                ]
                entries[digest] = (hits, entry["raw_output"], True, entry["execution_backend"])
            except (OSError, ValueError, KeyError, TypeError):
                missing.append(group[0])
        failed = set()
        if missing:
            batch = Path(workdir) / "cache-misses.faa"
            write_queries(missing, batch)
            hits = self.adapter.search(batch, workdir, threads)
            self.warnings.extend(getattr(self.adapter, "warnings", []))
            failed = set(getattr(hits, "failed_queries", ()))
            by_query = {}
            for hit in hits:
                by_query.setdefault(hit.query_id, []).append(hit)
            raw = None
            try:
                check_cancelled()
                _cache_directory(self.directory)
                raw = Path(tempfile.mkdtemp(prefix="raw-", dir=self.directory))
                control = current_control()
                if control:
                    control.on_rollback(
                        lambda raw=raw: shutil.rmtree(raw) if raw.exists() else None
                    )
                for asset in sorted(Path(workdir).rglob("*")):
                    check_cancelled()
                    if asset.is_symlink():
                        raise OSError("Cannot retain symbolic raw search evidence")
                    if asset.is_file():
                        target = raw / asset.relative_to(workdir)
                        target.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(asset, target)
                (raw / "provenance.json").write_text(
                    json.dumps(
                        {
                            "identity": self.identity,
                            "threads": threads,
                            "execution_backend": getattr(
                                self.adapter, "last_backend", "canonical-cli"
                            ),
                            "queries": {p.id: sequence_hash(p.sequence) for p in missing},
                        },
                        indent=2,
                    )
                )
                raw_files = _raw_manifest(raw)
            except OSError:
                raw = None
                self.warnings.append("Search succeeded; local cache could not retain raw output.")
            for protein in missing:
                check_cancelled()
                digest = sequence_hash(protein.sequence)
                rows = by_query.get(protein.id, [])
                backend = getattr(self.adapter, "last_backend", "canonical-cli")
                entries[digest] = (rows, str(raw) if raw else None, False, backend)
                if raw and protein.id not in failed:
                    entry = {
                        "identity": self.identity,
                        "query_sha256": digest,
                        "hits": [asdict(h) for h in rows],
                        "raw_output": str(raw),
                        "raw_files": raw_files,
                        "execution_backend": backend,
                    }
                    entry["hits_sha256"] = fingerprint(entry["hits"])
                    try:
                        _write_entry(self.directory / f"{digest}.json", entry)
                    except OSError:
                        self.warnings.append(
                            "Search succeeded; local cache could not save results."
                        )
        result = SearchHits()
        failed_ids = []
        for digest, group in groups.items():
            check_cancelled()
            hits, raw, cached, backend = entries[digest]
            self.provenance[digest] = {
                "query_sha256": digest,
                "cached": cached,
                "raw_output": raw,
                "identity": self.identity,
                "execution_backend": backend,
            }
            for protein in group:
                result.extend(replace(h, query_id=protein.id) for h in hits)
                if group[0].id in failed:
                    failed_ids.append(protein.id)
        result.failed_queries = tuple(failed_ids)
        return result
