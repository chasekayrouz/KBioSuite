"""Exclusive local checkpoints at the existing canonical file-group boundaries."""

import fcntl
import importlib
import json
import marshal
import os
import platform
import shutil
import tempfile
from contextlib import contextmanager
from dataclasses import asdict
from importlib.resources import files
from pathlib import Path

import Bio

from . import __version__
from .cancellation import check_cancelled, current_control
from .models import Alignment, Annotation, Hit, Protein, RunResult
from .search.identity import digest, file_identity, fingerprint, reference_identity, tool_version

SCHEMA = 1
BOUNDARIES = "iter-genbank-512-cds-128-records/iter-fasta-512/engine-512-v1"


def implementation_identity():
    modules = (
        "engine",
        "batch",
        "checkpoints",
        "normalization",
        "hits",
        "scoring",
        "selection",
        "architectures",
        "descriptions",
        "products",
        "models",
        "io.fasta",
        "io.genbank",
        "io.csv",
        "io.tabular",
        "search.hmmer",
        "search.rpsblast",
        "search.scheduling",
        "search.resident",
        "search.identity",
        "search.common",
        "cancellation",
    )
    code = {}
    for name in modules:
        module = importlib.import_module(f"kannotate.{name}")
        # Works for source and PyInstaller's frozen modules. Different Python builds,
        # source locations or code bytes conservatively require a fresh checkpoint.
        value = module.__loader__.get_code(module.__name__)
        code[name] = fingerprint(marshal.dumps(value).hex())
    resources = files("kannotate.resources")
    return {
        "version": __version__,
        "python": platform.python_version(),
        "biopython": Bio.__version__,
        "platform": platform.platform(),
        "code": code,
        "rules": {
            name: fingerprint(resources.joinpath(name).read_bytes().hex())
            for name in ("scoring.json", "search.json", "architectures.json", "products.json")
        },
    }


def engine_identity(engine):
    """Hash operative references/tools once per invocation, including every resume."""
    if not hasattr(engine, "config"):
        raise ValueError("Checkpoints require an engine with a bound runtime identity")
    result = {"implementation": implementation_identity(), "searches": {}}
    for source, adapter in (("pfam", engine.hmm_adapter), ("cdd", engine.rps_adapter)):
        if hasattr(adapter, "adapter"):
            raise ValueError("Batch checkpoints require an uncached engine")
        reference = reference_identity(engine.config["databases"][source])
        if getattr(adapter, "reference", reference) != reference:
            raise ValueError("Resident reference changed; construct a fresh engine")
        description = getattr(adapter, "description_identity", None)
        if description is not None and file_identity(description[0]) != description:
            raise ValueError("CDD descriptions changed; construct a fresh engine")
        result["searches"][source] = {
            **reference,
            "contents": {item[0]: digest(item[0]) for item in reference["assets"]},
            "executable": file_identity(adapter.executable),
            "executable_sha256": digest(adapter.executable),
            "executable_version": tool_version(tuple(file_identity(adapter.executable))),
            "backend": getattr(adapter, "backend_identity", {"backend": "canonical-cli"}),
            "runtime_options": getattr(adapter, "runtime_options", {}),
            "output_format": getattr(adapter, "output_format", "hmmer-domtbl"),
        }
    return result


def _sync_directory(path):
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _atomic_json(path, value):
    descriptor, temporary = tempfile.mkstemp(prefix=".json-", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w") as handle:
            json.dump(value, handle, sort_keys=True, allow_nan=False)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        _sync_directory(path.parent)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _manifest(directory):
    result = {}
    for path in sorted(directory.rglob("*")):
        check_cancelled()
        if path.is_symlink():
            raise ValueError("Checkpoint evidence contains a symbolic link")
        if path.is_file() and path.name != "receipt.json":
            result[str(path.relative_to(directory))] = [path.stat().st_size, digest(path)]
    return result


def _restore_run(value):
    def hit(row):
        return Hit(
            **{
                **row,
                "spans": tuple(map(tuple, row["spans"])),
                "alignments": tuple(Alignment(**a) for a in row["alignments"]),
            }
        )

    return RunResult(
        **{
            **value,
            "annotations": [
                Annotation(
                    **{
                        **row,
                        "protein": Protein(**row["protein"]),
                        "hmm_hits": tuple(hit(h) for h in row["hmm_hits"]),
                        "rps_hits": tuple(hit(h) for h in row["rps_hits"]),
                    }
                )
                for row in value["annotations"]
            ],
        }
    )


class CheckpointBatch:
    def __init__(self, directory, plan, engine):
        self.directory, self.plan, self.engine = directory, plan, engine
        self.plan_sha256 = fingerprint(plan)

    @classmethod
    @contextmanager
    def open(cls, directory, inputs, output_dir, engine, *, resume, parameters):
        directory = Path(directory).expanduser().resolve()
        inputs = [Path(p).expanduser().absolute() for p in inputs]
        output_dir = Path(output_dir or inputs[0].parent).expanduser().resolve()
        if any(p == directory or p.is_relative_to(directory) for p in inputs):
            raise ValueError("Checkpoint directory cannot contain an input")
        if output_dir == directory or output_dir.is_relative_to(directory):
            raise ValueError("Outputs must be outside the checkpoint directory")
        # Lock before any costly reference verification. flock is released by the OS
        # after a crash, so no unsafe stale-claim recovery or PID guessing is needed.
        if resume:
            if not any((directory / name).is_file() for name in ("plan.json", "initializing.json")):
                raise ValueError("No valid checkpoint plan to resume; choose a new directory")
        else:
            directory.mkdir(parents=True, exist_ok=False, mode=0o700)
        with (directory / ".lock").open("a+") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise ValueError("Checkpoint batch is already in use") from exc
            plan_path = directory / "plan.json"
            initialization = {
                "schema": SCHEMA,
                "sources": list(map(str, inputs)),
                "output_dir": str(output_dir),
                "parameters": parameters,
            }
            if not plan_path.is_file():
                marker = directory / "initializing.json"
                if resume:
                    if json.loads(marker.read_text()) != initialization or any(
                        p.name not in {".lock", "initializing.json"}
                        and not p.name.startswith(".json-")
                        for p in directory.iterdir()
                    ):
                        raise ValueError("Incomplete checkpoint initialization identity changed")
                else:
                    _atomic_json(marker, initialization)
            plan = {
                "schema": SCHEMA,
                "boundaries": BOUNDARIES,
                "parameters": parameters,
                "output_dir": str(output_dir),
                "runtime": engine_identity(engine),
                "sources": [
                    {
                        "path": str(p),
                        "filename": p.name,
                        "identity": file_identity(p),
                        "sha256": digest(p),
                    }
                    for p in inputs
                ],
            }
            if plan_path.is_file():
                stored = json.loads(plan_path.read_text())
                if stored != {"plan": plan, "sha256": fingerprint(plan)}:
                    raise ValueError(
                        "Checkpoint identity changed (input, settings, code, tool or reference); "
                        "use the original runtime/settings or a new checkpoint directory"
                    )
            else:
                _atomic_json(plan_path, {"plan": plan, "sha256": fingerprint(plan)})
            batch = cls(directory, plan, engine)
            yield batch

    def file(self, path):
        index = [s["path"] for s in self.plan["sources"]].index(str(path.expanduser().absolute()))
        directory = self.directory / f"file-{index:05d}"
        directory.mkdir(exist_ok=True, mode=0o700)
        return CheckpointFile(self, directory, self.plan["sources"][index])

    def check_runtime(self):
        for search in self.plan["runtime"]["searches"].values():
            for identity in search["assets"] + [
                search["executable"],
                search["installation_identity"],
            ]:
                if file_identity(identity[0]) != identity:
                    raise ValueError("Checkpoint tool/reference changed during the batch")


class CheckpointFile:
    def __init__(self, batch, directory, source):
        self.batch, self.directory, self.source = batch, directory, source

    def check_source(self, *, full=False):
        path = Path(self.source["path"])
        if file_identity(path) != self.source["identity"] or (
            full and digest(path) != self.source["sha256"]
        ):
            raise ValueError("Input changed during checkpointed annotation")
        self.batch.check_runtime()

    def group(self, index, proteins, warnings, *, threads, min_score):
        self.check_source()
        identity = {
            "plan": self.batch.plan_sha256,
            "index": index,
            "proteins_sha256": fingerprint([asdict(p) for p in proteins]),
            "parser_warnings": warnings,
        }
        destination = self.directory / f"group-{index:08d}"
        if destination.exists():
            receipt = json.loads((destination / "receipt.json").read_text())
            if receipt != {"identity": identity, "files": _manifest(destination)}:
                raise ValueError("Checkpoint evidence is incomplete, changed or out of order")
            run = _restore_run(json.loads((destination / "run.json").read_text()))
        else:
            with tempfile.TemporaryDirectory(prefix=".stage-", dir=self.directory) as temporary:
                stage = Path(temporary)
                run = self.batch.engine.annotate(
                    proteins, threads=threads, min_score=min_score, evidence_dir=stage / "evidence"
                )
                if run.warnings:
                    raise ValueError(
                        "Search warnings prevent a successful checkpoint: "
                        + "; ".join(run.warnings)
                    )
                self.check_source()
                _atomic_json(stage / "run.json", asdict(run))
                manifest = _manifest(stage)
                # fsync raw outputs as well as the commit receipt before directory promotion.
                for name in manifest:
                    with (stage / name).open("rb") as handle:
                        os.fsync(handle.fileno())
                _atomic_json(stage / "receipt.json", {"identity": identity, "files": manifest})
                with current_control().publication():
                    os.rename(stage, destination)
                    _sync_directory(self.directory)
        if [a.protein for a in run.annotations] != proteins or run.protein_count != len(proteins):
            raise ValueError("Checkpoint annotations differ from the canonical source group")
        return run

    def _publication(self):
        path = self.directory / "publication.json"
        if not path.exists():
            return None
        receipt = json.loads(path.read_text())
        if receipt["plan"] != self.batch.plan_sha256:
            raise ValueError("Checkpoint output receipt identity changed")
        return receipt

    def allows_existing(self, destination):
        receipt = self._publication()
        return bool(
            receipt
            and str(destination.resolve()) in receipt["outputs"]
            and digest(destination) == receipt["outputs"][str(destination.resolve())]
        )

    def publish(self, staged, summary, *, publish):
        self.check_source(full=True)
        outputs = {str(d.resolve()): digest(t) for t, d in staged}
        previous = self._publication()
        if previous and previous["outputs"] != outputs:
            raise ValueError("Checkpoint assembly differs from its prepared output receipt")
        # Write the intended bytes BEFORE publication; resume can recognize a file
        # linked just before a crash, including the first member of a Both pair.
        _atomic_json(
            self.directory / "publication.json",
            {
                "plan": self.batch.plan_sha256,
                "outputs": outputs,
                "summary": asdict(summary),
            },
        )
        missing = []
        for temporary, destination in staged:
            if destination.exists():
                if (
                    destination.is_symlink()
                    or digest(destination) != outputs[str(destination.resolve())]
                ):
                    raise FileExistsError(f"Output exists with different bytes: {destination.name}")
            else:
                missing.append((temporary, destination))
        publish(missing, force=False)

    def discard_abandoned_staging(self):
        # Exclusive batch lock proves no worker still owns these unpublished directories.
        for path in self.directory.glob(".stage-*"):
            if path.is_dir() and not path.is_symlink():
                shutil.rmtree(path)
