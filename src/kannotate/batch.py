import os
import tempfile
from contextlib import ExitStack
from dataclasses import dataclass
from itertools import batched
from pathlib import Path

from .cancellation import RunCancelled, run_scope
from .config import default_threads
from .io.csv import CsvWriter, write_csv
from .io.detection import detect_format
from .io.fasta import iter_fasta, read_fasta
from .io.genbank import iter_genbank, read_genbank, write_genbank, write_genbank_handle
from .io.xlsx import write_xlsx
from .models import RunResult

GENBANK_OUTPUTS = ("genbank", "csv", "both", "xlsx")


@dataclass
class FileResult:
    input: Path
    output: Path | None = None
    error: str | None = None
    run: object = None
    warnings: list | None = None
    outputs: tuple[Path, ...] = ()

    def __post_init__(self):
        if self.output is not None and not self.outputs:
            self.outputs = (self.output,)


def annotate_batch(
    inputs,
    output_dir: Path,
    engine,
    *,
    threads=None,
    min_score=8.0,
    force=False,
    control=None,
    on_result=None,
    genbank_output="genbank",
    table_format="csv",
    retain_annotations=True,
    on_progress=None,
    checkpoint_dir=None,
    resume=False,
):
    if genbank_output not in GENBANK_OUTPUTS:
        raise ValueError("GenBank output must be genbank, csv, both, or xlsx")
    if table_format not in ("csv", "xlsx"):
        raise ValueError("Table format must be csv or xlsx")
    inputs = list(inputs)
    if resume and checkpoint_dir is None:
        raise ValueError("Resume requires a checkpoint directory")
    if checkpoint_dir is not None:
        if not inputs:
            raise ValueError("Choose at least one protein FASTA or GenBank file")
        if retain_annotations or table_format == "xlsx" or genbank_output == "xlsx" or force:
            raise ValueError(
                "Checkpoints require streamed CSV/GenBank outputs without force replacement"
            )
    with run_scope(control) as active, ExitStack() as stack:
        checkpoints = None
        if checkpoint_dir is not None:
            from .checkpoints import CheckpointBatch

            checkpoints = stack.enter_context(
                CheckpointBatch.open(
                    checkpoint_dir,
                    inputs,
                    output_dir,
                    engine,
                    resume=resume,
                    parameters={
                        "threads": default_threads() if threads is None else threads,
                        "min_score": min_score,
                        "genbank_output": genbank_output,
                        "table_format": table_format,
                    },
                )
            )
        return _annotate_batch(
            inputs,
            output_dir,
            engine,
            threads=threads,
            min_score=min_score,
            force=force,
            control=active,
            on_result=on_result,
            genbank_output=genbank_output,
            table_format=table_format,
            retain_annotations=retain_annotations,
            on_progress=on_progress,
            checkpoints=checkpoints,
        )


def _temporary(output_dir, suffix):
    fd, name = tempfile.mkstemp(prefix=".kannotate-", suffix=suffix, dir=output_dir)
    os.close(fd)
    return Path(name)


def _publish(staged, *, force):
    """Publish a file's outputs together under the caller's cancellation lock.

    Each file is atomic. Roll back a partially published pair on a handled error;
    a process/OS crash cannot be made atomic across two separate filesystem names.
    """
    backups = {}
    published = []
    keep_backups = set()
    try:
        if force:
            for _, destination in staged:
                if destination.exists():
                    backup = _temporary(destination.parent, ".backup")
                    backups[destination] = backup
                    backup.unlink()
                    os.link(destination, backup)
        for temporary, destination in staged:
            identity = temporary.stat()
            if force:
                os.replace(temporary, destination)
            else:
                # A file created during searching must never be silently replaced.
                os.link(temporary, destination)
            published.append((destination, (identity.st_dev, identity.st_ino)))
    except BaseException as error:
        rollback_errors = []
        for destination, identity in reversed(published):
            try:
                current = destination.lstat()
                if (current.st_dev, current.st_ino) != identity:
                    raise OSError("Output changed during publication; left untouched")
                if destination in backups:
                    os.replace(backups[destination], destination)
                else:
                    destination.unlink()
            except OSError as rollback_error:
                if destination in backups:
                    keep_backups.add(backups[destination])
                rollback_errors.append(f"{destination}: {rollback_error}")
        if rollback_errors:
            raise OSError(
                f"{error}; output rollback needs attention: {'; '.join(rollback_errors)}; "
                f"recovery copies retained: {', '.join(map(str, keep_backups))}"
            ) from error
        raise
    finally:
        for backup in backups.values():
            if backup not in keep_backups:
                backup.unlink(missing_ok=True)


def _stream_file(
    path, fmt, staged, engine, *, threads, min_score, control, on_progress, checkpoint=None
):
    """Keep at most one record group and its selected evidence alive between writes."""
    summary = RunResult(
        [],
        {"hmmscan": 0.0, "rpsblast": 0.0, "total": 0.0},
        [],
        threads or default_threads(),
        getattr(engine, "database_versions", {}).copy(),
    )
    warning_counts = {}
    first = True
    with ExitStack() as stack:
        writers = []
        for temporary, destination in staged:
            csv = destination.suffix == ".csv"
            handle = stack.enter_context(
                temporary.open(
                    "w",
                    encoding="utf-8-sig" if csv else "utf-8",
                    newline="" if csv else None,
                    buffering=1024 * 1024,
                )
            )
            writers.append(
                (
                    CsvWriter(
                        handle,
                        genbank=fmt == "genbank",
                        database_versions=summary.database_versions,
                    )
                    if csv
                    else handle,
                    csv,
                )
            )
        groups = (
            iter_genbank(path)
            if fmt == "genbank"
            else batched(iter_fasta(path), 512, strict=False)
        )
        for index, group in enumerate(groups):
            control.check()
            data = group if fmt == "genbank" else None
            proteins = data.proteins if data is not None else list(group)
            if on_progress:
                on_progress(path, summary.protein_count)
            run = (
                checkpoint.group(
                    index,
                    proteins,
                    data.warnings if data else [],
                    threads=threads,
                    min_score=min_score,
                )
                if checkpoint is not None
                else engine.annotate(proteins, threads=threads, min_score=min_score)
            )
            for writer, csv in writers:
                control.check()
                if csv:
                    writer.append(run.annotations, genbank=data, input_name=path.name)
                else:
                    write_genbank_handle(data, run.annotations, writer, provenance=first)
            first = False
            summary.protein_count += run.protein_count
            summary.database_versions = run.database_versions
            summary.threads = run.threads
            for key, elapsed in run.timings.items():
                summary.timings[key] = summary.timings.get(key, 0.0) + elapsed
            for warning in run.warnings + (data.warnings if data else []):
                warning_counts[warning] = warning_counts.get(warning, 0) + 1
            # No completed-window RunResult or scientific evidence is retained by the UI.
            del run
        if on_progress:
            on_progress(path, summary.protein_count)
        control.check()
    summary.warnings = [
        message if count == 1 else f"{message} ({count} batches)"
        for message, count in warning_counts.items()
    ]
    return summary


def _annotate_batch(
    inputs,
    output_dir,
    engine,
    *,
    threads,
    min_score,
    force,
    control,
    on_result,
    genbank_output,
    table_format,
    retain_annotations,
    on_progress,
    checkpoints=None,
):
    control.check()
    inputs = [Path(p).expanduser() for p in inputs]
    if not inputs:
        raise ValueError("Choose at least one protein FASTA or GenBank file")
    input_paths = {p.resolve() for p in inputs}
    output_dir = (
        Path(output_dir).expanduser() if output_dir is not None else inputs[0].resolve().parent
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    used = set()
    for input_path in inputs:
        control.check()
        path = Path(input_path)
        staged = []
        try:
            checkpoint = checkpoints.file(path) if checkpoints is not None else None
            if checkpoint is not None:
                checkpoint.discard_abandoned_staging()
            fmt = detect_format(path)
            extensions = (
                [f".{table_format}"]
                if fmt == "fasta"
                else [f".{genbank_output}"]
                if genbank_output in ("csv", "xlsx")
                else [".gbk", f".{table_format}"]
                if genbank_output == "both"
                else [".gbk"]
            )
            destinations = [output_dir / f"{path.stem}.kannotate{ext}" for ext in extensions]
            # Two same-stem inputs never silently overwrite each other, even with --force.
            collisions = [p for p in destinations if p.resolve() in used]
            used.update(p.resolve() for p in destinations)
            if collisions:
                raise ValueError(f"Output name collision: {collisions[0].name}")
            for destination in destinations:
                if destination.resolve() in input_paths:
                    raise ValueError("Output cannot replace an input")
                if destination.is_symlink():
                    raise ValueError(f"Output is a symbolic link: {destination.name}")
                if (
                    destination.exists()
                    and not force
                    and not (checkpoint is not None and checkpoint.allows_existing(destination))
                ):
                    raise FileExistsError(
                        f"Output exists: {destination.name}; use --force to replace it"
                    )
            for destination in destinations:
                temporary = _temporary(output_dir, destination.suffix)
                staged.append((temporary, destination))
            data = None
            if not retain_annotations and ".xlsx" not in extensions:
                run = _stream_file(
                    path,
                    fmt,
                    staged,
                    engine,
                    threads=threads,
                    min_score=min_score,
                    control=control,
                    on_progress=on_progress,
                    checkpoint=checkpoint,
                )
            else:
                data = read_genbank(path) if fmt == "genbank" else None
                control.check()
                proteins = data.proteins if data else read_fasta(path)
                run = engine.annotate(proteins, threads=threads, min_score=min_score)
                for temporary, destination in staged:
                    control.check()
                    if destination.suffix == ".gbk":
                        write_genbank(data, run.annotations, temporary)
                    else:
                        writer = write_csv if destination.suffix == ".csv" else write_xlsx
                        writer(
                            run.annotations,
                            temporary,
                            database_versions=run.database_versions,
                            genbank=data,
                            input_name=path.name,
                        )
            with control.publication():
                if checkpoint is None:
                    _publish(staged, force=force)
                else:
                    checkpoint.publish(staged, run, publish=_publish)
                results.append(
                    FileResult(
                        path,
                        destinations[0],
                        run=run,
                        warnings=run.warnings + (data.warnings if data else []),
                        outputs=tuple(destinations),
                    )
                )
                if on_result:
                    on_result(results[-1])
        except RunCancelled:
            raise
        except Exception as exc:  # noqa: BLE001 -- a broken input must not abort the batch
            control.check()
            results.append(FileResult(path, error=str(exc)))
            if on_result:
                on_result(results[-1])
        finally:
            for temporary, _ in staged:
                temporary.unlink(missing_ok=True)
    return results
