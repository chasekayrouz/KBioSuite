import argparse
import math
import sys
from contextlib import closing
from pathlib import Path

from . import __version__
from .batch import GENBANK_OUTPUTS, annotate_batch
from .config import default_threads
from .databases import setup, status
from .engine import AnnotationEngine
from .scoring import DEFAULT_MIN_SCORE


def positive_int(value):
    value = int(value)
    if value < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return value


def score_value(value):
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise argparse.ArgumentTypeError("must be finite and nonnegative")
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="kannotate", description="Local bacterial protein reannotation"
    )
    parser.add_argument("--version", action="version", version=f"Kannotate {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("setup", help="install and validate the local Pfam and CDD databases")
    commands.add_parser("status", help="show dependency and database readiness")
    annotate = commands.add_parser("annotate", help="annotate protein FASTA and GenBank files")
    annotate.add_argument("inputs", nargs="+", type=Path)
    annotate.add_argument("-o", "--output-dir", type=Path, help="default: first input's folder")
    annotate.add_argument(
        "--genbank-output",
        choices=GENBANK_OUTPUTS,
        default="genbank",
        help="GenBank input output format (default: genbank); both uses the selected table format",
    )
    annotate.add_argument(
        "--table-format",
        choices=("csv", "xlsx"),
        default="csv",
        help="FASTA and combined-output table format (default: csv)",
    )
    annotate.add_argument("--threads", type=positive_int, default=default_threads())
    annotate.add_argument("--min-score", type=score_value, default=DEFAULT_MIN_SCORE)
    annotate.add_argument("--force", action="store_true")
    annotate.add_argument(
        "--checkpoint-dir",
        type=Path,
        help="keep canonical groups/evidence for streamed CSV/GenBank batches in a new directory",
    )
    annotate.add_argument("--resume", action="store_true", help="resume the exact checkpoint plan")
    annotate.add_argument(
        "--pfam-backend",
        choices=("cli", "auto"),
        default="cli",
        help="auto uses qualified resident models where available (default: cli)",
    )
    args = parser.parse_args(argv)
    try:
        if args.command == "setup":
            setup()
            return 0
        if args.command == "status":
            current = status()
            for name, path in current["executables"].items():
                print(f"{name}: {path or 'not found'}")
            for name in ("pfam", "cdd"):
                db = current["databases"].get(name, {})
                print(
                    f"{name}: {db.get('version', 'not installed')} "
                    f"({'ready' if db.get('ready') else 'not ready'})"
                )
            print(f"Database storage: {current['disk_bytes'] / 1024**3:.2f} GiB")
            return (
                0
                if all(current["databases"].get(k, {}).get("ready") for k in ("pfam", "cdd"))
                else 1
            )
        engine = AnnotationEngine.from_config(pfam_backend=args.pfam_backend)
        with closing(engine):
            results = annotate_batch(
                args.inputs,
                args.output_dir,
                engine,
                threads=args.threads,
                min_score=args.min_score,
                force=args.force,
                genbank_output=args.genbank_output,
                table_format=args.table_format,
                retain_annotations=False,
                checkpoint_dir=args.checkpoint_dir,
                resume=args.resume,
                on_progress=lambda path, count: print(
                    f"{path.name}: {count} proteins completed", file=sys.stderr, flush=True
                ),
            )
        succeeded = 0
        for result in results:
            if result.error:
                print(f"FAILED {result.input.name}: {result.error}", file=sys.stderr)
                continue
            succeeded += 1
            run = result.run
            print(
                f"OK {', '.join(p.name for p in result.outputs)}: {run.protein_count} proteins; "
                f"hmmscan {run.timings['hmmscan']:.2f}s; rpsblast {run.timings['rpsblast']:.2f}s; "
                f"total {run.timings['total']:.2f}s; {run.threads} threads"
            )
            for warning in result.warnings:
                print(f"  Warning: {warning}", file=sys.stderr)
        versions = ", ".join(f"{k} {v}" for k, v in engine.database_versions.items())
        print(f"Completed: {succeeded} succeeded, {len(results) - succeeded} failed; {versions}")
        return 0 if succeeded == len(results) else 1
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"Kannotate: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
