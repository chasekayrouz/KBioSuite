"""File-based API owner; never reconstructs environments during ordinary launch."""

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

from . import __version__
from .config import Config
from .database import database_status
from .pipeline import Pipeline
from .predictors import BACKENDS, AF3Backend, BoltzBackend
from .process import Cancellation
from .search import FoldseekBackend
from .setup import (
    check_command_line_tools,
    install_af3,
    install_boltz,
    install_foldseek,
    install_pdb,
    install_uv,
)


def status(config):
    return [asdict(item.status()) for item in (BoltzBackend(config), AF3Backend(config))] + [
        asdict(FoldseekBackend(config).executable_status()),
        asdict(database_status(config.pdb_prefix)),
    ]


def setup_component(component, config, cancel, progress):
    if component in ("boltz", "af3-code"):
        # Before anything is downloaded: the engine builds need Apple's developer tools.
        check_command_line_tools("Boltz" if component == "boltz" else "AlphaFold 3")
        config.home.mkdir(parents=True, exist_ok=True)
        # Kmatch's pinned uv (<home>/tools/uv), downloaded here on first engine setup only.
        uv = str(install_uv(config, cancel, progress))
        if component == "boltz":
            return install_boltz(config, cancel, progress, uv)
        return install_af3(config, cancel, progress, uv)
    config.home.mkdir(parents=True, exist_ok=True)
    if component == "foldseek":
        return install_foldseek(config, cancel, progress)
    return install_pdb(config, cancel, progress)


def build_parser():
    parser = argparse.ArgumentParser(
        description=f"Kmatch {__version__} · local single-protein folding and structural search"
    )
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("status")
    setup = sub.add_parser("setup")
    setup.add_argument("component", choices=["boltz", "foldseek", "pdb100", "af3-code"])
    run = sub.add_parser("run")
    run.add_argument("sequence_file", type=Path)
    run.add_argument(
        "--backend", choices=[key for key, definition in BACKENDS.items() if definition.production]
    )
    run.add_argument("--msa", type=Path)
    run.add_argument(
        "--network-msa",
        action="store_true",
        help=(
            "Explicitly send the sequence to api.colabfold.com for an MSA; inference remains local."
        ),
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        config = Config.load()
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"Kmatch could not load its configuration: {exc}\n")
    if args.command == "status":
        print(json.dumps(status(config), indent=2))
        return
    cancel = Cancellation()

    def progress(event):
        print(event.message, flush=True)

    try:
        if args.command == "setup":
            setup_component(args.component, config, cancel, progress)
        else:
            if args.backend:
                from dataclasses import replace

                config = replace(config, predictor=args.backend)
            result = Pipeline(config).run(
                args.sequence_file.read_text(),
                cancel=cancel,
                progress=progress,
                network_msa=args.network_msa,
                msa_path=args.msa,
            )
            print(result.directory)
    except KeyboardInterrupt:
        cancel.cancel()
        print("Cancelled.", file=sys.stderr)
        raise SystemExit(130) from None
    except (OSError, ValueError, RuntimeError) as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
