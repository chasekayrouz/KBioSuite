"""Resolve replaceable development or bundled dependencies, never install them."""
from __future__ import annotations

import os
import shutil
from pathlib import Path
from importlib.resources import files


def resource_root() -> Path:
    return Path(os.environ["KSTACK_RESOURCES"]) if os.environ.get("KSTACK_RESOURCES") else Path(str(files("kstack").joinpath("resources")))


def resolve_mmseqs() -> str:
    candidates = [os.environ.get("KSTACK_MMSEQS"),
                  str(resource_root() / "runtime/mmseqs/bin/mmseqs"),
                  shutil.which("mmseqs"),
                  str(Path(__file__).resolve().parents[2] / ".tools/mmseqs/bin/mmseqs")]
    for candidate in candidates:
        if candidate and Path(candidate).is_file() and os.access(candidate, os.X_OK):
            return candidate
    raise RuntimeError("MMseqs2 is unavailable. Set KSTACK_MMSEQS or run scripts/kstack/setup-mmseqs.sh.")
