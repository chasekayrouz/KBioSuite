"""Resolve existing installations; never download or vendor executables."""

import os
import shutil
from functools import lru_cache
from pathlib import Path

from .cancellation import check_cancelled
from .search.common import SearchError, run_search

REQUIRED = ("hmmscan", "hmmpress", "rpsblast", "blastdbcmd")
OPTIONAL = ("makeprofiledb",)


def _executable(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)


def _usable(path: Path) -> bool:
    if not _executable(path):
        return False
    try:
        _probe(_identity(path))
        return True
    except (OSError, SearchError):
        return False


def _identity(path):
    check_cancelled()
    path = Path(path).resolve()
    stat = path.stat()
    return str(path), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


@lru_cache(maxsize=32)
def _probe(identity):
    # Cache only successful probes; every call still verifies current file identity.
    path = identity[0]
    option = "-h" if Path(path).name.startswith("hmm") else "-version"
    return run_search([path, option], timeout=3)


def resolve_executable(name: str, configured: str | None = None) -> str | None:
    if configured:
        path = Path(configured).expanduser()
        return str(path.resolve()) if _usable(path) else None
    found = shutil.which(name)
    if found and _usable(Path(found)):
        return str(Path(found).resolve())
    for directory in ("/opt/homebrew/bin", "/usr/local/bin", "/usr/bin"):
        path = Path(directory) / name
        if _usable(path):
            return str(path.resolve())
    return None


def resolve_dependencies(configured=None):
    configured = configured or {}
    return {name: resolve_executable(name, configured.get(name)) for name in REQUIRED + OPTIONAL}


def executable_version(path: str) -> str:
    output = _probe(_identity(path))
    for line in output.splitlines():
        if (
            "HMMER " in line
            or "blast:" in line
            or "blastdbcmd:" in line
            or "makeprofiledb:" in line
        ):
            return line.lstrip("# ").strip()
    return output.splitlines()[0] if output else "unknown"
