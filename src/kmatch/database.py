"""Database readiness never triggers a download."""

import fcntl
import json
import struct
import time
from contextlib import contextmanager
from pathlib import Path

from .models import DependencyStatus

REQUIRED = (
    "",
    ".dbtype",
    ".index",
    ".lookup",
    "_ss",
    "_ss.dbtype",
    "_ss.index",
    "_ca",
    "_ca.dbtype",
    "_ca.index",
    "_h",
    "_h.dbtype",
    "_h.index",
)


@contextmanager
def database_lease(config, cancel, *, exclusive=False):
    """Coordinate managed searches/legacy exports with activation and retirement."""
    root = config.home / "databases/pdb100"
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".usage.lock").open("a") as lock:
        mode = fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH
        while True:
            cancel.check()
            try:
                fcntl.flock(lock, mode | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                time.sleep(0.1)
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def database_status(prefix: Path):
    missing = [
        suffix or "(amino-acid data)"
        for suffix in REQUIRED
        if not Path(str(prefix) + suffix).is_file()
        or Path(str(prefix) + suffix).stat().st_size == 0
    ]
    if missing:
        return DependencyStatus(
            "PDB100",
            False,
            "PDB100 is missing or incomplete. Install/update it with Setup. "
            "Missing: " + ", ".join(missing),
        )
    for suffix, expected in (
        (".dbtype", 0),
        ("_ss.dbtype", 0),
        ("_ca.dbtype", 101),
        ("_h.dbtype", 12),
    ):
        data = Path(str(prefix) + suffix).read_bytes()
        if len(data) != 4 or (struct.unpack("<I", data)[0] & 0xFFFF) != expected:
            return DependencyStatus(
                "PDB100", False, f"Invalid database type file: {suffix}. Reinstall PDB100."
            )
    manifest = prefix.parent / "installation.json"
    version = None
    if manifest.exists():
        try:
            data = json.loads(manifest.read_text())
            version = data.get("pdb_date")
            if data.get("status") != "complete":
                return DependencyStatus(
                    "PDB100", False, "Database installation was interrupted; resume Setup."
                )
        except (ValueError, OSError):
            return DependencyStatus(
                "PDB100", False, "Database installation manifest is unreadable; resume Setup."
            )
    return DependencyStatus(
        "PDB100", True, "Local PDB-derived database ready with Cα coordinates.", version
    )


def inspect_database(prefix, cancel):
    """Full index bounds/key agreement check, used on explicit setup only."""
    state = database_status(prefix)
    if not state.ready:
        raise RuntimeError(state.message)
    reference = None
    for suffix in ("", "_ss", "_ca", "_h"):
        size = Path(str(prefix) + suffix).stat().st_size
        keys = set()
        with Path(str(prefix) + suffix + ".index").open() as stream:
            for i, line in enumerate(stream):
                if i % 10000 == 0:
                    cancel.check()
                key, offset, length = map(int, line.split())
                if key in keys or offset < 0 or length <= 0 or offset + length > size:
                    raise RuntimeError(f"PDB100 {suffix or 'sequence'} index is inconsistent.")
                keys.add(key)
        if not keys or (reference is not None and keys != reference):
            raise RuntimeError(
                "PDB100 coordinate/sequence/header indices do not cover the same hits."
            )
        reference = keys
    return {"entries": len(reference), "sample_key": min(reference)}
