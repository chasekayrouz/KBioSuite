"""Atomic checkpoints and a deliberately run-local response cache."""

import json
import os
import re
from pathlib import Path

from .models import Run, utc_now


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", value)


def atomic_write(path: Path, data: str | bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    mode = "wb" if isinstance(data, bytes) else "w"
    with temporary.open(mode) as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


class RunStore:
    def __init__(self, root: Path):
        self.root = root
        root.mkdir(parents=True, exist_ok=True, mode=0o700)

    def directory(self, run_id: str) -> Path:
        if not re.fullmatch(r"[0-9a-f]{32}", run_id):
            raise ValueError("Invalid run identifier.")
        return self.root / run_id

    def save(self, run: Run) -> None:
        run.updated_at = utc_now()
        atomic_write(self.directory(run.id) / "run.json", json.dumps(run.to_dict(), indent=2))

    def load(self, run_id: str) -> Run:
        return Run.from_dict(json.loads((self.directory(run_id) / "run.json").read_text()))

    def list_runs(self) -> list[dict]:
        results = []
        for path in self.root.glob("*/run.json"):
            try:
                run = Run.from_dict(json.loads(path.read_text()))
                results.append(
                    {
                        "id": run.id,
                        "created_at": run.created_at,
                        "status": run.status,
                        "query": run.query.description,
                        "counts": run.to_dict()["counts"],
                    }
                )
            except (ValueError, KeyError, TypeError):
                continue
        return sorted(results, key=lambda r: r["created_at"], reverse=True)
