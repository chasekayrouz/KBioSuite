"""Normalized contracts. Engine-specific outputs remain available in provenance."""

import json
from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Protocol


@dataclass(frozen=True)
class Progress:
    phase: str
    message: str
    elapsed_seconds: float = 0


ProgressCallback = Callable[[Progress], None]


@dataclass(frozen=True)
class PredictionResult:
    structure_path: Path
    predictor: str
    predictor_version: str
    sequence_length: int
    confidence: dict[str, Any]
    per_residue_confidence: list[float]
    runtime_seconds: float
    provenance: dict[str, Any] = field(default_factory=dict)
    confidence_scale: str = "pLDDT 0–100"
    warnings: list[dict[str, Any]] = field(default_factory=list)


@dataclass(frozen=True)
class SearchResult:
    query: str
    target: str
    evalue: float
    probability: float
    bits: float
    sequence_identity: float
    alignment_length: int
    query_length: int
    target_length: int
    query_coverage: float
    target_coverage: float
    query_tm: float
    target_tm: float
    alignment_tm: float
    lddt: float | None
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    query_alignment: str
    target_alignment: str
    rotation: list[list[float]]
    translation: list[float]
    raw: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class DependencyStatus:
    name: str
    ready: bool
    message: str
    version: str | None = None


class PredictionBackend(Protocol):
    identity: str

    def status(self) -> DependencyStatus: ...
    def predict(
        self,
        sequence,
        output: Path,
        cancel,
        progress: ProgressCallback,
        *,
        msa_path: Path | None = None,
    ) -> PredictionResult: ...


class StructuralSearchBackend(Protocol):
    def status(self) -> DependencyStatus: ...
    def search(
        self, prediction: PredictionResult, output: Path, cancel, progress: ProgressCallback
    ) -> list[SearchResult]: ...


def write_json(path: Path, value):
    """Atomic, local run snapshots. Never serialize environment variables."""
    payload = asdict(value) if hasattr(value, "__dataclass_fields__") else value
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, default=str, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)
