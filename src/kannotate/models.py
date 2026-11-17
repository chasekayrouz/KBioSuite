from dataclasses import dataclass, field
from typing import Literal


@dataclass(frozen=True)
class Protein:
    id: str
    header: str
    sequence: str


@dataclass(frozen=True)
class Alignment:
    """One engine-reported query/model HSP; coordinates are 1-based, inclusive."""

    query_start: int
    query_end: int
    model_start: int
    model_end: int
    query: str
    model: str
    evalue: float
    bitscore: float


@dataclass(frozen=True)
class Hit:
    """A profile-level candidate. All coordinates are 1-based, inclusive.

    spans retain the actual aligned intervals; start/end show their outer envelope.
    Raw bit scores are diagnostic only and never compared between search engines.
    """

    engine: Literal["hmm", "rps"]
    query_id: str
    profile_id: str
    name: str
    description: str
    evalue: float
    start: int
    end: int
    query_length: int
    profile_length: int = 0
    bitscore: float = 0.0
    spans: tuple[tuple[int, int], ...] = ()
    score: float = 0.0
    alignments: tuple[Alignment, ...] = ()


@dataclass(frozen=True)
class Annotation:
    protein: Protein
    product: str
    score: float
    hmm_hits: tuple[Hit, ...] = ()
    rps_hits: tuple[Hit, ...] = ()
    selection: str = "single"


@dataclass
class RunResult:
    annotations: list[Annotation]
    timings: dict[str, float] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    threads: int = 1
    database_versions: dict[str, str] = field(default_factory=dict)
    search_provenance: dict[str, list[dict]] = field(default_factory=dict)
    protein_count: int | None = None

    def __post_init__(self):
        if self.protein_count is None:
            self.protein_count = len(self.annotations)
