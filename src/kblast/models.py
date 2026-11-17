"""Small serializable scientific model, independent of the GUI."""

import hashlib
import re
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from typing import Any

from . import __version__


def utc_now() -> str:
    return datetime.now(UTC).isoformat(timespec="microseconds")


def digest(data: str | bytes) -> str:
    return hashlib.sha256(data.encode() if isinstance(data, str) else data).hexdigest()


BLAST_SUMMARY_DEFINITIONS = {
    "query_coverage": "union of inclusive query positions across all HSPs / query length × 100",
    "percent_identity": "maximum HSP identity / alignment length × 100",
    "bit_score": "maximum HSP bit score",
    "evalue": "minimum HSP E-value",
}
NO_MAPPING_MESSAGE = "No nucleotide occurrence returned by IPG or NCBI protein-to-nucleotide links."
TARGET_CDS_ERROR = (
    "No target CDS matches the resolved coordinates and strand in the GenBank context."
)


def usable_nucleotide_location(accession: str, start: int | None, end: int | None) -> bool:
    return bool(
        re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*(?:\.\d+)?", accession)
        and type(start) is int
        and type(end) is int
        and min(start, end) > 0
    )


def summarize_hsps(hsps: list[dict[str, Any]], query_length: int) -> dict[str, float]:
    """Subject summaries only; original HSPs stay separate and unchanged."""
    if not hsps or query_length < 1:
        raise ValueError("BLAST subject needs valid HSPs and a positive query length.")
    intervals = []
    for hsp in hsps:
        low, high = sorted((hsp["query-from"], hsp["query-to"]))
        if not 1 <= low <= high <= query_length or hsp["align-len"] <= 0:
            raise ValueError("BLAST HSP has invalid query coordinates or alignment length.")
        intervals.append((low, high))
    covered, previous_end = 0, 0
    for low, high in sorted(intervals):
        covered += max(0, high - max(low, previous_end + 1) + 1)
        previous_end = max(previous_end, high)
    return {
        "query_coverage": 100 * covered / query_length,
        "percent_identity": max(100 * h["identity"] / h["align-len"] for h in hsps),
        "bit_score": max(h["bit-score"] for h in hsps),
        "evalue": min(h["evalue"] for h in hsps),
    }


@dataclass
class Query:
    identifier: str
    description: str
    sequence: str
    warnings: list[str] = field(default_factory=list)


@dataclass
class Settings:
    top_n: int = 10
    flank: int = 20_000
    database: str = "nr"

    def __post_init__(self) -> None:
        if type(self.top_n) is not int or not 1 <= self.top_n <= 5000:
            raise ValueError("BLAST hit count must be an integer from 1 to 5,000.")
        if type(self.flank) is not int or not 0 <= self.flank <= 10_000_000:
            raise ValueError("Flank must be an integer from 0 to 10,000,000 bp.")
        if self.database not in {"nr", "nr_cluster_seq"}:
            raise ValueError("Unsupported protein BLAST database.")

    def blast_parameters(self) -> dict[str, str | int]:
        return {
            "PROGRAM": "blastp",
            "DATABASE": self.database,
            "HITLIST_SIZE": self.top_n,
            # Intentional Kblast settings, aligned with the standard web form
            # observed 2026-10-07. These are overrides, not URL API defaults.
            "EXPECT": 0.05,
            "GAPCOSTS": "11 1",
            "WORD_SIZE": 5,
            "SHORT_QUERY_ADJUST": "true",
        }


@dataclass
class Hit:
    rank: int
    accession: str
    subject_id: str
    title: str
    organism: str
    taxid: str
    length: int
    percent_identity: float
    query_coverage: float
    evalue: float
    bit_score: float
    hsps: list[dict[str, Any]]
    descriptions: list[dict[str, Any]] = field(default_factory=list)
    mapping_status: str = "pending"
    mapping_error: str = ""
    sequence_status: str = "pending"
    sequence_error: str = ""
    sequence: str = ""
    sequence_accession: str = ""
    mapping_resolution_version: int = 0


@dataclass
class Occurrence:
    id: str
    hit_rank: int
    hit_accession: str
    ipg_row: int
    ipg_id: str = ""
    source: str = ""
    protein_accession: str = ""
    protein_name: str = ""
    nucleotide_accession: str = ""
    cds_start: int | None = None
    cds_end: int | None = None
    strand: str = ""
    organism: str = ""
    strain: str = ""
    assembly: str = ""
    raw_mapping: dict[str, str] = field(default_factory=dict)
    status: str = "pending"
    error: str = ""
    requested_start: int | None = None
    requested_end: int | None = None
    context_start: int | None = None
    context_end: int | None = None
    record_length: int | None = None
    record_accession: str = ""
    clipped_left: bool | None = None
    clipped_right: bool | None = None
    context_file: str = ""
    context_sha256: str = ""
    retrieved_at: str = ""
    feature_count: int = 0
    target_cds_found: bool | None = None
    mapping_method: str = "IPG"


@dataclass
class Run:
    id: str
    query: Query
    settings: Settings
    created_at: str = field(default_factory=utc_now)
    updated_at: str = field(default_factory=utc_now)
    version: str = __version__
    schema_version: int = 4
    status: str = "pending"
    stage: str = "blast"
    message: str = "Queued"
    error: str = ""
    rid: str = ""
    original_rid: str = ""
    rid_source: str = ""
    rid_adopted_at: str = ""
    rid_history: list[dict[str, Any]] = field(default_factory=list)
    rid_validation: dict[str, Any] = field(default_factory=dict)
    rtoe: int = 0
    submitted_at: str = ""
    submitted_blast_parameters: dict[str, Any] = field(default_factory=dict)
    blast_completed_at: str = ""
    submission_attempted: bool = False
    last_blast_get_at: str = ""
    blast_phase: str = ""
    ncbi_status: str = ""
    blast_events: list[dict[str, Any]] = field(default_factory=list)
    blast_metadata: dict[str, Any] = field(default_factory=dict)
    hits: list[Hit] = field(default_factory=list)
    occurrences: list[Occurrence] = field(default_factory=list)
    timings: dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["query_sha256"] = digest(self.query.sequence)
        result["planned_blast_parameters"] = self.settings.blast_parameters()
        result["blast_parameters"] = (
            dict(self.blast_metadata.get("parameters", {}))
            if self.rid_source in {"external", "browser"}
            else dict(self.submitted_blast_parameters) or self.settings.blast_parameters()
        )
        result["blast_summary_definitions"] = dict(BLAST_SUMMARY_DEFINITIONS)
        result["coordinate_convention"] = (
            "1-based inclusive, native orientation, no origin wrapping"
        )
        result["counts"] = {
            "hits": len(self.hits),
            "occurrences": len(self.occurrences),
            "contexts": sum(o.status == "success" for o in self.occurrences),
            "unresolved": sum(o.status == "unresolved_mapping" for o in self.occurrences),
            "failed": sum(o.status == "retrieval_failure" for o in self.occurrences),
            "integrity_failures": sum(o.status == "integrity_failure" for o in self.occurrences),
        }
        return result

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Run:
        fields = cls.__dataclass_fields__
        values = {k: v for k, v in data.items() if k in fields}
        values["query"] = Query(**data["query"])
        saved_parameters = data.get("submitted_blast_parameters") or data.get(
            "blast_parameters", {}
        )
        settings = dict(data["settings"])
        # Early checkpoints predating ClusteredNR must retain their actual nr search.
        settings.setdefault("database", saved_parameters.get("DATABASE", "nr"))
        values["settings"] = Settings(**settings)
        if data.get("schema_version", 1) < 4 and (
            data.get("submission_attempted") or data.get("rid")
        ):
            values["submitted_blast_parameters"] = dict(saved_parameters)
        values["hits"] = [Hit(**h) for h in data["hits"]]
        # Normalize older checkpoints in memory; raw IPG and GenBank evidence is untouched.
        values["occurrences"] = [
            Occurrence(**o)
            for o in data["occurrences"]
            if usable_nucleotide_location(
                o.get("nucleotide_accession", ""), o.get("cds_start"), o.get("cds_end")
            )
        ]
        run = cls(**values)
        for o in run.occurrences:
            if o.status == "success" and o.target_cds_found is not True:
                o.status, o.error = "integrity_failure", TARGET_CDS_ERROR
        if data.get("schema_version", 1) < 2:
            run.schema_version = 2
            for hit in run.hits:
                for name, value in summarize_hsps(hit.hsps, len(run.query.sequence)).items():
                    setattr(hit, name, value)
                if hit.mapping_status in {"success", "unresolved_mapping"} and not any(
                    o.hit_rank == hit.rank for o in run.occurrences
                ):
                    hit.mapping_status, hit.mapping_error = "no_mapping", NO_MAPPING_MESSAGE
            if run.status in {"success", "partial_success"}:
                run.finish()
        run.schema_version = 4
        if run.rid and not run.rid_source:
            # Legacy checkpoints cannot prove a manually recovered RID's origin.
            recovered = any(e["event"] == "rid_received" for e in run.blast_events)
            run.rid_source = (
                "kblast_url_api" if recovered or run.submitted_blast_parameters else "external"
            )
        if run.rid_source == "kblast_url_api" and not run.original_rid:
            run.original_rid = run.rid
        for hit in run.hits:
            if hit.mapping_status == "no_mapping" and hit.mapping_resolution_version < 3:
                hit.mapping_status, hit.mapping_error = "pending", ""
                if run.status in {"success", "partial_success"}:
                    run.status = "interrupted"
                    run.message = "Resume to verify saved no-map hits with NCBI nucleotide links."
        return run

    def finish(self) -> None:
        self.stage = "complete"
        if not self.hits:
            self.status = "no_hits"
            self.message = "BLAST completed normally; no protein subject hits were found."
            return
        has_failures = any(
            h.mapping_status not in {"success", "no_mapping"} or h.sequence_status != "success"
            for h in self.hits
        ) or any(o.status != "success" for o in self.occurrences)
        self.status = "partial_success" if has_failures else "success"
        n = sum(o.status == "success" for o in self.occurrences)
        self.message = (
            f"Finished · {len(self.hits)} hits → {len(self.occurrences)} genomic occurrences "
            f"→ {n} contexts"
        )
