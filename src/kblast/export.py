"""Reviewable export bundles; source GenBank bytes and all scientific rows survive."""

import csv
import io
import json
import zipfile
from pathlib import Path

from . import __version__
from .models import Run, digest, utc_now
from .storage import safe_name

HIT_COLUMNS = [
    "hit_rank",
    "hit_accession",
    "hit_subject_id",
    "hit_title",
    "hit_organism",
    "hit_taxid",
    "hit_length",
    "percent_identity",
    "query_coverage",
    "evalue",
    "bit_score",
    "hsp_count",
    "mapping_status",
    "mapping_error",
    "protein_sequence_status",
    "protein_sequence_error",
    "protein_fasta_exported",
]
OCCURRENCE_COLUMNS = [
    "occurrence_id",
    "ipg_row",
    "ipg_id",
    "mapping_method",
    "source",
    "protein_accession",
    "protein_name",
    "nucleotide_accession",
    "cds_start",
    "cds_end",
    "strand",
    "organism",
    "strain",
    "assembly",
    "requested_start",
    "requested_end",
    "record_accession",
    "record_length",
    "context_start",
    "context_end",
    "clipped_left",
    "clipped_right",
    "status",
    "error",
    "retrieved_at",
    "feature_count",
    "target_cds_found",
    "context_sha256",
    "context_exported",
    "genbank_file",
]
COLUMNS = [
    "run_id",
    "rid",
    "original_rid",
    "rid_source",
    "rid_adopted_at",
    "query_id",
    "query_sha256",
    "row_type",
    *HIT_COLUMNS,
    *OCCURRENCE_COLUMNS,
]


def selected_occurrences(run: Run, selected: list[str] | None):
    successful = {
        o.id: o for o in run.occurrences if o.status == "success" and o.target_cds_found is True
    }
    if selected is not None:
        if len(set(selected)) != len(selected) or any(i not in successful for i in selected):
            raise ValueError(
                "Export selection must contain unique, successful context identifiers."
            )
        chosen = set(selected)
    else:
        chosen = set(successful)
    return [o for o in run.occurrences if o.id in chosen]


def context_filename(o) -> str:
    return f"contexts/{o.id}_{safe_name(o.record_accession)}_{o.context_start}_{o.context_end}.gb"


def metadata_tsv(run: Run, selected: list[str]) -> str:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    chosen = set(selected)
    for hit in run.hits:
        common = {
            "run_id": run.id,
            "rid": run.rid,
            "original_rid": run.original_rid,
            "rid_source": run.rid_source,
            "rid_adopted_at": run.rid_adopted_at,
            "query_id": run.query.identifier,
            "query_sha256": digest(run.query.sequence),
            "hit_rank": hit.rank,
            "hit_accession": hit.accession,
            "hit_subject_id": hit.subject_id,
            "hit_title": hit.title,
            "hit_organism": hit.organism,
            "hit_taxid": hit.taxid,
            "hit_length": hit.length,
            "percent_identity": hit.percent_identity,
            "query_coverage": hit.query_coverage,
            "evalue": hit.evalue,
            "bit_score": hit.bit_score,
            "hsp_count": len(hit.hsps),
            "mapping_status": hit.mapping_status,
            "mapping_error": hit.mapping_error,
            "protein_sequence_status": hit.sequence_status,
            "protein_sequence_error": hit.sequence_error,
            "protein_fasta_exported": hit.sequence_status == "success",
        }
        # One hit row also preserves subjects that have no occurrence rows at all.
        writer.writerow({**common, "row_type": "hit", "status": hit.mapping_status})
        for o in run.occurrences:
            if o.hit_rank != hit.rank:
                continue
            values = {
                column: getattr(o, column) for column in OCCURRENCE_COLUMNS if hasattr(o, column)
            }
            values.update(
                occurrence_id=o.id,
                context_exported=o.id in chosen,
                genbank_file=context_filename(o) if o.id in chosen else "",
            )
            writer.writerow({**common, "row_type": "occurrence", **values})
    return output.getvalue()


def protein_fasta(run: Run) -> str:
    records = []
    for hit in run.hits:
        if hit.sequence_status != "success" or not hit.sequence:
            continue
        title = " ".join(hit.title.split())
        header = f">{hit.accession} Kblast_rank={hit.rank} {title}\n"
        records.append(
            header
            + "\n".join(hit.sequence[i : i + 70] for i in range(0, len(hit.sequence), 70))
            + "\n"
        )
    return "".join(records)


def export_bundle(directory: Path, run: Run, selected: list[str] | None = None) -> bytes:
    chosen = selected_occurrences(run, selected)
    selected_ids = [o.id for o in chosen]
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for o in chosen:
            raw = (directory / o.context_file).read_bytes()
            if digest(raw) != o.context_sha256:
                raise ValueError(
                    f"Stored context {o.id} failed its integrity check. Resume to retrieve it."
                )
            archive.writestr(context_filename(o), raw)
        with archive.open("contexts.gb", "w") as combined:
            for o in chosen:
                raw = (directory / o.context_file).read_bytes()
                if digest(raw) != o.context_sha256:
                    raise ValueError(f"Stored context {o.id} failed its integrity check.")
                combined.write(raw.rstrip(b"\r\n") + b"\n")
        archive.writestr("blast_hits.fasta", protein_fasta(run))
        archive.writestr("metadata.tsv", metadata_tsv(run, selected_ids))
        archive.writestr("query.fasta", f">{run.query.description}\n{run.query.sequence}\n")
        provenance = run.to_dict()
        provenance["export"] = {
            "at": utc_now(),
            "context_ids": selected_ids,
            "protein_hits_exported": sum(h.sequence_status == "success" for h in run.hits),
            "protein_hits_missing": [h.rank for h in run.hits if h.sequence_status != "success"],
            "note": "Metadata contains every hit and qualifying genomic occurrence, including "
            "unselected rows. Complete raw IPG responses retain excluded non-genomic rows.",
        }
        archive.writestr("provenance/run.json", json.dumps(provenance, indent=2))
        # Explicit allowlist keeps configuration and user originals out of the bundle.
        for path in sorted((directory / "raw").glob("*")):
            if (
                path.name == "blast.xml"
                or (path.name.startswith("blast-") and path.name.endswith(".xml"))
                or path.name == "blast-submission.txt"
                or path.name.endswith(".ipg.tsv")
                or path.name.endswith(".ipg-error.txt")
                or path.name.endswith(".ipg.xml")
                or path.name.endswith(".links.xml")
                or path.name.endswith(".linked-search.xml")
                or path.name.endswith(".linked.gb")
            ):
                archive.writestr(f"provenance/{path.name}", path.read_bytes())
        for path in sorted((directory / "cache").glob("*.record.json")):
            archive.writestr(f"provenance/nucleotide_records/{path.name}", path.read_bytes())
        log = directory / "requests.jsonl"
        if log.is_file():
            archive.writestr("provenance/requests.jsonl", log.read_bytes())
        lifecycle = directory / "blast-lifecycle.jsonl"
        if lifecycle.is_file():
            archive.writestr("provenance/blast-lifecycle.jsonl", lifecycle.read_bytes())
        archive.writestr(
            "READ_ME.txt",
            (
                f"Kblast {__version__} gene neighborhoods\n\n"
                "contexts/: one unmodified NCBI GenBank region per selected successful IPG row.\n"
                "contexts.gb: those same records concatenated in BLAST rank / IPG row order.\n"
                "blast_hits.fasta: full sequences of every successfully fetched BLAST subject.\n"
                "metadata.tsv: all hit/occurrence rows, statuses, errors, export selection.\n"
                "provenance/: search settings, all HSPs, raw BLAST and IPG, redacted request log.\n"
                "Coordinates are 1-based inclusive in the native nucleotide orientation.\n"
                "GenBank feature coordinates are relative to the extracted sequence.\n"
                "BLAST subject coverage: union of query positions across all HSPs.\n"
                "Identity and bit score: maximum HSP values. E-value: minimum HSP value.\n"
                "Protein-only IPG rows stay in raw provenance, not genomic occurrences.\n"
                "Contexts lacking a verified target CDS are excluded from successful exports.\n"
                "No biological deduplication: repeated accessions and regions are intentional.\n"
                "Missing protein sequences are listed in provenance/run.json and metadata.tsv.\n"
            ),
        )
    return output.getvalue()
