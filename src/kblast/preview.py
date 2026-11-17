"""Lightweight annotated-feature preview in native genomic coordinates."""

import io
from pathlib import Path

from Bio import SeqIO

from .models import Occurrence, digest


def context_preview(directory: Path, occurrence: Occurrence) -> dict:
    if occurrence.status != "success" or occurrence.target_cds_found is not True:
        raise ValueError(
            "This context lacks a verified target CDS and is not a successful context."
        )
    raw = (directory / occurrence.context_file).read_bytes()
    if digest(raw) != occurrence.context_sha256:
        raise ValueError(
            "Stored context failed its integrity check. Resume this run to retrieve it."
        )
    record = SeqIO.read(io.StringIO(raw.decode()), "genbank")
    features = []
    omitted_features = 0
    has_cds = any(f.type == "CDS" for f in record.features)
    wanted = {"CDS", "tRNA", "rRNA", "ncRNA", "misc_RNA"}
    if not has_cds:
        wanted.add("gene")
    for feature in record.features:
        if feature.type not in wanted or feature.location is None:
            continue
        location = feature.location
        try:
            local_start, local_end = int(location.start), int(location.end)
            if any(part.ref for part in location.parts) or not 0 <= local_start < local_end <= len(
                record
            ):
                omitted_features += 1
                continue
        except (TypeError, ValueError):
            omitted_features += 1
            continue
        start = local_start + occurrence.context_start
        end = local_end + occurrence.context_start - 1
        qualifiers = feature.qualifiers
        strand = location.strand
        target = (
            feature.type == "CDS"
            and start == occurrence.cds_start
            and end == occurrence.cds_end
            and ({"+": 1, "-": -1}.get(occurrence.strand, strand) == strand)
        )
        features.append(
            {
                "index": len(features),
                "type": feature.type,
                "start": start,
                "end": end,
                "strand": strand,
                "target": target,
                "label": next(
                    iter(qualifiers.get("gene", qualifiers.get("locus_tag", [feature.type])))
                ),
                "gene": "; ".join(qualifiers.get("gene", [])),
                "locus_tag": "; ".join(qualifiers.get("locus_tag", [])),
                "protein_id": "; ".join(qualifiers.get("protein_id", [])),
                "product": "; ".join(qualifiers.get("product", [])),
                "location": str(location),
                "partial": "<" in str(location) or ">" in str(location),
                "parts": [
                    {
                        "start": int(p.start) + occurrence.context_start,
                        "end": int(p.end) + occurrence.context_start - 1,
                        "strand": p.strand,
                    }
                    for p in location.parts
                ],
            }
        )
    return {
        "occurrence_id": occurrence.id,
        "record_id": record.id,
        "length": len(record),
        "start": occurrence.context_start,
        "end": occurrence.context_end,
        "organism": record.annotations.get("organism", occurrence.organism),
        "molecule_type": record.annotations.get("molecule_type", ""),
        "features": features,
        "target_cds_found": occurrence.target_cds_found,
        "omitted_features": omitted_features,
    }
