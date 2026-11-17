"""GenBank import. Nucleotides never cross this module's output boundary."""
from __future__ import annotations

import hashlib
import re
import warnings
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import ExactPosition


class ImportError(ValueError):
    pass


def parse_genbank(path: str | Path, display_name: str | None = None, slot: int = 0) -> dict:
    path = Path(path)
    notices: list[str] = []
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with path.open() as handle:
                records = list(SeqIO.parse(handle, "genbank"))
            notices.extend(str(w.message) for w in caught)
    except Exception as exc:
        raise ImportError(f"Cannot read {path.name}: {exc}") from exc
    if len(records) != 1:
        raise ImportError(f"Exactly one GenBank record is required; found {len(records)}.")
    record = records[0]
    if not 0 < len(record) <= 100_000:
        raise ImportError(f"Record length {len(record):,} bp is outside 1–100,000 bp.")
    # Slot disambiguates identical inputs. IDs do not depend on editable display names.
    digest = hashlib.sha256(path.read_bytes()).hexdigest()[:12]
    cid = f"c{slot:02d}_{digest}"
    cluster = {
        "id": cid, "name": display_name or path.stem, "source_name": path.name,
        "accession": record.id, "description": record.description,
        "length": len(record), "topology": record.annotations.get("topology", "linear").lower(),
        "annotations": {k: v for k, v in record.annotations.items()
                        if k in ("organism", "taxonomy", "molecule_type", "date", "accessions", "keywords", "source")},
        "cds": [], "warnings": notices,
    }
    for ordinal, feature in enumerate(f for f in record.features if f.type == "CDS"):
        qualifiers = {str(k): [str(v) for v in values] for k, values in feature.qualifiers.items()}
        problems: list[str] = []
        segments: list[dict] = []
        partial = False
        if feature.location is None:
            problems.append("Uninterpretable CDS location; retained without drawable coordinates.")
        else:
            for part in feature.location.parts:
                partial |= not isinstance(part.start, ExactPosition) or not isinstance(part.end, ExactPosition)
                try:
                    start, end = int(part.start), int(part.end)
                    if part.ref or start < 0 or end > len(record) or end <= start:
                        raise ValueError("remote, empty, or out-of-bounds segment")
                    segments.append({"start": start, "end": end, "strand": part.strand or 0})
                except (ValueError, TypeError) as exc:
                    problems.append(f"Problematic CDS segment: {exc}.")
        if partial:
            problems.append("Partial CDS; fuzzy boundary retained at its stated coordinate.")
        strand = feature.location.strand if feature.location else None
        if strand not in (-1, 1):
            problems.append("Unknown or mixed strand.")
        protein = ""
        origin = "unavailable"
        if qualifiers.get("translation"):
            origin = "translation"
            protein = re.sub(r"\s+", "", "".join(qualifiers["translation"])).upper()
        elif feature.location is not None:
            origin = "calculated"
            try:
                nucleotide = feature.extract(record.seq)
                codon_start = int(qualifiers.get("codon_start", ["1"])[0])
                if codon_start not in (1, 2, 3):
                    raise ValueError("codon_start must be 1, 2, or 3")
                nucleotide = nucleotide[codon_start - 1:]
                table = int(qualifiers.get("transl_table", ["1"])[0])
                if len(nucleotide) % 3:
                    problems.append("Incomplete terminal codon omitted from calculated translation.")
                nucleotide = nucleotide[:len(nucleotide) // 3 * 3]
                # cds=True correctly translates alternative start codons for complete CDSs.
                try:
                    protein = str(nucleotide.translate(table=table, cds=not partial))
                except Exception:
                    protein = str(nucleotide.translate(table=table))
                    if not partial:
                        problems.append("Complete-CDS translation checks failed; permissive translation retained.")
                if qualifiers.get("transl_except"):
                    problems.append("transl_except present; calculated translation does not apply exceptional codons.")
            except Exception as exc:
                problems.append(f"Protein extraction failed: {exc}.")
        if protein.endswith("*"):
            protein = protein[:-1]
        if "*" in protein:
            problems.append("Internal stop codon; preserved in the protein, masked as X for search.")
        if re.search(r"[^A-Z*]", protein):
            problems.append("Invalid protein characters; masked as X.")
            protein = re.sub(r"[^A-Z*]", "X", protein)
        if not protein:
            problems.append("No interpretable protein; CDS retained but excluded from homology/FASTA.")
        elif len(protein) < 20:
            problems.append("Very short protein; interpret homology evidence cautiously.")
        gene = {
            "id": f"{cid}_g{ordinal + 1:04d}", "cluster_id": cid,
            "gene": qualifiers.get("gene", [""])[0],
            "locus_tag": qualifiers.get("locus_tag", [""])[0],
            "product": qualifiers.get("product", [""])[0],
            "segments": segments, "location": str(feature.location),
            "strand": strand or 0, "cds_length": sum(s["end"] - s["start"] for s in segments),
            "protein": protein, "protein_length": len(protein), "translation_source": origin,
            "partial": partial, "warnings": problems, "qualifiers": qualifiers, "group_id": None,
            "origin_spanning": bool(segments and len(segments) > 1
                                    and any(s["start"] == 0 for s in segments)
                                    and any(s["end"] == len(record) for s in segments)),
        }
        cluster["cds"].append(gene)
    if not cluster["cds"]:
        cluster["warnings"].append("No CDS features found.")
    # Preserve GenBank feature order; rendering positions remain genomic.
    return cluster


def inspect_inputs(paths: list[str]) -> list[dict]:
    if not 1 <= len(paths) <= 50:
        raise ImportError("Select between 1 and 50 GenBank files.")
    results = []
    for slot, path in enumerate(paths):
        try:
            cluster = parse_genbank(path, slot=slot)
            results.append({"path": path, "cluster": cluster, "error": None})
        except ImportError as exc:
            results.append({"path": path, "cluster": None, "error": str(exc)})
    return results
