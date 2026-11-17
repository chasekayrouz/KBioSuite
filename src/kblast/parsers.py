"""Interpret official NCBI output without regrouping subjects or IPG rows."""

import csv
import io
import re
import xml.etree.ElementTree as ET
from typing import Any

from .models import Hit, Occurrence, summarize_hsps, usable_nucleotide_location

IPG_COLUMNS = [
    "Id",
    "Source",
    "Nucleotide Accession",
    "Start",
    "Stop",
    "Strand",
    "Protein",
    "Protein Name",
    "Organism",
    "Strain",
    "Assembly",
]


class ResponseParseError(ValueError):
    """A response does not establish a valid, complete scientific mapping."""


def ipg_xml_reports(text: str) -> list[tuple[str, set[str], str]]:
    """Keep report boundaries (also when two requested proteins share one IPG)."""
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise ResponseParseError("Malformed NCBI IPG XML.") from exc
    if root.tag != "IPGReportSet" or root.find(".//ERROR") is not None:
        raise ResponseParseError("NCBI did not return an IPGReportSet.")
    reports = []
    for report in root.findall("IPGReport"):
        proteins = report.findall("ProteinList/Protein")
        stats = report.find("Statistics")
        if stats is None or not report.get("ipg"):
            raise ResponseParseError("IPG report lacks identity or completeness statistics.")
        try:
            if int(stats.attrib["prot_count"]) != len(proteins) or int(
                stats.attrib["nuc_count"]
            ) != len(report.findall("ProteinList/Protein/CDSList/CDS")):
                raise ResponseParseError("IPG report counts disagree with its contents.")
        except (KeyError, ValueError) as exc:
            raise ResponseParseError("Invalid IPG completeness statistics.") from exc
        output = io.StringIO()
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(IPG_COLUMNS)
        members = set()
        for protein in proteins:
            accession = protein.get("accver", "")
            if not accession:
                raise ResponseParseError("IPG protein has no accession.")
            members.add(accession)
            for cds in protein.findall("CDSList/CDS") or [None]:
                attrs = cds.attrib if cds is not None else {}
                writer.writerow(
                    [
                        report.get("ipg"),
                        protein.get("source", ""),
                        attrs.get("accver", "N/A"),
                        attrs.get("start", ""),
                        attrs.get("stop", ""),
                        attrs.get("strand", ""),
                        accession,
                        protein.get("name", ""),
                        attrs.get("org", protein.get("org", "")),
                        attrs.get("strain", ""),
                        attrs.get("assembly", ""),
                    ]
                )
        reports.append((report.get("product_acc", ""), members, output.getvalue()))
    return reports


def parse_blast_xml(
    text: str,
    top_n: int,
    expected_query_length: int,
    expected_database: str = "nr",
    *,
    expected_query_sequence: str | None = None,
) -> tuple[list[Hit], dict]:
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise ResponseParseError("Malformed NCBI BLAST XML2 result.") from exc
    for element in root.iter():
        element.tag = element.tag.rsplit("}", 1)[-1]
    if root.tag != "BlastXML2" and root.tag != "BlastOutput2":
        raise ValueError("Expected single-query NCBI BLAST XML2 output.")
    report = root.find(".//Report")
    search = root.find(".//Search")
    if report is None or search is None:
        raise ValueError("BLAST XML2 contains no search report.")
    if len(root.findall(".//Search")) != 1 or len(root.findall(".//Report")) != 1:
        raise ValueError("Expected exactly one BLAST query result.")
    program = report.findtext("program")
    database = report.findtext("search-target/Target/db")
    if program != "blastp" or database != expected_database:
        raise ValueError(
            f"BLAST program/database mismatch: expected blastp/{expected_database}, "
            f"NCBI returned {program}/{database}."
        )
    query_length = int(search.findtext("query-len", "0"))
    if query_length != expected_query_length:
        raise ValueError("The BLAST report query length does not match this run.")
    if expected_query_sequence is not None and len(expected_query_sequence) != query_length:
        raise ValueError("The expected query sequence length is inconsistent.")
    if search.find("hits") is None:
        raise ValueError("BLAST XML2 lacks its subject hits element.")
    hits = []
    checked_positions = set()
    # Validate every returned subject, including those beyond the consumed Top N.
    for rank, element in enumerate(search.findall("hits/Hit"), 1):
        descriptions = [
            {child.tag: child.text or "" for child in desc}
            for desc in element.findall("description/HitDescr")
        ]
        if not descriptions:
            raise ValueError(f"BLAST hit {rank} has no subject description.")
        desc = descriptions[0]
        hsps: list[dict[str, Any]] = []
        integer_fields = {
            "num",
            "query-from",
            "query-to",
            "hit-from",
            "hit-to",
            "identity",
            "positive",
            "gaps",
            "align-len",
            "query-frame",
            "hit-frame",
        }
        float_fields = {"bit-score", "score", "evalue"}
        for hsp in element.findall("hsps/Hsp"):
            values: dict[str, Any] = {}
            for child in hsp:
                value = child.text or ""
                values[child.tag] = (
                    int(value)
                    if child.tag in integer_fields
                    else float(value)
                    if child.tag in float_fields
                    else value
                )
            hsps.append(values)
        summary = summarize_hsps(hsps, query_length)
        if expected_query_sequence is not None:
            for hsp in hsps:
                aligned = hsp.get("qseq", "").replace("-", "").upper()
                if not aligned:
                    continue
                low, high = hsp["query-from"], hsp["query-to"]
                if low > high or len(aligned) != high - low + 1:
                    raise ValueError("BLAST HSP query sequence disagrees with its coordinates.")
                for position, returned in enumerate(aligned, low - 1):
                    expected = expected_query_sequence[position].upper()
                    # Ambiguous/rare residues and masked X cannot reliably prove identity.
                    canonical = "ACDEFGHIKLMNPQRSTVWY"
                    if returned in canonical and expected in canonical:
                        if returned != expected:
                            raise ValueError(
                                "The BLAST report aligned query differs from this run."
                            )
                        checked_positions.add(position)
        organism = desc.get("sciname", "")
        if not organism:
            match = re.search(r"\[([^\[\]]+)\]\s*$", desc.get("title", ""))
            organism = match.group(1) if match else ""
        accession = desc.get("accession", "")
        if not accession:
            raise ValueError(f"BLAST hit {rank} has no accession.")
        if "." not in accession:
            for identifier in desc.get("id", "").split("|"):
                if re.fullmatch(re.escape(accession) + r"\.\d+", identifier):
                    accession = identifier
                    break
        hits.append(
            Hit(
                rank=rank,
                accession=accession,
                subject_id=desc.get("id", ""),
                title=desc.get("title", ""),
                organism=organism,
                taxid=desc.get("taxid", ""),
                length=int(element.findtext("len", "0")),
                **summary,
                hsps=hsps,
                descriptions=descriptions,
            )
        )
    metadata = {
        "program": report.findtext("program"),
        "version": report.findtext("version"),
        "database": report.findtext("search-target/Target/db"),
        "query_length": query_length,
        "parameters": {e.tag: e.text for e in report.findall("params/Parameters/*")},
        "statistics": {e.tag: e.text for e in search.findall("stat/Statistics/*")},
        "message": search.findtext("message", ""),
        "returned_hit_count": len(hits),
        "consumed_hit_count": min(top_n, len(hits)),
        "requested_top_n": top_n,
        "query_id": search.findtext("query-id", ""),
        "query_title": search.findtext("query-title", ""),
        "query_validation": {
            "method": "length_and_aligned_residues" if checked_positions else "length_only",
            "checked_positions": len(checked_positions),
            "query_length": query_length,
            "residual_risk": (
                "XML2 does not expose the complete original query. Equal length alone cannot "
                "prove identity; unaligned, ambiguous or masked residues remain unchecked."
                if len(checked_positions) < query_length
                else ""
            ),
        },
    }
    return hits[:top_n], metadata


def parse_ipg(text: str, hit: Hit, flank: int) -> list[Occurrence]:
    if not text.strip():
        raise ResponseParseError(
            "NCBI returned an empty IPG response, not a verified empty report."
        )
    # EFetch sometimes returns a plain error with HTTP 200. Never treat that as no mapping.
    if "\t" not in text.splitlines()[0]:
        raise ResponseParseError("NCBI did not return an IPG tabular report.")
    reader = csv.DictReader(io.StringIO(text.lstrip("\ufeff")), delimiter="\t")
    if not {"Id", "Nucleotide Accession", "Start", "Stop", "Protein"} <= set(
        reader.fieldnames or []
    ):
        raise ResponseParseError("Unrecognized NCBI IPG columns.")
    results = []
    for row_number, row in enumerate(reader, 1):
        if not any(row.values()):
            continue
        if None in row or any(v is None for v in row.values()):
            raise ResponseParseError(f"Truncated or overlong IPG row {row_number}.")
        raw = {str(k): v or "" for k, v in row.items() if k is not None}
        accession = raw.get("Nucleotide Accession", "").strip()
        if accession in {"", "-", "N/A"}:
            continue
        try:
            start, end = int(raw["Start"]), int(raw["Stop"])
        except (ValueError, KeyError) as exc:
            raise ResponseParseError(
                f"Invalid CDS coordinates in nucleotide-linked IPG row {row_number}."
            ) from exc
        if not usable_nucleotide_location(accession, start, end):
            raise ResponseParseError(f"Invalid nucleotide location in IPG row {row_number}.")
        low, high = min(start, end), max(start, end)
        o = Occurrence(
            id=f"h{hit.rank:04d}_o{row_number:06d}",
            hit_rank=hit.rank,
            hit_accession=hit.accession,
            ipg_row=row_number,
            ipg_id=raw.get("Id", ""),
            source=raw.get("Source", ""),
            protein_accession=raw.get("Protein", ""),
            protein_name=raw.get("Protein Name", ""),
            nucleotide_accession=accession,
            cds_start=low,
            cds_end=high,
            requested_start=low - flank,
            requested_end=high + flank,
            strand=raw.get("Strand", ""),
            organism=raw.get("Organism", ""),
            strain=raw.get("Strain", ""),
            assembly=raw.get("Assembly", ""),
            raw_mapping=raw,
        )
        results.append(o)
    return results
