import copy
import io
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

from Bio import BiopythonParserWarning, BiopythonWarning, SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition, SimpleLocation

from .. import __version__
from ..cancellation import CheckedText, check_cancelled
from ..models import Protein
from ..normalization import normalize_sequence

# Keep only structural/identity qualifiers. Functional legacy fields are removed by
# construction, including unknown annotation-tool-specific qualifiers.
PRESERVED_CDS_QUALIFIERS = {
    "locus_tag",
    "translation",
    "codon_start",
    "transl_table",
    "transl_except",
    "exception",
    "pseudo",
    "pseudogene",
    "partial",
    "ribosomal_slippage",
}


@dataclass
class GenBankInput:
    records: list
    proteins: list[Protein]
    cds_indices: list[tuple[int, int]]
    warnings: list[str]
    raw_locations: dict = field(default_factory=dict)


def raw_feature_locations(text):
    locations = {}
    ri, fi = 0, -1
    collecting = False
    for line in text.splitlines(keepends=True):
        check_cancelled()
        if line.startswith("//"):
            ri, fi, collecting = ri + 1, -1, False
        elif re.match(r"^     \S", line):
            fi += 1
            locations[(ri, fi)] = [line]
            collecting = True
        elif collecting and line.startswith("                     "):
            if line[21:].lstrip().startswith("/"):
                collecting = False
            else:
                locations[(ri, fi)].append(line)
        else:
            collecting = False
    return locations


def strip_functional_qualifiers(feature):
    feature.qualifiers = {
        key: value for key, value in feature.qualifiers.items() if key in PRESERVED_CDS_QUALIFIERS
    }


def derive_translation(record, feature) -> str:
    if (
        feature.location is None
        or "pseudo" in feature.qualifiers
        or "pseudogene" in feature.qualifiers
    ):
        return ""
    try:
        offset = int(feature.qualifiers.get("codon_start", ["1"])[0]) - 1
        table_id = int(feature.qualifiers.get("transl_table", ["1"])[0])
        if offset not in (0, 1, 2):
            return ""
        # Biopython extracts joins, complements and clipped boundary locations.
        nucleotide = str(feature.extract(record.seq))[offset:]
        nucleotide = nucleotide[: len(nucleotide) // 3 * 3]
        table = CodonTable.unambiguous_dna_by_id[table_id]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            amino = str(Seq(nucleotide).translate(table=table_id))
        start = feature.location.end if feature.location.strand == -1 else feature.location.start
        if (
            amino
            and offset == 0
            and isinstance(start, ExactPosition)
            and nucleotide[:3].upper() in table.start_codons
        ):
            amino = "M" + amino[1:]
        return amino.rstrip("*")
    except (ValueError, TypeError, KeyError, AttributeError):
        return ""


def _parse_genbank(text, protein_offset=0) -> GenBankInput:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", BiopythonParserWarning)
        records = list(SeqIO.parse(CheckedText(io.StringIO(text)), "genbank"))
    if not records:
        raise ValueError("No GenBank records found")
    proteins = []
    indices = []
    missing_location = 0
    for ri, record in enumerate(records):
        for fi, feature in enumerate(record.features):
            check_cancelled()
            if feature.type != "CDS":
                continue
            if "translation" in feature.qualifiers:
                # Retain the original qualifier verbatim; normalize only the search copy.
                raw = "".join(feature.qualifiers["translation"])
            else:
                raw = derive_translation(record, feature)
                feature.qualifiers["translation"] = [raw]
            if feature.location is None:
                missing_location += 1
            tag = feature.qualifiers.get(
                "locus_tag", [f"CDS_{protein_offset + len(proteins) + 1}"]
            )[0]
            proteins.append(
                Protein(
                    str(protein_offset + len(proteins)),
                    f"{record.id} {tag}",
                    normalize_sequence(raw).sequence,
                )
            )
            indices.append((ri, fi))
    messages = []
    if missing_location:
        messages.append(f"{missing_location} CDSs have unreadable locations and were retained")
    # Boundary clipping and imperfect sequences are expected in regional reannotation.
    significant = [
        w
        for w in caught
        if "location" not in str(w.message).lower()
        and "sequence length" not in str(w.message).lower()
    ]
    if significant:
        messages.append(f"GenBank parser tolerated {len(significant)} record-format irregularities")
    return GenBankInput(records, proteins, indices, messages, raw_feature_locations(text))


def _combine(parts):
    result = GenBankInput([], [], [], [], {})
    for part in parts:
        offset = len(result.records)
        result.records.extend(part.records)
        result.proteins.extend(part.proteins)
        result.cds_indices.extend((ri + offset, fi) for ri, fi in part.cds_indices)
        result.raw_locations.update(
            {(ri + offset, fi): lines for (ri, fi), lines in part.raw_locations.items()}
        )
        result.warnings.extend(part.warnings)
    return result


def iter_genbank(path: Path, *, batch_size=512):
    """Read each record once, bounding retained records except one oversized record."""
    parts, lines = [], []
    count = offset = records = 0
    with path.open() as handle:
        for line in handle:
            check_cancelled()
            lines.append(line)
            if not line.startswith("//"):
                continue
            part = _parse_genbank("".join(lines), offset)
            lines.clear()
            if parts and count + len(part.proteins) > batch_size:
                yield _combine(parts)
                parts, count = [], 0
            parts.append(part)
            count += len(part.proteins)
            offset += len(part.proteins)
            records += len(part.records)
            if count >= batch_size or len(parts) >= 128:
                yield _combine(parts)
                parts, count = [], 0
        text = "".join(lines)
        if text.strip():
            part = _parse_genbank(text, offset)
            parts.append(part)
            records += len(part.records)
    if parts:
        yield _combine(parts)
    if not records:
        raise ValueError("No GenBank records found")


def read_genbank(path: Path) -> GenBankInput:
    return _combine(iter_genbank(path))


def parse_genbank(text: str) -> GenBankInput:
    """Parse local in-memory input with the same CDS rules as file annotation."""
    return _parse_genbank(text)


def spreadsheet_context(data: GenBankInput, input_name: str):
    """Source fields and reported 1-based inclusive bounds, in CDS order.

    A joined CDS uses its outer bounds. Retain the original expression as a cell
    note when those bounds cannot express its full semantics. Remote coordinates
    are not local coordinates, and unknown positions/strands are never invented.
    """
    for ri, fi in data.cds_indices:
        check_cancelled()
        record = data.records[ri]
        feature = record.features[fi]
        location = feature.location
        expression = "".join(line[21:].strip() for line in data.raw_locations.get((ri, fi), []))
        start = end = strand = None
        note = None
        if location is not None:
            strand = {1: "+", -1: "-"}.get(location.strand)
            remote = any(part.ref or part.ref_db for part in location.parts)
            if not remote:
                try:
                    start, end = int(location.start) + 1, int(location.end)
                    if start > end:  # Between-base/empty features have no inclusive span.
                        start = end = None
                except (TypeError, ValueError):
                    start = end = None
            simple_exact = (
                isinstance(location, SimpleLocation)
                and isinstance(location.start, ExactPosition)
                and isinstance(location.end, ExactPosition)
                and start is not None
                and strand is not None
            )
            if not simple_exact:
                note = (
                    f"Input CDS location: {expression or str(location)}\n"
                    "From/To are the lower/upper reported bounds, 1-based inclusive. "
                    "Joins and partial/uncertain boundaries are not represented by the "
                    "two numbers. Remote or unknown bounds and unknown/mixed strands "
                    "are blank. Keep the source GenBank for the complete location."
                )
        else:
            note = f"Unreadable input CDS location: {expression or 'not available'}"

        def qualifier(name, feature=feature):
            values = feature.qualifiers.get(name, [])
            return "; ".join(values) if values else None

        yield (
            [
                input_name,
                record.id,
                qualifier("locus_tag"),
                start,
                end,
                strand,
                qualifier("product"),
            ],
            note,
        )


def annotated_records(data: GenBankInput, annotations, *, provenance=True):
    if len(annotations) != len(data.cds_indices):
        raise ValueError("CDS and annotation counts differ")
    if [a.protein.id for a in annotations] != [p.id for p in data.proteins]:
        raise ValueError("CDS and annotation order differs")
    records = []
    for record in data.records:
        check_cancelled()
        records.append(copy.deepcopy(record))
    for (ri, fi), annotation in zip(data.cds_indices, annotations):
        check_cancelled()
        feature = records[ri].features[fi]
        strip_functional_qualifiers(feature)
        feature.qualifiers["product"] = [annotation.product]
    for record in records:
        comment = record.annotations.get("comment", "")
        comment = re.sub(r"(?m)^Annotated with Kannotate v[^\n]*\n?", "", comment).strip()
        if comment:
            record.annotations["comment"] = comment
        else:
            record.annotations.pop("comment", None)
        record.annotations.setdefault("molecule_type", "DNA")
    if provenance:
        note = f"Annotated with Kannotate v{__version__}"
        existing = records[0].annotations.get("comment", "")
        records[0].annotations["comment"] = f"{existing}\n{note}".strip()
    return records


def write_genbank_handle(data: GenBankInput, annotations, handle, *, provenance=True):
    records = annotated_records(data, annotations, provenance=provenance)
    # Temporarily serialize unknown locations, then restore their exact input lines.
    # This retains a malformed feature without manufacturing coordinates.
    restore = {
        (ri, fi): data.raw_locations[(ri, fi)]
        for ri, record in enumerate(records)
        for fi, feature in enumerate(record.features)
        if feature.location is None and (ri, fi) in data.raw_locations
    }
    if any(f.location is None for r in records for f in r.features):
        for ri, record in enumerate(records):
            for fi, feature in enumerate(record.features):
                if feature.location is None:
                    if (ri, fi) not in restore:
                        raise ValueError(
                            "An unreadable feature has no original location text to preserve"
                        )
                    feature.location = SimpleLocation(0, 1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        output = io.StringIO() if restore else handle
        SeqIO.write(records, CheckedText(output), "genbank")
    if restore:
        lines = output.getvalue().splitlines(keepends=True)
        ri, fi = 0, -1
        rewritten = []
        skipping = False
        for line in lines:
            check_cancelled()
            if line.startswith("//"):
                ri, fi, skipping = ri + 1, -1, False
            elif re.match(r"^     \S", line):
                fi += 1
                skipping = (ri, fi) in restore
                if skipping:
                    rewritten.extend(restore[(ri, fi)])
                    continue
            elif skipping and line.startswith("                     "):
                if not line[21:].lstrip().startswith("/"):
                    continue
                skipping = False
            rewritten.append(line)
        handle.write("".join(rewritten))


def write_genbank(data: GenBankInput, annotations, path: Path):
    with path.open("w") as handle:
        write_genbank_handle(data, annotations, handle)
