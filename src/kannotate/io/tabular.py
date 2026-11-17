"""One ordered, source-preserving table schema for CSV and legacy workbooks."""

from ..cancellation import check_cancelled
from .genbank import spreadsheet_context

HEADERS = ["FASTA header", "Normalized/search sequence", "Kannotate annotation", "Kannotate score"]
for engine in ("HMM", "RPS"):
    for rank in (1, 2):
        HEADERS.extend(
            [
                f"{engine} hit {rank} profile",
                f"{engine} hit {rank} description",
                f"{engine} hit {rank} E-value",
                f"{engine} hit {rank} query span",
                f"{engine} hit {rank} Kannotate score",
            ]
        )

GENBANK_HEADERS = [
    "Input filename",
    "Record ID",
    "Locus tag",
    "From",
    "To",
    "Strand",
    "Original product",
] + HEADERS[1:]


def annotation_rows(annotations, *, genbank=None, input_name=None):
    """Yield values and a location note; never rewrite source fields for display."""
    contexts = None
    if genbank is not None:
        if not input_name:
            raise ValueError("GenBank tables require the input filename")
        if len(annotations) != len(genbank.cds_indices):
            raise ValueError("CDS and annotation counts differ")
        if [a.protein.id for a in annotations] != [p.id for p in genbank.proteins]:
            raise ValueError("CDS and annotation order differs")
        contexts = iter(spreadsheet_context(genbank, input_name))
    for annotation in annotations:
        check_cancelled()
        context, note = (
            next(contexts) if contexts is not None else ([annotation.protein.header], None)
        )
        values = context + [annotation.protein.sequence, annotation.product, annotation.score]
        for hits in (annotation.hmm_hits, annotation.rps_hits):
            for rank in range(2):
                hit = hits[rank] if len(hits) > rank else None
                values.extend(
                    [
                        f"{hit.name} ({hit.profile_id})",
                        hit.description,
                        hit.evalue,
                        f"{hit.start}-{hit.end}",
                        hit.score,
                    ]
                    if hit
                    else [None] * 5
                )
        yield values, note
