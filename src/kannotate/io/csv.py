"""Buffered, rectangular CSV output with lossless quoting."""

import csv
from pathlib import Path

from .. import __version__
from ..cancellation import check_cancelled
from ..scoring import SCORE_VERSION
from .tabular import GENBANK_HEADERS, HEADERS, annotation_rows

PROVENANCE_HEADERS = ["Kannotate version", "Scoring version", "Pfam version", "CDD version"]
GENBANK_CSV_HEADERS = GENBANK_HEADERS[:6] + ["Length"] + GENBANK_HEADERS[6:]


class CsvWriter:
    def __init__(self, handle, *, genbank=False, database_versions=None):
        self.genbank = genbank
        versions = database_versions or {}
        self.provenance = [__version__, SCORE_VERSION, versions.get("pfam"), versions.get("cdd")]
        self.writer = csv.writer(handle, lineterminator="\r\n")
        self.writer.writerow(GENBANK_CSV_HEADERS if genbank else HEADERS + PROVENANCE_HEADERS)

    def append(self, annotations, *, genbank=None, input_name=None):
        if self.genbank != (genbank is not None):
            raise ValueError("CSV input kind differs from its header")
        for values, _ in annotation_rows(annotations, genbank=genbank, input_name=input_name):
            check_cancelled()
            if self.genbank:
                # Count residues in the normalized protein, excluding stop symbols.
                values.insert(6, len(values[7].replace("*", "")))
            else:
                values += self.provenance
            self.writer.writerow(values)


def write_csv(annotations, path: Path, *, database_versions=None, genbank=None, input_name=None):
    with path.open("w", encoding="utf-8-sig", newline="", buffering=1024 * 1024) as handle:
        writer = CsvWriter(handle, genbank=genbank is not None, database_versions=database_versions)
        writer.append(annotations, genbank=genbank, input_name=input_name)
        check_cancelled()
