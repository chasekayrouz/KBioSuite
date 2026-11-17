"""One protein per run, with conservative input normalization."""

import re

from .models import Query

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def parse_query(text: str) -> Query:
    if not isinstance(text, str):
        raise ValueError("Protein input must be sequence text or single-record FASTA text.")
    text = text.lstrip("\ufeff").strip()
    if not text:
        raise ValueError("Paste a protein sequence or import a single-record FASTA.")
    identifier, description = "query", "Protein query"
    if text.startswith(">"):
        lines = text.splitlines()
        if any(line.lstrip().startswith(">") for line in lines[1:]):
            raise ValueError("Kblast supports one FASTA record per run.")
        header = lines[0][1:].strip()
        if not header:
            raise ValueError("The FASTA header needs an identifier.")
        identifier = header.split()[0]
        description = header
        text = "\n".join(lines[1:])
    elif ">" in text:
        raise ValueError("FASTA headers must be the first line of a single-record input.")
    sequence = re.sub(r"\s+", "", text).upper()
    warnings = []
    if sequence.endswith("*"):
        sequence = sequence[:-1]
        warnings.append("Removed one terminal stop symbol (*).")
    invalid = sorted(set(sequence) - AMINO_ACIDS)
    if invalid:
        raise ValueError("Invalid protein characters: " + ", ".join(repr(c) for c in invalid))
    if not sequence:
        raise ValueError("The protein sequence is empty.")
    if len(sequence) > 100_000:
        raise ValueError("Kblast accepts proteins up to 100,000 amino acids.")
    # ACGT-only sequences can also be valid proteins, so do not silently reject them.
    if set(sequence) <= set("ACGTUN"):
        warnings.append("This looks nucleotide-like; confirm it is an amino-acid sequence.")
    if len(sequence) < 30:
        warnings.append(
            "Short protein: NCBI may automatically adjust BLAST parameters, as on its web form."
        )
    return Query(identifier, description, sequence, warnings)
