"""Strict single-chain amino-acid input; no silent biological substitutions."""

from dataclasses import dataclass

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class ProteinSequence:
    identifier: str
    sequence: str

    @property
    def length(self):
        return len(self.sequence)


def parse_sequence(text: str) -> ProteinSequence:
    if not isinstance(text, str) or not text.strip():
        raise ValueError("Paste one protein sequence or a single FASTA record.")
    lines = text.strip().splitlines()
    identifier = "query"
    if lines[0].startswith(">"):
        header = lines.pop(0)[1:].strip()
        identifier = header.split(maxsplit=1)[0] if header else "query"
    if any(">" in line for line in lines):
        raise ValueError("Kmatch accepts exactly one FASTA record.")
    sequence = "".join("".join(lines).split()).upper()
    if not sequence:
        raise ValueError("The protein sequence is empty.")
    invalid = sorted(set(sequence) - AMINO_ACIDS)
    if invalid:
        raise ValueError(
            "Unsupported amino-acid characters: "
            + ", ".join(invalid)
            + ". Use the 20 standard amino acids; remove gaps or stop symbols explicitly."
        )
    if len(sequence) < 10:
        raise ValueError("Provide at least 10 residues for single-protein prediction.")
    return ProteinSequence(identifier, sequence)
