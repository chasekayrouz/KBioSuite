from dataclasses import dataclass

# The common unaligned protein alphabet, verified against the pinned local tools.
# The stop symbol is accepted without dropping coordinates by both search tools.
# Alignment gaps/punctuation are incompatible with common unaligned search input.
PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWYBJZOUX*")
CANONICAL_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWY")
MIN_SEARCH_LENGTH = 6


@dataclass(frozen=True)
class NormalizedSequence:
    sequence: str
    truncated: bool = False

    @property
    def searchable(self) -> bool:
        return is_searchable(self.sequence)


def normalize_sequence(raw: str) -> NormalizedSequence:
    sequence = []
    for char in raw:
        if char.isspace():
            continue
        char = char.upper()
        if char not in PROTEIN_ALPHABET:
            return NormalizedSequence("".join(sequence), True)
        sequence.append(char)
    return NormalizedSequence("".join(sequence))


def is_searchable(sequence: str) -> bool:
    return len(sequence) >= MIN_SEARCH_LENGTH and any(
        c in PROTEIN_ALPHABET - {"X", "*"} for c in sequence
    )
