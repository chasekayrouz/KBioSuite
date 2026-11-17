from pathlib import Path

from ..cancellation import check_cancelled
from ..models import Protein
from ..normalization import PROTEIN_ALPHABET, normalize_sequence


def read_fasta(path: Path) -> list[Protein]:
    with path.open(encoding="utf-8-sig", errors="replace") as handle:
        return list(_iter_lines(handle))


def iter_fasta(path: Path):
    with path.open(encoding="utf-8-sig", errors="replace") as handle:
        yield from _iter_lines(handle)


def parse_fasta(text: str, *, strict=False, allow_raw=False) -> list[Protein]:
    """Parse file or pasted protein input without changing file-job normalization.

    Quick Search opts into strict validation: reject incompatible characters rather
    than truncating them, and require nonempty records. Duplicate headers stay distinct.
    """
    if not isinstance(text, str):
        raise ValueError("Paste protein FASTA or an amino-acid sequence.")
    text = text.lstrip("\ufeff")
    if allow_raw and not any(line.lstrip().startswith(">") for line in text.splitlines()):
        text = ">query_1\n" + text
    return list(_iter_lines(text.splitlines(), strict=strict))


def _iter_lines(lines, *, strict=False):
    count = 0
    header = None
    sequence = []

    def finish():
        if header is not None:
            raw = "".join(sequence)
            if strict:
                invalid = next(
                    (
                        c
                        for c in raw
                        if not c.isspace()
                        and c not in PROTEIN_ALPHABET
                        and c not in "acdefghiklmnpqrstvwybjzoux"
                    ),
                    None,
                )
                if invalid is not None:
                    raise ValueError(
                        f"{header or 'Unnamed record'}: invalid protein character {invalid!r}."
                    )
                if not header or not raw.strip():
                    raise ValueError("Each FASTA record needs a name and a protein sequence.")
            return Protein(str(count), header, normalize_sequence(raw).sequence)

    for line in lines:
        check_cancelled()
        if line.lstrip().startswith(">"):
            if header is not None:
                yield finish()
                count += 1
            header = line.lstrip()[1:].strip()
            sequence = []
        elif line.strip() and not line.lstrip().startswith(";"):
            if header is None:
                raise ValueError("Protein FASTA must begin with a > header")
            sequence.append(line)
    if header is not None:
        yield finish()
        count += 1
    if not count:
        raise ValueError("No FASTA records found")
