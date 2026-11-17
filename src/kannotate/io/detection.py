from pathlib import Path


def detect_format(path: Path) -> str:
    with path.open(encoding="utf-8-sig", errors="replace") as handle:
        prefix = handle.read(8192).lstrip()
    if prefix.startswith("LOCUS"):
        return "genbank"
    if prefix.startswith(">") or prefix.startswith(";") and ">" in prefix:
        return "fasta"
    if path.suffix.lower() in {".gb", ".gbk", ".gbff", ".genbank"}:
        return "genbank"
    if path.suffix.lower() in {".fa", ".faa", ".fasta", ".fas", ".fsa"}:
        return "fasta"
    raise ValueError("Input is neither protein FASTA nor GenBank")
