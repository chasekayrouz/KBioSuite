import re


def normalize_description(description: str, fallback: str = "hypothetical protein") -> str:
    text = " ".join(description.split())
    text = re.sub(r"^(?:gnl\|CDD\|\d+\s+)?(?:[A-Za-z]+\d+(?:\.\d+)?,\s*)", "", text)
    # A leading biological word can be meaningful: "glycosyltransferase,
    # activator-dependent family" must retain its enzyme name. Remove a short
    # prefix only when the following description repeats that same word.
    text = re.sub(r"^(\S+),\s+(?=\1(?:\s|[,.;]))", "", text)
    text = re.sub(r"^\[?conserved domain\]?\s*[:\-]\s*", "", text, flags=re.IGNORECASE)
    # Long explanatory paragraphs are not suitable for a GenBank product name.
    text = re.split(r"\.\s+(?=[A-Z0-9])", text, maxsplit=1)[0]
    text = text.rstrip(".;, ")
    text = re.sub(r"\s*;\s*(?:Provisional|Reviewed|Validated)\s*$", "", text, flags=re.IGNORECASE)
    boilerplate = (
        "NCBI|Pfam|CDD|General function prediction only|Function unknown|"
        "Secondary metabolites biosynthesis|Coenzyme transport|Amino acid transport|"
        "Carbohydrate transport|Nucleotide transport|Lipid transport|"
        "Inorganic ion transport|Energy production|Translation|Transcription|"
        "Signal transduction|Defense mechanisms|Intracellular trafficking|"
        "Replication|Cell wall|Cell cycle|Posttranslational modification|Motility"
    )
    text = re.sub(rf"\s*\[(?:{boilerplate})[^\]]*\]\s*$", "", text)
    text = re.split(r"\s+(?:catalyzes|catalyses)\s+", text, maxsplit=1)[0]
    text = re.sub(r"[;,.\s]+$", "", text).strip()
    return text or fallback
