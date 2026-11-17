from .architectures import classify_architecture
from .descriptions import normalize_description
from .models import Annotation, Protein
from .products import profile_product
from .scoring import DEFAULT_MIN_SCORE, rank_hits


def select_annotation(protein: Protein, hmm_hits, rps_hits, min_score=DEFAULT_MIN_SCORE):
    hmm = rank_hits(hmm_hits)
    rps = rank_hits(rps_hits)
    all_hits = sorted(hmm + rps, key=lambda h: (-h.score, h.evalue, h.profile_id, h.engine))
    best = all_hits[0] if all_hits else None
    description = "hypothetical protein"
    mode = "single"
    if best and best.score >= min_score:
        architecture = classify_architecture(all_hits, len(protein.sequence), best)
        mapped = profile_product(best, len(protein.sequence))
        description = architecture or mapped or normalize_description(best.description)
        mode = "architecture" if architecture else "lookup" if mapped else "single"
    return Annotation(
        protein, description, best.score if best else 0.0, tuple(hmm), tuple(rps), mode
    )
