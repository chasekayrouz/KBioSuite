import json
import math
from dataclasses import replace
from importlib.resources import files

from .hits import merge_spans
from .models import Hit

PARAMETERS = json.loads(files("kannotate.resources").joinpath("scoring.json").read_text())
SCORE_VERSION = PARAMETERS["version"]
DEFAULT_MIN_SCORE = PARAMETERS["default_min_score"]


def coverage(hit: Hit) -> float:
    if hit.query_length <= 0:
        return 0.0
    spans = merge_spans(hit.spans or ((hit.start, hit.end),))
    length = sum(max(0, min(b, hit.query_length) - max(1, a) + 1) for a, b in spans)
    return min(1.0, length / hit.query_length)


def _score_hit(hit: Hit, cov: float) -> Hit:
    """v1: capped -log10(E) × (0.5 + 0.5√coverage) + 5×coverage + engine offset."""
    if not math.isfinite(hit.evalue) or hit.evalue < 0:
        raise ValueError("E-value must be finite and nonnegative")
    significance = min(PARAMETERS["significance_cap"], -math.log10(max(hit.evalue, 1e-300)))
    score = max(
        0.0,
        significance
        * (PARAMETERS["coverage_base"] + PARAMETERS["coverage_weight"] * math.sqrt(cov))
        + PARAMETERS["coverage_bonus"] * cov
        + PARAMETERS["engine_offsets"][hit.engine],
    )
    return replace(hit, score=round(score, 6))


def score_hit(hit: Hit) -> Hit:
    return _score_hit(hit, coverage(hit))


def rank_hits(hits: list[Hit]) -> list[Hit]:
    scored = [(_score_hit(hit, cov), cov) for hit in hits for cov in (coverage(hit),)]
    scored.sort(
        key=lambda pair: (
            -pair[0].score,
            pair[0].evalue,
            -pair[1],
            pair[0].profile_id,
            pair[0].engine,
        )
    )
    return [hit for hit, _ in scored]
