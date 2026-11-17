"""Small profile-supported product mappings; never accession/sequence-specific."""

import json
from importlib.resources import files

from .scoring import coverage

RULES = json.loads(files("kannotate.resources").joinpath("products.json").read_text())


def profile_product(hit, query_length):
    ids = {hit.profile_id.split(".")[0], hit.name}
    for rule in RULES["rules"]:
        spans = hit.spans or ((hit.start, hit.end),)
        if (
            ids.intersection(rule["profiles"])
            and query_length <= rule["max_query_length"]
            and hit.score >= rule["min_score"]
            and coverage(hit) >= rule["min_coverage"]
            and max(b - a + 1 for a, b in spans) >= rule["min_alignment_span"]
        ):
            return rule["product"]
    return None
