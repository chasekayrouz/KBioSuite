import json
from importlib.resources import files
from itertools import pairwise

from .cancellation import check_cancelled
from .models import Hit
from .scoring import coverage

RULES = json.loads(files("kannotate.resources").joinpath("architectures.json").read_text())


def classify_architecture(hits: list[Hit], query_length: int, best: Hit):
    if query_length < min(
        r.get("min_protein_length", RULES["min_protein_length"]) for r in RULES["rules"]
    ):
        return None
    check_cancelled()
    groups = {}
    for hit in hits:
        check_cancelled()
        if hit.score < RULES["min_domain_score"]:
            continue
        ids = {hit.profile_id.split(".")[0], hit.name}
        for group, members in RULES["groups"].items():
            if ids.intersection(members):
                for a, b in hit.spans or ((hit.start, hit.end),):
                    if b - a + 1 >= RULES.get("group_min_domain_span", {}).get(
                        group, RULES["min_domain_span"]
                    ):
                        groups.setdefault(group, set()).add((a, b))
    for rule in RULES["rules"]:
        if query_length < rule.get("min_protein_length", RULES["min_protein_length"]):
            continue
        required = rule["requires"]
        if not all(g in groups for g in required):
            continue
        order = tuple(pairwise(rule.get("ordered", [])))
        positions = {}

        def compatible(index, required=required, positions=positions, order=order):
            check_cancelled()
            if index == len(required):
                return True
            group = required[index]
            for a, b in sorted(groups[group]):
                check_cancelled()
                if any(
                    max(0, min(b, d) - max(a, c) + 1) / min(b - a + 1, d - c + 1)
                    > RULES["max_overlap_fraction"]
                    for c, d in positions.values()
                ):
                    continue
                positions[group] = (a, b)
                if not any(
                    x in positions and y in positions and positions[x][0] >= positions[y][0]
                    for x, y in order
                ) and compatible(index + 1):
                    return True
                del positions[group]
            return False

        if compatible(0):
            # A strongly supported full-length family is more informative.
            if coverage(best) >= 0.8 and best.score >= 40:
                desc = best.description.lower()
                if "synthetase" not in desc and "synthase" not in desc:
                    return None
            return rule["product"]
    return None
