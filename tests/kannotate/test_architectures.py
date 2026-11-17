"""Compare pruned architecture selection with an exhaustive, independent oracle."""

import random
from dataclasses import replace
from itertools import combinations, pairwise, product

import pytest

from kannotate import Hit
from kannotate import architectures as module
from kannotate.cancellation import RunCancelled, RunControl
from kannotate.scoring import coverage


def exhaustive(hits, length, best):
    rules = module.RULES
    candidates = {}
    for group, members in rules["groups"].items():
        candidates[group] = [
            span
            for hit in hits
            if hit.score >= rules["min_domain_score"]
            and {hit.profile_id.split(".")[0], hit.name}.intersection(members)
            for span in hit.spans or ((hit.start, hit.end),)
            if span[1] - span[0] + 1
            >= rules.get("group_min_domain_span", {}).get(group, rules["min_domain_span"])
        ]
    for rule in rules["rules"]:
        if length < rule.get("min_protein_length", rules["min_protein_length"]):
            continue
        for spans in product(*(candidates[group] for group in rule["requires"])):
            positions = dict(zip(rule["requires"], spans, strict=True))
            ordered = all(
                positions[a][0] < positions[b][0] for a, b in pairwise(rule.get("ordered", []))
            )
            overlaps = [
                max(0, min(b, d) - max(a, c) + 1) / min(b - a + 1, d - c + 1)
                for (a, b), (c, d) in combinations(spans, 2)
            ]
            if ordered and all(value <= rules["max_overlap_fraction"] for value in overlaps):
                if coverage(best) >= 0.8 and best.score >= 40:
                    if not any(
                        word in best.description.lower() for word in ("synthetase", "synthase")
                    ):
                        return None
                return rule["product"]
    return None


def domain(group, spans, score=20):
    profile = module.RULES["groups"][group][0]
    return Hit(
        "hmm",
        "q",
        profile + ".1",
        profile,
        "domain",
        1e-20,
        spans[0][0],
        spans[-1][1],
        3000,
        score=score,
        spans=tuple(spans),
    )


def test_pruned_selection_matches_exhaustive_rules_and_boundaries():
    rng = random.Random(28317)
    best = Hit("hmm", "q", "other", "other", "family", 1e-20, 1, 100, 3000, score=30)
    for rule in module.RULES["rules"]:
        hits = [
            domain(group, [(1 + i * 300, 100 + i * 300)])
            for i, group in enumerate(rule["requires"])
        ]
        assert module.classify_architecture(hits, 3000, best) == exhaustive(hits, 3000, best)
    for _ in range(400):
        hits = []
        for group in module.RULES["groups"]:
            spans = []
            for _ in range(rng.randrange(5)):
                start = rng.randrange(1, 31) * 80
                spans.append((start, start + rng.choice((49, 50, 79, 80, 100, 200)) - 1))
            if spans:
                hit = domain(group, spans, score=rng.choice((14.9, 15, 20)))
                hits.extend([hit] * rng.choice((1, 2)))
        length = rng.choice((699, 700, 1499, 1500, 3000))
        strongest = replace(best, end=rng.choice((100, 2399, 2400)), score=rng.choice((39.9, 40)))
        assert module.classify_architecture(hits, length, strongest) == exhaustive(
            hits, length, strongest
        )
        assert module.classify_architecture(list(reversed(hits)), length, strongest) == exhaustive(
            hits, length, strongest
        )


def test_architecture_stop_interrupts_combination_search(monkeypatch):
    hits = [domain(group, [(i, i + 100) for i in range(1, 33)]) for group in module.RULES["groups"]]
    control = RunControl()
    original = module.check_cancelled
    checks = 0

    def stop_inside_combinations():
        nonlocal checks
        checks += 1
        if checks == len(hits) + 5:
            control.request_stop()
        original()

    monkeypatch.setattr(module, "check_cancelled", stop_inside_combinations)
    with control.activate(), pytest.raises(RunCancelled):
        module.classify_architecture(hits, 3000, hits[0])
    assert checks == len(hits) + 5
