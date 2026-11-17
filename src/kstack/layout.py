"""Deterministic overlap ordering, two-state orientation DP, median alignment."""
from __future__ import annotations

import statistics


def centers(cluster: dict, orientation: int) -> dict[str, float]:
    result = {}
    for gene in cluster["cds"]:
        if not gene["segments"]:
            continue
        length = gene["cds_length"]
        center = sum((s["start"] + s["end"]) / 2 * (s["end"] - s["start"])
                     for s in gene["segments"]) / max(1, length)
        if orientation == -1:
            center = cluster["length"] - center
        result.setdefault(gene["group_id"], []).append(center)
    return {group: statistics.median(values) for group, values in result.items()}


def automatic_layout(clusters: list[dict], preserve: bool = False) -> dict:
    if preserve or len(clusters) < 2:
        return {"order": [c["id"] for c in clusters],
                "rows": {c["id"]: {"orientation": 1, "offset": 0.0} for c in clusters}}
    by_id = {c["id"]: c for c in clusters}
    ids = sorted(by_id)
    group_sets = {c["id"]: {g["group_id"] for g in c["cds"]} for c in clusters}
    similarity = {}
    for a in ids:
        for b in ids:
            union = group_sets[a] | group_sets[b]
            similarity[a, b] = len(group_sets[a] & group_sets[b]) / max(1, len(union))

    # Best nearest-neighbour path from every start, then deterministic best-improving 2-opt.
    paths = []
    for start in ids:
        path, remaining = [start], set(ids) - {start}
        while remaining:
            nxt = min(remaining, key=lambda b: (-similarity[path[-1], b], b))
            path.append(nxt)
            remaining.remove(nxt)
        paths.append(path)
    score = lambda path: sum(similarity[a, b] for a, b in zip(path, path[1:]))
    order = min(paths, key=lambda path: (-score(path), tuple(path)))
    while True:
        best, best_score = order, score(order)
        for i in range(len(order) - 1):
            for j in range(i + 1, len(order)):
                candidate = order[:i] + list(reversed(order[i:j + 1])) + order[j + 1:]
                candidate_score = score(candidate)
                if candidate_score > best_score + 1e-12:
                    best, best_score = candidate, candidate_score
        if best is order:
            break
        order = best
    order = min(order, list(reversed(order)))
    positions = {(cid, orientation): centers(by_id[cid], orientation)
                 for cid in ids for orientation in (1, -1)}

    def correspondence(a, oa, b, ob):
        ca, cb = positions[a, oa], positions[b, ob]
        shared = sorted(ca.keys() & cb.keys())
        if len(shared) < 2:
            return 0.0
        delta = statistics.median(ca[g] - cb[g] for g in shared)
        return sum(abs(ca[g] - cb[g] - delta) for g in shared) / len(shared)

    # Fix first row forward. DP minimizes median-centered adjacent correspondence error.
    dp = {1: (0.0, (1,))}
    for a, b in zip(order, order[1:]):
        next_dp = {}
        for ob in (1, -1):
            options = [(cost + correspondence(a, oa, b, ob), path + (ob,))
                       for oa, (cost, path) in dp.items()]
            next_dp[ob] = min(options, key=lambda item: (item[0], tuple(-x for x in item[1])))
        dp = next_dp
    orientations = min(dp.values(), key=lambda item: (item[0], tuple(-x for x in item[1])))[1]
    rows = {order[0]: {"orientation": orientations[0], "offset": 0.0}}
    for index, (a, b) in enumerate(zip(order, order[1:])):
        ca, cb = positions[a, orientations[index]], positions[b, orientations[index + 1]]
        shared = ca.keys() & cb.keys()
        delta = statistics.median(ca[g] - cb[g] for g in shared) if shared else 0.0
        rows[b] = {"orientation": orientations[index + 1], "offset": round(rows[a]["offset"] + delta, 6)}
    return {"order": order, "rows": rows}
