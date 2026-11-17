from dataclasses import replace

from .models import Hit


def merge_spans(spans):
    merged = []
    for start, end in sorted(spans):
        if start > end:
            start, end = end, start
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return tuple(merged)


def collapse_profiles(hits: list[Hit]) -> list[Hit]:
    groups = {}
    for hit in hits:
        groups.setdefault((hit.engine, hit.query_id, hit.profile_id), []).append(hit)
    collapsed = []
    for key in sorted(groups):
        rows = groups[key]
        best = min(rows, key=lambda h: (h.evalue, -h.bitscore, h.start, h.end))
        spans = merge_spans(s for h in rows for s in (h.spans or ((h.start, h.end),)))
        alignments = tuple(
            sorted(
                {a for h in rows for a in h.alignments},
                key=lambda a: (a.query_start, a.query_end, a.evalue, -a.bitscore),
            )
        )
        collapsed.append(
            replace(best, start=spans[0][0], end=spans[-1][1], spans=spans, alignments=alignments)
        )
    return collapsed
