"""Conservative presentation advice, separate from search acceptance/ranking."""

import math


def boltz_warnings(identity, confidence, provenance):
    if identity != "boltz2" or provenance.get("msa") != "single-sequence":
        return []
    thresholds = {"mean_plddt": 70.0, "ptm": 0.5}
    low = {}
    for key, threshold in thresholds.items():
        value = confidence.get(key)
        if isinstance(value, (int, float)) and math.isfinite(value) and value < threshold:
            low[key] = value
    if not low:
        return []
    return [
        {
            "code": "boltz_msa_free_low_confidence",
            "message": (
                "Low-confidence MSA-free prediction: structural-search results may be unreliable. "
                "Rerun with an MSA (a local A3M, or explicitly enable network-assisted MSA). "
                "This is a conservative warning, not a biological guarantee."
            ),
            "observed_below_threshold": low,
            "thresholds": thresholds,
            "rule": "MSA-free Boltz-2: mean pLDDT <70 or pTM <0.5",
        }
    ]
