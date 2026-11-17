"""Deterministic farthest-point palette in OKLab, gamut-mapped to sRGB."""
from __future__ import annotations

import math


def rgb(l, a, b):
    ll = (l + .3963377774 * a + .2158037573 * b) ** 3
    mm = (l - .1055613458 * a - .0638541728 * b) ** 3
    ss = (l - .0894841775 * a - 1.291485548 * b) ** 3
    linear = (4.0767416621 * ll - 3.3077115913 * mm + .2309699292 * ss,
              -1.2684380046 * ll + 2.6097574011 * mm - .3413193965 * ss,
              -.0041960863 * ll - .7034186147 * mm + 1.707614701 * ss)
    if any(v < 0 or v > 1 for v in linear):
        return None
    srgb = [12.92 * v if v <= .0031308 else 1.055 * v ** (1 / 2.4) - .055 for v in linear]
    return "#" + "".join(f"{round(v * 255):02x}" for v in srgb)


def generate_colors(groups: list[dict]) -> None:
    count = sum(len(g["members"]) > 1 for g in groups)
    # Increase candidate density with figure size, never repeat a fixed categorical set.
    hue_step = max(1, min(5, 720 // max(1, count)))
    levels = max(4, math.ceil(count / (360 / hue_step)))
    candidates = []
    lightnesses = [.65, .76, .54, .86] + [.50 + .38 * i / (levels - 3) for i in range(1, levels - 3)]
    for lightness in lightnesses:
        for hue in range(0, 360, hue_step):
            angle = math.radians(hue)
            chroma = .20
            color = None
            while color is None and chroma > .035:
                lab = (lightness, chroma * math.cos(angle), chroma * math.sin(angle))
                color = rgb(*lab)
                chroma -= .008
            if color:
                candidates.append((lab, color))
    chosen = []
    distances = [float("inf")] * len(candidates)
    for group in sorted(groups, key=lambda g: g["id"]):
        if len(group["members"]) == 1:
            group["generated_color"] = "#a3aab4"
            continue
        if not chosen:
            index = 8  # warm coral anchor, followed by maximum perceptual distances
        else:
            index = max(range(len(candidates)), key=lambda i: (distances[i], -i))
        lab, color = candidates[index]
        group["generated_color"] = color
        chosen.append(lab)
        for i, (other, _) in enumerate(candidates):
            distance = sum((a - b) ** 2 for a, b in zip(lab, other))
            distances[i] = min(distances[i], distance)
