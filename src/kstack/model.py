"""Project schema ownership and validation at every persistence/backend boundary."""
from __future__ import annotations

import copy
import math
import re
from datetime import datetime, timezone

from . import SCHEMA_VERSION, __version__
from .colors import generate_colors
from .homology import connected_components, search, shorter_protein_coverage
from .layout import automatic_layout


def timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def validate_settings(settings: dict) -> None:
    for key in ("identity", "coverage"):
        value = settings.get(key)
        if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"{key} must be between 0 and 1.")


def validate(project: dict) -> dict:
    if not isinstance(project, dict) or project.get("format") != "Kstack":
        raise ValueError("This is not a Kstack project.")
    version = project.get("schema_version")
    if type(version) is not int or version != SCHEMA_VERSION:
        raise ValueError(f"Unsupported Kstack schema {version!r}; this app reads schema {SCHEMA_VERSION}.")
    validate_settings(project["analysis"])
    if type(project["analysis"]["preserve_input"]) is not bool:
        raise ValueError("preserve_input must be boolean.")
    if not all(isinstance(project["metadata"][key], str) for key in ("title", "created", "modified")):
        raise ValueError("Invalid project metadata.")
    valid_id = lambda value: isinstance(value, str) and re.fullmatch(r"[A-Za-z0-9_-]+", value)
    clusters = project["clusters"]
    if not 1 <= len(clusters) <= 50:
        raise ValueError("A project must have 1–50 clusters.")
    cluster_ids, genes = set(), {}
    for cluster in clusters:
        cid = cluster["id"]
        if not valid_id(cid) or cid in cluster_ids or not isinstance(cluster["name"], str):
            raise ValueError("Invalid/duplicate cluster ID or name.")
        cluster_ids.add(cid)
        if type(cluster["length"]) is not int or not 0 < cluster["length"] <= 100_000:
            raise ValueError("Invalid cluster length.")
        if "sequence" in cluster or "nucleotide" in cluster:
            raise ValueError("Project cannot retain nucleotide sequences.")
        for gene in cluster["cds"]:
            if not valid_id(gene["id"]) or gene["id"] in genes or gene["cluster_id"] != cid:
                raise ValueError("Invalid/duplicate gene reference.")
            genes[gene["id"]] = gene
            if not isinstance(gene["protein"], str) or re.search(r"[^A-Z*]", gene["protein"]):
                raise ValueError("Invalid protein sequence.")
            if gene["protein_length"] != len(gene["protein"]):
                raise ValueError("Inconsistent protein length.")
            if gene["strand"] not in (-1, 0, 1):
                raise ValueError("Invalid strand.")
            for segment in gene["segments"]:
                if not 0 <= segment["start"] < segment["end"] <= cluster["length"]:
                    raise ValueError("Invalid CDS segment coordinates.")
            if gene["cds_length"] != sum(s["end"] - s["start"] for s in gene["segments"]):
                raise ValueError("Inconsistent CDS length.")
    layout = project["layout"]
    if set(layout["order"]) != cluster_ids or len(layout["order"]) != len(cluster_ids):
        raise ValueError("Cluster order must be a permutation of cluster IDs.")
    if set(layout["rows"]) != cluster_ids:
        raise ValueError("Layout rows do not match clusters.")
    for row in layout["rows"].values():
        if row["orientation"] not in (-1, 1) or not math.isfinite(row["offset"]):
            raise ValueError("Invalid orientation or offset.")
    group_ids, members = set(), set()
    for group in project["homolog_groups"]:
        if not valid_id(group["id"]) or group["id"] in group_ids or not group["members"]:
            raise ValueError("Invalid homolog group.")
        group_ids.add(group["id"])
        for gid in group["members"]:
            if gid not in genes or gid in members or genes[gid]["group_id"] != group["id"]:
                raise ValueError("Inconsistent homolog membership.")
            members.add(gid)
        if not re.fullmatch(r"#[0-9a-fA-F]{6}", group["generated_color"]):
            raise ValueError("Invalid generated color.")
    if members != set(genes):
        raise ValueError("Every CDS must belong to one homolog component.")
    seen_hits = set()
    for hit in project["pairwise_hits"]:
        q, t = hit["query"], hit["target"]
        if q not in genes or t not in genes or q >= t or (q, t) in seen_hits:
            raise ValueError("Invalid direct protein relationship.")
        seen_hits.add((q, t))
        if not 0 <= hit["identity"] <= 1 or type(hit["alignment_length"]) is not int or hit["alignment_length"] < 1:
            raise ValueError("Invalid alignment values.")
        if not math.isfinite(hit["bits"]) or not math.isfinite(hit["evalue"]) or hit["evalue"] < 0:
            raise ValueError("Invalid alignment score/e-value.")
        if hit["query_length"] != len(genes[q]["protein"]) or hit["target_length"] != len(genes[t]["protein"]):
            raise ValueError("Alignment protein lengths do not match.")
        coverage = hit.get("shorter_coverage")
        if type(coverage) not in (int, float) or not math.isfinite(coverage) or not 0 <= coverage <= 1:
            raise ValueError("Shorter-protein coverage must be between 0 and 1.")
        expected = shorter_protein_coverage(hit["query_length"], hit["target_length"],
                                            hit.get("query_coverage"), hit.get("target_coverage"))
        if not math.isclose(hit["shorter_coverage"], expected, rel_tol=1e-6):
            raise ValueError("Inconsistent shorter-protein coverage.")
    visual = project["visualization"]
    if not 0 <= visual["connector_identity"] <= 1:
        raise ValueError("Invalid connector display identity.")
    for key in ("connectors", "labels"):
        if type(visual[key]) is not bool:
            raise ValueError(f"{key} must be boolean.")
    for key in ("width_scale", "height_scale", "row_spacing_scale"):
        value = visual.get(key, 1)
        if type(value) not in (int, float) or not math.isfinite(value) or not .25 <= value <= 4:
            raise ValueError(f"{key} must be between 0.25 and 4.")
    for namespace, valid_ids in (("genes", genes), ("groups", group_ids)):
        for key, color in project["color_overrides"][namespace].items():
            if key not in valid_ids or not re.fullmatch(r"#[0-9a-fA-F]{6}", color):
                raise ValueError("Invalid color override.")
    if len(set(project["fasta_selection"])) != len(project["fasta_selection"]) or not set(project["fasta_selection"]) <= genes.keys():
        raise ValueError("Invalid FASTA selection.")
    return project


def create_project(clusters: list[dict], settings: dict | None = None,
                   preserve: bool = False, connectors: bool = True) -> dict:
    settings = settings or {"identity": .30, "coverage": .70}
    validate_settings(settings)
    if not 1 <= len(clusters) <= 50:
        raise ValueError("Generate requires 1–50 clusters.")
    clusters = copy.deepcopy(clusters)
    if len({c["id"] for c in clusters}) != len(clusters):
        raise ValueError("Input cluster IDs must be unique.")
    hits, provenance = search(clusters)
    groups = connected_components(clusters, hits, settings)
    generate_colors(groups)
    now = timestamp()
    return validate({
        "format": "Kstack", "schema_version": SCHEMA_VERSION, "app_version": __version__,
        "metadata": {"created": now, "modified": now, "title": "Kstack", "analysis_provenance": provenance},
        "analysis": {**settings, "preserve_input": preserve},
        "visualization": {"connectors": connectors, "connector_identity": settings["identity"], "labels": False},
        "clusters": clusters, "homolog_groups": groups, "pairwise_hits": hits,
        "layout": automatic_layout(clusters, preserve), "color_overrides": {"genes": {}, "groups": {}},
        "fasta_selection": [],
    })


def recompute(project: dict, settings: dict, relayout: bool = False) -> dict:
    validate(project)
    validate_settings(settings)
    result = copy.deepcopy(project)
    result["analysis"].update(settings)
    hits, provenance = search(result["clusters"])
    groups = connected_components(result["clusters"], hits, result["analysis"])
    generate_colors(groups)
    # Group overrides survive only when the component's exact membership survives.
    old = {tuple(g["members"]): g["id"] for g in project["homolog_groups"]}
    overrides = {}
    for group in groups:
        previous = old.get(tuple(group["members"]))
        if previous in project["color_overrides"]["groups"]:
            overrides[group["id"]] = project["color_overrides"]["groups"][previous]
    result["color_overrides"]["groups"] = overrides
    result["pairwise_hits"], result["homolog_groups"] = hits, groups
    result["metadata"]["analysis_provenance"] = provenance
    result["metadata"]["modified"] = timestamp()
    if relayout:
        result["layout"] = automatic_layout(result["clusters"], result["analysis"]["preserve_input"])
    return validate(result)
