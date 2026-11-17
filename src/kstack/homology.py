"""Sensitive pairwise search, explicit threshold graph, separate direct relationships."""
from __future__ import annotations

import os
import math
import re
import subprocess
import tempfile
import time
from pathlib import Path

from .dependencies import resolve_mmseqs


def passing(hit: dict, settings: dict) -> bool:
    return (hit["identity"] + 1e-9 >= settings["identity"]
            and 0 <= hit["shorter_coverage"] <= 1
            and hit["shorter_coverage"] >= settings["coverage"])


def shorter_protein_coverage(qlen: int, tlen: int, qcov: float, tcov: float) -> float:
    """Sequence coverage, never the number of gapped alignment columns."""
    if qlen <= 0 or tlen <= 0:
        raise ValueError("Alignment protein lengths must be positive.")
    for value in (qcov, tcov):
        if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("Query/target sequence coverage must be between 0 and 1.")
    return qcov if qlen < tlen else tcov if tlen < qlen else min(qcov, tcov)


def parse_alignment(line: str) -> dict | None:
    q, t, identity, length, qlen, tlen, qcov, tcov, score, evalue = line.split("\t")
    if q == t:
        return None
    qlen, tlen, qcov, tcov = int(qlen), int(tlen), float(qcov), float(tcov)
    if q > t:
        q, t, qlen, tlen, qcov, tcov = t, q, tlen, qlen, tcov, qcov
    return {"query": q, "target": t, "identity": float(identity), "alignment_length": int(length),
            "query_length": qlen, "target_length": tlen,
            "query_coverage": qcov, "target_coverage": tcov,
            "shorter_coverage": shorter_protein_coverage(qlen, tlen, qcov, tcov),
            "bits": float(score), "evalue": float(evalue)}


def connected_components(clusters: list[dict], hits: list[dict], settings: dict) -> list[dict]:
    ids = sorted(g["id"] for c in clusters for g in c["cds"])
    parent = {gid: gid for gid in ids}

    def find(gid):
        while parent[gid] != gid:
            parent[gid] = parent[parent[gid]]
            gid = parent[gid]
        return gid

    for hit in hits:
        if passing(hit, settings):
            a, b = find(hit["query"]), find(hit["target"])
            if a != b:
                parent[max(a, b)] = min(a, b)
    components = {}
    for gid in ids:
        components.setdefault(find(gid), []).append(gid)
    groups = [{"id": f"h_{root}", "members": members, "generated_color": "#a3aab4"}
              for root, members in sorted(components.items())]
    assignments = {gid: group["id"] for group in groups for gid in group["members"]}
    for cluster in clusters:
        for gene in cluster["cds"]:
            gene["group_id"] = assignments[gene["id"]]
    return groups


def search(clusters: list[dict], threads: int | None = None) -> tuple[list[dict], dict]:
    proteins = {g["id"]: g["protein"] for c in clusters for g in c["cds"] if g["protein"]}
    mmseqs = resolve_mmseqs()
    version = subprocess.check_output([mmseqs, "version"], text=True).strip()
    cores = max(1, min(threads or os.cpu_count() or 1, os.cpu_count() or 1))
    short = {gid: seq for gid, seq in proteins.items() if len(seq) < 20}
    prefilter = 2 if short and len(short) == len(proteins) else 0
    settings = {"engine": "MMseqs2", "version": version, "sensitivity": 7.5,
                "evalue": 10.0, "alignment_mode": 3, "min_seq_id": 0.0, "coverage_prefilter": 0.0,
                "max_seqs": max(1, len(proteins)), "max_accept": 2147483647,
                "max_rejected": 2147483647, "threads": cores, "iterations": 1,
                "prefilter_mode": prefilter, "short_protein_supplement": bool(short and prefilter == 0),
                "mask": 0, "comp_bias_corr": 1}
    start = time.perf_counter()
    if len(proteins) < 2:
        return [], {**settings, "elapsed_seconds": time.perf_counter() - start, "protein_count": len(proteins)}
    with tempfile.TemporaryDirectory(prefix="kstack-search-") as temp:
        temp = Path(temp)
        fasta, tsv = temp / "proteins.faa", temp / "alignments.tsv"
        fasta.write_text("".join(f">{gid}\n{re.sub('[^ACDEFGHIKLMNPQRSTVWYBXZJUO]', 'X', seq)}\n"
                                  for gid, seq in sorted(proteins.items())))
        command = [mmseqs, "easy-search", str(fasta), str(fasta), str(tsv), str(temp / "work"),
                   "--search-type", "1", "-s", "7.5", "-e", "10", "--min-seq-id", "0", "-c", "0",
                   "--alignment-mode", "3", "--seq-id-mode", "0", "--num-iterations", "1", "--max-seqs", str(len(proteins)),
                   "--max-accept", "2147483647", "--max-rejected", "2147483647",
                   "--threads", str(cores), "--mask", "0", "--comp-bias-corr", "1",
                   "--prefilter-mode", str(prefilter),
                   "--format-output", "query,target,fident,alnlen,qlen,tlen,qcov,tcov,bits,evalue",
                   "--remove-tmp-files", "1", "-v", "1"]
        result = subprocess.run(command, capture_output=True, text=True, timeout=110)
        if result.returncode:
            raise RuntimeError(f"MMseqs2 failed ({result.returncode}): {(result.stderr or result.stdout)[-4000:]}")
        if short and prefilter == 0:
            # Seedless supplementary search prevents partial short CDS proteins from
            # disappearing solely because they cannot supply a spaced k-mer seed.
            short_fasta, short_tsv = temp / "short.faa", temp / "short.tsv"
            short_fasta.write_text("".join(f">{gid}\n{re.sub('[^ACDEFGHIKLMNPQRSTVWYBXZJUO]', 'X', seq)}\n"
                                         for gid, seq in sorted(short.items())))
            supplement = list(command)
            supplement[2], supplement[4], supplement[5] = str(short_fasta), str(short_tsv), str(temp / "short-work")
            supplement[supplement.index("--prefilter-mode") + 1] = "2"
            result = subprocess.run(supplement, capture_output=True, text=True, timeout=110)
            if result.returncode:
                raise RuntimeError(f"Short-protein MMseqs2 search failed: {(result.stderr or result.stdout)[-4000:]}")
            tsv.write_text(tsv.read_text() + short_tsv.read_text())
        pairs = {}
        for line in tsv.read_text().splitlines():
            hit = parse_alignment(line)
            if hit is None:
                continue
            # Reciprocal local alignments can differ; keep one deterministic best-scoring alignment.
            key = (hit["query"], hit["target"])
            rank = lambda h: (h["bits"], -h["evalue"], h["identity"], h["alignment_length"])
            if key not in pairs or rank(hit) > rank(pairs[key]):
                pairs[key] = hit
    return [pairs[key] for key in sorted(pairs)], {
        **settings, "elapsed_seconds": time.perf_counter() - start, "protein_count": len(proteins)}
