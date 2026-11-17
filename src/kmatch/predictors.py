"""Exchangeable external folding engines with normalized confidence/provenance."""

import json
import math
import shutil
import subprocess
import zipfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import gemmi

from .confidence import boltz_warnings
from .models import DependencyStatus, PredictionBackend, PredictionResult, Progress, write_json
from .process import environment, run_process

RUNNER = Path(__file__).with_name("engine_runner.py")


def package_version(python, package):
    if not python.is_file():
        return "not installed"
    try:
        return subprocess.check_output(
            [
                str(python),
                "-c",
                f"import importlib.metadata;print(importlib.metadata.version({package!r}))",
            ],
            env=environment(),
            text=True,
            stderr=subprocess.DEVNULL,
            timeout=5,
        ).strip()
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def validate_a3m(path, sequence):
    text = path.read_text()
    rows = []
    current = ""
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        if line.startswith(">"):
            if current:
                rows.append(current)
            current = ""
        else:
            current += line.strip()
    if current:
        rows.append(current)
    if not rows:
        raise ValueError("The A3M has no aligned sequences.")
    query = "".join(c for c in rows[0] if not c.islower()).replace("-", "")
    if query != sequence:
        raise ValueError("The first A3M sequence must match the normalized query exactly.")
    columns = len("".join(c for c in rows[0] if not c.islower()))
    for row in rows:
        if len("".join(c for c in row if not c.islower())) != columns:
            raise ValueError("The A3M contains inconsistent alignment lengths.")
        if set(row.upper()) - set("ACDEFGHIKLMNPQRSTVWYX-"):
            raise ValueError("The A3M contains unsupported alignment characters.")
    return text


def residue_confidence(structure, length):
    model = gemmi.read_structure(str(structure))[0]
    chains = [chain for chain in model if any(res.find_atom("CA", "*") for res in chain)]
    if len(chains) != 1:
        raise ValueError("Prediction must contain one protein chain.")
    values = [
        float(res.find_atom("CA", "*").b_iso) for res in chains[0] if res.find_atom("CA", "*")
    ]
    if len(values) != length or any(not math.isfinite(v) or not 0 <= v <= 100 for v in values):
        raise ValueError("Prediction has missing residues or invalid residue confidence.")
    return values


def normalize_boltz(raw: Path, output: Path, identity, version, length, runtime, provenance):
    candidates = []
    for path in raw.glob("boltz_results_*/predictions/*/*_model_*.pdb"):
        confidence = json.loads(path.with_name(f"confidence_{path.stem}.json").read_text())
        score = float(confidence["confidence_score"])
        if not math.isfinite(score):
            raise ValueError("Boltz returned non-finite confidence.")
        candidates.append((score, path, confidence))
    if not candidates:
        raise ValueError("Boltz completed without a usable prediction. See predictor.log.")
    _, source, confidence = max(candidates, key=lambda x: x[0])
    destination = output / "predicted.pdb"
    shutil.copy2(source, destination)
    per_residue = residue_confidence(destination, length)
    confidence["mean_plddt"] = sum(per_residue) / len(per_residue)
    result = PredictionResult(
        destination,
        identity,
        version,
        length,
        confidence,
        per_residue,
        runtime,
        {**provenance, "selected_model": source.name, "samples_produced": len(candidates)},
        warnings=boltz_warnings(identity, confidence, provenance),
    )
    write_json(output / "prediction.json", result)
    return result


def normalize_af3(raw: Path, output: Path, version, length, runtime, provenance):
    structures = list(raw.glob("*/*_model.cif"))
    if len(structures) != 1:
        raise ValueError("AF3 completed without its ranked model mmCIF.")
    source = structures[0]
    summary = json.loads(
        source.with_name(source.name.replace("_model.cif", "_summary_confidences.json")).read_text()
    )
    destination = output / "predicted.cif"
    shutil.copy2(source, destination)
    per_residue = residue_confidence(destination, length)
    summary["mean_plddt"] = sum(per_residue) / len(per_residue)
    result = PredictionResult(
        destination,
        "af3",
        version,
        length,
        summary,
        per_residue,
        runtime,
        {
            **provenance,
            "experimental_mps": True,
            "output_terms": "https://github.com/google-deepmind/alphafold3/blob/main/OUTPUT_TERMS_OF_USE.md",
        },
    )
    write_json(output / "prediction.json", result)
    return result


class BoltzBackend:
    def __init__(self, config, model="boltz2"):
        self.config, self.identity = config, model

    def status(self):
        cfg = self.config
        weights = (
            ["boltz2_conf.ckpt", "boltz2_aff.ckpt"]
            if self.identity == "boltz2"
            else ["boltz1_conf.ckpt", "ccd.pkl"]
        )
        ready = cfg.boltz_python.is_file() and all(
            (cfg.boltz_models / name).is_file() and (cfg.boltz_models / name).stat().st_size > 0
            for name in weights
        )
        for name in weights:
            if name.endswith(".ckpt") and (cfg.boltz_models / name).is_file():
                ready = ready and zipfile.is_zipfile(cfg.boltz_models / name)
        if self.identity == "boltz2":
            ready = ready and (cfg.boltz_models / "mols/.kmatch-complete").is_file()
        version = package_version(cfg.boltz_python, "boltz-community")
        ready = ready and version not in ("unknown", "not installed")
        return DependencyStatus(
            self.identity,
            ready,
            "Local MPS engine and assets ready."
            if ready
            else "Install the Boltz environment and model assets with Setup.",
            version,
        )

    def predict(self, sequence, output, cancel, progress, *, msa_path=None):
        status = self.status()
        if not status.ready:
            raise RuntimeError(status.message)
        output.mkdir(parents=True, exist_ok=True)
        if msa_path:
            validate_a3m(msa_path, sequence.sequence)
            copied = output / "query.a3m"
            shutil.copy2(msa_path, copied)
            msa = json.dumps(str(copied))
        else:
            msa = "empty"
        input_path = output / "input.yaml"
        input_path.write_text(
            f"version: 1\nsequences:\n  - protein:\n      id: A\n"
            f"      sequence: {sequence.sequence}\n      msa: {msa}\n"
        )
        raw = output / "raw"
        cfg = self.config
        command = [
            cfg.boltz_python,
            RUNNER,
            "boltz",
            "predict",
            input_path,
            "--out_dir",
            raw,
            "--cache",
            cfg.boltz_models,
            "--accelerator",
            "mps",
            "--model",
            self.identity,
            "--seed",
            str(cfg.seed),
            "--diffusion_samples",
            str(cfg.diffusion_samples),
            "--recycling_steps",
            "3",
            "--sampling_steps",
            "200",
            "--output_format",
            "pdb",
            "--preprocessing-threads",
            "2",
            "--num_workers",
            "0",
        ]
        progress(
            Progress(
                "prediction", f"Folding {sequence.length} residues locally with {self.identity}…"
            )
        )
        runtime = run_process(
            command,
            output / "predictor.log",
            cancel,
            progress,
            "prediction",
            extra_env={"OMP_NUM_THREADS": str(cfg.threads)},
            scratch_root=cfg.cache / "process",
        )
        return normalize_boltz(
            raw,
            output,
            self.identity,
            status.version,
            sequence.length,
            runtime,
            {
                "command": list(map(str, command)),
                "msa": "custom-a3m" if msa_path else "single-sequence",
                "network_during_inference": False,
                "seed": cfg.seed,
                "diffusion_samples": cfg.diffusion_samples,
            },
        )


class AF3Backend:
    identity = "af3"

    def __init__(self, config):
        self.config = config

    def status(self):
        cfg = self.config
        code = cfg.af3_python.is_file() and (cfg.af3_source / "run_alphafold.py").is_file()
        weights = any(
            p.is_file() and p.stat().st_size > 0
            for pattern in ("*.bin", "*.bin.zst")
            for p in cfg.af3_models.glob(pattern)
        )
        message = (
            "External AF3 ready; MPS folding is experimental and unqualified."
            if code and weights
            else (
                "AF3 code installed; supply authorized af3.bin.zst in "
                f"{cfg.af3_models}. See AF3 setup documentation."
            )
            if code
            else (
                "Install optional external AF3 code, then supply authorized weights. "
                "See af3.md in the KBioSuite documentation (docs/kmatch/af3.md)."
            )
        )
        return DependencyStatus(
            "af3", code and weights, message, package_version(cfg.af3_python, "alphafold3")
        )

    def predict(self, sequence, output, cancel, progress, *, msa_path=None):
        status = self.status()
        if not status.ready:
            raise RuntimeError(status.message)
        output.mkdir(parents=True, exist_ok=True)
        msa = validate_a3m(msa_path, sequence.sequence) if msa_path else ""
        input_path = output / "input.json"
        write_json(
            input_path,
            {
                "name": "query",
                "modelSeeds": [self.config.seed],
                "dialect": "alphafold3",
                "version": 3,
                "sequences": [
                    {
                        "protein": {
                            "id": "A",
                            "sequence": sequence.sequence,
                            "unpairedMsa": msa,
                            "pairedMsa": "",
                            "templates": [],
                        }
                    }
                ],
            },
        )
        cfg, raw = self.config, output / "raw"
        command = [
            cfg.af3_python,
            RUNNER,
            "af3",
            cfg.af3_source,
            f"--json_path={input_path}",
            f"--output_dir={raw}",
            f"--model_dir={cfg.af3_models}",
            "--run_data_pipeline=false",
            "--jax_backend=mps",
            "--flash_attention_implementation=xla",
            f"--num_diffusion_samples={cfg.diffusion_samples}",
            "--num_recycles=10",
            f"--jax_compilation_cache_dir={cfg.cache / 'jax'}",
        ]
        progress(Progress("prediction", "Running external AF3 locally on the Apple GPU…"))
        runtime = run_process(
            command,
            output / "predictor.log",
            cancel,
            progress,
            "prediction",
            cwd=cfg.af3_source,
            scratch_root=cfg.cache / "process",
        )
        return normalize_af3(
            raw,
            output,
            status.version,
            sequence.length,
            runtime,
            {
                "command": list(map(str, command)),
                "msa": "custom-a3m" if msa_path else "single-sequence",
                "network_during_inference": False,
                "seed": cfg.seed,
            },
        )


@dataclass(frozen=True)
class PredictorDefinition:
    name: str
    label: str
    note: str
    factory: Callable[..., PredictionBackend]
    production: bool = True


# New engines register here; the GUI consumes descriptors and normalized results.
BACKENDS = {
    "boltz2": PredictorDefinition(
        "Boltz-2",
        "Boltz-2 · Apple GPU",
        "{samples} samples; highest model-confidence prediction selected.",
        BoltzBackend,
    ),
    "af3": PredictorDefinition(
        "AlphaFold 3",
        "AlphaFold 3 · deferred",
        "Deferred: authorized external weights and separate qualification are still required.",
        AF3Backend,
        production=False,
    ),
    "boltz1": PredictorDefinition(
        "Boltz-1",
        "Boltz-1 · benchmark control",
        "Benchmark control only.",
        lambda config: BoltzBackend(config, "boltz1"),
        production=False,
    ),
}


def predictor(config):
    if config.predictor not in BACKENDS:
        raise ValueError(f"Unknown prediction backend: {config.predictor}.")
    return BACKENDS[config.predictor].factory(config)
