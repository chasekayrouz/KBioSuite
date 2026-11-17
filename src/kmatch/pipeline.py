"""One local single-query run, independently callable by the CLI, desktop app or a library."""

import hashlib
import shutil
import time
import uuid
from dataclasses import dataclass
from pathlib import Path

from .models import PredictionResult, Progress, SearchResult, write_json
from .predictors import RUNNER, predictor, validate_a3m
from .process import Cancellation, Cancelled, run_process
from .search import FoldseekBackend, hit_rankings, retained_hit_indices
from .sequence import parse_sequence


@dataclass(frozen=True)
class RunResult:
    directory: Path
    prediction: PredictionResult
    hits: list[SearchResult]


class MSACache:
    def __init__(self, config):
        self.config = config

    def get(self, sequence, cancel, progress):
        key = hashlib.sha256(("colabfold-v1:" + sequence.sequence).encode()).hexdigest()
        root = self.config.cache / "msa" / key
        root.mkdir(parents=True, exist_ok=True, mode=0o700)
        alignment = root / "query.a3m"
        if alignment.exists():
            validate_a3m(alignment, sequence.sequence)
            progress(
                Progress("msa", "Using cached network-assisted MSA; no new sequence submission.")
            )
            return alignment
        sequence_file = root / "query.txt"
        sequence_file.write_text(sequence.sequence)
        progress(
            Progress(
                "msa",
                (
                    "Network-assisted MSA: submitting sequence to api.colabfold.com. "
                    "Folding stays local."
                ),
            )
        )
        run_process(
            [self.config.boltz_python, RUNNER, "msa", sequence_file, root],
            root / "msa.log",
            cancel,
            progress,
            "MSA retrieval",
            timeout=900,
            scratch_root=self.config.cache / "process",
        )
        validate_a3m(alignment, sequence.sequence)
        return alignment


class Pipeline:
    def __init__(self, config, prediction_backend=None, search_backend=None):
        self.config = config
        self.prediction_backend = prediction_backend or predictor(config)
        self.search_backend = search_backend or FoldseekBackend(config)

    def run(
        self,
        text,
        *,
        cancel=None,
        progress=lambda event: None,
        network_msa=False,
        msa_path=None,
        on_prediction=lambda result: None,
        on_directory=lambda path: None,
    ):
        sequence = parse_sequence(text)
        cancel = cancel or Cancellation()
        cancel.check()
        for backend in (self.prediction_backend, self.search_backend):
            status = backend.status()
            if not status.ready:
                raise RuntimeError(status.message)
        if network_msa and msa_path:
            raise ValueError("Choose network-assisted MSA or a custom A3M, not both.")
        run_id = time.strftime("%Y%m%d-%H%M%S") + "-" + uuid.uuid4().hex[:8]
        directory = self.config.runs / run_id
        directory.mkdir(parents=True, mode=0o700)
        on_directory(directory)
        manifest = {
            "schema_version": 1,
            "run_id": run_id,
            "status": "running",
            "sequence_length": sequence.length,
            "sequence_identifier": sequence.identifier,
            "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "msa_mode": "network-assisted"
            if network_msa
            else "custom-a3m"
            if msa_path
            else "single-sequence",
            "structure_inference": "local",
            "search": "local",
        }
        write_json(directory / "run.json", manifest)
        (directory / "sequence.fasta").write_text(f">{sequence.identifier}\n{sequence.sequence}\n")
        started = time.monotonic()
        try:
            if network_msa:
                msa_path = MSACache(self.config).get(sequence, cancel, progress)
            if msa_path:
                validate_a3m(Path(msa_path), sequence.sequence)
            prediction = self.prediction_backend.predict(
                sequence,
                directory / "prediction",
                cancel,
                progress,
                msa_path=Path(msa_path) if msa_path else None,
            )
            cancel.check()
            on_prediction(prediction)
            for warning in prediction.warnings:
                progress(Progress("warning", warning["message"]))
            hits = self.search_backend.search(prediction, directory / "search", cancel, progress)
            rankings = hit_rankings(hits)
            manifest.update(
                status="complete",
                runtime_seconds=time.monotonic() - started,
                displayed_hits=len(rankings["qtm"]),
                returned_hits=len(hits),
                display_candidate_hits=len(retained_hit_indices(rankings)),
                prediction=str(prediction.structure_path),
                prediction_warnings=prediction.warnings,
            )
            write_json(directory / "run.json", manifest)
            progress(
                Progress(
                    "complete",
                    f"Complete · {len(hits)} homologs · up to 50 per ranking",
                    manifest["runtime_seconds"],
                )
            )
            return RunResult(directory, prediction, hits)
        except BaseException as exc:
            manifest.update(
                status="cancelled" if isinstance(exc, (Cancelled, KeyboardInterrupt)) else "failed",
                error=str(exc),
                runtime_seconds=time.monotonic() - started,
            )
            write_json(directory / "run.json", manifest)
            raise
        finally:
            # Only Kmatch-owned transient search files are removed. Scientific outputs/logs survive.
            shutil.rmtree(directory / "search/work", ignore_errors=True)
