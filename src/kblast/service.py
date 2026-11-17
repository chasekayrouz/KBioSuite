"""One cancellable worker; the GUI consumes immutable on-disk checkpoints."""

import threading
import uuid
from dataclasses import replace

from .config import Config
from .lifecycle import record, remember_rid, validate_rid
from .models import Run, Settings, utc_now
from .ncbi import NcbiClient
from .network import RateLimiter, Transport
from .pipeline import Pipeline
from .sequence import parse_query
from .storage import RunStore


class RunManager:
    def __init__(self, config: Config):
        self.config = config
        self.store = RunStore(config.data_dir)
        self.lock = threading.Lock()
        self.worker: threading.Thread | None = None
        self.active_id: str | None = None
        self.cancel_event = threading.Event()
        self.eutils_limiter = RateLimiter(0.11 if config.api_key else 0.36)
        self.blast_limiter = RateLimiter(10)
        for entry in self.store.list_runs():
            run = self.store.load(entry["id"])
            if run.status in {"running", "pending"}:
                record(run, self.store.directory(run.id), "interrupted", phase="interrupted")
                run.status = "interrupted"
                run.message = (
                    "Application stopped before completion. Resume using the saved checkpoint."
                )
                for hit in run.hits:
                    if hit.mapping_status == "running":
                        hit.mapping_status = "pending"
                    if hit.sequence_status == "running":
                        hit.sequence_status = "pending"
                for occurrence in run.occurrences:
                    if occurrence.status == "running":
                        occurrence.status = "pending"
                self.store.save(run)

    def is_active(self, run_id: str | None = None) -> bool:
        return bool(
            self.worker and self.worker.is_alive() and (run_id is None or self.active_id == run_id)
        )

    def start(
        self,
        text: str,
        top_n: int,
        flank: int,
        email: str = "",
        rid: str = "",
        rid_source: str = "external",
    ) -> Run:
        query, settings = parse_query(text), Settings(top_n, flank)
        if rid != "":
            validate_rid(rid, rid_source)
        with self.lock:
            if self.is_active():
                raise ValueError("A run is already active. Cancel it before starting another.")
            run = Run(uuid.uuid4().hex, query, settings)
            if rid:
                self.adopt(run, rid, rid_source)
            self.store.save(run)
            self.launch(run, email)
            return run

    def resume(
        self,
        run_id: str,
        email: str = "",
        rid: str = "",
        rid_source: str = "external",
    ) -> Run:
        if rid != "":
            validate_rid(rid, rid_source)
        with self.lock:
            if self.is_active():
                raise ValueError("A run is already active.")
            run = self.store.load(run_id)
            if rid:
                self.adopt(run, rid, rid_source)
                self.store.save(run)
            self.launch(run, email)
            return run

    def adopt(self, run: Run, rid: str, source: str) -> None:
        """Called only before launch or after the previous worker has stopped."""
        validate_rid(rid, source)
        if rid == run.rid:
            return  # Resuming the same job does not reset its pacing or provenance.
        if run.hits or run.occurrences or run.blast_metadata:
            raise ValueError(
                "This run already has parsed results. Use a fresh run for another RID."
            )
        if run.settings.database != "nr":
            raise ValueError("Existing RIDs require a fresh true-nr run.")
        self.store.directory(run.id).mkdir(parents=True, exist_ok=True)
        remember_rid(run)
        at = utc_now()
        old = run.rid
        previous = next((e for e in reversed(run.rid_history) if e["rid"] == rid), None)
        for entry in reversed(run.rid_history):
            if entry["rid"] == old:
                entry["superseded_at"] = at
                entry["superseded_by"] = rid
                break
        # Retain rejected/unfinished evidence without reusing it for the new job.
        xml = self.store.directory(run.id) / "raw" / "blast.xml"
        if xml.is_file():
            from .storage import safe_name

            xml.replace(xml.with_name(f"blast-{safe_name(old)}-{len(run.rid_history)}.xml"))
        run.rid, run.rid_source, run.rid_adopted_at = (
            rid,
            previous["source"] if previous else source,
            at,
        )
        run.rtoe = (previous.get("rtoe_seconds") or 0) if previous else 0
        run.last_blast_get_at, run.ncbi_status = "", ""
        run.blast_completed_at = ""
        run.rid_validation = {"status": "pending"}
        remember_rid(run)
        record(
            run,
            self.store.directory(run.id),
            "rid_adopted",
            phase="waiting",
            source=run.rid_source,
            superseded_rid=old,
        )

    def launch(self, run: Run, email: str) -> None:
        self.active_id = run.id
        self.cancel_event = threading.Event()
        config = replace(self.config, email=email.strip() or self.config.email)
        http = Transport(
            config,
            self.cancel_event,
            self.store.directory(run.id) / "requests.jsonl",
            self.eutils_limiter,
            self.blast_limiter,
        )
        pipeline = Pipeline(self.store, run, NcbiClient(http, self.store.directory(run.id)))
        self.worker = threading.Thread(target=pipeline.execute, name="kblast-worker", daemon=True)
        self.worker.start()

    def cancel(self, run_id: str) -> None:
        with self.lock:
            if self.is_active(run_id):
                self.cancel_event.set()

    def snapshot(self, run_id: str) -> dict:
        result = self.store.load(run_id).to_dict()
        result["active"] = self.is_active(run_id)
        result["cancelling"] = result["active"] and self.cancel_event.is_set()
        result["data_directory"] = str(self.store.directory(run_id))
        # Full sequences and all HSP alignments are saved, but unnecessary for GUI polling.
        for hit in result["hits"]:
            hit["hsp_count"] = len(hit["hsps"])
            del hit["sequence"], hit["hsps"], hit["descriptions"]
        return result

    def shutdown(self) -> None:
        self.cancel_event.set()
        if self.worker:
            self.worker.join(timeout=2)
