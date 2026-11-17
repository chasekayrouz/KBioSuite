"""Manual, one-query classic BLASTP/true-nr canary. Never part of offline tests."""

import argparse
import json
import threading
import uuid
from dataclasses import replace
from datetime import datetime
from pathlib import Path

from kblast.app import InstanceLock
from kblast.config import Config
from kblast.models import Run, Settings
from kblast.ncbi import NcbiClient
from kblast.network import RateLimiter, Transport
from kblast.pipeline import Pipeline
from kblast.sequence import parse_query
from kblast.storage import RunStore, atomic_write

ROOT = Path(__file__).resolve().parents[2]


def qualification_details(run, directory):
    """Separate observed protocol behavior from stages a queued job never reached."""
    requests = []
    log = directory / "requests.jsonl"
    if log.exists():
        requests = [json.loads(line) for line in log.read_text().splitlines()]
    started = [
        e for e in requests if e.get("service") == "blast" and e.get("outcome") == "started"
    ]
    polls = [e for e in started if e.get("poll")]
    intervals = [
        (datetime.fromisoformat(b["at"]) - datetime.fromisoformat(a["at"])).total_seconds()
        for a, b in zip(polls, polls[1:], strict=False)
    ]
    statuses = [e for e in run.blast_events if e["event"] == "status_received"]
    protocol = "not run"
    if polls or statuses:
        protocol = (
            "passed"
            if statuses
            and statuses[-1]["status"] in {"WAITING", "READY", "FAILED", "UNKNOWN"}
            and all(seconds >= 60 for seconds in intervals)
            else "failed"
        )

    def gate(passed, attempted):
        return "passed" if passed else "failed" if attempted else "not run"

    mapping_states = [hit.mapping_status for hit in run.hits]
    return {
        "gates": {
            "submission_rid": "passed" if run.rid else "failed",
            "poll_protocol": protocol,
            "xml_parsing": gate(
                bool(run.blast_metadata),
                any(e["event"] == "parser_started" for e in run.blast_events),
            ),
            "ipg_mapping": gate(
                bool(mapping_states)
                and all(s in {"success", "no_mapping"} for s in mapping_states),
                bool(mapping_states) and run.stage in {"mapping", "contexts", "complete"},
            ),
            "context_retrieval": gate(
                any(o.status == "success" for o in run.occurrences),
                bool(run.occurrences) and run.stage in {"contexts", "complete"},
            ),
        },
        "measurements": {
            "submission_requests": sum(e["parameters"].get("CMD") == "Put" for e in started),
            "poll_requests": len(polls),
            "minimum_poll_interval_seconds": min(intervals) if intervals else None,
            "protocol_warnings": sum(e["event"] == "protocol_warning" for e in run.blast_events),
            "status_transitions": [
                {"at": e["at"], "status": e["status"]}
                for e in run.blast_events
                if e["event"] == "status_transition"
            ],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--live",
        action="store_true",
        required=True,
        help="Explicit consent to send the fixed public fixture to NCBI",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Evidence directory; reuse resumes its saved RID, never resubmits",
    )
    parser.add_argument("--max-wait-seconds", type=int, default=900)
    args = parser.parse_args()
    if args.max_wait_seconds < 60:
        parser.error("Use at least 60 seconds; each RID is polled at most once/minute.")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    config = replace(Config.discover(output / "runs"), blast_timeout_seconds=args.max_wait_seconds)
    lock = InstanceLock(config.data_dir)
    try:
        store = RunStore(config.data_dir)
        pointer = output / "canary-run.json"
        if pointer.exists():
            run = store.load(json.loads(pointer.read_text())["run_id"])
        else:
            run = Run(
                uuid.uuid4().hex,
                parse_query((ROOT / "examples/kblast/actinorhodin_cyclase.fasta").read_text()),
                Settings(top_n=1, flank=20_000, database="nr"),
            )
            store.save(run)
            atomic_write(pointer, json.dumps({"run_id": run.id}) + "\n")
        transport = Transport(
            config,
            threading.Event(),
            store.directory(run.id) / "requests.jsonl",
            RateLimiter(0.11 if config.api_key else 0.36),
            RateLimiter(10),
        )
        report = {
            "status": "not run",
            "scope": "fresh BLASTP/nr → XML → IPG → contexts",
            "public_query": "CAC44203.1",
            "run_id": run.id,
            "parameters": run.settings.blast_parameters(),
        }

        def save():
            report.update(
                rid=run.rid,
                rtoe_seconds=run.rtoe,
                run_status=run.status,
                latest_ncbi_status=run.ncbi_status,
                events=run.blast_events,
                stage_seconds=run.timings,
                error=run.error,
            )
            atomic_write(output / "report.json", json.dumps(report, indent=2) + "\n")

        class ReportingPipeline(Pipeline):
            def checkpoint(self, message=None):
                super().checkpoint(message)
                save()
                if message and not message.startswith("Retrieving context ·"):
                    print(message, flush=True)

        ReportingPipeline(store, run, NcbiClient(transport, store.directory(run.id))).execute()
        report.update(qualification_details(run, store.directory(run.id)))
        if run.blast_metadata:
            try:
                from validate_live import audit_outputs

                audit = audit_outputs(store, run)
                report["audit"] = audit
                report["status"] = (
                    "passed"
                    if audit["passed"] and run.status == "success" and audit["contexts"] > 0
                    else "failed"
                )
            except Exception as exc:
                report.update(status="failed", audit_error=transport.safe_message(str(exc)))
        elif run.ncbi_status == "WAITING" and "BLAST polling timed out" in run.error:
            report.update(
                status="external/variable",
                reason=(
                    "NCBI remained WAITING within the observation window; saved RID is resumable"
                ),
            )
        else:
            report["status"] = "failed"
        if report["gates"]["poll_protocol"] == "failed":
            report["status"] = "failed"
        save()
        print(json.dumps({k: report[k] for k in ["status", "rid", "latest_ncbi_status", "gates"]}))
        if report["status"] == "failed":
            raise SystemExit(1)
    finally:
        lock.close()


if __name__ == "__main__":
    main()
