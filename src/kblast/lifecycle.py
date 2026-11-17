"""Credential-free BLAST events and an injectable clock for remote wait logic."""

import json
import re
import time
from datetime import UTC, datetime


class Clock:
    def monotonic(self):
        return time.monotonic()

    def now(self):
        return datetime.now(UTC)

    def wait(self, cancel, seconds):
        return cancel.wait(max(0, seconds))


def validate_rid(rid, source="external"):
    if not isinstance(rid, str) or not re.fullmatch(r"[A-Z0-9][A-Z0-9-]{5,63}", rid):
        raise ValueError("Enter a valid NCBI RID (6–64 uppercase letters, digits or hyphens).")
    if source not in {"external", "browser"}:
        raise ValueError("RID source must be external or browser.")
    return rid


def remember_rid(run):
    """Capture each job independently; submitted settings always belong to its RID."""
    if run.rid and not any(
        e["rid"] == run.rid and e["adopted_at"] == run.rid_adopted_at for e in run.rid_history
    ):
        run.rid_history.append(
            {
                "rid": run.rid,
                "source": run.rid_source or "kblast_url_api",
                "acquired_at": run.rid_adopted_at or run.submitted_at,
                "adopted_at": run.rid_adopted_at,
                "submitted_at": run.submitted_at if run.rid_source == "kblast_url_api" else "",
                "submitted_parameters": (
                    dict(run.submitted_blast_parameters)
                    if run.rid_source == "kblast_url_api"
                    else {}
                ),
                "rtoe_seconds": run.rtoe if run.rid_source == "kblast_url_api" else None,
                "returned_metadata": dict(run.blast_metadata),
                "validation": dict(run.rid_validation),
            }
        )


def update_rid_result(run):
    remember_rid(run)
    for entry in reversed(run.rid_history):
        if entry["rid"] == run.rid:
            entry["returned_metadata"] = dict(run.blast_metadata)
            entry["validation"] = dict(run.rid_validation)
            break


def record(run, directory, event, *, at=None, phase=None, **details):
    if phase:
        run.blast_phase = phase
    entry = {"at": (at or datetime.now(UTC)).isoformat(), "event": event, "rid": run.rid, **details}
    run.blast_events.append(entry)
    with (directory / "blast-lifecycle.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(entry) + "\n")
    return entry
