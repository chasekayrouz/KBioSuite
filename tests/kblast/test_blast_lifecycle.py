"""Hour-scale URL API lifecycles in milliseconds, with actual paced Transport."""

import json
import threading
from datetime import UTC, datetime, timedelta

import pytest
import requests
from conftest import PROTEIN, blast_xml

from kblast.config import Config
from kblast.models import Run, Settings
from kblast.ncbi import NcbiClient
from kblast.network import Cancelled, RateLimiter, RemoteError, Transport
from kblast.pipeline import Pipeline
from kblast.sequence import parse_query
from kblast.storage import RunStore, atomic_write


class FakeClock:
    def __init__(self):
        self.seconds = 100.0
        self.base = datetime(2026, 10, 7, tzinfo=UTC)
        self.cancel_at = None

    def monotonic(self):
        return self.seconds

    def now(self):
        return self.base + timedelta(seconds=self.seconds)

    def wait(self, cancel, delay):
        if cancel.is_set():
            return True
        if self.cancel_at is not None and self.seconds + delay >= self.cancel_at:
            self.seconds = self.cancel_at
            cancel.set()
            return True
        self.seconds += max(0, delay)
        return False


def rig(tmp_path, monkeypatch, responses, *, clock=None, run=None):
    clock = clock or FakeClock()
    store = RunStore(tmp_path / "runs")
    run = run or Run(
        "a" * 32,
        parse_query(PROTEIN),
        Settings(2, 100),
        rid="TEST-RID",
        submitted_at=clock.now().isoformat(),
    )
    store.save(run)
    http = Transport(
        Config(api_key="test_secret"),
        threading.Event(),
        store.directory(run.id) / "requests.jsonl",
        RateLimiter(0.11, clock=clock),
        RateLimiter(10, clock=clock),
        clock=clock,
    )
    sent = []

    def post(url, data, timeout):
        sent.append((clock.seconds, dict(data)))
        reply = responses.pop(0)
        if isinstance(reply, Exception):
            raise reply
        response = requests.Response()
        response.status_code = 200
        response._content = reply.encode()
        response.request = requests.Request("POST", url, data=data).prepare()
        return response

    monkeypatch.setattr(http.session, "post", post)
    client = NcbiClient(http, store.directory(run.id), clock=clock)
    return clock, store, run, client, sent


@pytest.mark.parametrize(
    "responses",
    [
        ["Status=WAITING", "Status=READY", blast_xml()],
        ["Status=READY", blast_xml()],
        [blast_xml()],
        ["<html>temporary gateway page</html>", "Status=WAITING", blast_xml()],
        [requests.Timeout("test_secret"), "Status=WAITING", blast_xml()],
    ],
)
def test_protocol_transitions_and_retry(tmp_path, monkeypatch, responses):
    clock, store, run, client, sent = rig(tmp_path, monkeypatch, list(responses))
    xml = client.wait_for_blast(run, lambda store=store, run=run: store.save(run))
    assert "BlastXML2" in xml and run.ncbi_status == "READY"
    polls = [t for t, _ in sent]
    # Retrieval immediately after READY is a separate request under the 10-second gate.
    assert all(b - a >= 10 for a, b in zip(polls, polls[1:], strict=False))
    logged = [
        json.loads(line)
        for line in (store.directory(run.id) / "requests.jsonl").read_text().splitlines()
    ]
    status_polls = [e for e in logged if e.get("poll") and e["outcome"] == "started"]
    assert all(
        (datetime.fromisoformat(b["at"]) - datetime.fromisoformat(a["at"])).total_seconds() >= 60
        for a, b in zip(status_polls, status_polls[1:], strict=False)
    )
    assert all(p["CMD"] == "Get" and p["RID"] == "TEST-RID" and "api_key" not in p for _, p in sent)
    events = [e["event"] for e in run.blast_events]
    assert "ready_detected" in events and "result_received" in events
    assert store.load(run.id).blast_events == run.blast_events
    assert "test_secret" not in (store.directory(run.id) / "blast-lifecycle.jsonl").read_text()


@pytest.mark.parametrize(
    "reply,match",
    [
        ("Status=FAILED", "failed"),
        ("Status=UNKNOWN", "unknown"),
        ("<bad>no status", "neither a BLAST status"),
    ],
)
def test_terminal_and_malformed_are_bounded(tmp_path, monkeypatch, reply, match):
    _, store, run, client, sent = rig(tmp_path, monkeypatch, [reply] * 3)
    with pytest.raises(RemoteError, match=match):
        client.wait_for_blast(run, lambda store=store, run=run: store.save(run))
    assert len(sent) == (3 if reply.startswith("<bad>") else 1)
    assert all(p["CMD"] == "Get" for _, p in sent)
    assert list((store.directory(run.id) / "raw").glob("blast-status-*.txt"))


def test_hours_of_waiting_and_elapsed_timeout(tmp_path, monkeypatch):
    clock, store, run, client, sent = rig(tmp_path, monkeypatch, ["Status=WAITING"] * 400)
    with pytest.raises(RemoteError, match="polling timed out"):
        client.wait_for_blast(run, lambda store=store, run=run: store.save(run))
    assert clock.seconds >= 6 * 3600 + 100
    assert 350 <= len(sent) <= 360
    assert [e["status"] for e in run.blast_events if e["event"] == "status_transition"] == [
        "WAITING"
    ]
    assert len([e for e in run.blast_events if e["event"] == "status_received"]) == len(sent)


def test_repeated_cancel_resume_preserves_rid_and_minute_gate(tmp_path, monkeypatch):
    clock = FakeClock()
    run = None
    all_sent = []
    for _ in range(12):
        clock, store, run, client, sent = rig(
            tmp_path, monkeypatch, ["Status=WAITING"], clock=clock, run=run
        )
        clock.cancel_at = clock.seconds + 60
        with pytest.raises(Cancelled):
            client.wait_for_blast(run, lambda store=store, run=run: store.save(run))
        all_sent.extend(sent)
        run = store.load(run.id)
    clock.cancel_at = None
    clock, store, run, client, sent = rig(
        tmp_path, monkeypatch, [blast_xml()], clock=clock, run=run
    )
    assert "BlastXML2" in client.wait_for_blast(run, lambda store=store, run=run: store.save(run))
    all_sent.extend(sent)
    assert all(b[0] - a[0] >= 60 for a, b in zip(all_sent, all_sent[1:], strict=False))
    assert all(p["RID"] == "TEST-RID" and p["CMD"] == "Get" for _, p in all_sent)


def test_pipeline_reuses_saved_rid_and_local_xml(tmp_path, monkeypatch):
    clock, store, run, client, sent = rig(tmp_path, monkeypatch, [blast_xml()])
    pipeline = Pipeline(store, run, client)
    pipeline.blast()
    assert len(run.hits) == 2 and run.blast_phase == "complete"
    restored = store.load(run.id)
    Pipeline(store, restored, client).blast()
    assert len(sent) == 1 and sent[0][1]["CMD"] == "Get"
    events = [e["event"] for e in restored.blast_events]
    assert (
        "parser_started" in events
        and "parser_completed" in events
        and "local_result_reused" in events
    )


def test_submission_receipt_recovers_after_interruption_without_put(tmp_path, monkeypatch):
    clock, store, run, client, sent = rig(tmp_path, monkeypatch, [])
    run.rid = ""
    run.submission_attempted = True
    atomic_write(store.directory(run.id) / "raw/blast-submission.txt", "RID=RECOVERED\nRTOE=99\n")
    Pipeline(store, run, client).ensure_blast_submission()
    assert run.rid == "RECOVERED" and run.rtoe == 99 and not sent
    assert run.blast_events[-1]["event"] == "rid_recovered"


def test_waiting_cancellation_is_checkpointed_and_resume_finishes_same_job(tmp_path, monkeypatch):
    clock, store, run, client, first_sent = rig(tmp_path, monkeypatch, ["Status=WAITING"])
    clock.cancel_at = clock.seconds + 20
    Pipeline(store, run, client).execute()
    saved = store.load(run.id)
    assert saved.status == "cancelled" and saved.rid == "TEST-RID"
    assert saved.blast_events[-1]["event"] == "cancelled"
    assert saved.ncbi_status == "WAITING"
    clock.cancel_at = None
    xml = blast_xml()
    xml = xml.replace(xml.split("<hits>")[1].split("</hits>")[0], "")
    clock, store, saved, client, second_sent = rig(
        tmp_path, monkeypatch, [xml], clock=clock, run=saved
    )
    Pipeline(store, saved, client).execute()
    assert saved.status == "no_hits" and saved.blast_phase == "complete"
    assert all(params["CMD"] == "Get" for _, params in first_sent + second_sent)
    assert second_sent[0][0] - first_sent[0][0] >= 60
    assert "resumed" in [e["event"] for e in saved.blast_events]
    assert "parser_completed" in [e["event"] for e in saved.blast_events]
