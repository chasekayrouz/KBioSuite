"""A queued or malformed live canary must never claim downstream qualification."""

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/kblast/live_canary.py"
spec = importlib.util.spec_from_file_location("live_canary", SCRIPT)
canary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(canary)


@pytest.mark.parametrize(
    "statuses,spacing,expected",
    [
        (["WAITING", "WAITING"], 60, "passed"),
        (["WAITING", "UNRECOGNIZED"], 60, "failed"),
        (["UNRECOGNIZED", "READY"], 60, "passed"),
        (["WAITING", "READY"], 59, "failed"),
        ([], 60, "not run"),
    ],
)
def test_canary_reports_only_observed_gates(tmp_path, statuses, spacing, expected):
    events = [{"event": "status_received", "status": status} for status in statuses]
    run = SimpleNamespace(
        rid="PUBLIC-TEST-RID", blast_events=events, blast_metadata={},
        hits=[], occurrences=[], stage="blast",
    )
    entries = [
        {
            "service": "blast",
            "outcome": "started",
            "poll": True,
            "parameters": {"CMD": "Get", "RID": run.rid},
            "at": f"2026-10-07T00:{i * spacing // 60:02d}:{i * spacing % 60:02d}+00:00",
        }
        for i in range(len(statuses))
    ]
    entries.append({
        "operation": "ipg_batch", "outcome": "individual_fallback",
        "error": "Offline fallback fixture",
    })
    (tmp_path / "requests.jsonl").write_text("".join(json.dumps(e) + "\n" for e in entries))
    details = canary.qualification_details(run, tmp_path)
    assert details["gates"]["poll_protocol"] == expected
    for stage in ["xml_parsing", "ipg_mapping", "context_retrieval"]:
        assert details["gates"][stage] == "not run"
    assert details["measurements"]["submission_requests"] == 0
    assert details["measurements"]["poll_requests"] == len(statuses)


def test_canary_distinguishes_failed_stages_from_unexercised_stages(tmp_path):
    run = SimpleNamespace(
        rid="PUBLIC-TEST-RID", blast_events=[{"event": "parser_started"}],
        blast_metadata={}, hits=[], occurrences=[], stage="blast",
    )
    assert canary.qualification_details(run, tmp_path)["gates"]["xml_parsing"] == "failed"
    run.blast_metadata = {"program": "blastp"}
    run.hits = [SimpleNamespace(mapping_status="retrieval_failure")]
    run.stage = "mapping"
    gates = canary.qualification_details(run, tmp_path)["gates"]
    assert gates["xml_parsing"] == "passed" and gates["ipg_mapping"] == "failed"
    assert gates["context_retrieval"] == "not run"
    run.hits[0].mapping_status = "success"
    run.occurrences = [SimpleNamespace(status="retrieval_failure")]
    run.stage = "contexts"
    gates = canary.qualification_details(run, tmp_path)["gates"]
    assert gates["ipg_mapping"] == "passed" and gates["context_retrieval"] == "failed"
