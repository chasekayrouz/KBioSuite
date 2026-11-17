"""Use disposable inputs and real local child processes to qualify Stop ownership."""

import json
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import psutil
import pytest

from kannotate import AnnotationEngine
from kannotate.batch import annotate_batch
from kannotate.cancellation import (
    RunCancelled,
    RunControl,
    check_cancelled,
    current_control,
    run_scope,
)
from kannotate.gui import AnnotationJob, DesktopApi, create_app
from kannotate.quick_search import parse_queries
from kannotate.search.cache import CachedAdapter, sequence_hash
from kannotate.search.common import run_search
from kannotate.search.hmmer import HmmerAdapter
from kannotate.search.rpsblast import RpsBlastAdapter


def wait_for(check, timeout=5):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        value = check()
        if value:
            return value
        time.sleep(0.01)
    raise AssertionError("Test condition did not complete")


def still_running(pid):
    try:
        return psutil.Process(pid).status() != psutil.STATUS_ZOMBIE
    except psutil.NoSuchProcess:
        return False


def test_search_wait_intervals_are_not_an_overall_deadline(monkeypatch):
    calls = []

    class Process:
        returncode = 0

        def communicate(self, timeout):
            calls.append(timeout)
            if len(calls) < 4:
                raise subprocess.TimeoutExpired("synthetic search", timeout=1801)
            return "", ""

    options = {}

    def popen(command, **kwargs):
        options.update(kwargs)
        return Process()

    monkeypatch.setattr("kannotate.search.common.subprocess.Popen", popen)
    run_search(["synthetic search"])
    assert calls == [0.1] * 4
    assert "timeout" not in options
    assert options["start_new_session"]


def test_stop_kills_both_process_groups_and_removes_workdirs_and_cdd_links(tmp_path, monkeypatch):
    workspace = tmp_path / "temporary search data"
    workspace.mkdir()
    monkeypatch.setattr("tempfile.tempdir", str(workspace))
    database = tmp_path / "installed database"
    database.mkdir()
    asset = database / "Cdd.pal"
    asset.write_bytes(b"synthetic installed database; preserve these bytes")
    executable = tmp_path / "search"
    # Both leaders and workers ignore TERM, exercising bounded KILL escalation.
    executable.write_text(
        f"#!{sys.executable}\n"
        "import json, os, pathlib, signal, subprocess, sys, time\n"
        "signal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
        "args=sys.argv[1:]\n"
        "hmm='--domtblout' in args\n"
        "query=pathlib.Path(args[-1] if hmm else args[args.index('-query')+1])\n"
        "out=pathlib.Path(args[args.index('--domtblout' if hmm else '-out')+1])\n"
        "out.write_text('partial synthetic search output')\n"
        "child=subprocess.Popen([sys.executable,'-c',"
        "'import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); time.sleep(60)'])\n"
        "info={'pid':os.getpid(),'child':child.pid,'workdir':str(query.parent.parent)}\n"
        "if not hmm: info['cdd_view']=str(pathlib.Path(args[args.index('-db')+1]).parent)\n"
        f"marker=pathlib.Path({str(tmp_path)!r})/('hmm.json' if hmm else 'rps.json')\n"
        "marker.write_text(json.dumps(info))\n"
        "time.sleep(60)\n"
    )
    executable.chmod(0o755)
    engine = AnnotationEngine(
        HmmerAdapter(str(executable), str(database / "Pfam")),
        RpsBlastAdapter(str(executable), str(database / "Cdd")),
    )
    control = RunControl()
    markers = [tmp_path / "hmm.json", tmp_path / "rps.json"]
    with ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(engine.annotate, parse_queries("MACDEFG"), threads=4, control=control)
        try:
            wait_for(lambda: all(p.exists() and p.stat().st_size for p in markers))
            evidence = [json.loads(p.read_text()) for p in markers]
            assert all(still_running(e["pid"]) and still_running(e["child"]) for e in evidence)
            assert control.request_stop()
            with pytest.raises(RunCancelled):
                future.result(timeout=5)
        finally:
            control.request_stop()
        wait_for(lambda: all(not still_running(e[k]) for e in evidence for k in ("pid", "child")))
    assert not list(workspace.iterdir())
    assert all(not Path(e["workdir"]).exists() for e in evidence)
    assert not Path(evidence[1]["cdd_view"]).exists()
    assert asset.read_bytes() == b"synthetic installed database; preserve these bytes"


def test_stop_during_process_launch_cannot_leave_a_child_running(monkeypatch):
    control = RunControl()
    original_popen = subprocess.Popen
    children = []

    def popen(*args, **kwargs):
        process = original_popen(*args, **kwargs)
        children.append(process)
        control.request_stop()
        return process

    monkeypatch.setattr("kannotate.search.common.subprocess.Popen", popen)
    with pytest.raises(RunCancelled), run_scope(control):
        run_search([sys.executable, "-c", "import time; time.sleep(60)"])
    assert len(children) == 1 and not still_running(children[0].pid)


class NoHits:
    def search(self, queries, workdir, threads):
        (workdir / "raw.tsv").write_text("synthetic zero-hit search")
        return []


def file_bytes(directory):
    return {
        str(p.relative_to(directory)): p.read_bytes() for p in directory.rglob("*") if p.is_file()
    }


def test_cancel_rolls_back_completed_engine_cache_and_restores_preexisting_bytes(tmp_path):
    control = RunControl()
    hmm = CachedAdapter(NoHits(), tmp_path / "hmm", {"version": 1})
    rps = CachedAdapter(NoHits(), tmp_path / "rps", {"version": 1})
    engine = AnnotationEngine(hmm, rps)
    original = parse_queries("MACDEFG")
    engine.annotate(original, threads=2)
    # A corrupt pre-existing entry is recomputed during the next run, then restored on Stop.
    old = hmm.directory / f"{sequence_hash(original[0].sequence)}.json"
    old.write_bytes(b"pre-existing corrupt cache; preserve on rollback")
    before = file_bytes(tmp_path)
    entered = threading.Event()

    class Blocked:
        def search(self, *args):
            entered.set()
            while True:
                check_cancelled()
                time.sleep(0.01)

    rps.adapter = Blocked()
    proteins = parse_queries(">old\nMACDEFG\n>new\nACDEFGH")
    with ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(engine.annotate, proteins, threads=2, control=control)
        try:
            assert entered.wait(2)
            new = hmm.directory / f"{sequence_hash(proteins[1].sequence)}.json"
            wait_for(new.exists)
            control.request_stop()
            with pytest.raises(RunCancelled):
                future.result(timeout=3)
        finally:
            control.request_stop()
    assert file_bytes(tmp_path) == before


@pytest.mark.parametrize("force", [False, True])
def test_stop_during_output_write_keeps_completed_files_inputs_and_existing_output(
    tmp_path, monkeypatch, force
):
    from kannotate.io.xlsx import write_xlsx

    inputs = [tmp_path / f"{name}.faa" for name in ("first", "second", "third")]
    for path in inputs:
        path.write_text(f">{path.stem}\nMACDEFG\n")
    originals = [p.read_bytes() for p in inputs]
    output = tmp_path / "outputs"
    output.mkdir()
    existing = output / "second.kannotate.xlsx"
    if force:
        existing.write_bytes(b"pre-existing completed output")
    entered = threading.Event()
    written = []

    def writer(annotations, path, **kwargs):
        written.append(annotations[0].protein.header)
        if written[-1] == "second":
            path.write_bytes(b"unfinished workbook")
            entered.set()
            while True:
                check_cancelled()
                time.sleep(0.01)
        return write_xlsx(annotations, path, **kwargs)

    monkeypatch.setattr("kannotate.batch.write_xlsx", writer)
    control = RunControl()
    reports = []
    with ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(
            annotate_batch,
            inputs,
            output,
            AnnotationEngine(NoHits(), NoHits()),
            threads=2,
            force=force,
            control=control,
            on_result=reports.append,
            table_format="xlsx",
        )
        try:
            assert entered.wait(3)
            control.request_stop()
            with pytest.raises(RunCancelled):
                future.result(timeout=3)
        finally:
            control.request_stop()
    assert written == ["first", "second"]
    assert len(reports) == 1 and reports[0].output.name == "first.kannotate.xlsx"
    assert (output / "first.kannotate.xlsx").is_file()
    assert not (output / "third.kannotate.xlsx").exists()
    assert not list(output.glob(".kannotate-*"))
    assert [p.read_bytes() for p in inputs] == originals
    if force:
        assert existing.read_bytes() == b"pre-existing completed output"
    else:
        assert not existing.exists()


def test_stop_immediately_before_publication_discards_finished_temp_output(tmp_path, monkeypatch):
    source = tmp_path / "input.faa"
    source.write_text(">synthetic\nMACDEFG\n")
    control = RunControl()

    def writer(annotations, path, **kwargs):
        path.write_bytes(b"finished but not published")
        control.request_stop()

    monkeypatch.setattr("kannotate.batch.write_xlsx", writer)
    output = tmp_path / "output"
    with pytest.raises(RunCancelled):
        annotate_batch(
            [source],
            output,
            AnnotationEngine(NoHits(), NoHits()),
            control=control,
            table_format="xlsx",
        )
    assert not list(output.iterdir())


def test_stop_after_publication_keeps_completed_file(tmp_path):
    source = tmp_path / "input.faa"
    source.write_text(">synthetic\nMACDEFG\n")
    control = RunControl()
    output = tmp_path / "output"
    with pytest.raises(RunCancelled):
        annotate_batch(
            [source],
            output,
            AnnotationEngine(NoHits(), NoHits()),
            control=control,
            on_result=lambda result: control.request_stop(),
            table_format="xlsx",
        )
    assert [p.name for p in output.iterdir()] == ["input.kannotate.xlsx"]


def test_confirmed_stop_is_scoped_and_stays_busy_until_cleanup_finishes(monkeypatch):
    entered, cleaning, release = threading.Event(), threading.Event(), threading.Event()
    job = AnnotationJob()
    job.quick = {"earlier": "completed search"}
    job.quick_revision = 7
    job.update(results=[{"output": "earlier-completed.xlsx"}])

    def search(proteins, **kwargs):
        def cleanup():
            cleaning.set()
            assert release.wait(3)

        current_control().on_rollback(cleanup)
        entered.set()
        while True:
            check_cancelled()
            time.sleep(0.01)

    monkeypatch.setattr("kannotate.gui.search", search)
    client = create_app(job).test_client()
    token = re.search(r'name="desktop-token" content="([^"]+)"', client.get("/").data.decode())[1]
    headers = {"X-Desktop-Token": token}
    job.quick_search(parse_queries("MACDEFG"))
    assert entered.wait(2)
    job_id = job.snapshot()["job_id"]
    confirmed = {"job_id": job_id, "confirmed": True}
    try:
        assert client.post("/api/stop", json=confirmed).status_code == 403
        assert (
            client.post(
                "/api/stop", headers={**headers, "Origin": "https://example.com"}, json=confirmed
            ).status_code
            == 403
        )
        for data in ({"job_id": job_id}, {"job_id": job_id, "confirmed": "true"}, []):
            assert client.post("/api/stop", headers=headers, json=data).status_code == 400
        assert (
            client.post(
                "/api/stop", headers=headers, json={"job_id": "earlier-run", "confirmed": True}
            ).status_code
            == 400
        )
        assert not current_control()
        assert client.post("/api/stop", headers=headers, json=confirmed).status_code == 202
        assert cleaning.wait(2)
        state = job.snapshot()
        assert state["busy"] and state["stopping"] and not state["can_stop"]
        assert not DesktopApi(job).can_close()
        with pytest.raises(ValueError, match="current operation"):
            job.start(lambda: None)
    finally:
        release.set()
    wait_for(lambda: not job.snapshot()["busy"])
    state = job.snapshot()
    assert state["stopped"] and not state["error"] and not state["stopping"]
    assert state["results"] == [{"output": "earlier-completed.xlsx"}]
    assert job.quick == {"earlier": "completed search"} and job.quick_revision == 7
    assert DesktopApi(job).can_close()
    assert client.post("/api/stop", headers=headers, json=confirmed).status_code == 400
    job.start(lambda: None)
    wait_for(lambda: not job.snapshot()["busy"])
    assert not job.snapshot()["stopped"]


def test_database_setup_does_not_offer_stop():
    entered, release = threading.Event(), threading.Event()
    job = AnnotationJob()
    job.start(lambda: (entered.set(), release.wait(3)))
    try:
        assert entered.wait(1)
        assert not job.snapshot()["can_stop"]
        with pytest.raises(ValueError, match="cannot be stopped"):
            job.stop(job.snapshot()["job_id"])
    finally:
        release.set()


def test_stop_reports_cleanup_failure_instead_of_claiming_data_removed():
    entered = threading.Event()
    job = AnnotationJob()

    def operation():
        def cleanup():
            raise OSError("synthetic cleanup failure")

        current_control().on_rollback(cleanup)
        entered.set()
        while True:
            check_cancelled()
            time.sleep(0.01)

    job.start(operation, kind="files")
    assert entered.wait(1)
    job.stop(job.snapshot()["job_id"])
    wait_for(lambda: not job.snapshot()["busy"])
    state = job.snapshot()
    assert state["stopped"] and "cleanup needs attention" in state["message"]
    assert "synthetic cleanup failure" in state["error"]
