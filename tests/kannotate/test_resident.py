"""Real disposable worker processes; no databases or PyHMMER models in unit tests."""

import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import pytest

from kannotate import AnnotationEngine
from kannotate.cancellation import RunCancelled, RunControl
from kannotate.quick_search import parse_queries
from kannotate.search.cache import CachedAdapter
from kannotate.search.resident import ResidentHmmer

WORKER = """
import json, os, pathlib, socket, sys, time
connection=socket.socket(fileno=int(sys.argv[1]))
stream=connection.makefile('rwb',buffering=0)
for line in stream:
    request=json.loads(line)
    workdir=pathlib.Path(request['workdir'])
    (workdir/'started').write_text(str(os.getpid()))
    if request['identity'].get('hang'): time.sleep(60)
    (workdir/'resident.domtbl').write_text('# synthetic zero-hit table\\n')
    stream.write(json.dumps({'ok':True,'pid':os.getpid()}).encode()+b'\\n')
"""


class Fallback:
    executable = sys.executable
    database = "synthetic"

    def __init__(self):
        self.calls = 0

    def search(self, queries, workdir, threads):
        self.calls += 1
        (workdir / "canonical.domtbl").write_text("# canonical zero-hit table\n")
        return []


def fixture_worker(monkeypatch):
    original = subprocess.Popen
    processes = []

    def popen(command, **kwargs):
        process = original([sys.executable, "-c", WORKER, command[-1]], **kwargs)
        processes.append(process)
        return process

    monkeypatch.setattr("kannotate.search.resident.subprocess.Popen", popen)
    return processes


def test_cache_hits_never_load_models_and_new_sequences_reuse_one_worker(tmp_path, monkeypatch):
    processes = fixture_worker(monkeypatch)
    adapter = ResidentHmmer(Fallback(), {"fixture": 1}, minimum_queries=1)
    cache = CachedAdapter(adapter, tmp_path / "cache", {"backend": "fixture"})
    engine = AnnotationEngine(cache, Fallback())
    try:
        engine.annotate(parse_queries("MACDEFG"), threads=3)
        pid = adapter.process.pid
        engine.annotate(parse_queries("ACDEFGH"), threads=3)
        assert adapter.process.pid == pid and len(processes) == 1
        engine.close()
        assert processes[0].poll() is not None
        engine.annotate(parse_queries("MACDEFG"), threads=3)
        assert len(processes) == 1 and adapter.process is None
        assert all(p["cached"] for p in cache.provenance.values())
    finally:
        engine.close()


def test_stop_reaps_resident_and_rolls_back_only_interrupted_cache(tmp_path, monkeypatch):
    processes = fixture_worker(monkeypatch)
    workspace = tmp_path / "temp"
    workspace.mkdir()
    monkeypatch.setattr("tempfile.tempdir", str(workspace))
    adapter = ResidentHmmer(Fallback(), {"hang": True}, minimum_queries=1)
    engine = AnnotationEngine(CachedAdapter(adapter, tmp_path / "cache", {}), Fallback())
    control = RunControl()
    try:
        with ThreadPoolExecutor(max_workers=1) as pool:
            future = pool.submit(
                engine.annotate, parse_queries("MACDEFG"), threads=3, control=control
            )
            deadline = time.monotonic() + 3
            while not list(workspace.rglob("started")) and time.monotonic() < deadline:
                time.sleep(0.01)
            assert list(workspace.rglob("started"))
            control.request_stop()
            with pytest.raises(RunCancelled):
                future.result(timeout=3)
        assert all(p.poll() is not None for p in processes)
        assert not list(workspace.iterdir()) and not list((tmp_path / "cache").rglob("*.json"))
    finally:
        engine.close()


def test_idle_expiry_and_memory_failure_release_native_process(tmp_path, monkeypatch):
    processes = fixture_worker(monkeypatch)
    monkeypatch.setattr("kannotate.search.resident.IDLE_SECONDS", 0.1)
    adapter = ResidentHmmer(Fallback(), {}, minimum_queries=1)
    query = tmp_path / "queries.faa"
    query.write_text(">q\nMACDEFG\n")
    first = tmp_path / "first"
    first.mkdir()
    adapter.search(query, first, 3)
    deadline = time.monotonic() + 3
    while adapter.process is not None and time.monotonic() < deadline:
        time.sleep(0.02)
    assert adapter.process is None and processes[0].poll() is not None
    monkeypatch.setattr("kannotate.search.resident.MAX_RSS", 1)
    second = tmp_path / "second"
    second.mkdir()
    adapter.search(query, second, 3)
    assert adapter.disabled and adapter.process is None
    assert adapter.fallback.calls == 1
    assert adapter.last_backend == "canonical-cli" and adapter.warnings
    assert all(p.poll() is not None for p in processes)


def test_stale_idle_callback_cannot_close_new_lease(tmp_path, monkeypatch):
    fixture_worker(monkeypatch)
    adapter = ResidentHmmer(Fallback(), {}, minimum_queries=1)
    try:
        adapter._start()
        lease = adapter.lease
        adapter._start()
        adapter._expire(lease)
        assert adapter.process.poll() is None
    finally:
        adapter.close()


def test_unknown_platform_and_changed_backend_options_stay_on_cli(monkeypatch):
    import kannotate.search.resident as resident

    monkeypatch.setattr(resident.sys, "platform", "linux")
    monkeypatch.setattr(resident, "package_identity", lambda: pytest.fail("must not import models"))
    assert resident.qualified_identity(Fallback(), {}) is None


def test_worker_exit_is_fallback_and_never_an_empty_resident_success(tmp_path, monkeypatch):
    original = subprocess.Popen

    def popen(command, **kwargs):
        return original([sys.executable, "-c", "pass"], **kwargs)

    monkeypatch.setattr("kannotate.search.resident.subprocess.Popen", popen)
    adapter = ResidentHmmer(Fallback(), {}, minimum_queries=1)
    query = tmp_path / "query.faa"
    query.write_text(">q\nMACDEFG\n")
    engine = AnnotationEngine(CachedAdapter(adapter, tmp_path / "cache", {}), Fallback())
    result = engine.annotate(parse_queries("MACDEFG"), threads=3)
    assert result.warnings == ["Resident Pfam unavailable; used the canonical CLI."]
    assert adapter.disabled and adapter.fallback.calls == 1
    entry = json.loads(next((tmp_path / "cache").rglob("*.json")).read_text())
    if "execution_backend" not in entry:
        entry = json.loads(next((tmp_path / "cache").glob("*/*.json")).read_text())
    assert entry["execution_backend"] == "canonical-cli"
    assert adapter.process is None


def test_measured_mac_policy_keeps_small_misses_on_cli_and_resides_for_twenty(
    tmp_path, monkeypatch
):
    processes = fixture_worker(monkeypatch)
    adapter = ResidentHmmer(Fallback(), {})
    query = tmp_path / "queries.faa"
    try:
        for count in [1, 5, 20, 5, 20]:
            query.write_text("".join(f">q{i}\nMACDEFG\n" for i in range(count)))
            workdir = tmp_path / f"call-{count}-{adapter.fallback.calls}"
            workdir.mkdir()
            adapter.search(query, workdir, 6)
            assert adapter.last_backend == ("resident" if count == 20 else "canonical-cli")
            assert len(processes) == (1 if count == 20 or adapter.fallback.calls == 3 else 0)
        assert adapter.fallback.calls == 3 and len(processes) == 1
    finally:
        adapter.close()
