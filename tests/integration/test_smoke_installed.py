"""Installed-smoke harness units with fake running applications; no app is ever launched."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import psutil
import pytest

pytest.importorskip("AppKit")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))  # the harness imports its sibling scripts
try:
    spec = importlib.util.spec_from_file_location(
        "smoke_installed", ROOT / "scripts/smoke-installed.py"
    )
    smoke = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(smoke)
finally:
    sys.path.remove(str(ROOT / "scripts"))

KMATCH = smoke.APPS["Kmatch"]


class FakeApp:
    """The NSRunningApplication surface the harness uses."""

    def __init__(self, running, identity, path, *, pid=None, quits=True):
        self.running, self.identity, self.path = running, identity, path
        self.pid, self.quits = pid or 2_000_000_000, quits
        self.terminated = self.forced = 0
        running.append(self)

    def bundleIdentifier(self):
        return self.identity

    def bundleURL(self):
        return None if self.path is None else SimpleNamespace(path=lambda: str(self.path))

    def processIdentifier(self):
        return self.pid

    def isFinishedLaunching(self):
        return True

    def terminate(self):
        self.terminated += 1
        if self.quits and self in self.running:
            self.running.remove(self)
        return True

    def forceTerminate(self):
        self.forced += 1
        if self in self.running:
            self.running.remove(self)
        return True


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    running = []
    monkeypatch.setattr(smoke, "running_applications", lambda: list(running))
    monkeypatch.setattr(smoke, "NORMAL_QUIT_SECONDS", 0.3)
    monkeypatch.setattr(smoke, "FORCED_QUIT_SECONDS", 0.3)
    monkeypatch.setattr(smoke, "CHILD_EXIT_SECONDS", 0.3)
    applications = tmp_path / "Applications"
    directory = applications / "KBioSuite"
    for name in smoke.APPS:
        (directory / f"{name}.app").mkdir(parents=True)
    return SimpleNamespace(running=running, applications=applications, directory=directory)


def test_detection_counts_only_apps_inside_the_smoked_installation(workspace):
    running, applications, directory = (
        workspace.running,
        workspace.applications,
        workspace.directory,
    )
    suite = FakeApp(running, KMATCH, directory / "Kmatch.app")
    standalone = FakeApp(running, KMATCH, applications / "Kmatch.app")
    previous = FakeApp(running, KMATCH, applications / "Kmatch.previous-20261001.app")
    sibling = FakeApp(running, KMATCH, applications / "KBioSuite-old/Kmatch.app")
    unknown = FakeApp(running, smoke.APPS["Kstack"], None)
    FakeApp(running, "com.apple.finder", Path("/System/Library/CoreServices/Finder.app"))
    assert smoke.select(running, KMATCH, directory) == [suite]
    assert smoke.apps(KMATCH, directory) == [suite]
    assert smoke.foreign(running, directory) == [standalone, previous, sibling, unknown]
    assert smoke.ready(KMATCH, directory) == [suite]
    assert smoke.ready(smoke.APPS["Kblast"], directory) is None


def test_close_never_quits_a_standalone_app_with_the_same_identifier(workspace):
    suite = FakeApp(workspace.running, KMATCH, workspace.directory / "Kmatch.app")
    standalone = FakeApp(workspace.running, KMATCH, workspace.applications / "Kmatch.app")
    smoke.close(KMATCH, workspace.directory)
    assert suite.terminated == 1 and suite not in workspace.running
    assert standalone.terminated == standalone.forced == 0 and standalone in workspace.running


def test_quit_kmatch_is_detected_even_while_a_standalone_one_runs(workspace):
    FakeApp(workspace.running, KMATCH, workspace.applications / "Kmatch.app")
    with pytest.raises(RuntimeError, match="quit before writing report.json"):
        smoke.running_report(workspace.directory / "report.json", KMATCH, workspace.directory)


def test_refuses_to_start_while_a_non_suite_app_shares_a_suite_identifier(
    workspace, tmp_path, monkeypatch
):
    standalone = FakeApp(workspace.running, KMATCH, workspace.applications / "Kmatch.app")
    prepare(workspace, monkeypatch)
    monkeypatch.setattr(smoke, "launch", lambda *a: pytest.fail("nothing may be launched"))
    with pytest.raises(RuntimeError, match="non-suite app with a suite bundle identifier") as error:
        smoke.smoke(workspace.directory, tmp_path / "evidence")
    assert str(workspace.applications / "Kmatch.app") in str(error.value)
    assert "never quits or cancels others" in str(error.value)
    assert standalone.terminated == standalone.forced == 0


def prepare(workspace, monkeypatch):
    source = {"commit": "fixture", "dirty": False}
    (workspace.directory / ".kbiosuite-install.json").write_text(json.dumps({"source": source}))
    monkeypatch.setattr(smoke, "validate", lambda directory: None)
    monkeypatch.setattr(smoke, "check", lambda directory: None)
    monkeypatch.setattr(smoke, "source_identity", lambda: dict(source))


@pytest.mark.parametrize("failure", ["assertion", "timeout", "exception"])
def test_every_failure_path_quits_the_apps_the_smoke_launched(
    workspace, tmp_path, monkeypatch, failure
):
    prepare(workspace, monkeypatch)
    output = tmp_path / "evidence"
    real_wait = smoke.wait
    started = []

    def launch(directory, name, *args):
        started.append(FakeApp(workspace.running, smoke.APPS[name], directory / f"{name}.app"))
        if failure == "assertion":
            (output / "kblast").mkdir(parents=True)
            (output / "kblast/report.json").write_text(json.dumps({"ok": False}))
        elif failure == "exception":
            raise OSError("LaunchServices failure after the app started")

    def wait(check, timeout=150, what="the installed app"):
        if failure == "timeout" and what == "the Kblast report":
            raise RuntimeError(f"timed out waiting for {what}")
        return real_wait(check, timeout, what)

    monkeypatch.setattr(smoke, "launch", launch)
    monkeypatch.setattr(smoke, "wait", wait)
    with pytest.raises((AssertionError, RuntimeError, OSError)):
        smoke.smoke(workspace.directory, output)
    assert [app.terminated for app in started] == [1] and not workspace.running
    report = json.loads((output / "report.json").read_text())
    assert not report["ok"]
    assert report["cleanup"] == {
        "closed": ["Kblast"],
        "forced": [],
        "children_stopped": [],
        "errors": [],
    }


def test_cleanup_forces_a_stuck_app_and_stops_its_engine_but_not_others(workspace):
    # A stand-in "app" process whose child plays a Boltz engine that outlives the quit.
    app_process = subprocess.Popen(["/bin/sh", "-c", "sleep 60 & wait"])
    try:
        parent = psutil.Process(app_process.pid)
        engine = smoke_wait(lambda: parent.children(recursive=True))[0]
        stuck = FakeApp(
            workspace.running,
            KMATCH,
            workspace.directory / "Kmatch.app",
            pid=app_process.pid,
            quits=False,
        )
        standalone = FakeApp(workspace.running, KMATCH, workspace.applications / "Kmatch.app")
        outcome = smoke.stop_launched(workspace.directory, ["Kmatch", "Kblast", "Kmatch"])
        assert stuck.terminated == 1 and stuck.forced == 1
        assert outcome["forced"] == ["Kmatch"] and not outcome["errors"]
        assert outcome["children_stopped"] == [f"Kmatch child pid {engine.pid}"]
        assert not engine.is_running() or engine.status() == psutil.STATUS_ZOMBIE
        assert standalone.terminated == standalone.forced == 0
        assert workspace.running == [standalone]
    finally:
        app_process.kill()
        app_process.wait()


def test_cleanup_reports_an_app_that_cannot_be_quit(workspace, monkeypatch):
    stuck = FakeApp(workspace.running, KMATCH, workspace.directory / "Kmatch.app", quits=False)
    monkeypatch.setattr(stuck, "forceTerminate", lambda: True)
    outcome = smoke.stop_launched(workspace.directory, ["Kmatch"])
    assert outcome["errors"] == ["Kmatch is still running after a forced quit"]


def smoke_wait(check, timeout=5):
    return smoke.wait(check, timeout, "the test helper process")


def test_kmatch_report_wait_outlasts_the_apps_own_deadline():
    assert smoke.KMATCH_REPORT_TIMEOUT_SECONDS == 2700
    text = (ROOT / "scripts/smoke-installed.py").read_text()
    body = text[text.index('start("Kmatch", "--smoke-dir"') :]
    assert "KMATCH_REPORT_TIMEOUT_SECONDS" in body[: body.index("check_kmatch_job")]
    gui = importlib.import_module("kmatch.gui")
    if not hasattr(gui, "SMOKE_DEADLINE_SECONDS"):
        pytest.skip("Kmatch does not define SMOKE_DEADLINE_SECONDS in this checkout yet")
    assert smoke.KMATCH_REPORT_TIMEOUT_SECONDS > gui.SMOKE_DEADLINE_SECONDS == 2400
