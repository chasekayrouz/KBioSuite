"""Suite defaults, launch errors and installed-acceptance hooks, without engines or a window."""

import json
import re
import sys
import threading
import time
import types
from dataclasses import asdict, replace
from pathlib import Path

import gemmi
import numpy as np
import pytest
from test_alignment import make_structure

from kbiosuite.version import __version__
from kmatch import alignment, cli, gui
from kmatch.config import Config, default_cache, default_home
from kmatch.models import PredictionResult, Progress, SearchResult, write_json
from kmatch.pipeline import Pipeline, RunResult
from kmatch.process import Cancellation, Cancelled, environment, run_process
from kmatch.sequence import parse_sequence

FIXTURE = Path(__file__).with_name("fixtures") / "1csp-A.fasta"
SEQUENCE = "ACDEFGKACDEF"


def roots(monkeypatch, tmp_path):
    monkeypatch.setenv("KMATCH_HOME", str(tmp_path / "home"))
    monkeypatch.setenv("KMATCH_CACHE", str(tmp_path / "cache"))
    return tmp_path / "home", tmp_path / "cache"


def test_suite_default_home_and_cache_follow_overrides(tmp_path, monkeypatch):
    monkeypatch.delenv("KMATCH_HOME", raising=False)
    monkeypatch.delenv("KMATCH_CACHE", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    assert default_home() == tmp_path / "Library/Application Support/KBioSuite/Kmatch"
    assert default_cache() == tmp_path / "Library/Caches/KBioSuite/Kmatch"
    monkeypatch.setenv("KMATCH_HOME", "~/custom/home")
    monkeypatch.setenv("KMATCH_CACHE", str(tmp_path / "x/../cache"))
    assert default_home() == (tmp_path / "custom/home").resolve()
    assert default_cache() == (tmp_path / "cache").resolve()
    config = Config.load()
    assert (config.home, config.cache) == (default_home(), default_cache())
    assert config.runs == config.home / "runs" and config.pdb_cache == config.cache / "pdb"
    # Loading reads configuration only; it never creates the home or cache.
    assert not (tmp_path / "custom").exists() and not (tmp_path / "cache").exists()


def test_fixture_is_the_public_1csp_chain_a_sequence():
    protein = parse_sequence(FIXTURE.read_text())
    # Copied from the fetched public RCSB 1CSP entity 1 benchmark input (chain A).
    assert protein.identifier == "1CSP_A" and protein.length == 67
    assert protein.sequence.startswith("MLEGKVKWFNSEKGFGFIEV")


def test_caller_selected_runs_root_isolates_pipeline_runs(tmp_path):
    from test_kmatch_core import FakePrediction, FakeSearch

    config = Config(tmp_path / "home", tmp_path / "cache")
    isolated = replace(config, runs_root=tmp_path / "acceptance/runs")
    assert isolated.runs == tmp_path / "acceptance/runs" and config.runs == config.home / "runs"
    result = Pipeline(isolated, FakePrediction(), FakeSearch()).run(SEQUENCE + "GHIK")
    assert result.directory.parent == tmp_path / "acceptance/runs"
    assert not (tmp_path / "home/runs").exists()


def test_external_interpreters_never_write_bytecode(tmp_path):
    assert environment()["PYTHONDONTWRITEBYTECODE"] == "1"
    assert environment({"OMP_NUM_THREADS": "2"})["PYTHONDONTWRITEBYTECODE"] == "1"
    log = tmp_path / "child.log"
    run_process(
        [sys.executable, "-c", "import sys;print('dont_write', sys.dont_write_bytecode)"],
        log,
        Cancellation(),
        lambda event: None,
        "probe",
    )
    assert "dont_write True" in log.read_text()


@pytest.mark.parametrize("module", [cli, gui])
def test_version_flag(module, capsys):
    with pytest.raises(SystemExit) as exit:
        module.main(["--version"])
    assert exit.value.code == 0 and capsys.readouterr().out.strip() == __version__


@pytest.mark.parametrize(
    "argv, message",
    [
        (["--smoke-input", "x.fasta"], "--smoke-input requires --smoke-dir"),
        (["--smoke-dir", "out"], "--smoke-dir requires --smoke-input"),
        (["--browser", "--bridge-check-dir", "out"], "omit --browser"),
        (["--browser", "--smoke-dir", "out", "--smoke-input", "x"], "omit --browser"),
    ],
)
def test_acceptance_arguments_are_validated(argv, message, capsys):
    with pytest.raises(SystemExit) as exit:
        gui.parse_args(argv)
    assert exit.value.code == 2 and message in capsys.readouterr().err


def test_acceptance_paths_are_resolved(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _, args = gui.parse_args(["--smoke-dir", "out", "--smoke-input", "in.fasta"])
    assert args.smoke_dir == tmp_path.resolve() / "out"
    assert args.smoke_input == tmp_path.resolve() / "in.fasta"
    _, args = gui.parse_args(["--bridge-check-dir", "~/bridge"])
    assert args.bridge_check_dir == Path.home() / "bridge"


class Event:
    def __init__(self):
        self.handlers = []

    def __iadd__(self, handler):
        self.handlers.append(handler)
        return self


def fake_webview(monkeypatch):
    """Stand-in pywebview: start() runs the startup callback, then closes the window."""
    module = types.ModuleType("webview")
    module.windows = []

    def create_window(title, url=None, **options):
        window = types.SimpleNamespace(
            title=title, url=url, options=options, events=types.SimpleNamespace(closed=Event())
        )
        module.windows.append(window)
        return window

    def start(func=None, debug=False):
        if func:
            func()
        for window in module.windows:
            for handler in window.events.closed.handlers:
                handler()

    module.create_window, module.start = create_window, start
    monkeypatch.setitem(sys.modules, "webview", module)
    return module


def broken_config(monkeypatch, tmp_path):
    home, _ = roots(monkeypatch, tmp_path)
    home.mkdir(parents=True)
    (home / "config.json").write_text(json.dumps({"<b>unknown</b>": 1}))
    return home


def test_cli_and_browser_startup_errors_exit_nonzero_with_message(tmp_path, monkeypatch, capsys):
    broken_config(monkeypatch, tmp_path)
    for main, argv in ((cli.main, ["status"]), (gui.main, ["--browser", "--no-open"])):
        with pytest.raises(SystemExit) as exit:
            main(argv)
        assert exit.value.code == 1
        assert "Kmatch could not load its configuration: Unknown Kmatch configuration" in (
            capsys.readouterr().err
        )


def test_desktop_startup_error_is_shown_and_recorded_for_checks(tmp_path, monkeypatch):
    home = broken_config(monkeypatch, tmp_path)
    webview = fake_webview(monkeypatch)
    with pytest.raises(SystemExit) as exit:
        gui.main([])
    assert exit.value.code == 1
    window = webview.windows[0]
    assert window.title == f"Kmatch · {__version__}" and window.url is None
    page = window.options["html"]
    assert "Kmatch could not start" in page and "&lt;b&gt;unknown&lt;/b&gt;" in page
    assert str(home / "config.json") in page and "<b>" not in page
    smoke_dir, bridge_dir = tmp_path / "smoke", tmp_path / "bridge"
    with pytest.raises(SystemExit):
        gui.main(["--bridge-check-dir", str(bridge_dir)])
    with pytest.raises(SystemExit):
        gui.main(["--smoke-dir", str(smoke_dir), "--smoke-input", str(FIXTURE)])
    for path in (bridge_dir / "bridge.json", smoke_dir / "bridge.json", smoke_dir / "report.json"):
        report = json.loads(path.read_text())
        assert report["ok"] is False and report["version"] == __version__
        assert report["home"] == str(home) and "<b>unknown</b>" in report["error"]
    assert not (home / "runs").exists() and not (smoke_dir / "runs").exists()


def test_desktop_main_titles_window_and_isolates_acceptance_runs(tmp_path, monkeypatch):
    home, _ = roots(monkeypatch, tmp_path)
    webview = fake_webview(monkeypatch)
    calls = []
    monkeypatch.setattr(gui, "bridge_contract", lambda api, path: calls.append(("bridge", path)))
    monkeypatch.setattr(
        gui, "smoke", lambda api, path, fasta, closed: calls.append((api._job.config, path, fasta))
    )
    gui.main(["--smoke-dir", str(tmp_path / "out"), "--smoke-input", str(FIXTURE)])
    window = webview.windows[0]
    assert window.title == f"Kmatch · {__version__}" and window.url.startswith("http://127.0.0.1:")
    assert isinstance(window.options["js_api"], gui.DesktopApi)
    config, path, fasta = calls[0]
    assert path == (tmp_path / "out").resolve() and fasta == FIXTURE.resolve()
    assert config.runs == path / "runs" and config.home == home.resolve()
    calls.clear()
    gui.main(["--bridge-check-dir", str(tmp_path / "bridge")])
    assert calls == [("bridge", (tmp_path / "bridge").resolve())]


def test_js_api_exposes_exactly_the_bridge_methods():
    public = {name for name in dir(gui.DesktopApi) if not name.startswith("_")}
    assert public == set(gui.BRIDGE_METHODS)
    bridge = (Path(gui.__file__).with_name("static") / "bridge.js").read_text()
    assert all(f"'{name}'" in bridge for name in gui.BRIDGE_METHODS)
    page = (Path(gui.__file__).with_name("static") / "index.html").read_text()
    assert page.index("/static/bridge.js") < page.index("/static/app.js")


READY = {
    "event_seen": True,
    "ready": True,
    "desktop": True,
    "methods": dict.fromkeys(gui.BRIDGE_METHODS, True),
    "error": None,
}


class BridgeWindow:
    def __init__(self, *states):
        self.states = list(states)

    def run_js(self, script):
        assert script == "JSON.stringify(window.kmatchBridge || null)"
        state = self.states.pop(0) if len(self.states) > 1 else self.states[0]
        if isinstance(state, Exception):
            raise state
        return json.dumps(state)


def bridge_api(tmp_path, window):
    api = gui.DesktopApi(gui.Job(Config(tmp_path / "home", tmp_path / "cache")))
    api._window = window
    return api


def test_bridge_check_waits_for_native_readiness_and_writes_contract(tmp_path):
    api = bridge_api(tmp_path, BridgeWindow(None, {**READY, "ready": False}, READY))
    report = gui.bridge_contract(api, tmp_path / "out", timeout=5, interval=0)
    assert report == json.loads((tmp_path / "out/bridge.json").read_text())
    assert report["ok"] is True and "error" not in report
    assert report["version"] == __version__ and report["home"] == str(tmp_path / "home")
    assert report["bridge"] == READY
    # Readiness reuses status checks; AF3 is deferred and never reported or required.
    assert report["dependencies_ready"] == {"boltz2": False, "foldseek": False, "pdb100": False}
    assert not (tmp_path / "home").exists() and not (tmp_path / "cache").exists()


@pytest.mark.parametrize(
    "state, error",
    [
        (
            {**READY, "ready": False, "error": "Native controls could not initialize."},
            "could not initialize",
        ),
        ({**READY, "methods": {**READY["methods"], "open_pdb": False}}, "did not become ready"),
        ({**READY, "desktop": False}, "did not become ready"),
        (None, "did not become ready"),
        (RuntimeError("window closed"), "Native bridge check failed: window closed"),
    ],
)
def test_bridge_check_failures_are_recorded(tmp_path, state, error):
    api = bridge_api(tmp_path, BridgeWindow(state))
    report = gui.bridge_contract(api, tmp_path, timeout=0.05, interval=0.01)
    assert report["ok"] is False and error in report["error"]
    assert json.loads((tmp_path / "bridge.json").read_text())["ok"] is False


def acceptance_run(tmp_path):
    """A completed run with two hits; the second has the higher query TM (top-ranked)."""
    xyz = np.array([[3.8 * i, (i % 3) * 1.5, (i % 2) * 2.0] for i in range(len(SEQUENCE))])
    query, target = make_structure(SEQUENCE, xyz, full=True), make_structure(SEQUENCE, xyz)
    template = tmp_path / "template"
    (template / "prediction").mkdir(parents=True)
    (template / "search").mkdir()
    source = template / "prediction/predicted.pdb"
    query.write_pdb(str(source))
    hits = [
        SearchResult(
            query="query",
            target=name,
            evalue=evalue,
            probability=probability,
            bits=80,
            sequence_identity=1,
            alignment_length=len(SEQUENCE),
            query_length=len(SEQUENCE),
            target_length=len(SEQUENCE),
            query_coverage=1,
            target_coverage=1,
            query_tm=query_tm,
            target_tm=query_tm,
            alignment_tm=query_tm,
            lddt=0.8,
            query_start=1,
            query_end=len(SEQUENCE),
            target_start=1,
            target_end=len(SEQUENCE),
            query_alignment=SEQUENCE,
            target_alignment=SEQUENCE,
            rotation=np.eye(3).tolist(),
            translation=[0.0, 0.0, 0.0],
        )
        for name, query_tm, evalue, probability in (
            ("2xyz_B", 0.5, 1e-3, 0.7),
            ("1abc-assembly1_A", 0.9, 1e-12, 1.0),
        )
    ]
    write_json(
        template / "search/search.json",
        {
            "database": str(tmp_path / "db/pdb"),
            "database_version": "fixture",
            "hits": [asdict(hit) for hit in hits],
            "query_structure_sha256": alignment.digest(source),
        },
    )
    return template, target, hits


class AcceptancePipeline:
    """Completes like Pipeline.run, writing only under the configured runs root."""

    template = None
    hits = None
    block = False

    def __init__(self, config):
        self.config = config

    def run(self, text, *, cancel, progress, network_msa, msa_path, on_directory, on_prediction):
        import shutil

        assert network_msa is False and msa_path is None
        sequence = parse_sequence(text)
        directory = self.config.runs / "20261116-120000-fixture"
        shutil.copytree(self.template, directory)
        on_directory(directory)
        write_json(
            directory / "run.json",
            {"sequence_identifier": sequence.identifier, "msa_mode": "single-sequence"},
        )
        while self.block:
            cancel.check()
            time.sleep(0.01)
        prediction = PredictionResult(
            directory / "prediction/predicted.pdb",
            "boltz2",
            "fixture",
            sequence.length,
            {"mean_plddt": 88.5, "ptm": 0.81, "confidence_score": 0.84},
            [88.5] * sequence.length,
            2.0,
            {"msa": "single-sequence", "network_during_inference": False},
            warnings=[{"code": "fixture", "message": "fixture warning"}],
        )
        write_json(directory / "prediction/prediction.json", prediction)
        on_prediction(prediction)
        progress(Progress("complete", "Complete · fixture"))
        return RunResult(directory, prediction, self.hits)


class Page:
    """The page steps the acceptance driver performs, answered from the real Job state."""

    def __init__(self, job, closed=None):
        self.job, self.closed, self.steps, self.sequence = job, closed, [], None

    def run_js(self, script):
        if script == "JSON.stringify(window.kmatchBridge || null)":
            return json.dumps(READY)
        step = re.match(r"/\*kmatch-smoke:([a-z-]+)\*/", script)[1]
        self.steps.append(step)
        state = self.job.snapshot()
        if step == "configure":
            self.sequence = json.loads(re.search(r"sequence\.value=(\".*?\");", script)[1])
            value = dict(
                backend="boltz2",
                msa_mode="single",
                msa_label="Local only · no MSA",
                pdb_network=False,
                full_atom=False,
            )
        elif step == "validated":
            value = {"text": "Valid single protein sequence", "error": False}
        elif step == "run":
            # What the Run button posts to /api/run with "Local only · no MSA".
            self.job.run(self.sequence, "boltz2", False, None)
            if self.closed:
                self.closed.set()
            value = True
        elif step == "rendered":
            value = dict(
                rows=len(state["hit_rankings"]["qtm"]),
                results_visible=True,
                viewer_empty_hidden=True,
                canvas=True,
                atoms=36,
                loaded_run=state["run_id"],
                error="",
                summary="",
                pdb_status="PDB details: 0/2 hits available.",
            )
        elif step == "select-top-hit":
            index = state["hit_rankings"]["qtm"][0]
            target = state["hits"][index]["target"]
            value = dict(
                hit_index=index,
                rank="1",
                target=target,
                selected=True,
                detail_visible=True,
                title=f"#1 · {target}",
                metrics=12,
                rcsb="https://www.rcsb.org/structure/" + target[:4].upper(),
                ranking="qtm",
                export_enabled=True,
                save_enabled=True,
            )
        elif step == "alignment-result":
            value = "1abc · CIF · Saved to " + state["alignment_export"]["saved_to"]
        else:
            value = {"ready": True, "page-error": "", "final-error": ""}[step]
        return json.dumps({"value": value})


@pytest.fixture
def acceptance(tmp_path, monkeypatch):
    template, target, hits = acceptance_run(tmp_path)
    monkeypatch.setattr(AcceptancePipeline, "template", template)
    monkeypatch.setattr(AcceptancePipeline, "hits", hits)
    monkeypatch.setattr(gui, "Pipeline", AcceptancePipeline)
    monkeypatch.setattr(alignment, "materialize_target", lambda *a: (target.clone(), "7"))
    monkeypatch.setattr(
        gui,
        "dependency_readiness",
        lambda config: dict.fromkeys(("boltz2", "foldseek", "pdb100"), True),
    )
    config = Config(tmp_path / "home", tmp_path / "cache", pymol_override=str(tmp_path / "none"))
    out = tmp_path / "acceptance"
    fasta = tmp_path / "public.fasta"
    fasta.write_text(">public_fixture\n" + SEQUENCE + "\n")
    job = gui.Job(replace(config, runs_root=out / "runs"))
    api = gui.DesktopApi(job)
    api._window = Page(job)
    return api, out, fasta


def test_smoke_drives_one_isolated_job_and_exports_without_dialogs(acceptance, tmp_path):
    api, out, fasta = acceptance
    report = gui.smoke(api, out, fasta, threading.Event())
    assert report == json.loads((out / "report.json").read_text())
    assert report["ok"] is True and "error" not in report, report
    assert report["checks"] and all(report["checks"].values())
    assert json.loads((out / "bridge.json").read_text())["ok"] is True
    assert api._window.steps[:4] == ["ready", "configure", "validated", "run"]
    run = Path(report["run_directory"])
    assert run.parent == out / "runs" and report["run_id"] == run.name
    assert not (tmp_path / "home/runs").exists()
    assert report["prediction"]["mean_plddt"] == 88.5 and report["prediction"]["ptm"] == 0.81
    assert report["warnings"] == [{"code": "fixture", "message": "fixture warning"}]
    assert report["hits"] == 2 and report["returned_hits"] == 2
    assert report["top_hit"]["target"] == "1abc-assembly1_A" and report["top_hit"]["hit_index"] == 1
    assert report["top_hit"]["rcsb_url"] == "https://www.rcsb.org/structure/1ABC"
    assert set(report["exports"]) == {"predicted_structure", "superposition_cif"}
    for path in report["exports"].values():
        path = Path(path)
        assert (
            path.parent == out / "exports" and gemmi.read_structure(str(path))[0].count_atom_sites()
        )
    superposition = gemmi.read_structure(report["exports"]["superposition_cif"])
    assert [chain.name for chain in superposition[0]] == ["Q", "T"]
    assert report["pymol_session"]["status"] == "skipped"
    # Coordinate-only export: Foldseek's saved transform, no full-atom (network) target.
    provenance = json.loads(next((run / "exports").glob("*/provenance.json")).read_text())
    assert "full_atom" not in provenance and provenance["refitted"] is False


def test_smoke_stops_before_running_when_dependencies_are_missing(acceptance, monkeypatch):
    api, out, fasta = acceptance
    monkeypatch.setattr(
        gui,
        "dependency_readiness",
        lambda config: {"boltz2": True, "foldseek": True, "pdb100": False},
    )
    report = gui.smoke(api, out, fasta)
    assert report["ok"] is False and report["error"].endswith("not ready: pdb100")
    assert api._window.steps == [] and not (out / "runs").exists()
    assert json.loads((out / "bridge.json").read_text())["dependencies_ready"]["pdb100"] is False


def test_smoke_records_bridge_failure_without_touching_the_page(acceptance):
    api, out, fasta = acceptance
    api._window = BridgeWindow({**READY, "ready": False, "error": "Native controls failed."})
    report = gui.smoke(api, out, fasta)
    assert report["ok"] is False and "Native controls failed." in report["error"]
    assert json.loads((out / "report.json").read_text())["ok"] is False
    assert not (out / "runs").exists()


def test_closing_the_window_cancels_the_acceptance_job(acceptance, monkeypatch):
    api, out, fasta = acceptance
    monkeypatch.setattr(AcceptancePipeline, "block", True)
    closed = threading.Event()
    api._window = Page(api._job, closed)
    report = gui.smoke(api, out, fasta, closed)
    assert report["ok"] is False and "window closed" in report["error"]
    api._job.thread.join(5)
    state = api._job.snapshot()
    assert not state["busy"] and state["phase"] == "cancelled"


def test_dependencies_neither_announce_nor_read_pre_suite_kmatch_data(tmp_path, monkeypatch):
    """The 0.5 standalone-data notice is retired: older data is not mentioned or touched."""
    user = tmp_path / "user"
    monkeypatch.setenv("HOME", str(user))
    monkeypatch.delenv("KMATCH_HOME", raising=False)
    legacy = user / "Library/Application Support/Kmatch"
    (legacy / "runs").mkdir(parents=True)
    (legacy / "config.json").write_text('{"threads": 8}')
    before = sorted((str(p), p.lstat().st_mtime_ns) for p in legacy.rglob("*"))
    assert default_home() == user / "Library/Application Support/KBioSuite/Kmatch"
    client = gui.create_app(gui.Job(Config(default_home(), tmp_path / "cache"))).test_client()
    page = client.get("/").text
    token = re.search(r'name="kmatch-token" content="([^"]+)"', page)[1]
    body = client.get("/api/dependencies", headers={"X-Kmatch-Token": token}).get_json()
    assert set(body) == {"items"} and len(body["items"]) == 4
    assert "standalone" not in page.lower() and "adopt" not in page.lower()
    assert sorted((str(p), p.lstat().st_mtime_ns) for p in legacy.rglob("*")) == before
    assert not default_home().exists()


class ConfiguredPage(Page):
    """The page reports a control the smoke could not set to its local-only value."""

    override = {}

    def run_js(self, script):
        raw = super().run_js(script)
        if "kmatch-smoke:configure" in script:
            return json.dumps({"value": {**json.loads(raw)["value"], **self.override}})
        return raw


@pytest.mark.parametrize(
    "override, unmet",
    [
        ({"backend": "af3"}, "boltz2_backend"),
        ({"msa_mode": "network"}, "local_only_no_msa"),
        ({"msa_label": "Network-assisted MSA · recommended for accuracy"}, "local_only_no_msa"),
        ({"pdb_network": True}, "pdb_network_disabled"),
        ({"full_atom": True}, "full_atom_disabled"),
    ],
)
def test_smoke_never_presses_run_unless_the_page_is_local_only(acceptance, override, unmet):
    api, out, fasta = acceptance
    api._window = ConfiguredPage(api._job)
    api._window.override = override
    report = gui.smoke(api, out, fasta)
    assert report["ok"] is False and report["checks"][unmet] is False
    assert unmet in report["error"] and "Run was not pressed" in report["error"]
    assert api._window.steps == ["ready", "configure"]
    assert not (out / "runs").exists() and api._job.thread is None


def test_reused_directories_never_show_a_stale_pass(acceptance):
    api, out, fasta = acceptance
    out.mkdir(parents=True)
    for name in ("bridge.json", "report.json"):
        (out / name).write_text('{"ok": true}\n')
    seen = []

    class Watching(Page):
        def run_js(self, script):
            seen.append(sorted(p.name for p in out.glob("*.json")))
            return super().run_js(script)

    api._window = Watching(api._job)
    report = gui.smoke(api, out, fasta)
    # Both stale files were gone before the first page call of the new check.
    assert seen[0] == [] and report["ok"] is True
    (out / "bridge.json").write_text('{"ok": true}\n')
    api._window = BridgeWindow(None)
    gui.bridge_contract(api, out, timeout=0.05, interval=0.01)
    assert json.loads((out / "bridge.json").read_text())["ok"] is False


class EndsBeforeDirectory(AcceptancePipeline):
    error = Cancelled("Run cancelled. Completed outputs and logs are retained.")

    def run(self, text, **options):
        raise self.error


@pytest.mark.parametrize(
    "error, phase",
    [
        (Cancelled("Run cancelled. Completed outputs and logs are retained."), "cancelled"),
        (RuntimeError("Boltz-2 could not start"), "failed"),
    ],
)
def test_a_run_cancelled_or_failed_before_its_folder_ends_the_smoke_promptly(
    acceptance, monkeypatch, error, phase
):
    api, out, fasta = acceptance
    monkeypatch.setattr(EndsBeforeDirectory, "error", error)
    monkeypatch.setattr(gui, "Pipeline", EndsBeforeDirectory)
    began = time.monotonic()
    report = gui.smoke(api, out, fasta)
    assert time.monotonic() - began < 10  # Not the 60 s launch wait.
    assert report["ok"] is False and report["error"].startswith("The run did not complete")
    assert report["final_job"] == {"busy": False, "phase": phase} and report["job_idle"] is True
    assert "rendered" not in api._window.steps


def test_overall_deadline_cancels_the_job_before_the_report(acceptance, monkeypatch):
    api, out, fasta = acceptance
    monkeypatch.setattr(AcceptancePipeline, "block", True)
    monkeypatch.setattr(gui, "SMOKE_DEADLINE_SECONDS", 1)
    began = time.monotonic()
    report = gui.smoke(api, out, fasta, threading.Event())
    assert time.monotonic() - began < 10  # The cancelled fixture job stops promptly.
    assert report["ok"] is False and report["deadline_seconds"] == 1
    assert report["error"] == "The 1 s smoke deadline expired while waiting for the run."
    # Kmatch is idle when report.json is written: the job was cancelled, not abandoned.
    assert report["job_idle"] is True and report["final_job"] == {
        "busy": False,
        "phase": "cancelled",
    }
    assert json.loads((out / "report.json").read_text())["final_job"]["busy"] is False


def test_optional_pymol_timeout_is_cancelled_and_recorded_without_failing(acceptance, monkeypatch):
    api, out, fasta = acceptance
    real = gui.export_alignment

    def export(*args, **options):
        if args[3] == "pse":  # A PyMOL session that never finishes.
            while True:
                options["cancel"].check()
                time.sleep(0.01)
        return real(*args, **options)

    monkeypatch.setattr(gui, "export_alignment", export)
    monkeypatch.setattr(gui, "export_support", lambda config: {"pse": True, "message": ""})
    monkeypatch.setattr(gui, "SMOKE_EXPORT_TIMEOUT_SECONDS", 2)
    report = gui.smoke(api, out, fasta)
    assert report["ok"] is True and "error" not in report, report
    assert report["pymol_session"] == {
        "status": "failed",
        "error": "Timed out after 2 s waiting for the 3D export.",
        "cancelled": True,
        "idle": True,
    }
    assert "pymol_session" not in report["exports"] and "pymol_session" not in report["checks"]
    assert report["final_job"] == {"busy": False, "phase": "cancelled"}


def test_optional_pymol_export_is_skipped_when_the_deadline_is_spent(acceptance, monkeypatch):
    api, out, fasta = acceptance
    monkeypatch.setattr(gui, "export_support", lambda config: {"pse": True, "message": ""})
    clock = types.SimpleNamespace(offset=0, sleep=time.sleep)
    clock.monotonic = lambda: time.monotonic() + clock.offset
    monkeypatch.setattr(gui, "time", clock)
    real = gui._structure_residues
    calls = []

    def spend_deadline(path, chain=None):
        calls.append(path)
        if len(calls) == 2:  # Just after the superposition export, no time remains.
            clock.offset = gui.SMOKE_DEADLINE_SECONDS + 1
        return real(path, chain)

    monkeypatch.setattr(gui, "_structure_residues", spend_deadline)
    report = gui.smoke(api, out, fasta)
    assert report["ok"] is True and "error" not in report, report
    assert report["pymol_session"] == {
        "status": "skipped",
        "reason": "The smoke deadline left no time for the optional session export.",
    }
    assert report["final_job"] == {"busy": False, "phase": "complete"}
