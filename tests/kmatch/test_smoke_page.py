"""Offline consistency of the installed smoke's page steps and of the setup page with app.js."""

import inspect
import json
import re
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from kmatch import cli, gui
from kmatch import setup as kmatch_setup
from kmatch.config import Config

STATIC = Path(gui.__file__).with_name("static")
INDEX = (STATIC / "index.html").read_text()
APP = (STATIC / "app.js").read_text()
BRIDGE = (STATIC / "bridge.js").read_text()
STATIC_IDS = set(re.findall(r'\sid="([^"]+)"', INDEX))
HARNESS = Path(__file__).with_name("app_page_harness.cjs")
NODE = shutil.which("node")
needs_node = pytest.mark.skipif(NODE is None, reason="node is not installed")

# Selectors the smoke reads that do not exist in index.html until app.js (or 3Dmol) creates
# them, each with the app.js source that creates it. A new dynamic selector must be added here.
DYNAMIC_SELECTORS = {
    "#hits tr": "const row = document.createElement('tr'); row.dataset.hitIndex = index;",
    "td:nth-child(2) button": "if (column === 1) {const button=document.createElement('button');",
    "a.pdb-link": "link.className='pdb-link'",
    "canvas": "viewer = $3Dmol.createViewer($('viewer')",
}


def snippet(step):
    body = gui.SMOKE_STEPS[step]
    return (
        body.replace("TEXT", json.dumps(">public\nMLEGKVKWFN\n")) if step == "configure" else body
    )


def test_smoke_steps_are_the_complete_set_of_page_calls():
    calls = re.findall(r'_page\(\s*window,\s*"([\w-]+)"\s*([,)])', inspect.getsource(gui.smoke))
    assert {step for step, _ in calls} == set(gui.SMOKE_STEPS)
    # Only configure substitutes its sequence; every other step evaluates SMOKE_STEPS verbatim.
    assert {step for step, end in calls if end == ","} == {"configure"}


def test_every_smoke_dom_reference_exists_in_the_page():
    ids, selectors = set(), set()
    for step in gui.SMOKE_STEPS:
        body = gui.SMOKE_STEPS[step]
        ids.update(re.findall(r"(?:\bel|getElementById)\('([^']+)'\)", body))
        selectors.update(re.findall(r"querySelector(?:All)?\('([^']+)'\)", body))
    assert ids and ids <= STATIC_IDS, ids - STATIC_IDS
    assert selectors == set(DYNAMIC_SELECTORS), selectors ^ set(DYNAMIC_SELECTORS)
    for selector, source in DYNAMIC_SELECTORS.items():
        assert source in APP, selector
        assert set(re.findall(r"#([\w-]+)", selector)) <= STATIC_IDS
    assert "$('hits').append(row)" in APP and '<tbody id="hits">' in INDEX
    # Page globals the steps read: app.js top-level bindings and bridge.js's contract.
    assert re.search(r"^let .*\bviewer = null\b.*\bloadedRun = null\b", APP, re.M)
    assert "window.kmatchBridge = {" in BRIDGE
    # Text the smoke compares against is what the page actually shows.
    assert '<option value="single">Local only · no MSA</option>' in INDEX
    assert '<option value="boltz2">' in INDEX
    assert "'Valid single protein sequence'" in APP and "` · Saved to ${result.saved_to}`" in APP
    assert "`#${hitOrder.indexOf(index)+1} · ${hit.target}`" in APP


@needs_node
@pytest.mark.parametrize("step", list(gui.SMOKE_STEPS))
def test_each_wrapped_smoke_step_is_valid_javascript(tmp_path, step):
    script = tmp_path / f"{step}.js"
    script.write_text(gui._page_script(step, snippet(step)))
    completed = subprocess.run(
        [NODE, "--check", str(script)], capture_output=True, text=True, timeout=60
    )
    assert completed.returncode == 0, completed.stderr


def page(tmp_path, dependencies, job=None):
    scenario = tmp_path / "scenario.json"
    scenario.write_text(
        json.dumps({"ids": sorted(STATIC_IDS), "dependencies": dependencies, "job": job or {}})
    )
    completed = subprocess.run(
        [NODE, str(HARNESS), str(STATIC / "app.js"), str(scenario)],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0, completed.stderr
    return json.loads(completed.stdout)


def setup_client(tmp_path, monkeypatch):
    """The real app for a default suite home; pre-suite Kmatch data is present but irrelevant."""
    user = tmp_path / "user"
    monkeypatch.setenv("HOME", str(user))
    monkeypatch.delenv("KMATCH_HOME", raising=False)
    (user / "Library/Application Support/Kmatch/runs").mkdir(parents=True)
    home = user / "Library/Application Support/KBioSuite/Kmatch"
    job = gui.Job(Config(home, tmp_path / "cache"))
    client = gui.create_app(job).test_client()
    token = re.search(r'name="kmatch-token" content="([^"]+)"', client.get("/").text)[1]
    return job, client, {"X-Kmatch-Token": token}


@needs_node
def test_setup_page_starts_setup_without_a_prompt(tmp_path, monkeypatch):
    _, client, headers = setup_client(tmp_path, monkeypatch)
    payload = client.get("/api/dependencies", headers=headers).get_json()
    assert set(payload) == {"items"}
    shown = page(tmp_path, payload)
    assert shown["summary"] == "Setup required" and shown["panel_open"] is True
    assert shown["confirms"] == [] and shown["setup"] == [{"component": "boltz"}]
    assert shown["error"] == ""
    assert shown["setup_error"] == ""


@needs_node
def test_setup_page_shows_the_command_line_tools_message(tmp_path, monkeypatch):
    """A Mac without Command Line Tools: Boltz setup fails before any download, visibly."""
    job, client, headers = setup_client(tmp_path, monkeypatch)
    ran = []

    def run(command, **kwargs):
        ran.append(command)
        return subprocess.CompletedProcess(
            command, 2, "", "xcode-select: error: no developer tools"
        )

    def offline(*args, **kwargs):
        pytest.fail("setup must not download or build without Command Line Tools")

    fake = SimpleNamespace(
        run=run,
        check_output=offline,
        DEVNULL=subprocess.DEVNULL,
        SubprocessError=subprocess.SubprocessError,
    )
    monkeypatch.setattr(kmatch_setup, "subprocess", fake)
    for name in ("download", "run_process"):
        monkeypatch.setattr(kmatch_setup, name, offline)
    for name in ("install_uv", "install_boltz"):
        monkeypatch.setattr(cli, name, offline)
    payload = client.get("/api/dependencies", headers=headers).get_json()
    assert (
        client.post("/api/setup", headers=headers, json={"component": "boltz"}).status_code == 200
    )
    job.thread.join(10)
    state = client.get("/api/status", headers=headers).get_json()["job"]
    message = (
        "Boltz setup needs Apple's Command Line Tools. Run `xcode-select --install` in Terminal, "
        "then try again. Nothing was downloaded."
    )
    assert state["error"] == message and state["phase"] == "failed" and not state["busy"]
    assert ran == [["/usr/bin/xcode-select", "-p"]]  # git/otool stubs were never run
    assert not job.config.home.exists()
    shown = page(tmp_path, payload, state)
    assert shown["error"] == message and shown["setup"] == [{"component": "boltz"}]
    assert shown["setup_error"] == message  # visible inside the setup panel too
