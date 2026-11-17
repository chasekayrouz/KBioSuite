"""Small macOS desktop host. Scientific jobs belong to the independent pipeline."""

import argparse
import copy
import html
import io
import json
import secrets
import shutil
import subprocess
import sys
import threading
import time
import zipfile
from dataclasses import asdict, replace
from pathlib import Path

from flask import Flask, jsonify, request, send_file

from kbiosuite.desktop import loopback_server

from . import __version__
from .alignment import FORMATS, export_alignment, export_support, safe_name
from .cli import setup_component, status
from .config import Config, default_home
from .database import database_status
from .models import SearchResult
from .pdb_data import load_details, pdb_url
from .pipeline import Pipeline
from .predictors import BACKENDS, BoltzBackend
from .process import Cancellation, Cancelled
from .search import FoldseekBackend, hit_rankings, retained_hit_indices
from .sequence import parse_sequence

STATIC = Path(__file__).parent / "static"
# Native js_api methods the page requires; bridge.js checks the same names.
BRIDGE_METHODS = ("choose_msa", "export_structure", "reveal_run", "export_alignment", "open_pdb")
BRIDGE_TIMEOUT_SECONDS = 15
# One overall budget from smoke start; every wait is also bounded by the time remaining.
# The installed harness waits 2700 s for report.json: deadline + cancel grace stays below it.
SMOKE_DEADLINE_SECONDS = 2400
SMOKE_RUN_TIMEOUT_SECONDS = 25 * 60
SMOKE_EXPORT_TIMEOUT_SECONDS = 5 * 60
SMOKE_CANCEL_GRACE_SECONDS = 60
# After the window closes, the app's own shutdown waits only briefly for this thread.
SMOKE_CLOSED_GRACE_SECONDS = 5
AF3_EXPORT_NOTICE = (
    "AlphaFold 3 output. This exported copy is unmodified. "
    "Use and redistribution are subject to "
    "https://github.com/google-deepmind/alphafold3/blob/main/OUTPUT_TERMS_OF_USE.md\n"
)


class Job:
    def __init__(self, config):
        self.config, self.lock = config, threading.Lock()
        self.cancel_token, self.thread = Cancellation(), None
        self.details_cancel, self.details_thread = Cancellation(), None
        self.state = dict(
            busy=False,
            phase="ready",
            message="Ready for one protein sequence",
            error=None,
            prediction=None,
            hits=[],
            hit_rankings={},
            directory=None,
            run_id=None,
            alignment_export=None,
            operation=None,
            pdb_details={},
            pdb_details_status="",
        )

    def snapshot(self):
        with self.lock:
            return copy.deepcopy(self.state)

    def update(self, **values):
        with self.lock:
            self.state.update(values)

    def progress(self, event):
        self.update(phase=event.phase, message=event.message)

    def start(self, operation, *, clear=True, expected_run=None, kind="run"):
        with self.lock:
            if self.state["busy"]:
                raise ValueError("A Kmatch operation is already running.")
            if expected_run is not None and self.state["run_id"] != expected_run:
                raise ValueError("The active run changed. Select the hit again before exporting.")
            self.cancel_token = Cancellation()
            self.state.update(
                busy=True, error=None, phase="starting", message="Starting…", operation=kind
            )
            if clear:
                self.details_cancel.cancel()
                self.state.update(pdb_details={}, pdb_details_status="")
                self.state.update(
                    prediction=None,
                    hits=[],
                    hit_rankings={},
                    directory=None,
                    run_id=None,
                    alignment_export=None,
                )

        def work():
            try:
                operation()
            except Cancelled as exc:
                self.update(phase="cancelled", error=None, message=str(exc))
            except Exception as exc:
                self.update(
                    phase="failed",
                    error=str(exc),
                    message="Operation failed; logs and completed outputs are retained",
                )
            finally:
                self.update(busy=False)

        self.thread = threading.Thread(target=work, daemon=False, name="Kmatch worker")
        self.thread.start()

    def run(self, text, backend, network_msa, msa_path=None):
        parse_sequence(text)
        if backend not in BACKENDS or not BACKENDS[backend].production:
            raise ValueError("Choose an available prediction backend.")
        config = replace(self.config, predictor=backend)

        def directory(path):
            self.update(directory=str(path), run_id=path.name)

        def operation():
            result = Pipeline(config).run(
                text,
                cancel=self.cancel_token,
                progress=self.progress,
                network_msa=network_msa,
                msa_path=msa_path,
                on_directory=directory,
                on_prediction=lambda result: self.update(prediction=asdict(result)),
            )
            self.update(
                prediction=asdict(result.prediction),
                hits=[asdict(hit) for hit in result.hits],
                hit_rankings=hit_rankings(result.hits),
            )

        self.start(operation)

    def setup(self, component):
        if component not in ("boltz", "foldseek", "pdb100", "af3-code"):
            raise ValueError("Unknown setup component.")

        def operation():
            setup_component(component, self.config, self.cancel_token, self.progress)
            self.update(phase="complete", message=f"{component} setup finished")

        self.start(operation, clear=False, kind="setup")

    def selected_hit(self, run_id, hit_index):
        state = self.snapshot()
        if (
            not run_id
            or state["run_id"] != run_id
            or not state["prediction"]
            or not state["directory"]
        ):
            raise ValueError("Select a hit from the current completed search.")
        if type(hit_index) is not int or not 0 <= hit_index < len(state["hits"]):
            raise ValueError("Select an existing structural-search hit.")
        return state, state["hits"][hit_index]

    def export_alignment(
        self, run_id, hit_index, format, destination=None, *, full_atom=False, allow_network=False
    ):
        state, hit = self.selected_hit(run_id, hit_index)

        def operation():
            self.update(alignment_export=None)
            result = export_alignment(
                self.config,
                Path(state["directory"]),
                hit_index,
                format,
                cancel=self.cancel_token,
                progress=self.progress,
                destination=destination,
                expected_target=hit["target"],
                full_atom=full_atom,
                allow_network=allow_network,
            )
            payload = asdict(result)
            payload["directory"], payload["artifact"] = str(result.directory), str(result.artifact)
            self.update(alignment_export=payload)

        self.start(operation, clear=False, expected_run=run_id, kind="alignment_export")

    def load_pdb_details(self, run_id, allow_network):
        if type(allow_network) is not bool:
            raise ValueError("Choose whether public PDB network access is enabled.")
        state, _ = self.selected_hit(run_id, 0)
        with self.lock:
            if self.state["run_id"] != run_id:
                raise ValueError("The active run changed. Reload PDB details for the current run.")
            self.details_cancel.cancel()
            cancel = self.details_cancel = Cancellation()
            self.state["pdb_details_status"] = (
                "Loading public PDB details…" if allow_network else "Reading cached PDB details…"
            )

        def work():
            try:
                rankings = hit_rankings([SearchResult(**h) for h in state["hits"]])
                targets = [state["hits"][i]["target"] for i in retained_hit_indices(rankings)]
                details = load_details(self.config, targets, allow_network, cancel)
                ready = sum(d["status"] == "ready" for d in details.values())
                with self.lock:
                    if self.state["run_id"] == run_id and self.details_cancel is cancel:
                        self.state.update(
                            pdb_details=details,
                            pdb_details_status=(
                                f"PDB details: {ready}/{len(details)} hits available (cached)."
                            ),
                        )
            except Cancelled:
                pass
            except Exception as exc:
                with self.lock:
                    if self.state["run_id"] == run_id and self.details_cancel is cancel:
                        self.state["pdb_details_status"] = "PDB details unavailable: " + str(exc)

        self.details_thread = threading.Thread(target=work, daemon=True, name="Kmatch PDB metadata")
        self.details_thread.start()

    def cancel(self):
        self.cancel_token.cancel()
        if self.snapshot()["busy"]:
            self.update(message="Cancelling; waiting for the engine to stop…")

    def close(self):
        self.details_cancel.cancel()
        self.cancel()
        if self.thread:
            self.thread.join(timeout=6)


class DesktopApi:
    """Methods exposed to the page through pywebview's js_api; helpers stay private."""

    def __init__(self, job):
        self._job, self._window = job, None

    def _require_window(self):
        if self._window is None:
            raise RuntimeError("Native window is unavailable. Close and reopen Kmatch.")
        return self._window

    def export_structure(self):
        state = self._job.snapshot()
        if not state["prediction"]:
            raise ValueError("There is no predicted structure to export.")
        import webview

        source = Path(state["prediction"]["structure_path"])
        destination = self._require_window().create_file_dialog(
            webview.FileDialog.SAVE,
            save_filename=source.name,
            file_types=("Structure (*.pdb;*.cif)",),
        )
        if not destination:
            return None
        return self._save_structure(
            state["prediction"], destination if isinstance(destination, str) else destination[0]
        )

    def _save_structure(self, prediction, destination):
        """Copy the predicted structure; the save dialog and installed acceptance share this."""
        source, target = Path(prediction["structure_path"]), Path(destination)
        if target.resolve() == source.resolve():
            return str(target)
        shutil.copy2(source, target)
        if prediction["predictor"] == "af3":
            notice = target.with_name(target.name + ".terms.txt")
            notice.write_text(AF3_EXPORT_NOTICE)
        return str(target)

    def reveal_run(self):
        directory = self._job.snapshot()["directory"]
        if not directory or not Path(directory).is_dir():
            raise ValueError("Run output is not available yet.")
        subprocess.run(["/usr/bin/open", directory], check=True)
        return directory

    def open_pdb(self, run_id, hit_index):
        _, hit = self._job.selected_hit(run_id, hit_index)
        subprocess.Popen(["/usr/bin/open", pdb_url(hit["target"])])

    def export_alignment(self, run_id, hit_index, format, full_atom=False, allow_network=False):
        _, hit = self._job.selected_hit(run_id, hit_index)
        if format not in FORMATS:
            raise ValueError("Choose a supported 3D alignment format.")
        import webview

        result = self._require_window().create_file_dialog(
            webview.FileDialog.SAVE,
            save_filename=f"Kmatch_{safe_name(hit['target'])}_alignment.{format}",
            file_types=(f"3D alignment (*.{format})",),
        )
        if not result:
            return False
        self._start_alignment_export(
            run_id,
            hit_index,
            format,
            result if isinstance(result, str) else result[0],
            full_atom=full_atom,
            allow_network=allow_network,
        )
        return True

    def _start_alignment_export(
        self, run_id, hit_index, format, destination, *, full_atom, allow_network
    ):
        """Begin a background export to a chosen path; shared by the dialog and acceptance."""
        if format not in FORMATS:
            raise ValueError("Choose a supported 3D alignment format.")
        destination = Path(destination)
        if not destination.suffix:
            destination = destination.with_suffix("." + format)
        if destination.suffix.lower() != "." + format:
            raise ValueError(f"Use a .{format} filename for this export format.")
        self._job.export_alignment(
            run_id, hit_index, format, destination, full_atom=full_atom, allow_network=allow_network
        )

    def choose_msa(self):
        import webview

        result = self._require_window().create_file_dialog(
            webview.FileDialog.OPEN, allow_multiple=False, file_types=("A3M alignment (*.a3m)",)
        )
        return result[0] if result else None


def create_app(job, desktop=False):
    app = Flask(__name__, static_folder=str(STATIC), static_url_path="/static")
    app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024
    token = secrets.token_urlsafe(32)

    @app.before_request
    def guard():
        if request.host.split(":")[0] not in ("127.0.0.1", "localhost"):
            return jsonify(error="Only loopback access is allowed."), 403
        if request.path.startswith("/api/"):
            if not secrets.compare_digest(request.headers.get("X-Kmatch-Token", ""), token):
                return jsonify(error="Reload Kmatch to reconnect."), 403
            origin = request.headers.get("Origin")
            if origin and origin != request.host_url.rstrip("/"):
                return jsonify(error="Cross-origin requests are disabled."), 403

    @app.after_request
    def headers(response):
        response.headers["Cache-Control"] = "no-store"
        response.headers["X-Content-Type-Options"] = "nosniff"
        script = "'self' 'unsafe-eval'" if desktop else "'self'"
        response.headers["Content-Security-Policy"] = (
            f"default-src 'self'; script-src {script}; style-src 'self' 'unsafe-inline'; "
            "connect-src 'self'; img-src 'self' data:; frame-ancestors 'none'; object-src 'none'"
        )
        return response

    @app.errorhandler(ValueError)
    @app.errorhandler(OSError)
    @app.errorhandler(RuntimeError)
    def error(exc):
        return jsonify(error=str(exc)), 400

    @app.get("/")
    def index():
        return (
            (STATIC / "index.html")
            .read_text()
            .replace("__TOKEN__", token)
            .replace("__VERSION__", __version__)
        )

    @app.get("/api/status")
    def current():
        state = job.snapshot()
        if state["prediction"]:
            state["prediction"]["structure_path"] = str(state["prediction"]["structure_path"])
        return jsonify(
            job=state,
            desktop=desktop,
            home=str(job.config.home),
            runs=str(job.config.runs),
            defaults={"predictor": job.config.predictor, "samples": job.config.diffusion_samples},
            backends=[
                {
                    "id": key,
                    "name": definition.name,
                    "label": definition.label,
                    "note": definition.note.format(samples=job.config.diffusion_samples),
                }
                for key, definition in BACKENDS.items()
                if definition.production
            ],
            alignment_support=export_support(job.config),
        )

    @app.get("/api/dependencies")
    def dependencies():
        return jsonify(items=status(job.config))

    @app.post("/api/validate")
    def validate():
        protein = parse_sequence(request.get_json()["sequence"])
        return jsonify(length=protein.length, identifier=protein.identifier)

    @app.post("/api/run")
    def run():
        body = request.get_json()
        if not isinstance(body.get("network_msa", False), bool):
            raise ValueError("network_msa must be an explicit boolean selection.")
        msa = body.get("msa_path")
        if msa and not desktop:
            raise ValueError("Custom A3M selection is available in the native desktop app.")
        job.run(
            body["sequence"],
            body.get("backend", "boltz2"),
            bool(body.get("network_msa", False)),
            msa,
        )
        return jsonify(started=True)

    @app.post("/api/cancel")
    def cancel():
        job.cancel()
        return jsonify(cancelled=True)

    @app.post("/api/setup")
    def setup():
        job.setup(request.get_json()["component"])
        return jsonify(started=True)

    @app.get("/api/structure")
    def structure():
        prediction = job.snapshot()["prediction"]
        if not prediction:
            raise ValueError("No predicted structure is available.")
        return send_file(prediction["structure_path"], mimetype="text/plain")

    @app.get("/api/export")
    def export():
        prediction = job.snapshot()["prediction"]
        if not prediction:
            raise ValueError("No predicted structure is available.")
        if prediction["predictor"] == "af3":
            source = Path(prediction["structure_path"])
            archive = io.BytesIO()
            with zipfile.ZipFile(archive, "w") as bundle:
                bundle.write(source, source.name)
                bundle.writestr(source.name + ".terms.txt", AF3_EXPORT_NOTICE)
            archive.seek(0)
            return send_file(
                archive,
                mimetype="application/zip",
                as_attachment=True,
                download_name=source.name + ".zip",
            )
        return send_file(prediction["structure_path"], as_attachment=True)

    @app.post("/api/pdb-details")
    def pdb_details():
        body = request.get_json()
        job.load_pdb_details(body.get("run_id"), body.get("allow_network", False))
        return jsonify(started=True)

    @app.post("/api/alignment-export")
    def begin_alignment_export():
        body = request.get_json()
        # HTTP clients may select a saved hit/format, never paths or arbitrary transforms.
        job.export_alignment(
            body.get("run_id"),
            body.get("hit_index"),
            body.get("format"),
            full_atom=body.get("full_atom", False),
            allow_network=body.get("allow_network", False),
        )
        return jsonify(started=True)

    @app.get("/api/alignment-export/<export_id>")
    def download_alignment(export_id):
        result = job.snapshot()["alignment_export"]
        if not result or export_id != result["export_id"]:
            raise ValueError("This alignment export is no longer the current result.")
        return send_file(
            result["artifact"],
            as_attachment=True,
            download_name=f"Kmatch_{safe_name(result['target'])}_alignment.{result['format']}",
        )

    @app.get("/api/af3-help")
    def af3_help():
        return jsonify(
            url="https://github.com/google-deepmind/alphafold3#obtaining-model-parameters",
            location=str(job.config.af3_models),
            weights="af3.bin.zst",
            note=(
                "AF3 weights and outputs have separate restricted terms. Obtain "
                "parameters directly from Google; Kmatch never bundles them."
            ),
        )

    return app


def write_report(path, report):
    """Atomic acceptance evidence; a polling harness never reads a partial file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(report, indent=2, default=str) + "\n")
    temporary.replace(path)


def dependency_readiness(config):
    """Production readiness from the ordinary status checks. AF3 is deferred, never required."""
    return {
        "boltz2": BoltzBackend(config).status().ready,
        "foldseek": FoldseekBackend(config).executable_status().ready,
        "pdb100": database_status(config.pdb_prefix).ready,
    }


def bridge_contract(api, directory, *, timeout=BRIDGE_TIMEOUT_SECONDS, interval=0.1):
    """Installed native/JS gate: no setup, network or scientific work is triggered."""
    path = Path(directory) / "bridge.json"
    path.unlink(missing_ok=True)  # A reused directory must never show a stale pass.
    report = {"ok": False, "version": __version__, "home": str(api._job.config.home), "bridge": {}}
    try:
        deadline = time.monotonic() + timeout
        while True:
            raw = api._window.run_js("JSON.stringify(window.kmatchBridge || null)")
            bridge = json.loads(raw) if raw else None
            if bridge:
                report["bridge"] = bridge
                if (bridge.get("ready") and bridge.get("desktop")) or bridge.get("error"):
                    break
            if time.monotonic() >= deadline:
                break
            time.sleep(interval)
        bridge = report["bridge"]
        methods = bridge.get("methods") or {}
        report["ok"] = bool(
            bridge.get("ready")
            and bridge.get("desktop")
            and all(methods.get(name) is True for name in BRIDGE_METHODS)
        )
        report["dependencies_ready"] = dependency_readiness(api._job.config)
        if not report["ok"]:
            report["error"] = bridge.get("error") or (
                f"Kmatch's native bridge did not become ready within {timeout:g} s."
            )
    except Exception as exc:
        report["ok"] = False
        report["error"] = f"Native bridge check failed: {exc}"
    write_report(path, report)
    return report


def _page_script(step, body):
    """The exact script one acceptance step evaluates; tests syntax-check these."""
    return (
        f"/*kmatch-smoke:{step}*/JSON.stringify((()=>{{try{{return {{value:(()=>{{{body}}})()}};}}"
        "catch(e){return {error:String(e&&e.message||e)};}})())"
    )


def _page(window, step, body=None):
    """Run one acceptance step in the real page; page exceptions become Python errors."""
    raw = window.run_js(_page_script(step, SMOKE_STEPS[step] if body is None else body))
    if not raw:
        raise RuntimeError(f"The Kmatch page did not answer ({step}).")
    result = json.loads(raw)
    if "error" in result:
        raise RuntimeError(f"Kmatch page step '{step}' failed: {result['error']}")
    return result.get("value")


class SmokeTimeout(RuntimeError):
    pass


def _wait(check, timeout, closed, what, deadline, interval=0.25):
    """Poll CHECK for at most TIMEOUT seconds and never past the smoke's overall DEADLINE."""
    limit = min(time.monotonic() + timeout, deadline)
    while True:
        if closed.is_set():
            raise RuntimeError("The Kmatch window closed during acceptance; its job was cancelled.")
        value = check()
        if value:
            return value
        now = time.monotonic()
        if now >= limit:
            if limit >= deadline:
                raise SmokeTimeout(
                    f"The {SMOKE_DEADLINE_SECONDS:g} s smoke deadline expired while waiting "
                    f"for {what}."
                )
            raise SmokeTimeout(f"Timed out after {timeout:g} s waiting for {what}.")
        time.sleep(min(interval, limit - now))


def _finished(job):
    state = job.snapshot()
    return None if state["busy"] else state


def _stop(job, until, closed):
    """Cancel any running job or export and wait, bounded by UNTIL, until Kmatch is idle."""
    if not job.snapshot()["busy"]:
        return True
    job.cancel()
    if closed.is_set():
        until = min(until, time.monotonic() + SMOKE_CLOSED_GRACE_SECONDS)
    while job.snapshot()["busy"] and time.monotonic() < until:
        time.sleep(0.05)
    return not job.snapshot()["busy"]


def _structure_residues(path, chain=None):
    import gemmi

    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        return None
    chains = {item.name: len(item) for item in structure[0]}
    return chains if chain is None else chains.get(chain)


SMOKE_READY = (
    "return window.kmatchBridge?.desktop===true&&window.kmatchBridge?.ready===true"
    "&&!document.getElementById('run').disabled;"
)
SMOKE_CONFIGURE = """
const el=id=>document.getElementById(id);
const sequence=el('sequence'); sequence.value=TEXT; sequence.dispatchEvent(new Event('input'));
const mode=el('msa-mode'); mode.value='single'; mode.dispatchEvent(new Event('change'));
if(el('pdb-network').checked) el('pdb-network').click();
if(el('full-atom').checked) el('full-atom').click();
return {backend:el('backend').value, msa_mode:mode.value,
  msa_label:mode.selectedOptions[0]?mode.selectedOptions[0].textContent:'',
  pdb_network:el('pdb-network').checked, full_atom:el('full-atom').checked};
"""
SMOKE_VALIDATED = (
    "const v=document.getElementById('validation');"
    "return {text:v.textContent,error:v.classList.contains('error')};"
)
SMOKE_RUN = "document.getElementById('run').click();return true;"
SMOKE_PAGE_ERROR = "const e=document.getElementById('error');return e.hidden?'':e.textContent;"
SMOKE_RENDERED = """
const el=id=>document.getElementById(id);
let atoms=0;
try{atoms=viewer?(viewer.getModel()?.selectedAtoms({}).length||0):0;}catch(e){atoms=-1;}
return {rows:document.querySelectorAll('#hits tr').length, results_visible:!el('results').hidden,
  viewer_empty_hidden:el('viewer-empty').hidden, canvas:!!el('viewer').querySelector('canvas'),
  atoms, loaded_run:typeof loadedRun==='undefined'?null:loadedRun,
  error:el('error').hidden?'':el('error').textContent,
  summary:el('prediction-summary').textContent, pdb_status:el('pdb-status').textContent};
"""
SMOKE_SELECT_TOP = """
const el=id=>document.getElementById(id);
const row=document.querySelector('#hits tr');
if(!row) throw new Error('The results table has no rows.');
const button=row.querySelector('td:nth-child(2) button'); button.click();
const link=row.querySelector('a.pdb-link');
return {hit_index:Number(row.dataset.hitIndex), rank:row.cells[0].textContent,
  target:button.textContent, selected:row.classList.contains('selected'),
  detail_visible:!el('hit-detail').hidden, title:el('hit-title').textContent,
  metrics:el('metrics').children.length, rcsb:link?link.href:null, ranking:el('hit-sort').value,
  export_enabled:!el('export-alignment').disabled, save_enabled:!el('export').disabled};
"""
SMOKE_ALIGNMENT_RESULT = (
    "const r=document.getElementById('alignment-result');"
    "return !r.hidden&&r.textContent.includes('Saved to')?r.textContent:'';"
)
# Every page step the smoke evaluates, by step name (TEXT in configure is the JSON sequence).
SMOKE_STEPS = {
    "ready": SMOKE_READY,
    "configure": SMOKE_CONFIGURE,
    "validated": SMOKE_VALIDATED,
    "run": SMOKE_RUN,
    "page-error": SMOKE_PAGE_ERROR,
    "rendered": SMOKE_RENDERED,
    "select-top-hit": SMOKE_SELECT_TOP,
    "alignment-result": SMOKE_ALIGNMENT_RESULT,
    "final-error": SMOKE_PAGE_ERROR,
}
# Page configuration the smoke requires before it may press Run.
SMOKE_PRECONDITIONS = (
    "boltz2_backend",
    "local_only_no_msa",
    "pdb_network_disabled",
    "full_atom_disabled",
)


def smoke(api, directory, fasta, closed=None):
    """Opt-in installed acceptance: one real local MSA-free job through the desktop page.

    Runs write under DIRECTORY/runs and exports under DIRECTORY/exports. Run is pressed only
    after the page shows Boltz-2, "Local only · no MSA", public PDB details off and full-atom
    export off, so the job and exports read only local/cached data. The whole smoke shares
    one SMOKE_DEADLINE_SECONDS budget; a job or export still running when it stops (timeout,
    failure or a closed window) is cancelled so Kmatch is idle when report.json is written.
    """
    directory = Path(directory)
    (directory / "report.json").unlink(missing_ok=True)  # Never leave a stale pass.
    closed = closed or threading.Event()
    job, window = api._job, api._window
    report = {
        "ok": False,
        "version": __version__,
        "deadline_seconds": SMOKE_DEADLINE_SECONDS,
        "run_id": None,
        "run_directory": None,
        "prediction": None,
        "warnings": [],
        "hits": 0,
        "top_hit": None,
        "exports": {},
        "checks": {},
    }
    checks = report["checks"]
    started = time.monotonic()
    deadline = started + SMOKE_DEADLINE_SECONDS
    stop_by = deadline + SMOKE_CANCEL_GRACE_SECONDS

    def wait(check, timeout, what):
        return _wait(check, timeout, closed, what, deadline)

    def stop():
        return _stop(job, min(time.monotonic() + SMOKE_CANCEL_GRACE_SECONDS, stop_by), closed)

    try:
        bridge = bridge_contract(
            api, directory, timeout=min(BRIDGE_TIMEOUT_SECONDS, deadline - time.monotonic())
        )
        if not bridge["ok"]:
            raise RuntimeError("Native bridge check failed: " + bridge.get("error", "not ready"))
        text = Path(fasta).read_text()
        protein = parse_sequence(text)
        report["input"] = {
            "path": str(fasta),
            "identifier": protein.identifier,
            "length": protein.length,
        }
        report["dependencies_ready"] = bridge["dependencies_ready"]
        checks["dependencies_ready"] = all(bridge["dependencies_ready"].values())
        if not checks["dependencies_ready"]:
            missing = [name for name, ready in bridge["dependencies_ready"].items() if not ready]
            raise RuntimeError("Kmatch dependencies are not ready: " + ", ".join(missing))
        user_runs = job.config.home / "runs"
        before = sorted(p.name for p in user_runs.iterdir()) if user_runs.is_dir() else []
        wait(lambda: _page(window, "ready"), 30, "the Kmatch page")
        controls = _page(
            window,
            "configure",
            SMOKE_CONFIGURE.replace("TEXT", json.dumps(text)),
        )
        report["controls"] = controls
        checks["boltz2_backend"] = controls["backend"] == "boltz2"
        checks["local_only_no_msa"] = (
            controls["msa_mode"] == "single" and controls["msa_label"] == "Local only · no MSA"
        )
        checks["pdb_network_disabled"] = controls["pdb_network"] is False
        checks["full_atom_disabled"] = controls["full_atom"] is False
        unmet = [name for name in SMOKE_PRECONDITIONS if not checks[name]]
        if unmet:
            raise RuntimeError(
                "The page did not reach the required local-only configuration ("
                + ", ".join(unmet)
                + "); Run was not pressed."
            )

        def validated():
            result = _page(window, "validated")
            done = result["text"] == "Valid single protein sequence" or result["error"]
            return result if done else None

        validation = wait(validated, 15, "sequence validation")
        if validation["error"]:
            raise RuntimeError("The page rejected the sequence: " + validation["text"])
        if job.snapshot()["busy"]:
            raise RuntimeError("Kmatch is already busy; the smoke did not press Run.")
        _page(window, "run")
        run_started = time.monotonic()

        def launched():
            state = job.snapshot()
            if state["run_id"] or state["phase"] in ("failed", "cancelled"):
                return state
            page_error = _page(window, "page-error")
            if page_error:
                raise RuntimeError("The page could not start the run: " + page_error)
            return None

        wait(launched, 60, "the run to start")
        state = wait(lambda: _finished(job), SMOKE_RUN_TIMEOUT_SECONDS, "the run")
        report["run_id"], report["run_directory"] = state["run_id"], state["directory"]
        report["timings"] = {"run_seconds": time.monotonic() - run_started}
        if state["phase"] != "complete" or not state["prediction"]:
            raise RuntimeError("The run did not complete: " + (state["error"] or state["message"]))
        prediction = state["prediction"]
        confidence, provenance = prediction["confidence"], prediction.get("provenance", {})
        report["prediction"] = {
            "predictor": prediction["predictor"],
            "predictor_version": prediction["predictor_version"],
            "sequence_length": prediction["sequence_length"],
            "mean_plddt": confidence.get("mean_plddt"),
            "ptm": confidence.get("ptm"),
            "confidence_score": confidence.get("confidence_score"),
            "runtime_seconds": prediction["runtime_seconds"],
            "structure_path": str(prediction["structure_path"]),
            "msa": provenance.get("msa"),
            "network_during_inference": provenance.get("network_during_inference"),
        }
        report["warnings"] = prediction.get("warnings", [])
        report["returned_hits"] = len(state["hits"])
        run_directory = Path(state["directory"]).resolve()
        manifest = json.loads((run_directory / "run.json").read_text())
        checks["run_isolated"] = run_directory.is_relative_to((directory / "runs").resolve())
        checks["msa_free_local_inference"] = (
            manifest.get("msa_mode") == "single-sequence"
            and provenance.get("msa") == "single-sequence"
            and provenance.get("network_during_inference") is False
        )
        checks["prediction_length"] = prediction["sequence_length"] == protein.length

        def rendered():
            view = _page(window, "rendered")
            done = view["rows"] and view["viewer_empty_hidden"] and view["loaded_run"]
            return view if done or view["error"] else None

        view = wait(rendered, 120, "the structure and results to render")
        report["hits"] = view["rows"]
        report["pdb_details_status"] = view["pdb_status"]
        checks["structure_rendered"] = bool(
            view["viewer_empty_hidden"]
            and view["canvas"]
            and view["atoms"] > 0
            and view["loaded_run"] == state["run_id"]
            and not view["error"]
        )
        checks["results_table"] = view["results_visible"] and 1 <= view["rows"] <= 50
        top = _page(window, "select-top-hit")
        index = top["hit_index"]
        hit = state["hits"][index]
        try:
            expected_link = pdb_url(hit["target"])
        except ValueError:
            expected_link = None
        checks["top_hit_selected"] = (
            top["selected"]
            and top["detail_visible"]
            and top["target"] == hit["target"]
            and top["title"] == f"#1 · {hit['target']}"
            and index == state["hit_rankings"][top["ranking"]][0]
            and top["metrics"] == 12
        )
        checks["rcsb_link"] = bool(expected_link) and top["rcsb"] == expected_link
        checks["native_export_controls"] = top["export_enabled"] and top["save_enabled"]
        report["top_hit"] = {
            "rank": 1,
            "ranking": top["ranking"],
            "hit_index": index,
            "rcsb_url": top["rcsb"],
            **{
                key: hit[key]
                for key in (
                    "target",
                    "query_tm",
                    "target_tm",
                    "alignment_tm",
                    "lddt",
                    "probability",
                    "evalue",
                    "bits",
                    "sequence_identity",
                    "query_coverage",
                    "target_coverage",
                    "alignment_length",
                )
            },
        }
        exports = directory / "exports"
        exports.mkdir(parents=True, exist_ok=True)
        source = Path(prediction["structure_path"])
        structure = exports / f"Kmatch_{safe_name(protein.identifier)}_predicted{source.suffix}"
        report["exports"]["predicted_structure"] = api._save_structure(prediction, structure)
        chains = _structure_residues(structure)
        checks["structure_export"] = chains is not None and [*chains.values()] == [protein.length]

        def alignment(format):
            destination = exports / f"Kmatch_{safe_name(hit['target'])}_alignment.{format}"
            api._start_alignment_export(
                state["run_id"], index, format, destination, full_atom=False, allow_network=False
            )
            result = wait(lambda: _finished(job), SMOKE_EXPORT_TIMEOUT_SECONDS, "the 3D export")
            exported = result["alignment_export"]
            if result["phase"] != "complete" or not exported:
                raise RuntimeError(result["error"] or result["message"])
            if exported["saved_to"] != str(destination) or not destination.is_file():
                raise RuntimeError(f"The 3D export was not saved to {destination}.")
            return destination

        superposition = alignment("cif")
        report["exports"]["superposition_cif"] = str(superposition)
        chains = _structure_residues(superposition)
        checks["superposition_export"] = (
            chains is not None and set(chains) == {"Q", "T"} and chains["Q"] == protein.length
        )
        shown = wait(lambda: _page(window, "alignment-result"), 15, "the export result to appear")
        checks["superposition_export_shown"] = str(superposition) in shown
        checks["no_page_error"] = not _page(window, "final-error")
        after = sorted(p.name for p in user_runs.iterdir()) if user_runs.is_dir() else []
        checks["user_runs_unchanged"] = before == after
        # Optional, never part of ok: licensed local PyMOL may be absent or unusable.
        if not export_support(job.config)["pse"]:
            report["pymol_session"] = {
                "status": "skipped",
                "reason": "PyMOL was not detected; session export is optional.",
            }
        elif time.monotonic() >= deadline:
            report["pymol_session"] = {
                "status": "skipped",
                "reason": "The smoke deadline left no time for the optional session export.",
            }
        else:
            try:
                session = alignment("pse")
                report["exports"]["pymol_session"] = str(session)
                report["pymol_session"] = {"status": "passed", "path": str(session)}
            except Exception as exc:
                # A timed-out optional export is cancelled; it must not leave Kmatch busy.
                cancelled = job.snapshot()["busy"]
                report["pymol_session"] = {
                    "status": "failed",
                    "error": str(exc),
                    "cancelled": cancelled,
                    "idle": stop(),
                }
    except Exception as exc:
        report["error"] = str(exc)
    # A job or export still running after a timeout/failure is cancelled, never left unattended.
    report["job_idle"] = stop()
    final = job.snapshot()
    report["final_job"] = {"busy": final["busy"], "phase": final["phase"]}
    if not report["job_idle"]:
        report.setdefault(
            "error",
            "A cancelled Kmatch job or export was still running when the report was written.",
        )
    report["timings"] = {**report.get("timings", {}), "total_seconds": time.monotonic() - started}
    report["ok"] = (
        "error" not in report and bool(checks) and all(checks.values()) and report["job_idle"]
    )
    write_report(directory / "report.json", report)
    return report


def startup_error_page(message, home):
    config = html.escape(str(Path(home) / "config.json"))
    return (
        "<!doctype html><html lang='en'><head><meta charset='utf-8'><title>Kmatch</title>"
        "<style>body{font:14px -apple-system,BlinkMacSystemFont,sans-serif;margin:28px 32px;"
        "color:#23352d;background:#f4f6f5}h1{font-size:19px;margin:0 0 12px}"
        "pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#fff;color:#9b3126;"
        "border:1px solid #dce4df;border-radius:6px;padding:12px}</style></head><body>"
        "<h1>Kmatch could not start</h1>"
        f"<pre>{html.escape(message)}</pre>"
        f"<p>Kmatch reads <code>{config}</code> and the shared-resource pins it names or, "
        "without a setting, the registered default pins. "
        "Correct the setting or repair the shared resources, then reopen Kmatch. "
        "Nothing was changed.</p></body></html>"
    )


def report_startup_error(args, message):
    """Configuration failures are shown in a window rather than ending a launch silently."""
    print(message, file=sys.stderr, flush=True)
    home = default_home()
    failure = {"ok": False, "version": __version__, "home": str(home), "error": message}
    for directory in filter(None, (args.bridge_check_dir, args.smoke_dir)):
        write_report(directory / "bridge.json", failure)
    if args.smoke_dir:
        write_report(args.smoke_dir / "report.json", failure)
    try:
        import webview
    except ImportError:
        raise SystemExit(1) from None
    webview.create_window(
        f"Kmatch · {__version__}",
        html=startup_error_page(message, home),
        width=640,
        height=360,
        background_color="#f4f6f5",
        text_select=True,
    )
    webview.start(debug=False)
    raise SystemExit(1)


def build_parser():
    parser = argparse.ArgumentParser(
        description=f"Kmatch {__version__} · local protein structure prediction and search"
    )
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument(
        "--browser", action="store_true", help="Serve a local browser page instead of a window"
    )
    parser.add_argument(
        "--no-open", action="store_true", help="With --browser, do not open a browser"
    )
    parser.add_argument("--port", type=int, default=0, help="Loopback port (default: free port)")
    parser.add_argument(
        "--bridge-check-dir",
        type=Path,
        help="Installed acceptance: write DIR/bridge.json once the native bridge initializes",
    )
    parser.add_argument(
        "--smoke-dir",
        type=Path,
        help="Installed acceptance: bridge check plus one real local job; runs/exports stay in DIR",
    )
    parser.add_argument(
        "--smoke-input", type=Path, help="Public single-protein FASTA for --smoke-dir"
    )
    return parser


def parse_args(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.smoke_input and not args.smoke_dir:
        parser.error("--smoke-input requires --smoke-dir")
    if args.smoke_dir and not args.smoke_input:
        parser.error("--smoke-dir requires --smoke-input FASTA")
    if args.browser and (args.smoke_dir or args.bridge_check_dir):
        parser.error("--bridge-check-dir and --smoke-dir check the desktop window; omit --browser")
    for name in ("bridge_check_dir", "smoke_dir", "smoke_input"):
        if getattr(args, name):
            setattr(args, name, getattr(args, name).expanduser().resolve())
    return parser, args


def main(argv=None):
    parser, args = parse_args(argv)
    try:
        config = Config.load()
    except (OSError, ValueError, RuntimeError) as exc:
        message = f"Kmatch could not load its configuration: {exc}"
        if args.browser:
            parser.exit(1, message + "\n")
        report_startup_error(args, message)
    if args.smoke_dir:
        # Acceptance runs never enter the user's run history.
        config = replace(config, runs_root=args.smoke_dir / "runs")
    if not args.browser:
        try:
            import webview
        except ImportError:
            parser.exit(1, "Desktop dependency missing. Run uv sync, or use --browser.\n")
    job = Job(config)
    app = create_app(job, desktop=not args.browser)
    server = loopback_server(app, args.port)
    server_thread = threading.Thread(target=server.serve_forever, daemon=True)
    server_thread.start()
    url = f"http://127.0.0.1:{server.server_port}"
    print(f"Kmatch: {url}", flush=True)
    checker = None
    try:
        if args.browser:
            if not args.no_open:
                import webbrowser

                webbrowser.open(url)
            server_thread.join()
        else:
            api = DesktopApi(job)
            closed = threading.Event()
            api._window = webview.create_window(
                f"Kmatch · {__version__}",
                url,
                js_api=api,
                width=1250,
                height=880,
                min_size=(1024, 720),
                background_color="#f4f6f5",
                text_select=True,
            )
            api._window.events.closed += closed.set
            api._window.events.closed += job.close

            def acceptance():
                if args.bridge_check_dir:
                    bridge_contract(api, args.bridge_check_dir)
                if args.smoke_dir:
                    smoke(api, args.smoke_dir, args.smoke_input, closed)

            def started():
                nonlocal checker
                if args.bridge_check_dir or args.smoke_dir:
                    # Daemon: a stalled page call must never prevent a normal quit.
                    checker = threading.Thread(target=acceptance, daemon=True, name="Kmatch check")
                    checker.start()

            webview.start(func=started, debug=False)
    except KeyboardInterrupt:
        pass
    finally:
        job.close()
        if checker:
            checker.join(timeout=10)
        server.shutdown()
        server.server_close()


if __name__ == "__main__":
    main()
