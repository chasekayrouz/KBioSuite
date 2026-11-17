"""Local file annotation and disposable protein search around the canonical API."""

import argparse
import copy
import json
import subprocess
import threading
import time
import uuid
from importlib.resources import files
from pathlib import Path

from flask import jsonify, request

from kbiosuite.desktop import host, local_app

from . import AnnotationEngine, __version__
from .batch import GENBANK_OUTPUTS, annotate_batch
from .cancellation import RunCancelled, RunControl, run_scope
from .config import default_root, default_threads
from .databases import setup, status
from .quick_search import QuickSearchSession, import_queries, parse_queries, search
from .scoring import DEFAULT_MIN_SCORE
from .search.identity import fingerprint

STATIC = Path(__file__).parent / "static"


class AnnotationJob:
    def __init__(self):
        self.lock = threading.Lock()
        self.state = {
            "busy": False,
            "message": "Ready",
            "results": [],
            "error": None,
            "job_id": None,
            "kind": None,
            "can_stop": False,
            "stopping": False,
            "stopped": False,
        }
        self._control = None
        self.quick = None
        self.quick_revision = 0
        self.search_session = QuickSearchSession()
        # None keeps the installation's persistent cache; acceptance runs isolate it.
        self.quick_cache_dir = None

    def snapshot(self, *, include_quick=False):
        with self.lock:
            result = copy.deepcopy(self.state)
            result["quick_revision"] = self.quick_revision
            if include_quick:
                result["quick"] = copy.deepcopy(self.quick)
            return result

    def update(self, **values):
        with self.lock:
            self.state.update(values)

    def clear_quick(self):
        with self.lock:
            if self.state["busy"]:
                raise ValueError("Wait for the current operation to finish.")
            self.quick = None
            self.quick_revision += 1
            self.state.update(message="Ready", error=None, stopped=False)
            return self.quick_revision

    def start(self, operation, *, preserve_results=False, kind=None, checkpoint_dir=None):
        control = RunControl()
        with self.lock:
            if self.state["busy"]:
                raise ValueError("Wait for the current operation to finish.")
            results = self.state["results"] if preserve_results else []
            self._control = control
            self.state = {
                "busy": True,
                "message": "Starting…",
                "results": results,
                "error": None,
                "job_id": uuid.uuid4().hex,
                "kind": kind,
                "can_stop": kind in {"files", "quick"},
                "stopping": False,
                "stopped": False,
                "checkpoint_dir": str(checkpoint_dir) if checkpoint_dir else None,
            }

        def work():
            try:
                with run_scope(control):
                    publish = operation()
                    control.finish(publish if callable(publish) else None)
            except RunCancelled:
                completed = sum(
                    len(r.get("outputs", [])) or bool(r.get("output"))
                    for r in self.snapshot()["results"]
                )
                self.update(
                    stopped=True,
                    message=(
                        "Stopped. Unfinished work discarded; temporary search data removed. "
                        f"{completed} completed output(s) kept."
                        + (" Completed checkpoints kept for resume." if checkpoint_dir else "")
                    ),
                )
            except Exception as exc:
                stopping = self.snapshot()["stopping"]
                self.update(
                    stopped=stopping,
                    error=str(exc),
                    message="Stopped; cleanup needs attention" if stopping else "Operation failed",
                )
            finally:
                self.update(busy=False, can_stop=False, stopping=False)

        threading.Thread(target=work, daemon=False).start()

    def stop(self, job_id):
        with self.lock:
            if self.state["job_id"] != job_id or not self.state["busy"]:
                raise ValueError("That search has already finished. Refresh before stopping a run.")
            if not self.state["can_stop"] and not self.state["stopping"]:
                raise ValueError("This operation cannot be stopped.")
            control = self._control

        def mark_stopping():
            with self.lock:
                if self._control is control and self.state["busy"]:
                    self.state.update(
                        can_stop=False,
                        stopping=True,
                        message="Stopping both searches and removing interrupted search data…",
                    )

        if not control.request_stop(mark_stopping):
            raise ValueError("That search has already finished.")

    def quick_search(self, proteins):
        def operation():
            self.update(message=f"Searching {len(proteins)} proteins locally · Pfam + CDD…")
            result = search(
                proteins, cache_dir=self.quick_cache_dir, session=self.search_session
            )

            def publish():
                with self.lock:
                    self.quick = result
                    self.quick_revision += 1
                    self.state.update(
                        message=f"Quick Search complete · {result['timings']['wall']:.2f} s"
                    )

            return publish

        self.start(operation, preserve_results=True, kind="quick")

    def annotate(self, inputs, destination, threads, genbank_output="genbank", checkpoint=False):
        checkpoint_dir = (
            default_root()
            / "checkpoints"
            / fingerprint(
                {
                    "inputs": [str(Path(p).resolve()) for p in inputs],
                    "destination": str(Path(destination).resolve()),
                }
            )
            if checkpoint
            else None
        )

        def operation():
            self.update(message="Loading local databases…")
            self.search_session.close()
            engine = AnnotationEngine.from_config()
            job = self

            def progress(path, count):
                index = [Path(p) for p in inputs].index(path) + 1
                job.update(
                    message=f"Searching file {index} of {len(inputs)} · {count} proteins completed"
                )

            def report(result):
                return {
                    "input": str(result.input),
                    "output": str(result.output) if result.output else None,
                    "outputs": [str(p) for p in result.outputs],
                    "error": result.error,
                    "warnings": result.warnings or [],
                    "proteins": result.run.protein_count if result.run else 0,
                    "timings": result.run.timings if result.run else {},
                }

            def file_completed(result):
                with self.lock:
                    self.state["results"].append(report(result))

            results = annotate_batch(
                inputs,
                destination,
                engine,
                threads=threads,
                min_score=DEFAULT_MIN_SCORE,
                force=False,
                on_result=file_completed,
                genbank_output=genbank_output,
                retain_annotations=False,
                on_progress=progress,
                checkpoint_dir=checkpoint_dir,
                resume=bool(checkpoint_dir and checkpoint_dir.exists()),
            )
            succeeded = sum(not r.error for r in results)
            return lambda: self.update(
                message=f"Completed: {succeeded} succeeded, {len(results) - succeeded} failed",
            )

        self.start(operation, kind="files", checkpoint_dir=checkpoint_dir)


class DesktopApi:
    def __init__(self, job):
        self._job, self._window = job, None

    def _require_window(self):
        if self._window is None:
            raise RuntimeError("Native window is unavailable. Close and reopen Kannotate.")
        return self._window

    def choose_inputs(self):
        import webview

        paths = self._require_window().create_file_dialog(
            webview.FileDialog.OPEN,
            allow_multiple=True,
            file_types=(
                "Protein FASTA and GenBank "
                "(*.fa;*.faa;*.fasta;*.fas;*.fsa;*.gb;*.gbk;*.gbff;*.genbank)",
                "All files (*.*)",
            ),
        )
        return list(paths or [])

    def choose_output(self):
        import webview

        paths = self._require_window().create_file_dialog(webview.FileDialog.FOLDER)
        return paths[0] if paths else None

    def reveal(self, path):
        allowed = [
            path
            for r in self._job.snapshot()["results"]
            for path in r.get("outputs", [r["output"]] if r.get("output") else [])
        ]
        if path not in allowed:
            raise ValueError("Choose an output from this run.")
        subprocess.run(["/usr/bin/open", "-R", path], check=True)

    def can_close(self):
        if self._job.snapshot()["busy"]:
            return False
        self._job.search_session.close()
        return True


def create_app(job, desktop=False):
    app, token = local_app(__name__, STATIC, desktop=desktop)

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
        return jsonify(
            runtime=str(default_root()),
            status=status(),
            job=job.snapshot(),
            desktop=desktop,
            threads=default_threads(),
        )

    @app.get("/api/job")
    def current_job():
        # Frequent progress polls must not launch executable/version probes.
        return jsonify(job.snapshot())

    @app.post("/api/stop")
    def stop():
        data = request.get_json() or {}
        if (
            not isinstance(data, dict)
            or data.get("confirmed") is not True
            or not isinstance(data.get("job_id"), str)
            or not data["job_id"]
        ):
            raise ValueError("Confirm stopping the current search and discarding unfinished work.")
        job.stop(data["job_id"])
        return jsonify(ok=True), 202

    @app.post("/api/setup")
    def configure():
        job.start(
            lambda: (
                setup(report=lambda message: job.update(message=message)),
                job.update(message="Database validation succeeded"),
            )
        )
        return jsonify(ok=True)

    @app.post("/api/annotate")
    def annotate():
        data = request.get_json() or {}
        if not isinstance(data, dict):
            raise ValueError("Expected annotation options.")
        if not isinstance(data.get("inputs"), list) or not data["inputs"]:
            raise ValueError("Choose at least one protein FASTA or GenBank file.")
        if not all(isinstance(p, str) and p for p in data["inputs"]):
            raise ValueError("Invalid input paths.")
        if not isinstance(data.get("output"), str) or not data["output"].strip():
            raise ValueError("Choose an output destination.")
        threads = data.get("threads", default_threads())
        if type(threads) is not int or threads < 1:
            raise ValueError("Threads must be a positive integer.")
        genbank_output = data.get("genbank_output", "genbank")
        if genbank_output not in GENBANK_OUTPUTS:
            raise ValueError("GenBank output must be genbank, csv, both, or xlsx.")
        checkpoint = data.get("checkpoint", False)
        if type(checkpoint) is not bool:
            raise ValueError("Checkpoint option must be true or false.")
        if checkpoint and genbank_output == "xlsx":
            raise ValueError("Checkpoints support CSV and GenBank outputs.")
        inputs = [Path(p).expanduser().resolve() for p in data["inputs"]]
        if not all(p.is_file() for p in inputs):
            raise ValueError("An input file is missing.")
        destination = Path(data["output"]).expanduser().resolve()
        job.annotate(
            inputs, destination, threads, genbank_output=genbank_output, checkpoint=checkpoint
        )
        return jsonify(ok=True), 202

    @app.post("/api/quick/import")
    def import_quick():
        data = request.get_json() or {}
        if not isinstance(data, dict):
            raise ValueError("Expected a protein FASTA or GenBank file.")
        proteins, warnings = import_queries(data.get("text"), data.get("filename", ""))
        return jsonify(
            fasta="".join(f">{p.header}\n{p.sequence}\n" for p in proteins),
            warnings=warnings,
        )

    @app.post("/api/quick/validate")
    def validate_quick():
        data = request.get_json() or {}
        if not isinstance(data, dict):
            raise ValueError("Expected protein input text.")
        proteins = parse_queries(data.get("text"))
        return jsonify(count=len(proteins), lengths=[len(p.sequence) for p in proteins])

    @app.post("/api/quick/search")
    def quick_search():
        data = request.get_json() or {}
        if not isinstance(data, dict):
            raise ValueError("Expected protein input text.")
        proteins = parse_queries(data.get("text"))
        job.quick_search(proteins)
        return jsonify(ok=True), 202

    @app.get("/api/quick/results")
    def quick_results():
        state = job.snapshot(include_quick=True)
        return jsonify(result=state["quick"], revision=state["quick_revision"])

    @app.post("/api/quick/clear")
    def clear_quick():
        return jsonify(revision=job.clear_quick())

    return app


def bridge_contract(api, directory):
    """Installed native/JS gate: no inputs, setup or scientific work is triggered."""
    directory.mkdir(parents=True, exist_ok=True)
    report = {"ok": False, "version": __version__}
    try:
        for _ in range(150):
            state = api._window.run_js("JSON.stringify(window.kannotateBridge || {})")
            if state:
                report["bridge"] = json.loads(state)
                if (report["bridge"].get("ready") and report["bridge"].get("desktop")) or report[
                    "bridge"
                ].get("error"):
                    break
            time.sleep(0.1)
        report["ok"] = bool(report.get("bridge", {}).get("ready"))
        runtime_status = status()
        report["runtime"] = str(default_root())
        report["databases_ready"] = {
            name: value.get("ready", False) for name, value in runtime_status["databases"].items()
        }
        report["controls"] = json.loads(
            api._window.run_js(
                "JSON.stringify(['choose-inputs','choose-output'].map(id=>"
                "({id,disabled:document.getElementById(id).disabled})))"
            )
        )
        report["ok"] = report["ok"] and all(not c["disabled"] for c in report["controls"])
    except Exception as exc:
        report["error"] = str(exc)
    (directory / "bridge.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def quick_smoke(api, text):
    """Opt-in installed-bundle gate; use public qualification input supplied by the harness."""
    from io import StringIO

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    proteins = parse_queries(text)
    record = SeqRecord(Seq("N" * sum(len(p.sequence) * 3 for p in proteins)), id="public")
    record.annotations["molecule_type"] = "DNA"
    offset = 0
    for i, protein in enumerate(proteins):
        end = offset + len(protein.sequence) * 3
        record.features.append(
            SeqFeature(
                SimpleLocation(offset, end),
                type="CDS",
                qualifiers={"locus_tag": [f"query_{i + 1}"], "translation": [protein.sequence]},
            )
        )
        offset = end
    genbank = StringIO()
    SeqIO.write(record, genbank, "genbank")
    # File-job completion can precede the next GUI state poll. A real user cannot
    # select a file until this control unlocks; synthetic selection must wait too.
    for _ in range(100):
        if api._window.run_js("!document.getElementById('quick-file').disabled"):
            break
        time.sleep(0.05)
    else:
        raise RuntimeError("Installed Quick Search file control is unavailable")
    started = time.monotonic()
    api._window.run_js(
        "document.getElementById('open-quick').click();"
        "(()=>{const files=new DataTransfer();"
        f"files.items.add(new File([{json.dumps(genbank.getvalue())}], 'public quick.gbk'));"
        "const input=document.getElementById('quick-file');input.files=files.files;"
        "input.dispatchEvent(new Event('change'));})(); true"
    )
    for _ in range(100):
        if api._window.run_js("!document.getElementById('quick-search').disabled"):
            break
        time.sleep(0.05)
    else:
        raise RuntimeError("Installed Quick Search control is unavailable")
    api._window.run_js("document.getElementById('quick-search').click(); true")
    for _ in range(1200):
        if api._window.run_js(
            f"document.querySelectorAll('#quick-queries button').length === {len(proteins)}"
            " && !document.getElementById('quick-results').hidden"
        ):
            break
        if api._job.snapshot()["error"]:
            raise RuntimeError(api._job.snapshot()["error"])
        time.sleep(0.05)
    else:
        raise RuntimeError("Installed Quick Search did not render results")
    wall = time.monotonic() - started
    checks = json.loads(
        api._window.run_js("""JSON.stringify((()=>{
      const el=id=>document.getElementById(id), checks={};
      const queries=[...el('quick-queries').children];
      checks.five_queries=queries.length===5;
      checks.architecture=!!el('domain-architecture').querySelector('[data-hit]');
      document.querySelector('[data-source="CDD"]').click();
      const rows=[...el('domain-rows').children];
      checks.cdd_filter=rows.length>0&&rows.every(r=>r.cells[0].textContent==='CDD');
      rows[0].click(); const hit=rows[0].dataset.hit;
      checks.row_to_graph=el('domain-architecture').querySelector('.selected')?.dataset.hit===hit;
      const block=el('domain-architecture').querySelector('[data-hit]');
      block.dispatchEvent(new MouseEvent('click',{bubbles:true}));
      checks.graph_to_row=el('domain-rows').querySelector('.selected')?.dataset.hit===block.dataset.hit;
      el('show-details').click();
      checks.cdd_alignment=!el('alignment-text').hidden&&el('alignment-text').textContent.includes('Query')
        &&el('alignment-text').textContent.includes('Model');
      queries[1].click(); const second=el('query-name').textContent;
      queries[0].click();
      checks.query_navigation=second!==el('query-name').textContent
        &&el('domain-rows').querySelector('.selected')?.dataset.hit===block.dataset.hit;
      queries[4].click(); el('show-hits').click();
      checks.zero_hits=el('domain-rows').children.length===0;
      document.querySelector('[data-source="All"]').click(); queries[2].click();
      checks.multidomain=el('domain-architecture').querySelectorAll('[data-hit]').length>2;
      return checks;
    })())""")
    )
    imported = parse_queries(api._window.run_js("document.getElementById('quick-input').value"))
    checks["genbank_file_converted"] = api._window.run_js(
        "document.getElementById('quick-import-status').textContent"
    ) == "Loaded public quick.gbk." and [p.sequence for p in imported] == [
        p.sequence for p in proteins
    ]
    state = api._job.snapshot(include_quick=True)
    return {
        "ok": all(checks.values()),
        "checks": checks,
        "visual_wall_seconds": wall,
        "timings": state["quick"]["timings"],
        "query_names": [q["name"] for q in state["quick"]["queries"]],
    }


def resident_smoke(api, text):
    """Verify measured-workload routing and worker reuse through the native DOM."""
    from .search.resident import ResidentHmmer

    session = api._job.search_session
    adapter = session.engine.hmm_adapter.adapter
    if not isinstance(adapter, ResidentHmmer):
        return {"ok": True, "status": "not qualified; canonical CLI active"}
    checks = {"small_search_did_not_load_models": adapter.process is None}
    base = parse_queries(text)
    originals = [
        type(base[0])(
            str(i),
            f"public synthetic control {i}",
            base[i % len(base)].sequence + "A" * (i // len(base) + 1),
        )
        for i in range(20)
    ]
    warm = [type(p)(p.id, p.header, p.sequence + "C") for p in originals]
    trials = []
    for label, proteins in (("first", originals), ("warm_new", warm), ("cached", originals)):
        revision = api._job.quick_revision
        fasta = "".join(f">{p.header}\n{p.sequence}\n" for p in proteins)
        api._window.run_js(
            "document.getElementById('open-quick').click();"
            "if(document.getElementById('quick-input-panel').hidden)"
            "document.getElementById('quick-edit').click();"
            f"document.getElementById('quick-input').value={json.dumps(fasta)};"
            "document.getElementById('quick-input').dispatchEvent(new Event('input'));true"
        )
        for _ in range(100):
            if api._window.run_js("!document.getElementById('quick-search').disabled"):
                break
            time.sleep(0.05)
        started = time.monotonic()
        api._window.run_js("document.getElementById('quick-search').click();true")
        for _ in range(1200):
            state = api._job.snapshot(include_quick=True)
            if state["error"]:
                raise RuntimeError(state["error"])
            if state["quick_revision"] > revision and not state["busy"]:
                break
            time.sleep(0.05)
        else:
            raise RuntimeError("Resident native Quick Search did not complete")
        adapter = session.engine.hmm_adapter.adapter
        checks[f"{label}_complete"] = (
            len(state["quick"]["queries"]) == 20 and not state["quick"]["warnings"]
        )
        if label == "cached":
            checks["cached_evidence_replayed"] = all(
                q["cached"] for p in state["quick"]["provenance"].values() for q in p["queries"]
            )
        trials.append(
            {
                "label": label,
                "pid": adapter.process.pid,
                "model_load_seconds": adapter.last_call["load_seconds"],
                "wall_seconds": time.monotonic() - started,
            }
        )
    checks["one_worker_reused"] = len({t["pid"] for t in trials}) == 1
    checks["warm_models_reused"] = trials[1]["model_load_seconds"] == 0
    return {"ok": all(checks.values()), "status": "passed", "checks": checks, "trials": trials}


def stop_smoke(api, directory, query):
    """Opt-in installed acceptance: actual engines, DOM confirmation, cleanup, rerun."""
    import psutil

    def wait_for(check):
        deadline = time.monotonic() + 20
        while time.monotonic() < deadline:
            value = check()
            if value:
                return value
            time.sleep(0.05)
        raise RuntimeError("Installed Stop acceptance timed out")

    def js(code):
        return api._window.run_js(code)

    def start(inputs, output):
        wait_for(lambda: js("!document.getElementById('run').disabled"))
        js(
            "document.getElementById('open-files').click(); "
            f"setInputs({json.dumps(list(map(str, inputs)))}); "
            f"document.getElementById('output').value={json.dumps(str(output))}; "
            "document.getElementById('run').click(); true"
        )

    def engines():
        found = {}
        for process in psutil.Process().children(recursive=True):
            try:
                if process.name() in {"hmmscan", "rpsblast"}:
                    found[process.name()] = process
            except psutil.NoSuchProcess:
                pass
        return found if len(found) == 2 else None

    slow, pending = directory / "public repeated control.faa", directory / "not started.faa"
    slow.write_text(query.read_text() * 1000)
    pending.write_text(query.read_text())
    original = slow.read_bytes()
    output = directory / "stopped outputs"
    checks = {}
    try:
        start([slow, pending], output)
        active = wait_for(engines)
        directories = []
        for name, process in active.items():
            args = process.cmdline()
            path = Path(args[-1] if name == "hmmscan" else args[args.index("-query") + 1])
            directories.append(path.parent.parent)
            if name == "rpsblast":
                view = Path(args[args.index("-db") + 1].strip('"')).parent
                if view.name.startswith("kannotate-cdd-"):
                    directories.append(view)
        wait_for(lambda: js("!document.getElementById('stop').disabled"))
        js("document.getElementById('stop').click(); true")
        checks["confirmation_visible"] = js("document.getElementById('stop-dialog').open")
        checks["preservation_warning"] = js(
            "document.getElementById('stop-warning').textContent.includes('original inputs')"
        )
        checks["safe_default"] = js("document.activeElement.value === 'continue'")
        js("document.querySelector('#stop-dialog button[value=continue]').click(); true")
        checks["decline_keeps_running"] = api._job.snapshot()["busy"] and all(
            p.is_running() for p in active.values()
        )
        js("document.getElementById('stop').click(); true")
        started = time.monotonic()
        js("document.querySelector('#stop-dialog button[value=stop]').click(); true")
        wait_for(lambda: not api._job.snapshot()["busy"])
        elapsed = time.monotonic() - started
        state = api._job.snapshot()
        checks["stopped_without_error"] = state["stopped"] and not state["error"]
        checks["both_engines_exited"] = all(not p.is_running() for p in active.values())
        checks["temporary_data_removed"] = all(not p.exists() for p in directories)
        checks["no_interrupted_or_pending_output"] = not list(output.iterdir())
        checks["input_preserved"] = slow.read_bytes() == original
        start([query], directory / "after stop outputs")
        wait_for(lambda: not api._job.snapshot()["busy"] and api._job.snapshot()["results"])
        state = api._job.snapshot()
        checks["successful_rerun"] = (
            not state["error"]
            and not state["stopped"]
            and len(state["results"]) == 1
            and not state["results"][0]["error"]
            and not state["results"][0]["warnings"]
            and Path(state["results"][0]["output"]).is_file()
        )
        return {"ok": all(checks.values()), "checks": checks, "stop_seconds": elapsed}
    finally:
        state = api._job.snapshot()
        if state["busy"] and state["can_stop"]:
            api._job.stop(state["job_id"])
            wait_for(lambda: not api._job.snapshot()["busy"])


def file_run_smoke(api, inputs, output, mode, *, checkpoint=False):
    """Exercise format submission/result publication through the installed DOM."""
    for _ in range(100):
        if api._window.run_js("!document.getElementById('run').disabled"):
            break
        time.sleep(0.1)
    else:
        raise RuntimeError("Installed Run annotation control did not become ready")
    previous = api._job.snapshot()["job_id"]
    api._window.run_js(
        "document.getElementById('open-files').click(); "
        f"setInputs({json.dumps([str(p) for p in inputs])}); "
        f"document.getElementById('output').value={json.dumps(str(output))}; "
        "document.getElementById('output').dispatchEvent(new Event('input')); "
        f"document.getElementById('genbank-output').value={json.dumps(mode)}; "
        f"document.getElementById('checkpoint').checked={json.dumps(checkpoint)}; "
        "document.getElementById('run').click(); true"
    )
    for _ in range(1200):
        state = api._job.snapshot()
        if state["job_id"] != previous and not state["busy"]:
            if state["error"] or len(state["results"]) != len(inputs):
                raise RuntimeError(f"Installed file run failed: {state}")
            for result in state["results"]:
                if result["error"] or result["warnings"] or not result["outputs"]:
                    raise RuntimeError(f"Installed file output failed: {result}")
                if not all(Path(p).is_file() for p in result["outputs"]):
                    raise RuntimeError("Installed file output is missing")
            api._window.run_js("refresh(); true")
            return state
        time.sleep(0.1)
    raise RuntimeError("Installed file run did not complete")


def checkpoint_smoke(api, directory, protein):
    """Actual checkbox/Stop/resume flow, with public/synthetic disposable proteins."""
    from .search.identity import digest

    source = directory / "public checkpoint controls.faa"
    source.write_text(
        f">public EF Tu\n{protein}\n"
        + ">empty\n\n" * 511
        + f">synthetic repeated EF Tu\n{protein * 20}\n"
        + ">empty\n\n" * 512
    )
    original = source.read_bytes()
    output = directory / "checkpoint outputs"
    previous = api._job.snapshot()["job_id"]
    for _ in range(100):
        if api._window.run_js("!document.getElementById('run').disabled"):
            break
        time.sleep(0.05)
    api._window.run_js(
        "document.getElementById('open-files').click();"
        f"setInputs({json.dumps([str(source)])});"
        f"document.getElementById('output').value={json.dumps(str(output))};"
        "document.getElementById('output').dispatchEvent(new Event('input'));"
        "document.getElementById('genbank-output').value='csv';"
        "document.getElementById('checkpoint').checked=true;"
        "document.getElementById('run').click();true"
    )
    for _ in range(1200):
        state = api._job.snapshot()
        checkpoint = Path(state["checkpoint_dir"]) if state.get("checkpoint_dir") else None
        receipts = list(checkpoint.glob("file-*/group-*/receipt.json")) if checkpoint else []
        if receipts and state["busy"]:
            break
        if state["job_id"] != previous and not state["busy"]:
            raise RuntimeError("Checkpoint control finished before the Stop test")
        time.sleep(0.05)
    else:
        raise RuntimeError("Checkpoint control did not complete its first group")
    saved = {str(p): digest(p) for p in receipts}
    api._window.run_js("document.getElementById('stop').click();true")
    warning = api._window.run_js("document.getElementById('stop-resume-warning').textContent")
    api._window.run_js("document.querySelector('#stop-dialog button[value=stop]').click();true")
    for _ in range(200):
        state = api._job.snapshot()
        if not state["busy"]:
            break
        time.sleep(0.05)
    else:
        raise RuntimeError("Checkpoint control did not stop")
    stopped = state["stopped"] and not state["error"] and not list(output.iterdir())
    resumed = file_run_smoke(api, [source], output, "csv", checkpoint=True)
    checks = {
        "resume_warning": "Completed checkpoints will be kept" in warning,
        "stopped_without_partial_output": stopped,
        "completed_groups_preserved": all(digest(p) == value for p, value in saved.items()),
        "source_preserved": source.read_bytes() == original,
        "all_rows_resumed": resumed["results"][0]["proteins"] == 1025,
    }
    return {"ok": all(checks.values()), "checks": checks}


def smoke(api, directory, quick_input=None):
    """Opt-in acceptance through actual desktop DOM controls and real engines."""
    directory.mkdir(parents=True, exist_ok=True)
    query = directory / "public EF Tu.faa"
    query.write_text(files("kannotate.resources").joinpath("smoke.faa").read_text())
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord

    from .io.fasta import read_fasta

    protein = read_fasta(query)[0].sequence
    record = SeqRecord(
        Seq("N" * (len(protein) * 3)),
        id="public_control",
        description="Synthetic CDS carrying public UniProt P0CE47 translation",
    )
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(
            SimpleLocation(0, len(record), strand=1),
            type="CDS",
            qualifiers={
                "translation": [protein],
                "locus_tag": ["control"],
                "product": ["old test product"],
            },
        )
    ]
    genbank = directory / "public control.gbk"
    SeqIO.write([record], genbank, "genbank")
    output = directory / "annotation outputs"
    report = {"ok": False}
    try:
        for _ in range(100):
            if api._window.run_js("document.readyState === 'complete' && !!window.setInputs"):
                break
            time.sleep(0.1)
        checks = {
            "genbank_default": api._window.run_js(
                "document.getElementById('genbank-output').value === 'genbank'"
            ),
            "first_input_folder": api._window.run_js(
                f"setInputs({json.dumps([str(query)])}); "
                f"document.getElementById('output').value === {json.dumps(str(directory))}"
            ),
            "folder_follows_inputs": api._window.run_js(
                f"setInputs({json.dumps([str(output / 'other.faa')])}); "
                f"document.getElementById('output').value === {json.dumps(str(output))}"
            ),
            "manual_folder_persists": api._window.run_js(
                f"document.getElementById('output').value={json.dumps(str(output))}; "
                "document.getElementById('output').dispatchEvent(new Event('input')); "
                f"setInputs({json.dumps([str(query), str(genbank)])}); "
                f"document.getElementById('output').value === {json.dumps(str(output))}"
            ),
        }
        state = file_run_smoke(api, [query, genbank], output, "both")
        report = {
            "ok": all(checks.values()),
            "state": state,
            "version": __version__,
            "format_checks": checks,
            "genbank_modes": {},
        }
        for mode, suffix in (("genbank", ".gbk"), ("csv", ".csv")):
            mode_state = file_run_smoke(api, [genbank], directory / f"{mode} only", mode)
            paths = mode_state["results"][0]["outputs"]
            report["genbank_modes"][mode] = mode_state
            report["ok"] = report["ok"] and len(paths) == 1 and Path(paths[0]).suffix == suffix
        if report["ok"] and quick_input:
            # A fresh cache per acceptance run: resident routing must load models for
            # unseen queries, and the user's persistent search cache stays untouched.
            api._job.quick_cache_dir = directory / "quick-search-cache"
            report["quick_search"] = quick_smoke(api, quick_input.read_text())
            report["ok"] = report["quick_search"]["ok"]
            if report["ok"]:
                report["resident"] = resident_smoke(api, quick_input.read_text())
                report["ok"] = report["resident"]["ok"]
        if report["ok"]:
            report["stop"] = stop_smoke(api, directory, query)
            report["ok"] = report["stop"]["ok"]
        if report["ok"]:
            report["checkpoints"] = checkpoint_smoke(api, directory, protein)
            report["ok"] = report["checkpoints"]["ok"]
    except Exception as exc:
        report["ok"] = False
        report["error"] = str(exc)
    (directory / "report.json").write_text(json.dumps(report, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Kannotate local annotation desktop")
    parser.add_argument("--browser", action="store_true")
    parser.add_argument("--no-open", action="store_true")
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--smoke-dir", type=Path)
    parser.add_argument("--bridge-check-dir", type=Path)
    parser.add_argument("--quick-smoke-input", type=Path)
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()
    if args.quick_smoke_input and not args.smoke_dir:
        parser.error("--quick-smoke-input requires --smoke-dir")
    job = AnnotationJob()
    api = DesktopApi(job)
    app = create_app(job, desktop=not args.browser)

    def started():
        api._window.events.closing += api.can_close
        if args.bridge_check_dir or args.smoke_dir:
            checked = bridge_contract(api, (args.bridge_check_dir or args.smoke_dir).resolve())
            if not checked["ok"]:
                return
        if args.smoke_dir:
            smoke(api, args.smoke_dir.resolve(), args.quick_smoke_input)

    try:
        host(
            app,
            f"Kannotate · {__version__}",
            api=api,
            exposed=(api.choose_inputs, api.choose_output, api.reveal),
            browser=args.browser,
            no_open=args.no_open,
            port=args.port,
            on_started=started,
        )
    finally:
        job.search_session.close()


if __name__ == "__main__":
    main()
