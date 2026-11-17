"""Loopback-only GUI, optionally hosted in a native pywebview desktop window."""

import argparse
import io
import json
import secrets
import threading
import time
import webbrowser
import zipfile
from datetime import datetime
from pathlib import Path

from flask import Flask, jsonify, request, send_file

from kbiosuite.desktop import loopback_server

from . import __version__
from .config import Config
from .export import export_bundle
from .lifecycle import validate_rid
from .preview import context_preview
from .service import RunManager

STATIC = Path(__file__).parent / "static"


def create_app(manager: RunManager, desktop_api=None) -> Flask:
    app = Flask(__name__, static_folder=str(STATIC), static_url_path="/static")
    app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024
    token = secrets.token_urlsafe(32)
    app.config["KBLAST_TOKEN"] = token

    @app.before_request
    def local_only():
        host = request.host.split(":", 1)[0]
        if host not in {"127.0.0.1", "localhost"}:
            return jsonify(error="Kblast accepts loopback hosts only."), 403
        if request.path.startswith("/api/"):
            if not secrets.compare_digest(request.headers.get("X-Kblast-Token", ""), token):
                return jsonify(error="Reload Kblast to reconnect to this local session."), 403
            origin = request.headers.get("Origin")
            if origin and origin != request.host_url.rstrip("/"):
                return jsonify(error="Cross-origin requests are disabled."), 403

    @app.after_request
    def headers(response):
        response.headers["Cache-Control"] = "no-store"
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["Referrer-Policy"] = "no-referrer"
        response.headers["Content-Security-Policy"] = (
            "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; "
            "img-src 'self' data:; connect-src 'self'; frame-ancestors 'none'; object-src 'none'"
        )
        return response

    @app.errorhandler(ValueError)
    def bad_input(exc):
        return jsonify(error=str(exc)), 400

    @app.errorhandler(FileNotFoundError)
    def missing(_):
        return jsonify(error="The saved run or context could not be found."), 404

    @app.errorhandler(KeyError)
    @app.errorhandler(TypeError)
    def malformed(_):
        return jsonify(error="The request is missing a required field or has an invalid type."), 400

    @app.get("/")
    def index():
        return (
            (STATIC / "index.html")
            .read_text()
            .replace("__KBLAST_TOKEN__", token)
            .replace("__VERSION__", __version__)
        )

    @app.get("/api/config")
    def configuration():
        return jsonify(
            api_key_configured=bool(manager.config.api_key),
            email=manager.config.email,
            data_directory=str(manager.config.data_dir),
            desktop=desktop_api is not None,
        )

    @app.get("/api/example")
    def example():
        return jsonify(sequence=(STATIC / "example.fasta").read_text(), accession="CAC44203.1")

    @app.post("/api/open-protein-blast")
    def open_protein_blast():
        # Fixed official page only: no unpublished query in a URL or automatic submission.
        webbrowser.open("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins")
        return jsonify(ok=True)

    @app.get("/api/runs")
    def runs():
        return jsonify(
            runs=manager.store.list_runs(),
            active_id=manager.active_id if manager.is_active() else None,
        )

    @app.post("/api/runs")
    def start():
        data = request.get_json()
        run = manager.start(
            data["sequence"],
            data["top_n"],
            data["flank"],
            data.get("email", ""),
            data.get("rid", ""),
            data.get("rid_source", "external"),
        )
        return jsonify(id=run.id), 201

    @app.get("/api/runs/<run_id>")
    def run_state(run_id):
        return jsonify(manager.snapshot(run_id))

    @app.post("/api/runs/<run_id>/cancel")
    def cancel(run_id):
        manager.cancel(run_id)
        return jsonify(ok=True)

    @app.post("/api/runs/<run_id>/resume")
    def resume(run_id):
        data = request.get_json() or {}
        manager.resume(
            run_id, data.get("email", ""), data.get("rid", ""), data.get("rid_source", "external")
        )
        return jsonify(ok=True)

    @app.get("/api/runs/<run_id>/contexts/<occurrence_id>")
    def preview(run_id, occurrence_id):
        run = manager.store.load(run_id)
        occurrence = next(
            (o for o in run.occurrences if o.id == occurrence_id and o.status == "success"), None
        )
        if occurrence is None:
            raise ValueError("This context is not available for preview.")
        return jsonify(context_preview(manager.store.directory(run_id), occurrence))

    @app.post("/api/runs/<run_id>/export")
    def export(run_id):
        run = manager.store.load(run_id)
        data = request.get_json() or {}
        selected = data.get("selected")
        if selected is not None and not isinstance(selected, list):
            raise ValueError("Export selection must be a list of context identifiers.")
        bundle = export_bundle(manager.store.directory(run_id), run, selected)
        return send_file(
            io.BytesIO(bundle),
            mimetype="application/zip",
            as_attachment=True,
            download_name=f"Kblast_{run.created_at[:10]}_{run.id[:8]}.zip",
        )

    @app.post("/api/runs/<run_id>/export-desktop")
    def export_desktop(run_id):
        if desktop_api is None:
            raise ValueError("Desktop export is available in the native application window.")
        data = request.get_json() or {}
        return jsonify(desktop_api.export(run_id, data.get("selected")))

    return app


class DesktopApi:
    def __init__(self, manager: RunManager):
        self._manager, self._window = manager, None

    def export(self, run_id, selected=None):
        import webview

        folders = self._window.create_file_dialog(webview.FileDialog.FOLDER)
        if not folders:
            return {"cancelled": True}
        try:
            run = self._manager.store.load(run_id)
            bundle = export_bundle(self._manager.store.directory(run_id), run, selected)
            stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            target = Path(folders[0]) / f"Kblast_{stamp}_{run.id[:8]}"
            target.mkdir(exist_ok=False)
            with zipfile.ZipFile(io.BytesIO(bundle)) as archive:
                archive.extractall(target)
            return {"path": str(target)}
        except (ValueError, OSError) as exc:
            return {"error": str(exc)}


class InstanceLock:
    """Avoid concurrent GUI processes mutating the same run store."""

    def __init__(self, root: Path):
        root.mkdir(parents=True, exist_ok=True, mode=0o700)
        self.handle = (root / ".instance.lock").open("a+b")
        try:
            if __import__("os").name == "nt":
                import msvcrt

                self.handle.write(b"0")
                self.handle.flush()
                self.handle.seek(0)
                msvcrt.locking(self.handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl

                fcntl.flock(self.handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            self.handle.close()
            raise ValueError("Kblast is already using this data directory.") from None

    def close(self):
        self.handle.close()


def qualify_rid(window, manager, rid, output):
    """Opt-in fixed public-query gate through the installed native GUI controls."""
    from .storage import atomic_write

    query = (STATIC / "example.fasta").read_text()
    before = {r["id"] for r in manager.store.list_runs()}
    window.run_js(
        f"document.getElementById('sequence').value={json.dumps(query)};"
        "document.getElementById('top-n').value='1';"
        "document.getElementById('flank').value='20000';"
        "document.getElementById('sequence').dispatchEvent(new Event('input'));"
        "document.getElementById('use-rid').click();"
        f"document.getElementById('existing-rid').value={json.dumps(rid)};"
        "document.getElementById('rid-source').value='browser';"
        "document.getElementById('rid-form').requestSubmit();"
    )
    deadline = time.monotonic() + 1200
    run = None
    while time.monotonic() < deadline:
        fresh = [r for r in manager.store.list_runs() if r["id"] not in before]
        if fresh:
            run = manager.store.load(fresh[0]["id"])
            if not manager.is_active(run.id):
                break
        time.sleep(0.2)
    if run is None:
        raise ValueError("Installed RID form did not create the public qualification run.")
    if manager.is_active(run.id):
        manager.cancel(run.id)
        raise ValueError("Installed RID qualification exceeded its bounded observation window.")
    directory = manager.store.directory(run.id)
    entries = [json.loads(line) for line in (directory / "requests.jsonl").read_text().splitlines()]
    puts = sum(
        e.get("service") == "blast"
        and e.get("parameters", {}).get("CMD") == "Put"
        and e.get("outcome") == "started"
        for e in entries
    )
    successful = [o for o in run.occurrences if o.status == "success"]
    evidence = {
        "status": "failed",
        "public_query": "CAC44203.1",
        "rid": run.rid,
        "run_id": run.id,
        "run_status": run.status,
        "source": run.rid_source,
        "counts": run.to_dict()["counts"],
        "blast_put_requests": puts,
        "metadata": run.blast_metadata,
        "validation": run.rid_validation,
        "error": run.error,
    }
    if run.status == "success" and successful and puts == 0 and run.rid == rid:
        # Exercise the same exporter used by the GUI; destination is isolated evidence.
        bundle = export_bundle(directory, run)
        atomic_write(output / "rid-import-export.zip", bundle)
        with zipfile.ZipFile(io.BytesIO(bundle)) as archive:
            for o in successful:
                from .export import context_filename

                if archive.read(context_filename(o)) != (directory / o.context_file).read_bytes():
                    raise ValueError("Installed RID export changed source GenBank bytes.")
            provenance = json.loads(archive.read("provenance/run.json"))
            if provenance["rid_source"] != "browser":
                raise ValueError("Installed export lost browser RID provenance.")
        evidence["preview_features"] = len(context_preview(directory, successful[0])["features"])
        evidence["export_source_bytes_verified"] = len(successful)
        evidence["status"] = "passed"
        window.run_js("document.getElementById('tab-contexts').click();")
    return evidence


def main():
    parser = argparse.ArgumentParser(
        description=f"Kblast {__version__} · protein-to-genomic neighborhoods"
    )
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument(
        "--browser", action="store_true", help="Use a local browser instead of a desktop window"
    )
    parser.add_argument(
        "--no-open", action="store_true", help="Start browser server without opening a browser"
    )
    parser.add_argument(
        "--port", type=int,
        help="Loopback port (default: free port for desktop, 8765 for --browser; 0 selects free)",
    )
    parser.add_argument(
        "--data-dir", type=Path, help="Local run store (default suite Kblast run store)"
    )
    parser.add_argument("--smoke-dir", type=Path)
    parser.add_argument(
        "--qualify-rid",
        help="Opt-in installed GUI qualification: browser RID for public CAC44203.1, no submission",
    )
    parser.add_argument(
        "--check-eutils", action="store_true", help="Opt-in harmless EInfo in desktop smoke"
    )
    args = parser.parse_args()
    if args.qualify_rid:
        if not args.smoke_dir or args.browser:
            parser.error("--qualify-rid requires native --smoke-dir and the public CAC44203.1 RID")
        validate_rid(args.qualify_rid, "browser")
    try:
        config = Config.discover(args.data_dir)
        lock = InstanceLock(config.data_dir)
    except (ValueError, OSError) as exc:
        parser.exit(1, f"{exc}\n")
    manager = RunManager(config)
    desktop_api = None
    if not args.browser:
        try:
            import webview
        except ImportError:
            lock.close()
            parser.exit(1, "Desktop dependency missing. Run uv sync, or use --browser.\n")
        desktop_api = DesktopApi(manager)
    port = args.port if args.port is not None else (8765 if args.browser else 0)
    server = loopback_server(create_app(manager, desktop_api), port)
    url = f"http://127.0.0.1:{server.server_port}"
    server_thread = threading.Thread(target=server.serve_forever, daemon=True)
    server_thread.start()
    print(f"Kblast: {url}", flush=True)
    try:
        if args.browser:
            if not args.no_open:
                webbrowser.open(url)
            server_thread.join()
        else:
            desktop_api._window = webview.create_window(
                f"Kblast {__version__} · Gene neighborhoods",
                url,
                width=1360,
                height=900,
                min_size=(960, 680),
                background_color="#f4f6f5",
                text_select=True,
            )
            desktop_api._window.events.closing += manager.shutdown

            def smoke():
                if not args.smoke_dir:
                    return
                report = {
                    "ok": False,
                    "version": __version__,
                    "eutils_key_configured": bool(config.api_key),
                }
                try:
                    for _ in range(100):
                        state = desktop_api._window.run_js(
                            "document.readyState === 'complete' && "
                            "!!document.querySelector('meta[name=\"kblast-token\"]')"
                        )
                        if state:
                            report["ok"] = True
                            report["title"] = desktop_api._window.run_js("document.title")
                            report["controls"] = json.loads(
                                desktop_api._window.run_js(
                                    "JSON.stringify(Array.from(document.querySelectorAll('button'))"
                                    ".map(b=>b.textContent.trim()))"
                                )
                            )
                            break
                        time.sleep(0.1)
                    if args.qualify_rid and report["ok"]:
                        report["rid_import"] = qualify_rid(
                            desktop_api._window, manager, args.qualify_rid, args.smoke_dir
                        )
                        report["ok"] = report["rid_import"]["status"] == "passed"
                    if args.check_eutils:
                        from .network import RateLimiter, Transport

                        args.smoke_dir.mkdir(parents=True, exist_ok=True)
                        log = args.smoke_dir / "eutils-requests.jsonl"
                        http = Transport(
                            config,
                            threading.Event(),
                            log,
                            RateLimiter(0.11 if config.api_key else 0.36),
                            RateLimiter(10),
                        )
                        try:
                            result = http.request("eutils", "einfo.fcgi", {"db": "protein"})
                        finally:
                            http.close()
                        received = [
                            json.loads(line)
                            for line in log.read_text().splitlines()
                            if json.loads(line).get("outcome") == "received"
                        ]
                        report["eutils_check"] = {
                            "ok": "<eInfoResult>" in result and "<ERROR>" not in result,
                            "key_used_on_wire": received[-1].get("api_key_used", False),
                        }
                        report["ok"] = report["ok"] and report["eutils_check"]["ok"]
                except Exception as exc:
                    report["error"] = str(exc)
                args.smoke_dir.mkdir(parents=True, exist_ok=True)
                (args.smoke_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")

            webview.start(func=smoke, debug=False)
    except KeyboardInterrupt:
        pass
    finally:
        manager.shutdown()
        server.shutdown()
        server.server_close()
        lock.close()


if __name__ == "__main__":
    main()
