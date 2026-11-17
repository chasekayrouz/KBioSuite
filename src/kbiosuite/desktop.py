"""Small loopback host for local desktop tools; no scientific or job ownership."""

import secrets
import socketserver
import threading
import webbrowser
from pathlib import Path

from flask import Flask, jsonify, request
from werkzeug.serving import ThreadedWSGIServer


class LoopbackServer(ThreadedWSGIServer):
    def server_bind(self):
        # HTTPServer normally resolves its hostname here. A local desktop window
        # must not depend on DNS being responsive before it can open or quit.
        socketserver.TCPServer.server_bind(self)
        self.server_name, self.server_port = self.server_address[:2]


def loopback_server(app, port=0):
    return LoopbackServer("127.0.0.1", port, app)


def local_app(name: str, static: Path, *, desktop=False):
    app = Flask(name, static_folder=str(static), static_url_path="/static")
    app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024
    token = secrets.token_urlsafe(32)

    @app.before_request
    def local_only():
        if request.host.split(":", 1)[0] not in {"127.0.0.1", "localhost"}:
            return jsonify(error="Only loopback hosts are accepted."), 403
        if request.path.startswith("/api/"):
            if not secrets.compare_digest(request.headers.get("X-Desktop-Token", ""), token):
                return jsonify(error="Reload this window to reconnect."), 403
            origin = request.headers.get("Origin")
            if origin and origin != request.host_url.rstrip("/"):
                return jsonify(error="Cross-origin requests are disabled."), 403

    @app.after_request
    def headers(response):
        response.headers["Cache-Control"] = "no-store"
        response.headers["X-Content-Type-Options"] = "nosniff"
        # pywebview 6.2.1 creates exposed API functions with new Function.
        # Permit that only in the native host; browser/debug pages keep strict CSP.
        scripts = "'self' 'unsafe-eval'" if desktop else "'self'"
        response.headers["Content-Security-Policy"] = (
            f"default-src 'self'; script-src {scripts}; style-src 'self'; "
            "connect-src 'self'; frame-ancestors 'none'; object-src 'none'"
        )
        return response

    @app.errorhandler(ValueError)
    @app.errorhandler(OSError)
    @app.errorhandler(RuntimeError)
    def error(exc):
        return jsonify(error=str(exc)), 400

    return app, token


def host(
    app, title, *, api=None, exposed=(), browser=False, no_open=False, port=0, on_started=None
):
    server = loopback_server(app, port)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    url = f"http://127.0.0.1:{server.server_port}"
    print(f"{title}: {url}", flush=True)
    try:
        if browser:
            if not no_open:
                webbrowser.open(url)
            thread.join()
        else:
            import webview

            api._window = webview.create_window(
                title,
                url,
                width=1000,
                height=780,
                min_size=(800, 650),
                background_color="#f4f6f5",
                text_select=True,
            )
            api._window.expose(*exposed)
            webview.start(func=on_started, debug=False)
    except KeyboardInterrupt:
        pass
    finally:
        server.shutdown()
