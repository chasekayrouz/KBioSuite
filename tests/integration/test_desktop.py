import http.client
import socket
import threading

from flask import Flask

from kbiosuite.desktop import loopback_server


def test_desktop_loopback_startup_does_not_require_dns(monkeypatch):
    def unavailable_dns(*args):
        raise AssertionError("Desktop startup tried to resolve a hostname")

    monkeypatch.setattr(socket, "getfqdn", unavailable_dns)
    app = Flask(__name__)
    app.add_url_rule("/", view_func=lambda: "ready")
    server = loopback_server(app)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        client = http.client.HTTPConnection("127.0.0.1", server.server_port, timeout=2)
        client.request("GET", "/")
        response = client.getresponse()
        assert response.status == 200
        assert response.read() == b"ready"
        assert server.server_name == "127.0.0.1"
        client.close()
    finally:
        server.shutdown()
        thread.join(timeout=2)
        server.server_close()
