"""macOS, not OpenSSL's certificate file, decides HTTPS trust for the urllib clients. Offline."""

import ast
import ssl
import subprocess
import threading
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import pytest

from kbiosuite import tls

ROOT = Path(__file__).resolve().parents[2]
OPENSSL = "/usr/bin/openssl"
CONFIG = """[req]
distinguished_name = name
prompt = no
[name]
CN = {name}
[ca]
basicConstraints = critical, CA:TRUE
keyUsage = critical, keyCertSign
subjectKeyIdentifier = hash
[leaf]
basicConstraints = critical, CA:FALSE
extendedKeyUsage = serverAuth
subjectAltName = IP:127.0.0.1
subjectKeyIdentifier = hash
authorityKeyIdentifier = keyid
"""


def test_context_returns_fresh_truststore_contexts_with_peer_and_hostname_checks():
    import truststore

    made = tls.context()
    assert isinstance(made, truststore.SSLContext)
    assert made.verify_mode == ssl.CERT_REQUIRED and made.check_hostname
    assert made.minimum_version >= ssl.TLSVersion.TLSv1_2
    assert tls.context() is not made  # A second factory call returns a distinct context.


@pytest.fixture
def private_authority(tmp_path):
    """A loopback HTTPS server whose chain ends at a root only a certificate file knows."""
    if not Path(OPENSSL).is_file():
        pytest.skip("macOS's openssl is needed to make the throwaway certificates")
    # LibreSSL ignores -subj with prompt = no, so each subject gets its own configuration.
    root, leaf = tmp_path / "root.cnf", tmp_path / "leaf.cnf"
    root.write_text(CONFIG.format(name="KBioSuite test root"))
    leaf.write_text(CONFIG.format(name="127.0.0.1"))

    def openssl(*arguments):
        subprocess.run([OPENSSL, *arguments], cwd=tmp_path, check=True, capture_output=True)

    openssl("req", "-x509", "-newkey", "rsa:2048", "-nodes", "-days", "2", "-config", root,
            "-extensions", "ca", "-keyout", "ca.key", "-out", "ca.pem")
    openssl("req", "-new", "-newkey", "rsa:2048", "-nodes", "-config", leaf,
            "-keyout", "leaf.key", "-out", "leaf.csr")
    openssl("x509", "-req", "-in", "leaf.csr", "-CA", "ca.pem", "-CAkey", "ca.key",
            "-CAcreateserial", "-days", "2", "-extfile", leaf, "-extensions", "leaf",
            "-out", "leaf.pem")

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.send_header("Content-Length", "2")
            self.end_headers()
            self.wfile.write(b"ok")

        def log_message(self, *arguments):
            pass

    class Server(HTTPServer):
        def handle_error(self, request, client_address):
            pass  # a client that rejects the certificate hangs up mid-handshake

    server = Server(("127.0.0.1", 0), Handler)
    serving = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
    serving.load_cert_chain(tmp_path / "leaf.pem", tmp_path / "leaf.key")
    server.socket = serving.wrap_socket(server.socket, server_side=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"https://127.0.0.1:{server.server_port}/", tmp_path / "ca.pem"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(5)


def test_ssl_cert_file_does_not_make_a_file_only_root_trusted_by_the_suite(
    private_authority, monkeypatch
):
    """With SSL_CERT_FILE set to the fixture CA, default urllib accepts the server and both
    suite requests fail."""
    url, authority = private_authority
    monkeypatch.setenv("SSL_CERT_FILE", str(authority))
    by_file = urllib.request.build_opener(
        urllib.request.ProxyHandler({}),
        urllib.request.HTTPSHandler(context=ssl.create_default_context()),
    )
    with by_file.open(url, timeout=10) as response:
        assert response.read() == b"ok"  # control: the server and its chain are sound
    with pytest.raises(urllib.error.URLError) as refused:
        tls.urlopen(url, timeout=10, direct=True)
    assert isinstance(refused.value.reason, ssl.SSLCertVerificationError)
    with pytest.raises(urllib.error.URLError):
        tls.head(url)


def test_qualification_requests_refuse_redirects():
    """A direct HEAD succeeds on the control path and surfaces the redirect as HTTP 302."""

    class Handler(BaseHTTPRequestHandler):
        def do_HEAD(self):
            if self.path == "/moved":
                self.send_response(302)
                self.send_header("Location", "/elsewhere")
            else:
                self.send_response(200)
            self.send_header("Content-Length", "0")
            self.end_headers()

        def log_message(self, *arguments):
            pass

    server = HTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        base = f"http://127.0.0.1:{server.server_port}"
        assert tls.head(base + "/elsewhere") == 200
        with pytest.raises(urllib.error.HTTPError) as refused:
            tls.head(base + "/moved")
        assert refused.value.code == 302
    finally:
        server.shutdown()
        server.server_close()
        thread.join(5)


def test_exported_python_sources_use_the_expected_urllib_and_ssl_call_forms():
    """Flag direct urllib opener or default-context names outside the TLS helper."""
    # Development scripts are not part of every checkout of this repository.
    sources = [*sorted((ROOT / "src").rglob("*.py")), *sorted((ROOT / "scripts").rglob("*.py"))]
    assert len(sources) > 50
    offenders = []
    for path in sources:
        if path == ROOT / "src/kbiosuite/tls.py":
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            names = []
            if isinstance(node, ast.Attribute):
                names = [node.attr]
            elif isinstance(node, ast.ImportFrom) and node.module in ("urllib.request", "ssl"):
                names = [alias.name for alias in node.names]
            if {"urlopen", "build_opener", "create_default_context"} & set(names) and not (
                isinstance(node, ast.Attribute)
                and isinstance(node.value, ast.Name)
                and node.value.id == "tls"
            ):
                offenders.append(f"{path.relative_to(ROOT)}:{node.lineno}")
    assert not offenders, offenders
