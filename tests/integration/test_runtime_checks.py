"""Opt-in frozen probes: offline controls, never real network or installed state."""

import json
import sys
from pathlib import Path

import pytest

from kbiosuite import runtime_checks as checks
from kmatch import setup


def test_normal_launch_does_not_probe(monkeypatch):
    def never(*args):
        pytest.fail("normal launch attempted qualification")

    monkeypatch.setattr(checks, "network_check", never)
    monkeypatch.setattr(setup, "install_uv", never)
    assert checks.dispatch_network_check("Kblast", "requests", []) is None
    assert setup.dispatch_uv_check([]) is None


def test_network_report_records_failure_and_preserves_previous_evidence(tmp_path, monkeypatch):
    def failed(stack):
        raise RuntimeError("certificate validation failed")

    monkeypatch.setattr(checks, "network_check", failed)
    report = tmp_path / "report.json"
    args = [checks.NETWORK_FLAG, str(report)]
    assert checks.dispatch_network_check("Kannotate", "urllib", args) == 1
    assert json.loads(report.read_text())["passed"] is False
    before = report.read_bytes()
    with pytest.raises(FileExistsError):
        checks.dispatch_network_check("Kannotate", "urllib", args)
    assert report.read_bytes() == before


@pytest.mark.parametrize("variable", checks.CERTIFICATE_OVERRIDES)
def test_network_qualification_rejects_certificate_overrides(monkeypatch, variable):
    monkeypatch.setenv(variable, "/example/certificates")
    with pytest.raises(RuntimeError, match="certificate override"):
        checks.network_check("urllib")


def test_uv_check_refuses_existing_home_before_network(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    original = home / "keep"
    original.write_text("preserved")
    with pytest.raises(RuntimeError, match="must be new"):
        setup.dispatch_uv_check(["--kmatch-uv-check", str(home), str(tmp_path / "report.json")])
    assert original.read_text() == "preserved"
    assert not (tmp_path / "report.json").exists()


def test_uv_check_runs_only_uv_and_records_exact_version(tmp_path, monkeypatch):
    calls = []

    def install(config, cancel, progress):
        calls.append(config)
        return config.home / "tools/uv"

    monkeypatch.setattr(setup, "install_uv", install)
    monkeypatch.setattr(setup, "managed_uv", lambda home: {"version": setup.UV_VERSION})
    monkeypatch.setattr(setup.subprocess, "check_output", lambda *a, **kw: f"uv {setup.UV_VERSION}")
    home, report = tmp_path / "new-home", tmp_path / "report.json"
    assert setup.dispatch_uv_check(["--kmatch-uv-check", str(home), str(report)]) == 0
    result = json.loads(report.read_text())
    assert result["passed"] and result["engine_setup_attempted"] is False
    assert len(calls) == 1 and calls[0].home == home


def test_third_party_stack_requires_the_entry_points_probe(monkeypatch):
    with pytest.raises(ValueError, match="Unsupported"):
        checks.network_check("requests")
    seen = []

    def head(url):
        seen.append(url)
        return 200

    assert checks.network_check("requests", head) == {
        "stack": "requests", "url": checks.PUBLIC_URL, "http_status": 200, "tls_verified": True,
    }
    assert seen == [checks.PUBLIC_URL]
    with pytest.raises(RuntimeError, match="HTTP 302"):
        checks.network_check("requests", lambda url: 302)


def test_dispatch_passes_the_probe_and_records_its_stack(tmp_path):
    report = tmp_path / "report.json"
    status = checks.dispatch_network_check(
        "Kblast", "requests", [checks.NETWORK_FLAG, str(report)], head=lambda url: 200
    )
    assert status == 0
    assert json.loads(report.read_text())["stack"] == "requests"


def test_requests_probe_ignores_environment_and_verifies(monkeypatch):
    from kbiosuite import requests_probe

    calls = []

    class Response:
        status_code = 200

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def raise_for_status(self):
            calls.append("checked")

    class Session:
        trust_env = True

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def head(self, url, **options):
            calls.append((url, self.trust_env, options))
            return Response()

    monkeypatch.setattr(requests_probe.requests, "Session", Session)
    assert requests_probe.head(checks.PUBLIC_URL) == 200
    assert calls == [
        (checks.PUBLIC_URL, False, {"timeout": 20, "verify": True, "allow_redirects": False}),
        "checked",
    ]


def test_kmatch_and_kannotate_probe_their_stdlib_client():
    packaging = Path(__file__).resolve().parents[2] / "packaging"
    assert 'dispatch_network_check("Kmatch", "urllib")' in (packaging / "Kmatch.py").read_text()
    assert 'dispatch_network_check("Kannotate", "urllib")' in (
        packaging / "Kannotate.py"
    ).read_text()


def test_only_kblast_entry_reaches_requests():
    """PyInstaller's own import graph (function-level imports included), without freezing."""
    pytest.importorskip("PyInstaller")
    from PyInstaller.lib.modulegraph.modulegraph import ModuleGraph

    root = Path(__file__).resolve().parents[2]

    def graph(name):
        found = ModuleGraph(path=[str(root / "src"), *sys.path], excludes=["tkinter", "pytest"])
        found.add_script(str(root / "packaging" / f"{name}.py"))
        return found

    for name in ("Kannotate", "Kmatch"):
        reached = graph(name)
        assert not [m for m in ("requests", "urllib3", "certifi") if reached.find_node(m)], name
    assert graph("Kblast").find_node("kbiosuite.requests_probe") is not None
