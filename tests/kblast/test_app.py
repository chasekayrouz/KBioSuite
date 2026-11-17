import io
import sys
import zipfile
from types import SimpleNamespace
from urllib.request import urlopen

import pytest

from kblast.app import DesktopApi, create_app
from kblast.config import Config
from kblast.service import RunManager


@pytest.mark.parametrize(
    ("options", "expected_port"),
    [([], 0), (["--browser"], 8765), (["--port", "9123"], 9123),
     (["--browser", "--port", "9123"], 9123), (["--port", "0"], 0),
     (["--browser", "--port", "0"], 0)],
)
def test_launch_port_selection_and_local_page(tmp_path, monkeypatch, options, expected_port):
    from kblast import app

    monkeypatch.setattr(sys, "argv", ["kblast", *options])
    monkeypatch.setattr(
        app.Config, "discover", lambda _: Config(data_dir=tmp_path / "runs")
    )
    real_server = app.loopback_server
    servers = []

    def bind(application, port):
        assert port == expected_port
        # Keep fixed-port CLI cases independent of other local services.
        server = real_server(application, 0)
        servers.append(server)
        return server

    monkeypatch.setattr(app, "loopback_server", bind)
    opened = []

    def inspect_page(url):
        with urlopen(url, timeout=2) as response:
            assert response.status == 200
            assert b'kblast-token' in response.read()
        opened.append(url)

    class ClosingEvent:
        def __iadd__(self, callback):
            assert callable(callback)
            return self

    def create_window(title, url, **kwargs):
        inspect_page(url)
        return SimpleNamespace(events=SimpleNamespace(closing=ClosingEvent()))

    def browser_open(url):
        inspect_page(url)
        raise KeyboardInterrupt  # End the CLI's otherwise blocking server loop.

    monkeypatch.setitem(sys.modules, "webview", SimpleNamespace(
        create_window=create_window, start=lambda **kwargs: None,
    ))
    monkeypatch.setattr(app.webbrowser, "open", browser_open)
    app.main()
    assert len(opened) == 1
    assert servers[0].socket.fileno() == -1
    # A normal close must also release the run-store lock for the next launch.
    lock = app.InstanceLock(tmp_path / "runs")
    lock.close()


def test_local_api_token_host_checks_and_no_secret_in_gui(tmp_path):
    manager = RunManager(Config(api_key="private_fixture_key", data_dir=tmp_path / "runs"))
    app = create_app(manager)
    client = app.test_client()
    assert client.get("/", base_url="http://localhost").status_code == 200
    assert client.get("/api/config").status_code == 403
    headers = {"X-Kblast-Token": app.config["KBLAST_TOKEN"]}
    response = client.get("/api/config", headers=headers)
    assert response.status_code == 200
    assert response.json["api_key_configured"] is True
    assert b"private_fixture_key" not in response.data
    assert (
        client.get("/api/config", headers=headers, base_url="http://malicious.example").status_code
        == 403
    )
    assert (
        client.post(
            "/api/runs", headers={**headers, "Origin": "https://foreign.example"}, json={}
        ).status_code
        == 403
    )


def test_gui_saved_run_preview_and_selected_export(run_fixture):
    store, run, _, pipeline = run_fixture
    pipeline.execute()
    app = create_app(RunManager(Config(data_dir=store.root)))
    client = app.test_client()
    headers = {"X-Kblast-Token": app.config["KBLAST_TOKEN"]}
    data = client.get(f"/api/runs/{run.id}", headers=headers).json
    assert data["counts"]["contexts"] == 3
    assert "sequence" not in data["hits"][0] and data["hits"][0]["hsp_count"] == 2
    occurrence = run.occurrences[2]
    preview = client.get(f"/api/runs/{run.id}/contexts/{occurrence.id}", headers=headers)
    assert preview.status_code == 200
    assert preview.json["start"] == 800
    response = client.post(
        f"/api/runs/{run.id}/export", headers=headers, json={"selected": [occurrence.id]}
    )
    assert response.status_code == 200
    with zipfile.ZipFile(io.BytesIO(response.data)) as bundle:
        assert len([n for n in bundle.namelist() if n.startswith("contexts/")]) == 1
    invalid = client.post(f"/api/runs/{run.id}/export", headers=headers, json={"selected": "bad"})
    assert invalid.status_code == 400


def test_gui_rejects_invalid_input_before_remote_work(tmp_path):
    manager = RunManager(Config(data_dir=tmp_path / "runs"))
    app = create_app(manager)
    client = app.test_client()
    headers = {"X-Kblast-Token": app.config["KBLAST_TOKEN"]}
    response = client.post(
        "/api/runs",
        headers=headers,
        json={"sequence": ">a\nMA\n>b\nAA", "top_n": 10, "flank": 20_000},
    )
    assert response.status_code == 400
    assert not manager.is_active()
    assert manager.store.list_runs() == []


def test_restart_marks_active_checkpoint_interrupted(run_fixture):
    store, run, _, pipeline = run_fixture
    pipeline.execute()
    run.status = "running"
    run.rid = "SAVED-RID"
    run.hits[0].sequence_status = "running"
    run.occurrences[1].status = "running"
    store.save(run)
    manager = RunManager(Config(data_dir=store.root))
    saved = manager.store.load(run.id)
    assert saved.status == "interrupted"
    assert saved.rid == "SAVED-RID"
    assert saved.hits[0].sequence_status == "pending"
    assert saved.occurrences[1].status == "pending"
    assert saved.occurrences[0].status == "success"


def test_native_export_route_writes_selected_files_and_handles_cancel(
    run_fixture, tmp_path, monkeypatch
):
    store, run, _, pipeline = run_fixture
    pipeline.execute()
    manager = RunManager(Config(data_dir=store.root))
    api = DesktopApi(manager)
    destination = tmp_path / "exports"
    destination.mkdir()
    choices = [[str(destination)], None]
    api._window = SimpleNamespace(create_file_dialog=lambda _: choices.pop(0))
    monkeypatch.setitem(
        sys.modules, "webview", SimpleNamespace(FileDialog=SimpleNamespace(FOLDER=1))
    )
    app = create_app(manager, api)
    client = app.test_client()
    headers = {"X-Kblast-Token": app.config["KBLAST_TOKEN"]}
    assert client.get("/api/config", headers=headers).json["desktop"] is True
    url = f"/api/runs/{run.id}/export-desktop"
    response = client.post(url, headers=headers, json={"selected": [run.occurrences[2].id]})
    assert response.status_code == 200 and "path" in response.json
    from pathlib import Path

    target = Path(response.json["path"])
    assert len(list((target / "contexts").glob("*.gb"))) == 1
    assert (target / "metadata.tsv").is_file() and (target / "provenance/run.json").is_file()
    cancelled = client.post(url, headers=headers, json={})
    assert cancelled.json == {"cancelled": True}
    assert len(list(destination.iterdir())) == 1
