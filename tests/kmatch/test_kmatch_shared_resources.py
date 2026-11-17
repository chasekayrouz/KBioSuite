import json
import re
import shutil
from types import SimpleNamespace

import pytest

from kmatch import cli, gui
from kmatch.config import Config
from kmatch.process import Cancellation
from kmatch.setup import install_foldseek, install_pdb


def test_shared_pins_are_resolved_at_config_load_and_do_not_use_overrides(tmp_path, monkeypatch):
    home = tmp_path / "Kmatch"
    home.mkdir()
    (home / "config.json").write_text(
        json.dumps(
            {
                "shared_resources": "kmatch",
                "foldseek_override": "/old/tool",
                "pdb_prefix_override": "/old/pdb",
            }
        )
    )
    monkeypatch.setenv("KMATCH_HOME", str(home))
    monkeypatch.setenv("KMATCH_CACHE", str(tmp_path / "cache"))
    # Kmatch reads the suite's single shared-resource client, not a vendored copy.
    monkeypatch.setattr(
        "kbiosuite.resource_client.bindings",
        lambda _: {
            "foldseek": str(tmp_path / "central/foldseek"),
            "pdb100": str(tmp_path / "central/pdb"),
        },
    )
    config = Config.load()
    assert config.foldseek == tmp_path / "central/foldseek"
    assert config.pdb_prefix == tmp_path / "central/pdb"
    assert config.cache == tmp_path / "cache"


def test_unresolvable_shared_pins_fail_with_explicit_error(tmp_path, monkeypatch):
    home = tmp_path / "Kmatch"
    home.mkdir()
    (home / "config.json").write_text(json.dumps({"shared_resources": "kmatch"}))
    monkeypatch.setenv("KMATCH_HOME", str(home))
    monkeypatch.setenv("KMATCH_CACHE", str(tmp_path / "cache"))
    # An empty Resources root has no consumer record: the pin cannot silently fall back.
    monkeypatch.setenv("KBIOSUITE_RESOURCES", str(tmp_path / "Resources"))
    with pytest.raises(RuntimeError, match="shared resources .*kmatch.* need attention"):
        Config.load()
    monkeypatch.setattr("kbiosuite.resource_client.bindings", lambda _: {"foldseek": "/x"})
    with pytest.raises(RuntimeError, match="pdb100"):
        Config.load()


def test_shared_setup_only_checks_and_never_downloads_or_retires(tmp_path, monkeypatch):
    config = Config(tmp_path, tmp_path / "cache", shared_resources="kmatch")
    calls = []
    monkeypatch.setattr(
        "kmatch.setup.inspect_database", lambda *args: calls.append("inspect") or {}
    )
    monkeypatch.setattr("kmatch.setup.download", lambda *args, **kw: pytest.fail("download"))
    monkeypatch.setattr("kmatch.setup.activate_pdb_generation", lambda *args: pytest.fail("retire"))
    install_pdb(config, Cancellation(), lambda event: None)
    assert calls == ["inspect"]
    # Missing shared tool must fail; it cannot download a different executable.
    with pytest.raises(RuntimeError, match="missing"):
        install_foldseek(config, Cancellation(), lambda event: None)
    assert not (tmp_path / "databases").exists()


@pytest.mark.parametrize(
    "install, key", [(install_foldseek, "foldseek_override"), (install_pdb, "pdb_prefix_override")]
)
def test_managed_setup_refuses_an_external_override(tmp_path, monkeypatch, install, key):
    external = tmp_path / "external/tool"
    external.parent.mkdir()
    external.write_text("user-managed")
    config = Config(tmp_path / "home", tmp_path / "cache", **{key: str(external)})

    def offline(*args, **kwargs):
        pytest.fail("an override must not lead setup to download or query upstream")

    monkeypatch.setattr("kmatch.setup.subprocess", SimpleNamespace(check_output=offline))
    monkeypatch.setattr("kmatch.setup.download", offline)
    with pytest.raises(RuntimeError, match=f"external override.*remove {key}"):
        install(config, Cancellation(), lambda event: None)
    assert external.read_text() == "user-managed" and not (tmp_path / "home").exists()


def load(tmp_path, monkeypatch, settings=None):
    """Config.load() for a temporary home; without settings the home is never created."""
    home = tmp_path / "Kmatch"
    if settings is not None:
        home.mkdir(exist_ok=True)
        (home / "config.json").write_text(json.dumps(settings))
    monkeypatch.setenv("KMATCH_HOME", str(home))
    monkeypatch.setenv("KMATCH_CACHE", str(tmp_path / "cache"))
    return Config.load()


def managed_paths(config):
    return config.home / "tools/foldseek/bin/foldseek", config.home / "databases/pdb100/current/pdb"


@pytest.mark.parametrize("settings", [None, {}, {"threads": 4, "pymol_override": "/x/pymol"}])
def test_a_home_that_chooses_nothing_binds_the_default_kmatch_pins(
    tmp_path, monkeypatch, shared_library, settings
):
    config = load(tmp_path, monkeypatch, settings)
    assert config.shared_resources == "kmatch"
    assert (config.foldseek, config.pdb_prefix) == (
        shared_library.foldseek,
        shared_library.pdb_prefix,
    )
    assert config.foldseek.is_relative_to(shared_library.root)
    if settings is None:
        assert not config.home.exists()  # Loading reads only.


def test_default_shared_binding_reports_ready_and_setup_only_verifies(
    tmp_path, monkeypatch, shared_library
):
    config = load(tmp_path, monkeypatch)
    before = shared_library.snapshot()
    assert gui.dependency_readiness(config) == {"boltz2": False, "foldseek": True, "pdb100": True}
    job = gui.Job(config)
    client = gui.create_app(job).test_client()
    token = re.search(r'name="kmatch-token" content="([^"]+)"', client.get("/").text)[1]
    headers = {"X-Kmatch-Token": token}
    items = {
        item["name"]: item
        for item in client.get("/api/dependencies", headers=headers).get_json()["items"]
    }
    assert items["Foldseek"]["ready"] and items["Foldseek"]["version"] == "10.shared"
    assert items["PDB100"]["ready"]

    def offline(*args, **kwargs):
        pytest.fail("shared setup must not download, query upstream or activate")

    # Any managed-install step (curl, download, activation) fails the test instead.
    monkeypatch.setattr("kmatch.setup.subprocess", SimpleNamespace(check_output=offline))
    monkeypatch.setattr("kmatch.setup.download", offline)
    monkeypatch.setattr("kmatch.setup.activate_pdb_generation", offline)
    for component in ("foldseek", "pdb100"):
        response = client.post("/api/setup", headers=headers, json={"component": component})
        assert response.status_code == 200
        job.thread.join(10)
        state = job.snapshot()
        assert state["phase"] == "complete" and not state["error"], state
    assert shared_library.snapshot() == before
    assert not any(config.home.iterdir())  # No managed tools/ or databases/ copy.


def rewrite_consumer(library):
    path = library.root / "consumers/kmatch.json"
    record = json.loads(path.read_text())
    del record["bindings"]["pdb100"]
    path.write_text(json.dumps(record))


ABSENT = {
    "no library": lambda library: shutil.rmtree(library.root),
    "no kmatch consumer": lambda library: (library.root / "consumers/kmatch.json").unlink(),
}
BREAKAGES = {
    "no pdb100 role": rewrite_consumer,
    "consumer not JSON": lambda library: (library.root / "consumers/kmatch.json").write_text("{"),
    "catalog not an object": lambda library: (library.root / "catalog.json").write_text("[]"),
    "catalog not JSON": lambda library: (library.root / "catalog.json").write_text("{"),
    "changed pinned file": lambda library: library.pdb_prefix.write_text("changed"),
}
EXPLICIT_ERROR = r"^Kmatch shared resources \('kmatch'\) need attention"


@pytest.mark.parametrize("absent", ABSENT)
def test_without_a_library_or_kmatch_record_the_default_falls_back_silently(
    tmp_path, monkeypatch, shared_library, absent
):
    ABSENT[absent](shared_library)
    config = load(tmp_path, monkeypatch)
    assert config.shared_resources is None
    assert (config.foldseek, config.pdb_prefix) == managed_paths(config)
    # The same absence behind an explicit setting stays a startup error, never a fallback.
    with pytest.raises(RuntimeError, match=EXPLICIT_ERROR):
        load(tmp_path, monkeypatch, {"shared_resources": "kmatch"})


@pytest.mark.parametrize("breakage", BREAKAGES)
def test_a_broken_registered_kmatch_record_is_a_startup_error_by_default(
    tmp_path, monkeypatch, shared_library, breakage, capsys
):
    BREAKAGES[breakage](shared_library)
    before = shared_library.snapshot()
    # Never a silent fallback that Setup could fill with a multi-GB private download.
    with pytest.raises(
        RuntimeError, match=r"^Kmatch default shared resources \('kmatch'\)"
    ) as error:
        load(tmp_path, monkeypatch)
    message = str(error.value)
    assert "need attention" in message and "sets no shared_resources" in message
    assert '"shared_resources": null' in message
    with pytest.raises(SystemExit) as exit:
        cli.main(["status"])
    assert exit.value.code == 1
    assert "could not load its configuration: Kmatch default shared" in capsys.readouterr().err
    with pytest.raises(RuntimeError, match=EXPLICIT_ERROR):
        load(tmp_path, monkeypatch, {"shared_resources": "kmatch"})
    # The explicit opt-out the message suggests uses Kmatch-managed paths.
    config = load(tmp_path, monkeypatch, {"shared_resources": None})
    assert (config.foldseek, config.pdb_prefix) == managed_paths(config)
    assert shared_library.snapshot() == before


@pytest.mark.parametrize("kind", ["file", "dangling symlink", "unreadable consumers"])
def test_an_uninspectable_library_root_is_a_startup_error_by_default(
    tmp_path, monkeypatch, shared_library, kind
):
    # Only an absent root or record is a silent default; anything else might be a real library
    # (e.g. on an unmounted volume) that Setup must not replace with a private download.
    consumers = shared_library.root / "consumers"
    if kind == "file":
        root = tmp_path / "Resources-file"
        root.write_text("not a library")
        monkeypatch.setenv("KBIOSUITE_RESOURCES", str(root))
    elif kind == "dangling symlink":
        root = tmp_path / "Resources-link"
        root.symlink_to(tmp_path / "unmounted-volume/Resources")
        monkeypatch.setenv("KBIOSUITE_RESOURCES", str(root))
    else:
        consumers.chmod(0)
    try:
        with pytest.raises(RuntimeError, match=r"^Kmatch default shared resources \('kmatch'\)"):
            load(tmp_path, monkeypatch)
    finally:
        consumers.chmod(0o755)


@pytest.mark.parametrize("value", ["", "  ", 0, 1, True, [], {}, ["kmatch"]])
def test_shared_resources_must_be_a_consumer_name_or_an_opt_out(
    tmp_path, monkeypatch, shared_library, value
):
    with pytest.raises(ValueError, match="shared_resources must name a Resources consumer"):
        load(tmp_path, monkeypatch, {"shared_resources": value})


def test_null_path_overrides_still_use_the_default_pins(tmp_path, monkeypatch, shared_library):
    config = load(tmp_path, monkeypatch, {"foldseek_override": None, "pdb_prefix_override": None})
    assert config.shared_resources == "kmatch"
    assert (config.foldseek, config.pdb_prefix) == (
        shared_library.foldseek,
        shared_library.pdb_prefix,
    )


@pytest.mark.parametrize("value", [None, False])
def test_explicit_opt_out_keeps_home_paths_despite_resolvable_pins(
    tmp_path, monkeypatch, shared_library, value
):
    config = load(tmp_path, monkeypatch, {"shared_resources": value})
    assert not config.shared_resources
    assert (config.foldseek, config.pdb_prefix) == managed_paths(config)


@pytest.mark.parametrize("key", ["foldseek_override", "pdb_prefix_override"])
def test_an_explicit_path_override_disables_the_default_binding(
    tmp_path, monkeypatch, shared_library, key
):
    external = tmp_path / "external"
    config = load(tmp_path, monkeypatch, {key: str(external)})
    assert config.shared_resources is None
    foldseek, pdb_prefix = managed_paths(config)
    expected = (external, pdb_prefix) if key == "foldseek_override" else (foldseek, external)
    assert (config.foldseek, config.pdb_prefix) == expected
    # An invalid override is reported, not hidden behind the default binding.
    with pytest.raises(ValueError, match=f"{key} must be a nonempty path"):
        load(tmp_path, monkeypatch, {key: ""})
