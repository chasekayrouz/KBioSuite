import json
from pathlib import Path

import pytest

from kbiosuite.resource_client import bindings, resolve
from kbiosuite.resources import adopt, digest, refresh_view, register, verify


def recipe(source, version="1"):
    return {
        "name": "example",
        "version": version,
        "kind": "databases",
        "source": str(source),
        "files": {"db": digest(source / "db")},
        "entrypoints": {"database": "db"},
    }


def test_adoption_is_independent_idempotent_and_pins_survive_new_versions(tmp_path):
    source, root = tmp_path / "source", tmp_path / "store with spaces"
    source.mkdir()
    (source / "db").write_text("original")
    old = adopt(recipe(source), root)
    assert adopt(recipe(source), root) == old
    register(
        "test-app", {"bindings": {"database": {"resource": old, "entrypoint": "database"}}}, root
    )
    path = Path(bindings("test-app", root)["database"])
    (source / "db").write_text("new data")
    new = adopt(recipe(source, "2"), root)
    assert old != new and path.read_text() == "original"
    assert bindings("test-app", root)["database"] == str(path)
    assert verify(old, root)["passed"]
    assert verify(new, root)["passed"]
    assert (root / "index.html").is_file()


def test_changed_files_fail_verification_and_missing_pin_never_falls_back(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "db").write_text("original")
    root = tmp_path / "store"
    identity = adopt(recipe(source), root)
    path = resolve(identity, "database", root)
    path.write_text("modified")
    with pytest.raises(ValueError, match="Checksum"):
        verify(identity, root)
    path.unlink()
    with pytest.raises(ValueError, match="Missing"):
        resolve(identity, "database", root)
    with pytest.raises(KeyError):
        resolve("absent", root=root)


def test_rejects_external_links_paths_and_bad_source_hash_without_publication(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "db").write_text("original")
    root = tmp_path / "store"
    entry = recipe(source)
    entry["files"]["db"] = "0" * 64
    with pytest.raises(ValueError, match="Source checksum"):
        adopt(entry, root)
    assert not (root / "catalog.json").exists()
    outside = tmp_path / "outside"
    outside.write_text("external")
    (source / "db").unlink()
    (source / "db").symlink_to(outside)
    with pytest.raises(ValueError, match="escapes"):
        adopt(recipe(source), root)
    entry = recipe(source)
    entry["files"] = {"../outside": None}
    with pytest.raises(ValueError):
        adopt(entry, root)


def test_unregistered_staging_is_not_discoverable_and_bad_registration_is_atomic(tmp_path):
    root = tmp_path / "store"
    (root / "staging/unfinished").mkdir(parents=True)
    (root / "staging/unfinished/db").write_text("unfinished")
    source = tmp_path / "source"
    source.mkdir()
    (source / "db").write_text("complete")
    identity = adopt(recipe(source), root)
    good = {"bindings": {"database": {"resource": identity, "entrypoint": "database"}}}
    register("app", good, root)
    before = (root / "consumers/app.json").read_bytes()
    with pytest.raises(KeyError):
        register(
            "app",
            {"bindings": {"database": {"resource": "missing", "entrypoint": "database"}}},
            root,
        )
    assert (root / "consumers/app.json").read_bytes() == before
    assert len(json.loads((root / "catalog.json").read_text())["resources"]) == 1


def test_view_links_owner_notes_but_never_writes_them(tmp_path):
    source, root = tmp_path / "source", tmp_path / "store"
    source.mkdir()
    (source / "db").write_text("original")
    adopt(recipe(source), root)
    guide = (root / "README.md").read_text()
    assert "](LOCAL-NOTES.md)" not in guide
    assert not (root / "LOCAL-NOTES.md").exists()
    notes = root / "LOCAL-NOTES.md"
    notes.write_text("Owner notes stay as written.\n")
    before = notes.stat().st_mtime_ns
    for _ in range(2):
        refresh_view(root)
        readme = (root / "README.md").read_text()
        assert readme.startswith(guide.rstrip("\n"))
        assert readme.count("[LOCAL-NOTES.md](LOCAL-NOTES.md)") == 1
    assert notes.read_text() == "Owner notes stay as written.\n"
    assert notes.stat().st_mtime_ns == before
    notes.unlink()
    refresh_view(root)
    assert (root / "README.md").read_text() == guide
