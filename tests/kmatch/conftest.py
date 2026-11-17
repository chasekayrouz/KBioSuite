import struct
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.fixture(autouse=True)
def no_shared_resources(tmp_path, monkeypatch):
    """Never read the developer's shared Resources library.

    Config.load() binds the default ``kmatch`` pins when config.json chooses nothing, so a
    temporary home would otherwise resolve the real library. This root is never created;
    tests that need pins use ``shared_library`` or patch the client.
    """
    monkeypatch.setenv("KBIOSUITE_RESOURCES", str(tmp_path / "no-shared-resources"))


def snapshot(root):
    """Content, mode and mtime of every entry below root (detects any write or removal)."""
    return {
        str(path.relative_to(root)): (
            path.read_bytes() if path.is_file() else None,
            path.lstat().st_mode,
            path.lstat().st_mtime_ns,
        )
        for path in sorted(root.rglob("*"))
    }


@pytest.fixture
def shared_library(tmp_path, monkeypatch):
    """A real temporary Resources library whose ``kmatch`` pins are ready for Kmatch."""
    from kbiosuite.resources import adopt, register
    from kmatch.database import REQUIRED

    sources = tmp_path / "shared-sources"
    tool = sources / "foldseek/bin/foldseek"
    tool.parent.mkdir(parents=True)
    tool.write_text("#!/bin/sh\necho 10.shared\n")
    tool.chmod(0o755)
    dbtypes = {".dbtype": 0, "_ss.dbtype": 0, "_ca.dbtype": 101, "_h.dbtype": 12}
    for suffix in REQUIRED:
        path = sources / f"pdb100/pdb{suffix}"
        path.parent.mkdir(exist_ok=True)
        if suffix in dbtypes:
            path.write_bytes(struct.pack("<I", dbtypes[suffix]))
        else:
            path.write_text("1\t0\t5\n" if suffix.endswith(".index") else "data\n")
    root = tmp_path / "Resources"
    recipes = {
        "foldseek": ("tools", {"bin/foldseek": None}, "foldseek", "bin/foldseek"),
        "pdb100": ("databases", {f"pdb{s}": None for s in REQUIRED}, "database", "pdb"),
    }
    pins = {}
    for role, (kind, files, entrypoint, relative) in recipes.items():
        recipe = {
            "kind": kind,
            "name": role,
            "version": "test",
            "source": str(sources / role),
            "files": files,
            "entrypoints": {entrypoint: relative},
        }
        pins[role] = {"resource": adopt(recipe, root), "entrypoint": entrypoint}
    paths = register("kmatch", {"bindings": pins}, root)
    monkeypatch.setenv("KBIOSUITE_RESOURCES", str(root))
    return SimpleNamespace(
        root=root,
        foldseek=Path(paths["foldseek"]),
        pdb_prefix=Path(paths["pdb100"]),
        snapshot=lambda: snapshot(root),
    )
