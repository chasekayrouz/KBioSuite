"""Local, explicit adoption and discovery of immutable scientific resources.

No automatic deletion, downloading, shell changes, or dependency on a daemon.
"""

import argparse
import fcntl
import hashlib
import html
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from datetime import UTC, datetime
from importlib.resources import files
from pathlib import Path

from .resource_client import bindings, contained, read_json, resolve, resource_root


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def fingerprint(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as out:
        temp = Path(out.name)
        json.dump(value, out, indent=2)
        out.write("\n")
        out.flush()
        os.fsync(out.fileno())
    try:
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


@contextmanager
def writer(root):
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".writer.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def catalog(root):
    path = root / "catalog.json"
    return read_json(path) if path.exists() else {"schema_version": 1, "resources": {}}


def verify(identity, root=None):
    root = Path(root or resource_root())
    directory = resolve(identity, root=root)
    manifest = read_json(directory / "resource.json")
    for name, expected in manifest["files"].items():
        if digest(contained(directory, name)) != expected["sha256"]:
            raise ValueError(f"Checksum mismatch: {identity}/{name}")
    actual = fingerprint(
        {k: manifest[k] for k in ("kind", "name", "version", "entrypoints", "files")}
    )
    if actual != manifest["fingerprint"]:
        raise ValueError(f"Manifest fingerprint mismatch: {identity}")
    return {"id": identity, "files_verified": len(manifest["files"]), "passed": True}


def adopt(recipe, root=None):
    """Recipe explicitly maps public resource files to source paths and expected hashes.

    All source links must resolve within their declared source directory. Source
    files are never modified. A failed staging directory is retained for diagnosis.
    """
    root = Path(root or resource_root()).resolve()
    for field in ("name", "version"):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.+-]*", recipe[field]):
            raise ValueError(f"Invalid {field}")
    if recipe["kind"] not in {"databases", "tools", "models"}:
        raise ValueError("Invalid resource kind")
    source = Path(recipe["source"]).expanduser().resolve()
    requested = recipe["files"]
    if not requested:
        raise ValueError("No source files selected")
    files = {}
    for name, expected in sorted(requested.items()):
        if name == "resource.json" or any(p.startswith(".") for p in Path(name).parts):
            raise ValueError(f"Hidden/reserved file excluded: {name}")
        path = contained(source, name)
        if not path.is_file():
            raise ValueError(f"Missing source file: {path}")
        checksum = digest(path)
        if expected and checksum != expected:
            raise ValueError(f"Source checksum mismatch: {name}")
        files[name] = {
            "bytes": path.stat().st_size,
            "sha256": checksum,
            "executable": os.access(path, os.X_OK),
        }
    for relative in recipe["entrypoints"].values():
        contained(source, relative)
        if (
            relative != "."
            and relative not in files
            and not any(n.startswith((relative + ".", relative + "/")) for n in files)
        ):
            raise ValueError(f"Entrypoint has no selected files: {relative}")
    manifest = {
        "kind": recipe["kind"],
        "name": recipe["name"],
        "version": recipe["version"],
        "entrypoints": recipe["entrypoints"],
        "files": files,
    }
    identity_hash = fingerprint(manifest)
    identity = f"{recipe['name']}@{recipe['version']}:{identity_hash[:16]}"
    relative = f"{recipe['kind']}/{recipe['name']}/{recipe['version']}-{identity_hash[:16]}"
    manifest.update(
        schema_version=1,
        id=identity,
        fingerprint=identity_hash,
        provenance=recipe.get("provenance", {}),
        adopted_at=datetime.now(UTC).isoformat(),
        bytes=sum(x["bytes"] for x in files.values()),
    )
    with writer(root):
        current = catalog(root)
        if identity in current["resources"]:
            verify(identity, root)
            return identity
        destination = contained(root, relative)
        if destination.exists():
            raise ValueError(f"Unregistered installation exists; inspect: {destination}")
        staging = root / "staging"
        staging.mkdir(exist_ok=True)
        stage = Path(tempfile.mkdtemp(prefix=recipe["name"] + "-", dir=staging))
        for name, info in files.items():
            src = contained(source, name)
            target = stage / name
            target.parent.mkdir(parents=True, exist_ok=True)
            # A clone is independent of the source. Never create cross-install hardlinks.
            cloned = False
            if sys.platform == "darwin":
                cloned = (
                    subprocess.run(
                        ["/bin/cp", "-cp", str(src), str(target)], capture_output=True, check=False
                    ).returncode
                    == 0
                )
            if not cloned:
                if shutil.disk_usage(root).free < info["bytes"] + 1024**3:
                    raise RuntimeError(
                        "Insufficient space for an independent copy; staging retained"
                    )
                shutil.copy2(src, target)
            if digest(target) != info["sha256"]:
                raise ValueError(f"Copied file checksum mismatch: {name}")
        atomic_json(stage / "resource.json", manifest)
        destination.parent.mkdir(parents=True, exist_ok=True)
        stage.rename(destination)
        current["resources"][identity] = {
            "path": relative,
            "fingerprint": identity_hash,
            "name": recipe["name"],
            "version": recipe["version"],
            "kind": recipe["kind"],
            "bytes": manifest["bytes"],
        }
        atomic_json(root / "catalog.json", current)
        refresh_view(root)
    return identity


def register(name, record, root=None):
    root = Path(root or resource_root())
    if not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", name):
        raise ValueError("Invalid consumer name")
    record = dict(record, schema_version=1)
    with writer(root):
        for pin in record["bindings"].values():
            resolve(pin["resource"], pin["entrypoint"], root)
        for mapping in record.get("path_mappings", []):
            resolve(mapping["resource"], root=root)
            contained(Path("/"), mapping["prefix"])
        path = root / "consumers" / f"{name}.json"
        if path.exists():
            old = path.read_bytes()
            archive = root / "history" / f"{name}-{hashlib.sha256(old).hexdigest()[:16]}.json"
            if not archive.exists():
                atomic_json(archive, json.loads(old))
        atomic_json(path, record)
        refresh_view(root)
    return bindings(name, root)


LOCAL_NOTES = "LOCAL-NOTES.md"


def guide_text(root):
    """The packaged guide, plus a link to owner-maintained notes when they exist.

    The notes file belongs to the library's owner: it is only linked, never written.
    """
    guide = files("kbiosuite").joinpath("resource_guide.md").read_text()
    if (root / LOCAL_NOTES).is_file():
        guide = guide.rstrip("\n") + (
            f"\n\n## Local notes\n\nThis library has owner-maintained notes in "
            f"[{LOCAL_NOTES}]({LOCAL_NOTES}). KBioSuite links that file but never "
            "writes, replaces or removes it.\n"
        )
    return guide


def refresh_view(root):
    """Static, offline view; regenerated only on explicit management operations."""
    current = catalog(root)
    (root / "README.md").write_text(guide_text(root))
    users = {}
    for path in sorted((root / "consumers").glob("*.json")):
        record = read_json(path)
        for pin in record["bindings"].values():
            users.setdefault(pin["resource"], set()).add(path.stem)
    rows = []
    for identity, item in current["resources"].items():
        try:
            resolve(identity, root=root)
            state = "Files present"
        except (ValueError, OSError, KeyError):
            state = "Needs attention"
        cells = [
            item["name"],
            item["version"],
            f"{item['bytes'] / 1024**3:.2f} GiB",
            ", ".join(sorted(users.get(identity, []))) or "Unassigned",
            state,
        ]
        rows.append("<tr>" + "".join(f"<td>{html.escape(x)}</td>" for x in cells) + "</tr>")
    page = """<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>KBioSuite resources</title>
<style>body{font:16px/1.5 system-ui;margin:30px auto;padding:0 20px;
max-width:1000px;color:#18352c;background:#f6f8f6}
h1{font-size:28px}table{border-collapse:collapse;width:100%;background:white}
td,th{text-align:left;padding:12px;border-bottom:1px solid #d9e2da}
code{overflow-wrap:anywhere}.scroll{overflow-x:auto}small{color:#4d625a}</style>
<h1>Shared bioinformatics resources</h1><p>Installed versions and their registered applications.</p>
<div class="scroll"><table><thead><tr><th>Resource</th><th>Version</th>
<th>File size</th><th>Used by</th><th>Status</th></tr></thead><tbody>"""
    page += "".join(rows) + "</tbody></table></div>"
    page += f"<p><code>{html.escape(str(root))}</code></p>"
    page += (
        "<p>Versions remain available until explicitly retired. Updates and downloads never run "
        "when an app opens. This view lists installed releases; it does not check the internet "
        "for newer versions.</p><small>File presence is a quick check, "
        "not a checksum verification. "
        "Use the resource manager’s verify command for a full check. File sizes do not account "
        "for APFS shared blocks.</small></html>"
    )
    temporary = root / ".index-next.html"
    temporary.write_text(page)
    temporary.replace(root / "index.html")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=resource_root())
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list")
    sub.add_parser("view")
    p = sub.add_parser("resolve")
    p.add_argument("id")
    p.add_argument("--entrypoint")
    p = sub.add_parser("verify")
    p.add_argument("id", nargs="?")
    p = sub.add_parser("adopt")
    p.add_argument("recipe", type=Path)
    p = sub.add_parser("register")
    p.add_argument("name")
    p.add_argument("record", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "list":
            result = catalog(args.root)
        elif args.command == "resolve":
            result = str(resolve(args.id, args.entrypoint, args.root))
        elif args.command == "verify":
            ids = [args.id] if args.id else catalog(args.root)["resources"]
            result = [verify(identity, args.root) for identity in ids]
        elif args.command == "adopt":
            result = adopt(json.loads(args.recipe.read_text()), args.root)
        elif args.command == "register":
            result = register(args.name, json.loads(args.record.read_text()), args.root)
        else:
            with writer(args.root):
                refresh_view(args.root)
            result = str(args.root / "index.html")
        print(json.dumps(result, indent=2))
        return 0
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        print(f"Resources: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
