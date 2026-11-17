"""Read-only shared-resource protocol v1; portable stdlib client.

Canonical source: KBioSuite/src/kbiosuite/resource_client.py. Copies in independent
consumers implement this file protocol without importing an application checkout.
"""

import json
import os
from pathlib import Path


def resource_root():
    return (
        Path(
            os.environ.get(
                "KBIOSUITE_RESOURCES",
                Path.home() / "Library/Application Support/KBioSuite/Resources",
            )
        )
        .expanduser()
        .resolve()
    )


def contained(root, relative):
    relative = Path(relative)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Unsafe resource path: {relative}")
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError(f"Resource escapes its installation: {relative}")
    return path


def read_json(path):
    data = json.loads(path.read_text())
    if data.get("schema_version") != 1:
        raise ValueError(f"Unsupported resource schema: {path}")
    return data


def consumer(name, root=None):
    root = Path(root or resource_root())
    if not name or any(c not in "abcdefghijklmnopqrstuvwxyz0123456789-_" for c in name):
        raise ValueError("Invalid consumer name")
    return read_json(root / "consumers" / f"{name}.json")


def resolve(resource_id, entrypoint=None, root=None):
    root = Path(root or resource_root())
    catalog = read_json(root / "catalog.json")
    item = catalog["resources"][resource_id]
    directory = contained(root, item["path"])
    manifest = read_json(directory / "resource.json")
    if manifest["id"] != resource_id or manifest["fingerprint"] != item["fingerprint"]:
        raise ValueError(f"Resource identity mismatch: {resource_id}")
    # Cheap readiness at launch; explicit verification performs full hashes.
    for name, info in manifest["files"].items():
        path = contained(directory, name)
        if not path.is_file() or path.stat().st_size != info["bytes"]:
            raise ValueError(f"Missing or changed resource file: {resource_id}/{name}")
    if entrypoint is None:
        return directory
    return contained(directory, manifest["entrypoints"][entrypoint])


def bindings(name, root=None):
    root = Path(root or resource_root())
    record = consumer(name, root)
    # Resolve each immutable generation once per configuration load/job.
    directories = {}
    result = {}
    for role, pin in record["bindings"].items():
        identity = pin["resource"]
        if identity not in directories:
            directory = resolve(identity, root=root)
            directories[identity] = (directory, read_json(directory / "resource.json"))
        directory, manifest = directories[identity]
        result[role] = str(contained(directory, manifest["entrypoints"][pin["entrypoint"]]))
    return result
