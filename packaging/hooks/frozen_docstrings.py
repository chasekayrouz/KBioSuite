"""Narrow source-only docstring transforms for reviewed, locked upstream modules.

The pre-find hooks select generated source using PreFindModulePathAPI.search_dirs. Full
source hashes bind the reviewed wheel members; changes to versions, wheels or source require
review. Executable statements, assertions and all non-docstring constants remain untouched.
"""

import ast
import hashlib
import importlib.metadata
import tomllib
from pathlib import Path

from PyInstaller.config import CONF

PINS = {
    "numpy.lib._datasource": {
        "distribution": "numpy",
        "version": "2.5.3",
        "wheel_sha256": "adc1ada2662f8a5f960b8a10d9986897e7499ef07e06d4cfe7197f8cce923c07",
        "source_sha256": "9ead7eb599f6ed893aae85ddf9a8316502b4efd0595cea2b9ddf7af1ab097460",
        "docstrings": ("", "DataSource", "Repository"),
    },
    "urllib3.util.url": {
        "distribution": "urllib3",
        "version": "2.8.0",
        "wheel_sha256": "0cf3cae568d36aa9576b28dfb35f11328f1cb974ca7647d9475ebb86c75ac6e3",
        "source_sha256": "6a2814d91987dd9728278d413a154ee5277b15cb262d9945425377f2a693311f",
        "docstrings": ("Url.url",),
    },
}
POINTER = (
    "Documentation example omitted from frozen distribution; see upstream module documentation."
)
SCOPES = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)


def docstrings(tree):
    result = {}

    def visit(node, names):
        if isinstance(node, SCOPES):
            if not isinstance(node, ast.Module):
                names = (*names, node.name)
            if ast.get_docstring(node, clean=False) is not None:
                result[".".join(names)] = node.body[0]
        for child in ast.iter_child_nodes(node):
            visit(child, names)

    visit(tree, ())
    return result


def executable_ast(source):
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, SCOPES) and ast.get_docstring(node, clean=False) is not None:
            node.body.pop(0)
    return ast.dump(tree, include_attributes=False)


def transform(module, source):
    pin = PINS[module]  # no general-purpose allowlist or heuristic source rewriting
    if hashlib.sha256(source).hexdigest() != pin["source_sha256"]:
        raise RuntimeError(f"Unreviewed frozen documentation source: {module}")
    tree = ast.parse(source)
    docs = docstrings(tree)
    lines = source.splitlines(keepends=True)
    replacements = []
    for name in pin["docstrings"]:
        if name not in docs:
            raise RuntimeError(f"Reviewed docstring is missing: {module}")
        node = docs[name]
        start = sum(map(len, lines[: node.lineno - 1])) + node.col_offset
        end = sum(map(len, lines[: node.end_lineno - 1])) + node.end_col_offset
        replacement = repr(POINTER).encode() + b"\n" * source[start:end].count(b"\n")
        replacements.append((start, end, replacement))
    output = source
    for start, end, replacement in sorted(replacements, reverse=True):
        output = output[:start] + replacement + output[end:]
    if executable_ast(source) != executable_ast(output):
        raise RuntimeError(f"Frozen docstring transformation changed executable AST: {module}")
    compile(output, module.replace(".", "/") + ".py", "exec", optimize=0)
    return output


def pre_find_module_path(api):
    pin = PINS[api.module_name]
    root = Path(__file__).resolve().parents[2]
    lock = tomllib.loads((root / "uv.lock").read_text())
    package = next(p for p in lock["package"] if p["name"] == pin["distribution"])
    if package["version"] != pin["version"] or not any(
        wheel["hash"] == "sha256:" + pin["wheel_sha256"] for wheel in package["wheels"]
    ):
        raise RuntimeError(f"Unreviewed frozen documentation lock: {api.module_name}")
    distribution = importlib.metadata.distribution(pin["distribution"])
    if distribution.version != pin["version"]:
        raise RuntimeError(f"Unreviewed frozen documentation version: {api.module_name}")
    relative = Path(*api.module_name.split(".")).with_suffix(".py")
    original = Path(distribution.locate_file(relative))
    output = transform(api.module_name, original.read_bytes())
    destination = Path(CONF["workpath"]) / "kbiosuite-frozen-docstrings" / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(output)
    api.search_dirs = [str(destination.parent)]
