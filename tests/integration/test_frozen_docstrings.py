"""Pinned upstream documentation changes preserve every executable AST statement."""

import ast
import dis
import importlib
import importlib.util
import marshal
import re
import types
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[2]
HOOKS = ROOT / "packaging/hooks"
spec = importlib.util.spec_from_file_location("frozen_docstrings", HOOKS / "frozen_docstrings.py")
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
PERSONAL = re.compile(rb"/(?:Users|home)/[^/\s]+/")
AUTH = re.compile(rb"https?://[^/\s:@]+:[^/\s@]+@[A-Za-z0-9]")


def assertion_count(code):
    return sum(i.opname == "LOAD_ASSERTION_ERROR" for i in dis.get_instructions(code)) + sum(
        assertion_count(value) for value in code.co_consts if isinstance(value, types.CodeType)
    )


@pytest.mark.parametrize("module", helper.PINS)
def test_exact_source_preserves_executable_ast_assertions_and_other_docstrings(module):
    original = Path(importlib.import_module(module).__file__).read_bytes()
    transformed = helper.transform(module, original)
    assert helper.executable_ast(original) == helper.executable_ast(transformed)
    before, after = (
        helper.docstrings(ast.parse(original)),
        helper.docstrings(ast.parse(transformed)),
    )
    for name, node in before.items():
        expected = helper.POINTER if name in helper.PINS[module]["docstrings"] else node.value.value
        assert after[name].value.value == expected, name
    filename = module.replace(".", "/") + ".py"
    old, new = (compile(value, filename, "exec", optimize=0) for value in (original, transformed))
    assert assertion_count(old) == assertion_count(new)
    payload = marshal.dumps(new)  # Includes every nested code object and its docstrings/constants.
    assert PERSONAL.search(payload) is None
    assert AUTH.search(payload) is None
    assert Path(importlib.import_module(module).__file__).read_bytes() == original
    with pytest.raises(RuntimeError, match="Unreviewed frozen documentation source"):
        helper.transform(module, original + b"\n")


@pytest.mark.parametrize("module", helper.PINS)
def test_supported_pre_find_wrapper_selects_only_generated_pinned_source(
    module, tmp_path, monkeypatch
):
    monkeypatch.setitem(helper.CONF, "workpath", str(tmp_path))
    path = HOOKS / "pre_find_module_path" / ("hook-" + module + ".py")
    spec = importlib.util.spec_from_file_location("reviewed_hook", path)
    wrapper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(wrapper)
    api = SimpleNamespace(module_name=module, search_dirs=["original"])
    wrapper.pre_find_module_path(api)
    selected = Path(api.search_dirs[0]) / (module.rsplit(".", 1)[-1] + ".py")
    assert selected.is_file() and selected.is_relative_to(tmp_path)
    assert not PERSONAL.search(marshal.dumps(compile(selected.read_bytes(), "reviewed.py", "exec")))
    assert not AUTH.search(marshal.dumps(compile(selected.read_bytes(), "reviewed.py", "exec")))


def transformed_namespace(module):
    original = importlib.import_module(module)
    source = helper.transform(module, Path(original.__file__).read_bytes())
    namespace = {"__name__": module, "__package__": module.rpartition(".")[0]}
    exec(compile(source, module.replace(".", "/") + ".py", "exec", optimize=0), namespace)
    return original, namespace


def test_url_property_and_parsing_behavior_matches_original():
    original, transformed = transformed_namespace("urllib3.util.url")
    for value in ("https://example.org:8443/path?query#fragment", "http://[::1]/", "/relative"):
        expected = original.parse_url(value)
        actual = transformed["parse_url"](value)
        assert tuple(actual) == tuple(expected)
        assert actual.url == expected.url
        assert actual.request_uri == expected.request_uri


def test_numpy_datasource_local_file_behavior_matches_original(tmp_path):
    original, transformed = transformed_namespace("numpy.lib._datasource")
    path = tmp_path / "synthetic.txt"
    path.write_text("synthetic data\n")
    for cls in (original.DataSource, transformed["DataSource"]):
        source = cls(destpath=str(tmp_path))
        assert source.exists(str(path))
        with source.open(str(path), mode="rt") as handle:
            assert handle.read() == "synthetic data\n"


def test_unreviewed_wheel_or_distribution_upgrade_fails_closed(monkeypatch):
    module = "urllib3.util.url"
    api = SimpleNamespace(module_name=module)
    original = helper.PINS[module]
    monkeypatch.setitem(helper.PINS, module, {**original, "wheel_sha256": "0" * 64})
    with pytest.raises(RuntimeError, match="Unreviewed frozen documentation lock"):
        helper.pre_find_module_path(api)
    monkeypatch.setitem(helper.PINS, module, original)
    monkeypatch.setattr(
        helper.importlib.metadata, "distribution", lambda name: SimpleNamespace(version="999.0")
    )
    with pytest.raises(RuntimeError, match="Unreviewed frozen documentation version"):
        helper.pre_find_module_path(api)
