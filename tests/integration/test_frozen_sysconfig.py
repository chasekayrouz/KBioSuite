"""The selected interpreter's facts survive freezing without its installation location."""

import importlib
import importlib.util
import json
import marshal
import subprocess
import sys
import sysconfig
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[2]
HOOKS = ROOT / "packaging/hooks"
HOOK = HOOKS / "pre_find_module_path/hook-_sysconfigdata__darwin_darwin.py"
spec = importlib.util.spec_from_file_location("kbio_sysconfig_hook", HOOK)
hook = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hook)


def test_installation_rebasing_preserves_abi_and_compile_facts(monkeypatch):
    prefix = str(Path("/") / "Users" / "synthetic" / "managed-python")
    values = {
        "prefix": prefix,
        "exec_prefix": prefix,
        "LIBDEST": prefix + "/lib/python3.14",
        "INCLUDEPY": prefix + "/include/python3.14",
        "EXENAME": prefix + "/bin/python3.14",
        "BINDIR": prefix + "/bin",
        "DESTDIRS": prefix + " " + prefix + "/lib",
        "SOABI": "cpython-314-darwin",
        "EXT_SUFFIX": ".cpython-314-darwin.so",
        "SIZEOF_VOID_P": 8,
        "Py_GIL_DISABLED": 0,
        "CFLAGS": "-O3 -fPIC",
        "CC": "clang",
        "CONFIG_ARGS": "'--enable-shared' '--prefix=/install'",
    }
    source = hook.module_source(values, prefix)
    code = compile(source, "_sysconfigdata__darwin_darwin.py", "exec")
    assert prefix.encode() not in marshal.dumps(code)
    assert str(Path.home()).encode() not in marshal.dumps(code)
    monkeypatch.setattr(sys, "_MEIPASS", "/frozen/runtime", raising=False)
    monkeypatch.setattr(sys, "executable", "/frozen/bin/application")
    namespace = {}
    exec(code, namespace)
    actual = namespace["build_time_vars"]
    for key, value in values.items():
        if not isinstance(value, str) or prefix not in value:
            assert actual[key] == value, key
    assert actual["prefix"] == actual["exec_prefix"] == actual["LIBDEST"] == "/frozen/runtime"
    assert actual["INCLUDEPY"] == "/frozen/runtime/include/python3.14"
    assert actual["DESTDIRS"] == "/frozen/runtime /frozen/runtime/lib"
    assert actual["EXENAME"] == sys.executable
    assert actual["BINDIR"] == "/frozen/bin"


def test_unrecognized_personal_metadata_is_not_silently_removed():
    unknown = str(Path("/") / "home" / "synthetic" / "private-build")
    with pytest.raises(ValueError, match="Unrecognized personal path"):
        hook.module_source({"CFLAGS": "-I" + unknown}, "/managed/python")
    source = hook.module_source({"CC": "/managed/python-other/clang"}, "/managed/python")
    assert "/managed/python-other/clang" in source  # prefix lookalikes are not installation paths


def test_pre_find_hook_redirects_analysis_without_modifying_interpreter(tmp_path, monkeypatch):
    name = sysconfig._get_sysconfigdata_name()
    original = Path(importlib.import_module(name).__file__)
    before = original.read_bytes()
    monkeypatch.setitem(hook.CONF, "workpath", str(tmp_path))
    api = SimpleNamespace(module_name=name, search_dirs=["original"])
    hook.pre_find_module_path(api)
    generated = Path(api.search_dirs[0]) / (name + ".py")
    assert generated.is_file()
    assert str(Path(sys.base_prefix)).encode() not in generated.read_bytes()
    assert original.read_bytes() == before


@pytest.mark.skipif(sys.platform != "darwin", reason="macOS-specific sysconfig module")
def test_actual_pyinstaller_runtime_and_pyz_are_relocatable():
    from PyInstaller.archive.readers import CArchiveReader

    name = sysconfig._get_sysconfigdata_name()
    original = importlib.import_module(name).build_time_vars
    build_root = ROOT / "build.noindex"
    build_root.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="sysconfig-probe-", dir=build_root) as temporary:
        work = Path(temporary)
        probe = work / "probe.py"
        probe.write_text("""import importlib,json,os,sys,sysconfig
values=importlib.import_module(sysconfig._get_sysconfigdata_name()).build_time_vars
print(json.dumps({
 'facts': {key: values[key] for key in ('SOABI','EXT_SUFFIX','SIZEOF_VOID_P','CFLAGS','CONFIG_ARGS')},
 'root_matches': all(values[key] == sys._MEIPASS for key in ('prefix','exec_prefix','LIBDIR','LIBDEST','DESTSHARED')),
 'executable_matches': values['EXENAME'] == sys.executable,
 'bindir_matches': values['BINDIR'] == os.path.dirname(sys.executable),
 'suffix_matches': sysconfig.get_config_var('EXT_SUFFIX') == values['EXT_SUFFIX'],
 'scheme_rooted': sysconfig.get_path('stdlib').startswith(sys._MEIPASS + os.sep),
 'runtime_libdir_exists': os.path.isdir(sysconfig.get_config_var('LIBDEST')),
}))
""")
        command = [
            sys.executable,
            "-m",
            "PyInstaller",
            "--noconfirm",
            "--clean",
            "--onedir",
            "--additional-hooks-dir",
            str(HOOKS),
            "--name",
            "sysconfig-probe",
            "--distpath",
            str(work / "dist"),
            "--workpath",
            str(work / "work"),
            "--specpath",
            str(work),
            str(probe),
        ]
        # Build log stays in ignored temporary staging; no path-bearing log is printed.
        with (work / "build.log").open("w") as log:
            result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        assert result.returncode == 0, "sysconfig probe freeze failed (local build log)"
        binary = work / "dist/sysconfig-probe/sysconfig-probe"
        result = subprocess.run([str(binary)], capture_output=True, text=True)
        assert result.returncode == 0, "frozen sysconfig probe failed"
        report = json.loads(result.stdout)
        expected = {key: original[key] for key in report["facts"]}
        assert report.pop("facts") == expected, "frozen compiler/ABI facts changed"
        assert all(report.values()), "frozen runtime installation/sysconfig paths differ"
        archive = CArchiveReader(str(binary))
        pyz_name = next(key for key, entry in archive.toc.items() if entry[-1] == "z")
        code = archive.open_embedded_archive(pyz_name).extract(name)
        payload = marshal.dumps(code)
        assert str(Path(sys.base_prefix)).encode() not in payload, "build interpreter path in PYZ"
        assert str(Path.home()).encode() not in payload, "personal path in PYZ"
