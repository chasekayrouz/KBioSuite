"""Build retention and atomic install safety, isolated from actual installed apps."""

import hashlib
import importlib.util
import json
import plistlib
import runpy
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/build_apps.py"
spec = importlib.util.spec_from_file_location("build_apps", SCRIPT)
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


def test_pyhmmer_hook_excludes_upstream_test_corpus():
    pytest.importorskip("pyhmmer")
    hook = runpy.run_path(str(SCRIPT.parents[1] / "packaging/hooks/hook-pyhmmer.py"))
    assert {"pyhmmer.easel", "pyhmmer.plan7", "pyhmmer.hmmer"} <= set(hook["hiddenimports"])
    assert not any(name.startswith("pyhmmer.tests") for name in hook["hiddenimports"])
    assert hook["datas"]
    assert not any("pyhmmer/tests" in str(path) for entry in hook["datas"] for path in entry)


def products(path, version="new", apps=None):
    path.mkdir(parents=True)
    (path / ".kbiosuite-install.json").write_text(json.dumps({
        "apps": builder.APPS if apps is None else apps,
        "source": {"commit": "fixture", "dirty": False},
    }))
    (path / "payload").write_text(version)
    return path


def test_latest_retains_only_owned_completed_builds(tmp_path):
    root = tmp_path / "build.noindex"
    old = products(root / "old/apps")
    builder.retain_latest_build(old)
    unmarked = products(root / "unmarked/apps")
    failed = root / "failed"
    failed.mkdir()
    (failed / "diagnostic.log").write_text("failure")
    foreign = products(root / "foreign/apps")
    (foreign.parent / ".kbiosuite-build.json").write_text(json.dumps({"owner": "other"}))
    external = products(tmp_path / "external/apps")
    (root / "linked").symlink_to(external.parent, target_is_directory=True)
    new = products(root / "new/apps")
    builder.retain_latest_build(new)
    assert (root / "latest").resolve() == new
    assert not old.parent.exists()
    assert all(p.exists() for p in (unmarked, failed, foreign, external, new))
    assert (root / "linked").is_symlink()


def test_marked_build_from_smaller_suite_is_pruned_after_apps_grow(tmp_path):
    root = tmp_path / "build.noindex"
    earlier = {k: v for k, v in builder.APPS.items() if k != "Kmatch"}
    old = products(root / "old/apps", apps=earlier)
    builder.retain_latest_build(old)
    # Marked, but not a completed suite build: never pruned.
    foreign = products(root / "foreign/apps", apps={**earlier, "Other": "org.example.other"})
    builder.retain_latest_build(foreign)
    incomplete = products(root / "incomplete/apps")
    builder.retain_latest_build(incomplete)
    (incomplete / ".kbiosuite-install.json").unlink()
    empty = products(root / "empty/apps", apps={})
    builder.retain_latest_build(empty)
    new = products(root / "new/apps")
    builder.retain_latest_build(new)
    assert (root / "latest").resolve() == new
    assert not old.parent.exists()
    assert all(p.exists() for p in (foreign, incomplete, empty, new))


@pytest.fixture
def install_environment(tmp_path, monkeypatch):
    built = products(tmp_path / "build.noindex/current/apps")
    destination = products(tmp_path / "Applications/KBioSuite", "old")
    calls = []

    def run(*args):
        calls.append(args)
        if args[0] == "/usr/bin/ditto":
            shutil.copytree(args[1], args[2])

    monkeypatch.setattr(builder, "source_identity", lambda: {"commit": "fixture", "dirty": False})
    monkeypatch.setattr(builder, "run", run)
    monkeypatch.setattr(builder, "validate", lambda path: None)
    monkeypatch.setattr(builder.subprocess, "check_output", lambda *a, **kw: "")
    return built, destination, calls


def test_successful_install_cleans_nonindexed_staging_and_rollback(install_environment):
    built, destination, calls = install_environment
    builder.install(built, destination)
    assert (destination / "payload").read_text() == "new"
    staging = next(call[2] for call in calls if call[0] == "/usr/bin/ditto")
    assert staging.name.endswith(".noindex") and not staging.exists()
    assert not (destination.parent / ".KBioSuite-previous.noindex").exists()
    assert len([call for call in calls if call[0] == builder.LSREGISTER]) == len(builder.APPS)


def test_failed_promotion_restores_previous_install(install_environment, monkeypatch):
    built, destination, _ = install_environment

    def validate(path):
        if path == destination and (path / "payload").read_text() == "new":
            raise RuntimeError("signature failure")

    monkeypatch.setattr(builder, "validate", validate)
    with pytest.raises(RuntimeError, match="signature failure"):
        builder.install(built, destination)
    assert (destination / "payload").read_text() == "old"
    remaining = list(destination.parent.glob(".KBioSuite-install-*.noindex"))
    assert len(remaining) == 1 and (remaining[0] / "payload").read_text() == "new"


def test_running_install_is_never_replaced(install_environment, monkeypatch):
    built, destination, calls = install_environment
    monkeypatch.setattr(
        builder.subprocess, "check_output",
        lambda *a, **kw: str(destination / "Kannotate.app/Contents/MacOS/Kannotate"),
    )
    with pytest.raises(RuntimeError, match="Quit"):
        builder.install(built, destination)
    assert not calls
    assert (destination / "payload").read_text() == "old"


MMSEQS_BINARY = b"thin arm64 mmseqs fixture"


def mmseqs_setup(directory, binary=MMSEQS_BINARY):
    """What scripts/kstack/setup-mmseqs.sh leaves: thin binary, license and setup record."""
    (directory / "bin").mkdir(parents=True)
    (directory / "bin/mmseqs").write_bytes(binary)
    (directory / "LICENSE.md").write_text("MMseqs2 license fixture")
    (directory / builder.MMSEQS_PROVENANCE).write_text(json.dumps({
        "name": "MMseqs2", "version": "18-8cc5c", "architecture": "arm64",
        "archive_sha256": "a" * 64, "binary": "bin/mmseqs",
        "binary_sha256": hashlib.sha256(binary).hexdigest(),
    }))
    return directory


def bundle(directory):
    """Minimal unsigned suite products: frozen Kmatch layout, Kstack's MMseqs2 and notices."""
    for name, identity in builder.APPS.items():
        contents = directory / f"{name}.app/Contents"
        contents.mkdir(parents=True)
        (contents / "Info.plist").write_bytes(plistlib.dumps({
            "CFBundleIdentifier": identity, "CFBundleExecutable": name,
            "CFBundleShortVersionString": builder.__version__,
            "CFBundleVersion": builder.__version__,
            "LSMinimumSystemVersion": builder.MINIMUM_MACOS,
        }))
        (contents / "Resources/Notices").mkdir(parents=True)
        (contents / "Resources/Notices/LICENSE").write_text("MIT fixture")
    contents = directory / "Kmatch.app/Contents"
    for relative in builder.KMATCH_RESOURCES:
        path = contents / "Resources/kmatch" / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture")
    (contents / "Frameworks").mkdir()
    (contents / "Frameworks/kmatch").symlink_to("../Resources/kmatch")
    mmseqs_setup(directory / "Kstack.app" / builder.MMSEQS_RUNTIME)
    return directory


IDENTITY = {
    "built_utc": "20260101T000000000000Z",
    "source": {"commit": "a" * 40, "dirty": False, "files_sha256": {"src/x.py": "0" * 64}},
    "python": {"implementation": "CPython", "version": "3.14.10",
               "distribution": "cpython-3.14.10-macos-aarch64-none"},
}


def stamp_source(directory, source):
    for name in builder.APPS:
        plist = directory / f"{name}.app/Contents/Info.plist"
        info = plistlib.loads(plist.read_bytes())
        info.update(KBioSuiteSourceCommit=source["commit"], KBioSuiteSourceDirty=source["dirty"])
        plist.write_bytes(plistlib.dumps(info))


@pytest.mark.parametrize("problem, error", [
    (None, None),
    ("missing", "Missing Kmatch bundle resource: pymol_runner.py"),
    ("external", "external bundle resource"),
    ("manifest", "Bundled mmseqs does not match installation manifest"),
    ("license", "Bundled mmseqs license does not match installation manifest"),
    ("license-missing", "Bundled mmseqs license does not match installation manifest"),
    ("provenance", "Bundled mmseqs provenance does not match installation manifest"),
    ("absolute", "absolute paths"),
    ("floor", "minimum macOS: Kblast"),
    ("suite-license", "Missing suite license notice: KBioSuite"),
    ("source", "Build source does not match installation manifest"),
])
def test_validate_requires_frozen_resources_and_recorded_tools(
    tmp_path, monkeypatch, problem, error
):
    monkeypatch.setattr(builder, "run", lambda *args: None)
    # Native closure is exercised with real Mach-O files in test_notice_inventory.
    monkeypatch.setattr(builder, "validate_native_inventory", lambda app: None)
    directory = bundle(tmp_path / "apps")
    stamp_source(directory, IDENTITY["source"])
    builder.write_install_manifest(directory, **IDENTITY)
    kstack = directory / "Kstack.app" / builder.MMSEQS_RUNTIME
    manifest = directory / builder.MANIFEST
    package = directory / "Kmatch.app/Contents/Resources/kmatch"
    framework = directory / "Kmatch.app/Contents/Frameworks/kmatch"
    if problem == "missing":
        (package / "pymol_runner.py").unlink()
    elif problem == "external":
        shutil.move(package, tmp_path / "kmatch")
        framework.unlink()
        framework.symlink_to(tmp_path / "kmatch")
    elif problem == "manifest":
        (kstack / "bin/mmseqs").write_bytes(b"replaced")
    elif problem == "license":
        (kstack / "LICENSE.md").write_text("edited")
    elif problem == "license-missing":
        (kstack / "LICENSE.md").unlink()
    elif problem == "provenance":
        (kstack / builder.MMSEQS_PROVENANCE).write_text("{}")
    elif problem == "absolute":
        record = json.loads(manifest.read_text())
        record["bundled_tools"]["Kstack"]["mmseqs"]["source"] = str(Path("/") / "Users" / "synthetic" / ".tools/mmseqs")
        manifest.write_text(json.dumps(record))
    elif problem == "floor":
        plist = directory / "Kblast.app/Contents/Info.plist"
        plist.write_bytes(plistlib.dumps(
            {**plistlib.loads(plist.read_bytes()), "LSMinimumSystemVersion": "14.0"}
        ))
    elif problem == "suite-license":
        (directory / "KBioSuite.app/Contents/Resources/Notices/LICENSE").unlink()
    elif problem == "source":
        stamp_source(directory, {"commit": "b" * 40, "dirty": False})
    if problem is None:
        builder.validate(directory)
    else:
        with pytest.raises(RuntimeError, match=error):
            builder.validate(directory)


def test_manifest_is_relative_and_regenerates_after_resigning(tmp_path, monkeypatch):
    monkeypatch.setattr(builder, "run", lambda *args: None)
    # Native closure is exercised with real Mach-O files in test_notice_inventory.
    monkeypatch.setattr(builder, "validate_native_inventory", lambda app: None)
    directory = bundle(tmp_path / "apps")
    stamp_source(directory, IDENTITY["source"])
    first = builder.write_install_manifest(directory, **IDENTITY)
    assert not list(builder.absolute_paths(first))
    mmseqs = first["bundled_tools"]["Kstack"]["mmseqs"]
    thinned = hashlib.sha256(MMSEQS_BINARY).hexdigest()
    runtime = directory / "Kstack.app" / builder.MMSEQS_RUNTIME
    setup = hashlib.sha256((runtime / builder.MMSEQS_PROVENANCE).read_bytes()).hexdigest()
    assert mmseqs == {
        "path": f"{builder.MMSEQS_RUNTIME}/bin/mmseqs", "sha256": thinned,
        "version": "18-8cc5c", "architecture": "arm64", "archive_sha256": "a" * 64,
        "thinned_sha256": thinned,
        "provenance": {
            "path": f"{builder.MMSEQS_RUNTIME}/{builder.MMSEQS_PROVENANCE}",
            "sha256": setup,
        },
        "licenses": [{
            "path": f"{builder.MMSEQS_RUNTIME}/LICENSE.md",
            "sha256": hashlib.sha256(b"MMseqs2 license fixture").hexdigest(),
        }],
    }
    builder.validate(directory)
    # Release signing rewrites the binary; the stale manifest no longer validates ...
    binary = runtime / "bin/mmseqs"
    binary.write_bytes(MMSEQS_BINARY + b" + Developer ID signature")
    with pytest.raises(RuntimeError, match="Bundled mmseqs does not match"):
        builder.validate(directory)
    # ... until it is regenerated from the products, keeping the recorded build facts.
    regenerated = builder.write_install_manifest(directory)
    builder.validate(directory)
    for key in ("built_utc", "source", "python", "apps", "version"):
        assert regenerated[key] == first[key]
    resigned = regenerated["bundled_tools"]["Kstack"]["mmseqs"]
    assert resigned["sha256"] == hashlib.sha256(binary.read_bytes()).hexdigest() != thinned
    assert resigned["thinned_sha256"] == thinned


def test_manifest_refuses_absolute_paths_and_missing_build_facts(tmp_path):
    directory = bundle(tmp_path / "apps")
    with pytest.raises(RuntimeError, match="without its 'built_utc' record"):
        builder.write_install_manifest(directory)
    leaky = {**IDENTITY, "source": {**IDENTITY["source"], "checkout": str(Path("/") / "Users" / "synthetic" / "KBioSuite")}}
    with pytest.raises(RuntimeError, match="absolute paths"):
        builder.write_install_manifest(directory, **leaky)
    assert not (directory / builder.MANIFEST).exists()
    assert list(builder.absolute_paths({"~/x": ["a", {"b": "/c"}]})) == ["~/x", "/c"]


@pytest.mark.parametrize("required, running, managed, error", [
    ("3.14", "3.14.10", True, None),
    ("3.14", "3.14.11", True, None),
    ("3.14.10", "3.14.10", True, None),
    ("3.14.10", "3.14.11", True, "Build with Python 3.14.10 from .python-version"),
    ("3.14", "3.13.9", True, "Build with Python 3.14 from .python-version"),
    ("3.14", "3.14.10", False, "uv-managed CPython 3.14"),
])
def test_build_python_requires_uv_managed_selected_minor(
    tmp_path, monkeypatch, required, running, managed, error
):
    (tmp_path / ".python-version").write_text(required + "\n")
    uv_root = tmp_path / "uv/python"
    base = (uv_root if managed else tmp_path / "Cellar/python") / f"cpython-{running}-macos"
    base.mkdir(parents=True)
    (base / "BUILD").write_text("20261110\n")
    monkeypatch.setattr(builder, "ROOT", tmp_path)
    monkeypatch.setattr(builder, "uv_python_root", lambda: uv_root.resolve())
    monkeypatch.setattr(builder.platform, "python_version", lambda: running)
    monkeypatch.setattr(builder.sys, "base_prefix", str(base))
    if error:
        with pytest.raises(RuntimeError, match=error):
            builder.build_python()
    else:
        assert builder.build_python() == {
            "implementation": "CPython", "version": running, "distribution": base.name,
            "standalone_build": "20261110"
        }


def test_uv_python_root_asks_uv_then_falls_back(tmp_path, monkeypatch):
    uv = tmp_path / "bin/uv"
    uv.parent.mkdir()
    uv.write_text(f"#!/bin/sh\necho '{tmp_path / 'managed'}'\n")
    uv.chmod(0o755)
    monkeypatch.setenv("UV", str(uv))
    monkeypatch.delenv("UV_PYTHON_INSTALL_DIR", raising=False)
    assert builder.uv_python_root() == (tmp_path / "managed").resolve()
    uv.write_text("#!/bin/sh\nexit 2\n")
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    assert builder.uv_python_root() == (tmp_path / "home/.local/share/uv/python").resolve()
    monkeypatch.setenv("UV_PYTHON_INSTALL_DIR", str(tmp_path / "configured"))
    assert builder.uv_python_root() == (tmp_path / "configured").resolve()


def test_build_keeps_signed_source_provenance_and_removes_intermediate_apps(tmp_path, monkeypatch):
    root = tmp_path / "checkout"
    for name in (
        ".python-version", "LICENSE", "licenses/README.md", "docs/kstack/LICENSE.biopython.rst",
        "docs/kstack/LICENSE.mmseqs.md", "docs/kmatch/af3.md", "resources/kstack/Kstack.icns",
        "apps/kstack-macos/.build/release/Kstack", "apps/kbiosuite-macos/Launcher.swift",
        "apps/kbiosuite-macos/InstallLocation.swift",
        "src/kstack/resources/template.html", "dependency.dist-info/LICENSE",
    ):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture")
    mmseqs_setup(root / ".tools/mmseqs")
    source = {"commit": "a" * 40, "dirty": False}
    python = IDENTITY["python"]
    monkeypatch.setattr(builder, "ROOT", root)
    monkeypatch.setattr(builder, "BUILD_ROOT", root / "build.noindex")
    monkeypatch.setattr(builder.sys, "platform", "darwin")
    monkeypatch.setattr(builder, "source_identity", lambda: source.copy())
    monkeypatch.setattr(builder, "build_python", lambda: python)

    def validate(directory):
        # Validation runs on the finished manifest, so recorded tools are always checked.
        assert (directory / builder.MANIFEST).is_file()

    monkeypatch.setattr(builder, "validate", validate)
    # Exercise the real notice-copying loop, which must not overwrite build identity.
    distribution = SimpleNamespace(
        metadata={"Name": "dependency"}, version="1", files=[Path("dependency.dist-info/LICENSE")],
        locate_file=lambda item: root / item,
    )
    monkeypatch.setattr(builder.importlib.metadata, "distributions", lambda: [distribution])
    notice_path = "python/dependency/dependency.dist-info/LICENSE"
    notice = root / "licenses" / notice_path
    notice.parent.mkdir(parents=True)
    notice.write_text("fixture")
    inventory = {"schema": 1, "open_entries": [],
                 "runtime_identity": builder.collect_notices.__globals__["runtime_identity"](),
                 "components": [{
        "name": "dependency", "version": "1", "kind": "python-distribution",
        "sources": [{"url": "https://example.org/dependency.tar.gz", "sha256": "a" * 64}],
        "notices": [{"path": notice_path, "sha256": hashlib.sha256(b"fixture").hexdigest()}],
    }]}
    (root / "licenses/inventory.json").write_text(json.dumps(inventory))

    frozen_args, compiled = {}, []

    def run(*args):
        if "PyInstaller" in args:
            name = args[args.index("--name") + 1]
            frozen_args[name] = args
            stage = Path(args[args.index("--distpath") + 1]).parent
            frozen = stage / "frozen" / (name if name == "KBioRuntime" else f"{name}.app")
            frozen.mkdir(parents=True)
            (stage / "work").mkdir(exist_ok=True)
            (stage / "specs").mkdir(exist_ok=True)
            if name != "KBioRuntime":
                (frozen / "Contents").mkdir()
        elif args[0] == "/usr/bin/ditto":
            shutil.copytree(args[1], args[2])
        elif args[0] == "swiftc":
            compiled.append(args)
        elif args[0] == "swift" and str(args[1]).endswith("create-icons.swift"):
            args[2].mkdir()
            for name in builder.APPS:
                (args[2] / f"{name}.icns").write_bytes(b"icon")
        elif args[0] == "codesign":
            info = plistlib.loads((args[-1] / "Contents/Info.plist").read_bytes())
            assert info["KBioSuiteSourceCommit"] == source["commit"]
            assert info["KBioSuiteSourceDirty"] is False
            assert info["LSMinimumSystemVersion"] == builder.MINIMUM_MACOS
            assert info["NSHumanReadableCopyright"] == (
                f"{builder.COPYRIGHT}. {builder.LICENSE_NOTICE}."
            )

    monkeypatch.setattr(builder, "run", run)
    # This fixture mocks the compiler/freezer; separate native gate tests use real binaries.
    monkeypatch.setattr(builder, "native_source_index", lambda *args, **kwargs: {})
    monkeypatch.setattr(builder, "audit_native_inventory", lambda *args: None)
    output = builder.build()
    manifest = json.loads((output / builder.MANIFEST).read_text())
    assert manifest["source"] == source and manifest["apps"] == builder.APPS
    assert manifest["python"] == {**python, "redistribution_inventory": "reviewed"}
    assert not list(builder.absolute_paths(manifest))
    assert (root / "build.noindex/latest").resolve() == output
    assert {p.name for p in output.parent.iterdir()} == {"apps", ".kbiosuite-build.json"}
    assert set(frozen_args) == {"Kblast", "Kannotate", "Kmatch", "KBioRuntime"}
    kmatch = frozen_args["Kmatch"]
    assert "--windowed" in kmatch and "webview.platforms.cocoa" in kmatch
    assert kmatch[kmatch.index("--collect-data") + 1] == "kmatch"
    for runner in ("engine_runner.py", "pymol_runner.py"):
        assert f"{root / 'src/kmatch' / runner}:kmatch" in kmatch
    assert "tmtools" in kmatch and "tmtools" not in frozen_args["Kannotate"]
    # The launcher is compiled from every source file for the declared floor.
    (launcher,) = compiled
    assert launcher[launcher.index("-target") + 1] == f"arm64-apple-macos{builder.MINIMUM_MACOS}"
    assert {Path(a).name for a in launcher if str(a).endswith(".swift")} == {
        "Launcher.swift", "InstallLocation.swift"
    }
    # Kmatch no longer bundles uv; MMseqs2 is identified because source_identity() cannot see it.
    assert not (output / "Kmatch.app/Contents/Resources/kmatch/resources/tools").exists()
    runtime = output / "Kstack.app" / builder.MMSEQS_RUNTIME
    assert (runtime / "bin/mmseqs").read_bytes() == MMSEQS_BINARY
    assert (runtime / builder.MMSEQS_PROVENANCE).is_file()
    assert set(manifest["bundled_tools"]) == {"Kstack"}
    recorded = manifest["bundled_tools"]["Kstack"]["mmseqs"]
    assert recorded["thinned_sha256"] == recorded["sha256"]
    assert recorded["version"] == "18-8cc5c" and recorded["archive_sha256"] == "a" * 64
    # Every notice copied from the checkout comes from a fingerprinted source.
    for name in builder.APPS:
        notices = output / f"{name}.app/Contents/Resources/Notices"
        assert (notices / "LICENSE").read_text() == "fixture"
        assert not (notices / "uv").exists()
        copied = {
            p.relative_to(notices).as_posix() for p in notices.rglob("*")
            if p.is_file() and p.relative_to(notices).parts[0] != "dependencies"
        }
        fingerprinted = {
            Path(n).name for n in builder.NOTICE_SOURCES
            if n != "licenses" and (name == "Kmatch" or n not in builder.KMATCH_NOTICES)
        }
        fingerprinted |= {
            p.relative_to(root / "licenses").as_posix() for p in (root / "licenses").rglob("*")
            if p.is_file()
        }
        fingerprinted.add("notice-inventory.json")
        assert copied == fingerprinted
    notices = output / "Kmatch.app/Contents/Resources/Notices"
    assert (notices / "af3.md").read_text() == "fixture"
    for name in set(builder.APPS) - {"Kmatch"}:
        other = output / f"{name}.app/Contents/Resources"
        assert not (other / "Notices/af3.md").exists()
        assert not (other / "kmatch").exists()


class Frozen(Exception):
    """Raised by a mocked freezer: the build got past its runtime checks."""


@pytest.mark.parametrize("exact", [True, False])
def test_unreviewed_runtime_fails_release_before_freezing_but_not_local(
    tmp_path, monkeypatch, capsys, exact
):
    root = tmp_path / "checkout"
    (root / "licenses").mkdir(parents=True)
    identity = builder.collect_notices.__globals__["runtime_identity"]()
    inventory = {"schema": 1, "open_entries": [], "components": [],
                 "runtime_identity": {**identity, "build": "reviewed-build"}}
    (root / "licenses/inventory.json").write_text(json.dumps(inventory))
    monkeypatch.setattr(builder, "ROOT", root)
    monkeypatch.setattr(builder, "BUILD_ROOT", root / "build.noindex")
    monkeypatch.setattr(builder.sys, "platform", "darwin")
    monkeypatch.setattr(builder, "source_identity", lambda: {"commit": "a" * 40, "dirty": False})
    monkeypatch.setattr(builder, "build_python", lambda: IDENTITY["python"])

    def run(*args):
        raise Frozen

    monkeypatch.setattr(builder, "run", run)
    if exact:
        with pytest.raises(RuntimeError, match="Release builds require the reviewed standalone"):
            builder.build(exact_runtime=True)
        assert not (root / "build.noindex").exists()
    else:
        with pytest.raises(Frozen):
            builder.build()
        assert "not covered by the reviewed redistribution inventory" in capsys.readouterr().out


def test_build_refuses_mmseqs_without_matching_setup_record(tmp_path, monkeypatch):
    root = tmp_path / "checkout"
    mmseqs = mmseqs_setup(root / ".tools/mmseqs")
    monkeypatch.setattr(builder, "ROOT", root)
    monkeypatch.setattr(builder, "BUILD_ROOT", root / "build.noindex")
    monkeypatch.setattr(builder.sys, "platform", "darwin")
    monkeypatch.setattr(builder, "source_identity", lambda: {"commit": "a" * 40, "dirty": False})
    monkeypatch.setattr(builder, "build_python", lambda: IDENTITY["python"])

    def run(*args):
        if "PyInstaller" in args:
            name = args[args.index("--name") + 1]
            stage = Path(args[args.index("--distpath") + 1]).parent
            (stage / "frozen" / (name if name == "KBioRuntime" else f"{name}.app")).mkdir(
                parents=True
            )
        elif args[0] == "/usr/bin/ditto":
            shutil.copytree(args[1], args[2])

    monkeypatch.setattr(builder, "run", run)
    (root / "licenses").mkdir(parents=True)
    identity = builder.collect_notices.__globals__["runtime_identity"]()
    (root / "licenses/inventory.json").write_text(json.dumps({"runtime_identity": identity}))
    (root / "src/kstack/resources").mkdir(parents=True)
    (root / "apps/kstack-macos/.build/release").mkdir(parents=True)
    (root / "apps/kstack-macos/.build/release/Kstack").write_text("fixture")
    (mmseqs / "bin/mmseqs").write_bytes(b"universal, not the thinned binary")
    with pytest.raises(RuntimeError, match="differs from its setup record"):
        builder.build()
    shutil.rmtree(root / "build.noindex")  # failed builds stay for diagnosis; not reused here
    (mmseqs / builder.MMSEQS_PROVENANCE).unlink()
    with pytest.raises(RuntimeError, match="setup-mmseqs.sh"):
        builder.build()


def test_stale_build_never_replaces_install(install_environment, monkeypatch):
    built, destination, calls = install_environment
    monkeypatch.setattr(builder, "source_identity", lambda: {"commit": "new-tip", "dirty": False})
    with pytest.raises(RuntimeError, match="do not match current source"):
        builder.install(built, destination)
    assert (destination / "payload").read_text() == "old"
    assert not calls


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def test_source_identity_fingerprints_bundled_notice_sources(tmp_path, monkeypatch):
    root = tmp_path / "checkout"
    files = [
        "src/kbiosuite/version.py", "licenses/README.md", "licenses/LICENSE.boltz",
        "docs/kmatch/af3.md", "docs/kstack/LICENSE.biopython.rst",
        "docs/kstack/LICENSE.mmseqs.md", "docs/kmatch/architecture.md", "notes.md",
    ]
    for name in files:
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"{name} v1")
    git(root.parent, "init", "-q", "-b", "main", str(root))
    git(root, "config", "user.name", "Build test")
    git(root, "config", "user.email", "test@example.invalid")
    git(root, "add", ".")
    git(root, "commit", "-q", "-m", "fixture")
    monkeypatch.setattr(builder, "ROOT", root)
    clean = builder.source_identity()
    expected = set(files) - {"docs/kmatch/architecture.md", "notes.md"}
    assert set(clean["files_sha256"]) == expected and not clean["dirty"]
    for notice in ("docs/kmatch/af3.md", "licenses/README.md", "docs/kstack/LICENSE.mmseqs.md"):
        assert notice in clean["files_sha256"]

    # An uncommitted notice edit or new notice changes the fingerprint of a dirty build.
    (root / "docs/kmatch/af3.md").write_text("edited AlphaFold 3 terms")
    (root / "licenses/LICENSE.new").write_text("untracked notice")
    dirty = builder.source_identity()
    assert dirty["dirty"] and dirty["commit"] == clean["commit"]
    assert dirty["files_sha256"]["docs/kmatch/af3.md"] != clean["files_sha256"]["docs/kmatch/af3.md"]
    assert "licenses/LICENSE.new" in dirty["files_sha256"]


OTOOL = """\
apps/ok:
Mach header
      magic  cputype cpusubtype  caps    filetype ncmds sizeofcmds      flags
MH_MAGIC_64    ARM64        ALL  0x00     EXECUTE    16        744   NOUNDEFS DYLDLINK TWOLEVEL PIE
Load command 0
      cmd LC_SEGMENT_64
  segname __TEXT
Section
  sectname __text
Load command 1
      cmd LC_BUILD_VERSION
  cmdsize 32
 platform MACOS
    minos 26.0
      sdk 27.0
   ntools 1
     tool LD
  version 27037.1
Load command 2
          cmd LC_RPATH
      cmdsize 48
         path @loader_path/../lib with space (offset 12)
Load command 3
          cmd LC_LOAD_DYLIB
         name /usr/lib/libSystem.B.dylib (offset 24)
compatibility version 1.0.0
apps/fat (architecture x86_64):
Mach header
      magic  cputype cpusubtype  caps    filetype ncmds sizeofcmds      flags
MH_MAGIC_64   X86_64        ALL  0x00     EXECUTE    14        728   NOUNDEFS DYLDLINK TWOLEVEL PIE
Load command 0
      cmd LC_VERSION_MIN_MACOSX
  cmdsize 16
  version 10.13
      sdk n/a
Load command 1
      cmd LC_SOURCE_VERSION
  version 0.0
Load command 2
          cmd LC_LOAD_DYLIB
         name /opt/homebrew/opt/openssl@3/lib/libssl.3.dylib (offset 24)
apps/fat (architecture arm64e):
Mach header
      magic  cputype cpusubtype  caps    filetype ncmds sizeofcmds      flags
MH_MAGIC_64    ARM64          E USR00     EXECUTE    16        824   NOUNDEFS DYLDLINK TWOLEVEL PIE
Load command 0
      cmd LC_BUILD_VERSION
    minos 27.0
"""


def test_otool_parsing_covers_every_slice():
    assert builder.parse_otool(OTOOL, ["apps/ok", "apps/fat"]) == {
        "apps/ok": {
            "fat": False, "architectures": ["arm64"],
            "minimums": [("LC_BUILD_VERSION", "26.0")],
            "load_paths": ["@loader_path/../lib with space", "/usr/lib/libSystem.B.dylib"],
        },
        "apps/fat": {
            "fat": True, "architectures": ["x86_64", "arm64e"],
            "minimums": [("LC_VERSION_MIN_MACOSX", "10.13"), ("LC_BUILD_VERSION", "27.0")],
            "load_paths": ["/opt/homebrew/opt/openssl@3/lib/libssl.3.dylib"],
        },
    }


def test_openssl_directories_and_macho_detection(tmp_path):
    data = (
        b'\0OPENSSLDIR: "/opt/homebrew/etc/openssl@3"\0ENGINESDIR: "/tools/deps/engines-3"'
        b'\0/usr/local/Cellar/openssl@3/3.6.4/lib/ossl-modules\0/opt/homebrew/bin/python3\0'
    )
    assert builder.openssl_homebrew_directories(data) == [
        "/opt/homebrew/etc/openssl@3", "/usr/local/Cellar/openssl@3/3.6.4/lib/ossl-modules"
    ]
    assert builder.openssl_homebrew_directories(b'OPENSSLDIR: "/private/etc/ssl"') == []
    samples = {
        "thin": b"\xcf\xfa\xed\xfe\x0c\x00\x00\x01", "fat": b"\xca\xfe\xba\xbe\x00\x00\x00\x02",
        "fat-little": b"\xbe\xba\xfe\xca\x02\x00\x00\x00",
        "java": b"\xca\xfe\xba\xbe\x00\x00\x00\x41", "text": b"#!/bin/sh\n", "empty": b"",
    }
    for name, head in samples.items():
        (tmp_path / name).write_bytes(head)
    assert {n for n in samples if builder.is_macho(tmp_path / n)} == {"thin", "fat", "fat-little"}


TOOLS = ("cc", "otool", "lipo", "vtool")


@pytest.mark.skipif(not all(map(shutil.which, TOOLS)), reason="needs Xcode command line tools")
def test_macho_audit_on_real_binaries(tmp_path):
    source = tmp_path / "main.c"
    source.write_text("int main(void) { return 0; }\n")
    floor = builder.MINIMUM_MACOS
    newer = f"{int(floor.split('.')[0]) + 1}.0"
    contents = tmp_path / "Example.app/Contents"

    def compile_c(path, target, *extra):
        path.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["cc", "-target", target, *extra, str(source), "-o", str(path)],
            check=True, capture_output=True,
        )
        return path

    ok = compile_c(contents / "MacOS/Example", f"arm64-apple-macos{floor}")
    compile_c(contents / "MacOS/helper", "arm64-apple-macos11.0")
    compile_c(contents / "Frameworks/newer", f"arm64-apple-macos{newer}")
    compile_c(contents / "Frameworks/homebrew", f"arm64-apple-macos{floor}",
              "-Wl,-rpath,/opt/homebrew/lib")
    intel = compile_c(tmp_path / "intel", "x86_64-apple-macos10.13")
    subprocess.run(["vtool", "-set-version-min", "macos", newer, newer, "-replace",
                    "-output", str(tmp_path / "intel-new"), str(intel)], check=True)
    subprocess.run(["lipo", "-create", str(ok), str(tmp_path / "intel-new"),
                    "-output", str(contents / "Frameworks/universal")], check=True)
    (contents / "Frameworks/libcrypto.dylib").write_bytes(
        ok.read_bytes() + b'\0OPENSSLDIR: "/opt/homebrew/etc/openssl@3"\0'
    )
    (contents / "Frameworks/benign").write_bytes(
        ok.read_bytes() + b'\0OPENSSLDIR: "/private/etc/ssl"\0/opt/homebrew/bin/python3\0'
    )
    (contents / "Resources").mkdir()
    (contents / "Resources/Example.class").write_bytes(b"\xca\xfe\xba\xbe\x00\x00\x00\x41")
    (contents / "Resources/notes.txt").write_text('OPENSSLDIR: "/opt/homebrew/etc/openssl@3"')
    (contents / "Resources/link").symlink_to("../Frameworks/newer")

    found = {}
    for problem in builder.macho_problems(tmp_path / "Example.app"):
        name, detail = problem.split(": ", 1)
        found.setdefault(name, []).append(detail)
    assert found == {
        "Contents/Frameworks/newer": [f"LC_BUILD_VERSION minimum {newer} exceeds macOS {floor}"],
        "Contents/Frameworks/homebrew": ["Homebrew load path /opt/homebrew/lib"],
        "Contents/Frameworks/universal": [
            "universal x86_64+arm64, not thin arm64",
            f"LC_VERSION_MIN_MACOSX minimum {newer} exceeds macOS {floor}",
        ],
        "Contents/Frameworks/libcrypto.dylib": [
            "Homebrew OpenSSL directory /opt/homebrew/etc/openssl@3"
        ],
    }
