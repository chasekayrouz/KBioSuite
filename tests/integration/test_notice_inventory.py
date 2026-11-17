"""Offline adversarial fixtures for the redistribution gate."""
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/notice_inventory.py"
spec = importlib.util.spec_from_file_location("notice_inventory", SCRIPT)
notices = importlib.util.module_from_spec(spec)
spec.loader.exec_module(notices)


@pytest.fixture
def fixture(tmp_path):
    relative = "example-1.dist-info/licenses/nested/LICENCE.rst"
    original = tmp_path / "env" / relative
    original.parent.mkdir(parents=True)
    original.write_text("example permission\n")
    target = tmp_path / "licenses/python/example" / relative
    target.parent.mkdir(parents=True)
    target.write_bytes(original.read_bytes())
    data = {"schema": 1, "open_entries": [], "runtime_identity": notices.runtime_identity(), "components": [{
        "name": "example", "version": "1", "kind": "python-distribution",
        "sources": [{"url": "https://example.org/source.tar.gz", "sha256": "a" * 64}],
        "notices": [{"path": str(target.relative_to(tmp_path / "licenses")),
                     "sha256": hashlib.sha256(target.read_bytes()).hexdigest()}],
    }]}
    manifest = tmp_path / "licenses/inventory.json"
    manifest.write_text(json.dumps(data))
    dist = SimpleNamespace(metadata={"Name": "example"}, version="1",
                           files=[Path(relative)], locate_file=lambda p: tmp_path / "env" / p)
    return tmp_path, data, manifest, dist, target


def test_preserves_nested_distribution_paths(fixture):
    root, _, _, dist, _ = fixture
    notices.collect_notices(root, root / "result", [dist])
    path = root / "result/python/example/example-1.dist-info/licenses/nested/LICENCE.rst"
    assert path.read_text() == "example permission\n"
    public = (root / "result/notice-inventory.json").read_text()
    assert str(root) not in public


@pytest.mark.parametrize("failure", ["unknown", "version", "missing", "changed", "open", "escape"])
def test_fails_before_writing(fixture, failure):
    root, data, manifest, dist, target = fixture
    if failure == "unknown":
        dist.metadata["Name"] = "new-package"
    elif failure == "version":
        dist.version = "2"
    elif failure == "missing":
        target.unlink()
    elif failure == "changed":
        target.write_text("different terms")
    elif failure == "open":
        data["open_entries"] = ["pending audit"]
    else:
        data["components"][0]["notices"][0]["path"] = "../outside"
    manifest.write_text(json.dumps(data))
    with pytest.raises(RuntimeError):
        notices.collect_notices(root, root / "result", [dist])
    assert not (root / "result").exists()


@pytest.mark.parametrize("name", ["LICENCE.rst", "COPYING", "AUTHORS.txt", "NOTICE"])
def test_notice_spellings(name):
    assert notices.NOTICE_NAME.search(name)


def test_release_build_rejects_different_standalone_release(fixture):
    root, _, _, dist, _ = fixture
    runtime = notices.runtime_identity()
    runtime["build"] = "different-release"
    with pytest.raises(RuntimeError, match="runtime build"):
        notices.collect_notices(root, root / "result", [dist], runtime=runtime, exact_runtime=True)
    assert not (root / "result").exists()


def test_reviewed_runtime_is_recorded_without_extra_component(fixture):
    root, data, _, dist, _ = fixture
    published = notices.collect_notices(root, root / "result", [dist], exact_runtime=True)
    assert published["build_runtime"] == {
        "identity": notices.runtime_identity(), "status": notices.REVIEWED_RUNTIME
    }
    assert published["components"] == data["components"]
    assert not (root / "result" / notices.BUILD_RUNTIME_NOTICES).exists()


def fake_runtime(prefix):
    """A standalone-like prefix: own license, one native file, one site-packages native."""
    import shutil
    import sys
    stdlib = prefix / "lib" / f"python{sys.version_info.major}.{sys.version_info.minor}"
    (stdlib / "site-packages").mkdir(parents=True)
    (stdlib / "LICENSE.txt").write_text("local runtime terms\n")
    (stdlib / "os.py").write_text("# not a notice\n")
    (prefix / "bin").mkdir()
    shutil.copyfile("/bin/echo", prefix / "bin/python3")
    shutil.copyfile("/bin/cat", stdlib / "site-packages/extension.so")
    return stdlib


def test_local_build_ships_its_runtime_license_and_is_marked_uncovered(fixture, tmp_path):
    root, _, _, dist, _ = fixture
    stdlib = fake_runtime(tmp_path / "runtime")
    runtime = {**notices.runtime_identity(), "version": "3.14.99", "build": "local-build"}
    published = notices.collect_notices(root, root / "result", [dist], runtime=runtime,
                                        runtime_prefix=tmp_path / "runtime")
    assert published["build_runtime"] == {"identity": runtime, "status": notices.UNCOVERED_RUNTIME}
    component = published["components"][-1]
    relative = f"{notices.BUILD_RUNTIME_NOTICES}/{stdlib.relative_to(tmp_path / 'runtime')}"
    assert component["name"] == notices.BUILD_RUNTIME
    assert component["version"] == "3.14.99+local-build"
    assert component["review"] == notices.UNCOVERED_RUNTIME
    assert component["notices"] == [{"path": f"{relative}/LICENSE.txt",
                                     "sha256": notices.sha256(stdlib / "LICENSE.txt")}]
    assert (root / "result" / relative / "LICENSE.txt").read_text() == "local runtime terms\n"
    written = json.loads((root / "result/notice-inventory.json").read_text())
    assert written == published and str(tmp_path) not in json.dumps(written)


def test_ignores_legitimate_console_script_record(fixture):
    root, _, _, dist, _ = fixture
    dist.files.append(Path("../../../bin/example"))
    notices.collect_notices(root, root / "result", [dist])


def test_rejects_unreviewed_authors_file(fixture):
    root, _, _, dist, _ = fixture
    relative = Path("example-1.dist-info/AUTHORS.txt")
    path = root / "env" / relative
    path.write_text("new attribution")
    dist.files.append(relative)
    with pytest.raises(RuntimeError, match="Unreviewed distribution notice"):
        notices.collect_notices(root, root / "result", [dist])


@pytest.fixture
def native_app(tmp_path):
    import shutil
    import sys
    if sys.platform != "darwin":
        pytest.skip("Mach-O test requires macOS")
    app = tmp_path / "application"
    binary = app / "Contents/MacOS/example"
    binary.parent.mkdir(parents=True)
    shutil.copyfile("/bin/echo", binary)
    base = app / "Contents/Resources/Notices"
    base.mkdir(parents=True)
    (base / "LICENSE").write_text("suite permission")
    (base / "example-LICENSE").write_text("example permission")
    component = {"name": "example", "version": "1",
                 "notices": [{"path": "example-LICENSE",
                              "sha256": notices.sha256(base / "example-LICENSE")}]}
    inventory = {"schema": 1, "open_entries": [], "required_checks": [notices.NATIVE_CHECK],
                 "components": [component]}
    (base / "notice-inventory.json").write_text(json.dumps(inventory))
    owner = {"component": "example", "version": "1", "notices": ["example-LICENSE"]}
    index = {(arch, digest): owner for arch, digest in notices.native_identity(binary).items()}
    return app, binary, index


def test_native_unknown_component_fails_before_receipt(native_app):
    app, _, _ = native_app
    with pytest.raises(RuntimeError, match="Unmapped frozen native"):
        notices.audit_native_inventory(app, {})
    assert not (app / "Contents/Resources/Notices/native-inventory.json").exists()


def test_exported_and_resigned_native_inventory(native_app, tmp_path):
    import shutil
    import subprocess
    app, binary, index = native_app
    notices.audit_native_inventory(app, index)
    exported = tmp_path / "exported"
    shutil.copytree(app, exported)
    subprocess.run(["/usr/bin/codesign", "--force", "--sign", "-",
                    str(exported / binary.relative_to(app))], check=True, capture_output=True)
    notices.validate_native_inventory(exported)
    public = (exported / "Contents/Resources/Notices/native-inventory.json").read_text()
    assert str(tmp_path) not in public
    assert "passed" in public


def test_extra_native_library_rejected(native_app):
    import shutil
    app, _, index = native_app
    notices.audit_native_inventory(app, index)
    shutil.copyfile("/bin/cat", app / "Contents/MacOS/unknown.dylib")
    with pytest.raises(RuntimeError, match="file set differs"):
        notices.validate_native_inventory(app)
    with pytest.raises(RuntimeError, match="Unmapped frozen native"):
        notices.audit_native_inventory(app, index)


def test_changed_native_identity_rejected(native_app):
    import shutil
    app, binary, index = native_app
    notices.audit_native_inventory(app, index)
    shutil.copyfile("/bin/cat", binary)
    with pytest.raises(RuntimeError, match="identity changed"):
        notices.validate_native_inventory(app)


def test_changed_native_notice_rejected(native_app):
    app, _, index = native_app
    notices.audit_native_inventory(app, index)
    (app / "Contents/Resources/Notices/example-LICENSE").write_text("different")
    with pytest.raises(RuntimeError, match="Changed native notice"):
        notices.validate_native_inventory(app)


def test_native_record_hash_must_match(tmp_path):
    import shutil
    from importlib.metadata import PackagePath
    binary = tmp_path / "unreviewed.dylib"
    shutil.copyfile("/bin/echo", binary)
    relative = PackagePath(binary.name)
    relative.hash = SimpleNamespace(mode="sha256", value="changed-record")
    dist = SimpleNamespace(metadata={"Name": "example"}, version="1", files=[relative],
                           locate_file=lambda item: tmp_path / item)
    inventory = {"components": [{"name": "example", "version": "1", "notices": []}]}
    with pytest.raises(RuntimeError, match="differs from distribution RECORD"):
        notices.native_source_index(tmp_path, inventory, distributions=[dist])


def test_unpinned_runtime_file_cannot_gain_runtime_ownership(tmp_path, monkeypatch):
    import PyInstaller
    monkeypatch.setattr(PyInstaller, "__file__", str(tmp_path / "PyInstaller/__init__.py"))
    import shutil
    prefix = tmp_path / "runtime"
    (prefix / "bin").mkdir(parents=True)
    known = prefix / "bin/python"
    shutil.copyfile("/bin/echo", known)
    unknown = prefix / "bin/extra"
    shutil.copyfile("/bin/cat", unknown)
    mmseqs = tmp_path / ".tools/mmseqs"
    (mmseqs / "bin").mkdir(parents=True)
    shutil.copyfile(known, mmseqs / "bin/mmseqs")
    (mmseqs / "provenance.json").write_text(json.dumps({
        "version": "18-8cc5c", "binary_sha256": notices.sha256(known),
        "archive_sha256": "158102d4c3e3ce721424f455106fac348cd01d424ed5e66cb03eac3de5064eb6",
    }))
    inventory = {"runtime_identity": notices.runtime_identity(),
                 "runtime_native_sources": {"bin/python": notices.native_identity(known)},
                 "components": [{"name": name, "version": version, "notices": []}
                                for name, version in [("cpython", "3.14.10+20261110"),
                                                      ("pyinstaller-bootloader", "6.22.3"),
                                                      ("mmseqs2", "18-8cc5c")]]}
    index = notices.native_source_index(tmp_path, inventory, distributions=[], runtime_prefix=prefix)
    assert all((arch, digest) not in index
               for arch, digest in notices.native_identity(unknown).items())


def test_unrecorded_bootloader_native_rejected(tmp_path, monkeypatch):
    import shutil

    import PyInstaller
    package = tmp_path / "PyInstaller"
    bootloaders = package / "bootloader/Darwin-64bit"
    bootloaders.mkdir(parents=True)
    shutil.copyfile("/bin/echo", bootloaders / "extra-helper")
    monkeypatch.setattr(PyInstaller, "__file__", str(package / "__init__.py"))
    inventory = json.loads((SCRIPT.parents[1] / "licenses/inventory.json").read_text())
    with pytest.raises(RuntimeError, match="bootloader source outside RECORD"):
        notices.native_source_index(tmp_path, inventory, distributions=[])


def test_local_runtime_natives_map_to_uncovered_record_and_validate(tmp_path, monkeypatch):
    """Native index, published inventory and validation agree for an unreviewed local runtime."""
    import shutil
    import sys
    if sys.platform != "darwin":
        pytest.skip("Mach-O test requires macOS")
    import PyInstaller
    monkeypatch.setattr(PyInstaller, "__file__", str(tmp_path / "PyInstaller/__init__.py"))
    prefix = tmp_path / "runtime"
    stdlib = fake_runtime(prefix)
    root = tmp_path / "checkout"
    licenses = root / "licenses"
    entries = {}
    for name, version in [("cpython", "3.14.10+reviewed"), ("mmseqs2", "18-8cc5c")]:
        notice = licenses / name / "LICENSE"
        notice.parent.mkdir(parents=True)
        notice.write_text(f"{name} terms")
        entries[name] = {"name": name, "version": version, "kind": "native",
                         "sources": [{"url": "https://example.org/a.tar.gz", "sha256": "a" * 64}],
                         "notices": [{"path": f"{name}/LICENSE", "sha256": notices.sha256(notice)}]}
    mmseqs = root / ".tools/mmseqs"
    (mmseqs / "bin").mkdir(parents=True)
    shutil.copyfile("/usr/bin/true", mmseqs / "bin/mmseqs")
    (mmseqs / "provenance.json").write_text(json.dumps({
        "version": "18-8cc5c", "binary_sha256": notices.sha256(mmseqs / "bin/mmseqs"),
        "archive_sha256": "158102d4c3e3ce721424f455106fac348cd01d424ed5e66cb03eac3de5064eb6",
    }))
    inventory = {"schema": 1, "open_entries": [], "required_checks": [notices.NATIVE_CHECK],
                 "runtime_identity": {**notices.runtime_identity(), "build": "reviewed-build"},
                 "runtime_native_sources": {}, "components": list(entries.values())}
    (licenses / "inventory.json").write_text(json.dumps(inventory))
    with pytest.raises(RuntimeError, match="runtime build"):
        notices.native_source_index(root, inventory, distributions=[], runtime_prefix=prefix,
                                    exact_runtime=True)
    index = notices.native_source_index(root, inventory, distributions=[], runtime_prefix=prefix)
    owners = {index[key]["component"] for key in notices.native_identity(prefix / "bin/python3").items()}
    assert owners == {notices.BUILD_RUNTIME}
    assert all(key not in index
               for key in notices.native_identity(stdlib / "site-packages/extension.so").items())
    app = tmp_path / "Example.app"
    binary = app / "Contents/MacOS/Example"
    binary.parent.mkdir(parents=True)
    shutil.copyfile(prefix / "bin/python3", binary)
    base = app / "Contents/Resources/Notices"
    published = notices.collect_notices(root, base, [], runtime_prefix=prefix)
    (base / "LICENSE").write_text("suite permission")
    assert published["build_runtime"]["status"] == notices.UNCOVERED_RUNTIME
    receipt = notices.audit_native_inventory(app, index)
    (record,) = receipt["native_files"]
    assert record["owners"][0]["notices"][-1] == "cpython/LICENSE"
    notices.validate_native_inventory(app)
    (base / notices.BUILD_RUNTIME_NOTICES / stdlib.relative_to(prefix) / "LICENSE.txt").write_text(
        "changed"
    )
    with pytest.raises(RuntimeError, match="Changed native notice"):
        notices.validate_native_inventory(app)
