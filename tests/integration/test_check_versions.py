"""The version checker's macOS-floor rules: Swift targets, bundle metadata, Mach-O minimums."""

import importlib.util
import plistlib
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))  # the checker imports the builder's helpers
try:
    spec = importlib.util.spec_from_file_location("check_versions", ROOT / "scripts/check_versions.py")
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
finally:
    sys.path.remove(str(ROOT / "scripts"))

FLOOR = checker.MINIMUM_MACOS


@pytest.mark.parametrize("declaration, expected", [
    ("platforms: [.macOS(.v14)],", "14"),
    ("platforms: [.macOS(.v10_15), .iOS(.v17)],", "10.15"),
    ('platforms: [.macOS("26.0")],', "26.0"),
])
def test_swift_platform_forms(declaration, expected):
    assert checker.swift_platform(declaration) == expected


def test_source_floor_agrees_and_package_above_floor_is_rejected(tmp_path, monkeypatch):
    checker.check_floor()
    package = tmp_path / "Package.swift"
    package.write_text('platforms: [.macOS("99.0")],')
    monkeypatch.setattr(checker, "PACKAGE_SWIFT", package)
    with pytest.raises(AssertionError, match="above the suite floor"):
        checker.check_floor()
    package.write_text("platforms: [],")
    with pytest.raises(AssertionError, match="no macOS platform"):
        checker.check_floor()


def suite(directory):
    for name in checker.APPS:
        contents = directory / f"{name}.app/Contents"
        contents.mkdir(parents=True)
        (contents / "Info.plist").write_bytes(plistlib.dumps(checker.metadata(name)))
    return directory


def test_bundle_metadata_must_carry_floor_and_copyright(tmp_path):
    directory = suite(tmp_path / "apps")
    checker.check_floor(directory)
    plist = directory / "Kstack.app/Contents/Info.plist"
    info = plistlib.loads(plist.read_bytes())
    plist.write_bytes(plistlib.dumps({**info, "LSMinimumSystemVersion": "14.0"}))
    with pytest.raises(AssertionError, match="Kstack LSMinimumSystemVersion"):
        checker.check_floor(directory)
    plist.write_bytes(plistlib.dumps({k: v for k, v in info.items()
                                      if k != "NSHumanReadableCopyright"}))
    with pytest.raises(AssertionError, match="Kstack NSHumanReadableCopyright"):
        checker.check_floor(directory)
    assert checker.COPYRIGHT in info["NSHumanReadableCopyright"]
    assert checker.LICENSE_NOTICE in info["NSHumanReadableCopyright"]


@pytest.mark.skipif(not all(map(shutil.which, ("cc", "otool"))), reason="needs Xcode tools")
def test_bundled_binary_above_floor_is_rejected(tmp_path):
    directory = suite(tmp_path / "apps")
    source = tmp_path / "main.c"
    source.write_text("int main(void) { return 0; }\n")
    binary = directory / "Kmatch.app/Contents/MacOS/Kmatch"
    binary.parent.mkdir()
    for target, accepted in ((FLOOR, True), (f"{int(FLOOR.split('.')[0]) + 1}.0", False)):
        subprocess.run(["cc", "-target", f"arm64-apple-macos{target}", str(source),
                        "-o", str(binary)], check=True, capture_output=True)
        if accepted:
            checker.check_floor(directory)
        else:
            with pytest.raises(AssertionError, match="LC_BUILD_VERSION minimum"):
                checker.check_floor(directory)


def test_readme_label_survives_wrapping_and_updates_in_place():
    text = "`version.py` is the sole suite\nversion source: **0.0.1**.\n"
    assert checker.README_LABEL.search(text)[2] == "0.0.1"
    updated = checker.README_LABEL.sub(lambda label: f"{label[1]}**9.9.9**", text)
    assert updated == "`version.py` is the sole suite\nversion source: **9.9.9**.\n"
