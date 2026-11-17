"""The native launcher's install-location helper, compiled with a harness; no app is launched."""

import shutil
import subprocess
from pathlib import Path

import pytest

from kbiosuite.version import MINIMUM_MACOS

ROOT = Path(__file__).resolve().parents[2]
HARNESS = """
import Foundation
let cases: [(String, Bool?)] = [
    ("/Applications/KBioSuite/KBioSuite.app", false),
    ("/Applications/KBioSuite/KBioSuite.app", nil),
    ("/private/var/folders/xy/T/AppTranslocation/1F2E/d/KBioSuite.app", false),
    ("/Volumes/KBioSuite 0.6.0/KBioSuite/KBioSuite.app", true),
]
for (path, readOnly) in cases {
    print(InstallLocation.problem(bundlePath: path, volumeIsReadOnly: readOnly) ?? "ok")
}
for path in CommandLine.arguments.dropFirst() {
    print(InstallLocation.problem(bundleURL: URL(fileURLWithPath: path)) ?? "ok")
}
"""


@pytest.mark.skipif(shutil.which("swiftc") is None, reason="needs the Swift toolchain")
def test_install_location_flags_translocation_and_read_only_volumes(tmp_path):
    (tmp_path / "main.swift").write_text(HARNESS)
    harness = tmp_path / "harness"
    subprocess.run(
        ["swiftc", "-target", f"arm64-apple-macos{MINIMUM_MACOS}",
         str(ROOT / "apps/kbiosuite-macos/InstallLocation.swift"), str(tmp_path / "main.swift"),
         "-o", str(harness)],
        check=True, capture_output=True,
    )
    writable = tmp_path / "KBioSuite.app"
    writable.mkdir()
    lines = subprocess.run(
        [str(harness), str(writable)], check=True, capture_output=True, text=True
    ).stdout.splitlines()
    assert lines[:2] == ["ok", "ok"]
    assert "temporary, isolated location" in lines[2]
    assert "read-only volume" in lines[3]
    for line in lines[2:4]:
        assert "drag the whole KBioSuite folder into your Applications folder" in line
    assert lines[4:] == ["ok"]  # a writable install location resolves through the URL API
