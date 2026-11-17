"""Verify one application version and macOS floor across imports, documents and bundles."""

import argparse
import importlib
import importlib.metadata
import plistlib
import re
import tomllib
from pathlib import Path

from build_apps import APPS, inspect_macho, macho_files, metadata, version_key

from kbiosuite.version import COPYRIGHT, LICENSE_NOTICE, MINIMUM_MACOS, __version__

ROOT = Path(__file__).resolve().parents[1]
# The README label may wrap between any of its words.
README_LABEL = re.compile(r"(sole\s+suite\s+version\s+source:\s+)\*\*([^*]+)\*\*")
PACKAGE_SWIFT = ROOT / "apps/kstack-macos/Package.swift"


def swift_platform(text):
    """The macOS deployment target a Package.swift declares: .v14, .v10_15 or "26.0"."""
    match = re.search(r'\.macOS\(\s*(?:\.v(\d+(?:_\d+)*)|"(\d+(?:\.\d+)*)")\s*\)', text)
    assert match, "Package.swift declares no macOS platform"
    return match[1].replace("_", ".") if match[1] else match[2]


def check_floor(directory=None):
    """MINIMUM_MACOS bounds Swift targets, bundle metadata and every bundled Mach-O slice."""
    declared = swift_platform(PACKAGE_SWIFT.read_text())
    assert version_key(declared) <= version_key(MINIMUM_MACOS), (
        f"Package.swift targets macOS {declared}, above the suite floor {MINIMUM_MACOS}"
    )
    for name in APPS:
        info = metadata(name)
        assert info["LSMinimumSystemVersion"] == MINIMUM_MACOS, name
        notice = info["NSHumanReadableCopyright"]
        assert COPYRIGHT in notice and LICENSE_NOTICE in notice, name
        if directory:
            bundled = plistlib.loads((directory / f"{name}.app/Contents/Info.plist").read_bytes())
            for key in ("LSMinimumSystemVersion", "NSHumanReadableCopyright"):
                assert bundled.get(key) == info[key], f"{name} {key}"
    if directory:
        files = macho_files(directory)
        records = inspect_macho(files) if files else {}
        for path in files:
            assert str(path) in records, f"otool could not read {path}"
            for command, minimum in records[str(path)]["minimums"]:
                assert version_key(minimum) <= version_key(MINIMUM_MACOS), (
                    f"{path.relative_to(directory)}: {command} minimum {minimum} exceeds "
                    f"macOS {MINIMUM_MACOS}"
                )


def check(directory=None):
    manifest = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert manifest["project"]["dynamic"] == ["version"]
    assert "version" not in manifest["project"]
    assert manifest["tool"]["hatch"]["version"]["path"] == "src/kbiosuite/version.py"
    assert importlib.metadata.version("kbiosuite") == __version__
    for name in ["kbiosuite", "kblast", "kannotate", "kstack", "kmatch"]:
        assert importlib.import_module(name).__version__ == __version__, name
    from kstack.model import create_project

    # Writer metadata is checked with an empty-protein synthetic cluster, avoiding engines.
    project = create_project(
        [
            {
                "id": "check",
                "name": "check",
                "source_name": "check.gbk",
                "accession": "check",
                "description": "synthetic version probe",
                "topology": "linear",
                "annotations": {},
                "warnings": [],
                "length": 100,
                "cds": [],
            }
        ]
    )
    assert project["app_version"] == __version__
    from kstack.persistence import to_html

    assert f'content="Kstack; schema=1; app={__version__}"' in to_html(project)
    text = (ROOT / "README.md").read_text()
    label = README_LABEL.search(text)
    assert label and label[2] == __version__, "README version disagrees"
    if directory:
        for name in ["KBioSuite", "Kblast", "Kannotate", "Kstack", "Kmatch"]:
            info = plistlib.loads((directory / f"{name}.app/Contents/Info.plist").read_bytes())
            assert info["CFBundleShortVersionString"] == __version__, name
            assert info["CFBundleVersion"] == __version__, name
    check_floor(directory)
    print(f"All suite application versions agree: {__version__} (macOS {MINIMUM_MACOS}+)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apps", type=Path)
    parser.add_argument("--update-docs", action="store_true")
    args = parser.parse_args()
    if args.update_docs:
        path = ROOT / "README.md"
        path.write_text(
            README_LABEL.sub(lambda label: f"{label[1]}**{__version__}**", path.read_text())
        )
    check(args.apps)


if __name__ == "__main__":
    main()
