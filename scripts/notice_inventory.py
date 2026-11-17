#!/usr/bin/env python3
"""Verify pinned redistribution notices and write a path-safe public inventory.

The inventory is deliberately a superset of each app's imports. Python distributions
come from the build environment; a new or upgraded distribution needs explicit review.
Open source-audit entries block building rather than silently declaring compliance.

Release builds (exact_runtime=True) require the reviewed standalone CPython build. A local
build may use any uv-managed patch of the same series: it ships that runtime's own license
texts, records its identity and is marked as not covered by the reviewed inventory.
`python -I notice_inventory.py --runtime-identity` prints the running runtime's identity.
"""
from __future__ import annotations

import base64
import hashlib
import importlib.metadata
import json
import platform
import re
import shutil
import sys
from pathlib import Path, PurePosixPath

NOTICE_NAME = re.compile(r"licen[cs]e|copying|notice|authors", re.IGNORECASE)
BUILD_ONLY = {"kbiosuite", "pytest", "ruff", "iniconfig", "pluggy", "pygments"}
REVIEWED_RUNTIME = "reviewed"
UNCOVERED_RUNTIME = "not covered by the reviewed redistribution inventory"
# Component published only by local builds whose runtime differs from the reviewed one.
BUILD_RUNTIME = "cpython-build-runtime"
BUILD_RUNTIME_NOTICES = "cpython/build-runtime"


def normalize(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def safe_path(value: str) -> Path:
    path = PurePosixPath(value)
    if not value or path.is_absolute() or ".." in path.parts or "\\" in value:
        raise RuntimeError(f"Unsafe notice path: {value!r}")
    return Path(*path.parts)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def runtime_identity() -> dict:
    """Identify the managed standalone build without exposing its host path."""
    marker = Path(sys.base_prefix) / "BUILD"
    if not marker.is_file():
        raise RuntimeError("Standalone Python BUILD identity is missing")
    return {"implementation": sys.implementation.name, "version": platform.python_version(),
            "build": marker.read_text().strip(), "platform": sys.platform,
            "machine": platform.machine()}


def describe_runtime(identity: dict) -> str:
    return (f"{identity.get('implementation')} {identity.get('version')} standalone BUILD "
            f"{identity.get('build')} ({identity.get('platform')} {identity.get('machine')})")


def runtime_coverage(inventory: dict, runtime=None, *, exact_runtime=False) -> dict:
    """Whether the build runtime is the reviewed one. Release builds refuse any other runtime;
    a local build records it as not covered by the reviewed inventory instead of failing."""
    reviewed = inventory.get("runtime_identity")
    if not reviewed:
        raise RuntimeError("Notice inventory lacks the reviewed runtime_identity")
    actual = runtime_identity() if runtime is None else runtime
    if actual == reviewed:
        return {"identity": actual, "status": REVIEWED_RUNTIME}
    if exact_runtime:
        raise RuntimeError(
            f"Release builds require the reviewed standalone runtime build "
            f"{describe_runtime(reviewed)}, not {describe_runtime(actual)}. The release "
            "tooling provisions the reviewed runtime; another runtime needs a reviewed "
            "licenses/inventory.json update first."
        )
    return {"identity": actual, "status": UNCOVERED_RUNTIME, "reviewed_identity": reviewed}


def runtime_notice_files(prefix=None) -> dict:
    """The build runtime's own license texts by prefix-relative path (never site-packages)."""
    prefix = Path(prefix or sys.base_prefix).resolve()
    stdlib = prefix / "lib" / f"python{sys.version_info.major}.{sys.version_info.minor}"
    candidates = [p for p in sorted(stdlib.glob("*")) if NOTICE_NAME.search(p.name)]
    # Full standalone archives keep upstream library notices here; install-only ones do not.
    candidates += sorted((prefix / "licenses").rglob("*"))
    found = {p.relative_to(prefix).as_posix(): p for p in candidates
             if p.is_file() and not p.is_symlink()}
    if not found:
        raise RuntimeError("The build runtime has no license text to ship")
    return found


def build_runtime_component(coverage: dict, prefix=None) -> dict:
    """Published mapping for an unreviewed local runtime: its own license texts, plus the
    reviewed CPython notices for the libraries every standalone build links statically."""
    identity = coverage["identity"]
    return {
        "name": BUILD_RUNTIME, "kind": "native-runtime",
        "version": f"{identity['version']}+{identity['build']}",
        "review": UNCOVERED_RUNTIME,
        "notices": [{"path": f"{BUILD_RUNTIME_NOTICES}/{relative}", "sha256": sha256(path)}
                    for relative, path in runtime_notice_files(prefix).items()],
        "notice_dependencies": ["cpython"],
        "sources": [],
    }


def collect_notices(root: Path, destination: Path, distributions=None, *, runtime=None,
                    exact_runtime=False, runtime_prefix=None) -> dict:
    """Validate all mappings before copying any files; no network or host paths emitted.

    The published inventory records the build runtime and whether the reviewed inventory
    covers it; exact_runtime (release builds) refuses an unreviewed runtime instead.
    """
    source = root / "licenses"
    inventory = json.loads((source / "inventory.json").read_text())
    if inventory.get("schema") != 1:
        raise RuntimeError("Unsupported notice inventory schema")
    if inventory.get("open_entries"):
        raise RuntimeError("Unresolved notice inventory: " + "; ".join(inventory["open_entries"]))
    coverage = runtime_coverage(inventory, runtime, exact_runtime=exact_runtime)
    components = inventory["components"]
    mappings = {normalize(c["name"]): c for c in components
                if c.get("kind") == "python-distribution"}
    if distributions is None:
        distributions = importlib.metadata.distributions()
    copies = {}
    for component in components:
        if not component.get("notices") or not component.get("sources"):
            raise RuntimeError(f"Incomplete notice mapping: {component['name']}")
        for archive in component["sources"]:
            if not archive["url"].startswith("https://") or not re.fullmatch(
                r"[0-9a-f]{64}", archive["sha256"]
            ):
                raise RuntimeError(f"Invalid source identity: {component['name']}")
        for notice in component["notices"]:
            relative = safe_path(notice["path"])
            path = source / relative
            if (any(parent.is_symlink() for parent in [path, *path.parents]
                    if parent != source.parent and source.parent in parent.parents)
                    or not path.is_file() or sha256(path) != notice["sha256"]):
                raise RuntimeError(f"Missing or changed notice: {relative}")
            copies[relative] = path
    for distribution in distributions:
        name = normalize(distribution.metadata["Name"])
        if name in BUILD_ONLY:
            continue
        component = mappings.get(name)
        if component is None or component["version"] != distribution.version:
            raise RuntimeError(f"Unmapped distribution: {name}=={distribution.version}")
        # Preserve distribution-relative paths, including nested license directories.
        for item in distribution.files or []:
            if ".dist-info" not in str(item) or not NOTICE_NAME.search(Path(str(item)).name):
                continue
            relative = safe_path(str(item))
            path = Path(distribution.locate_file(item))
            target = Path("python") / name / relative
            if not path.is_file() or target not in copies or sha256(path) != sha256(copies[target]):
                raise RuntimeError(f"Unreviewed distribution notice: {name}/{relative}")
    published = {**inventory, "build_runtime": {key: coverage[key]
                                                for key in ("identity", "status")}}
    if coverage["status"] != REVIEWED_RUNTIME:
        published["components"] = [*components, build_runtime_component(coverage, runtime_prefix)]
        for relative, path in runtime_notice_files(runtime_prefix).items():
            copies[Path(BUILD_RUNTIME_NOTICES) / safe_path(relative)] = path
    destination.mkdir(parents=True, exist_ok=True)
    for relative, source_path in copies.items():
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, target)
    (destination / "notice-inventory.json").write_text(json.dumps(published, indent=2) + "\n")
    return published


MACHO_MAGICS = {
    b"\xfe\xed\xfa\xce", b"\xce\xfa\xed\xfe", b"\xfe\xed\xfa\xcf", b"\xcf\xfa\xed\xfe",
    b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca", b"\xca\xfe\xba\xbf", b"\xbf\xba\xfe\xca",
}
NATIVE_CHECK = "frozen-native-coverage-v1"


def is_native(path: Path) -> bool:
    if not path.is_file():
        return False
    with path.open("rb") as stream:
        return stream.read(4) in MACHO_MAGICS


def native_identity(path: Path) -> dict:
    """Hash actual Mach-O sections per architecture, excluding mutable signatures.

    PyInstaller rewrites load commands; release signing replaces __LINKEDIT code
    signatures. Neither changes section code/data. UUIDs and absolute filenames
    are deliberately absent; identities survive export, thinning and re-signing.
    """
    from macholib.MachO import MachO

    slices = {}
    with path.open("rb") as stream:
        for header in MachO(str(path)).headers:
            digest = hashlib.sha256()
            count = 0
            for _, _command, sections in header.commands:
                if not isinstance(sections, list):
                    continue
                for section in sections:
                    count += 1
                    description = [section.segname.decode().rstrip("\0"),
                                   section.sectname.decode().rstrip("\0"),
                                   section.size, section.flags]
                    digest.update(json.dumps(description).encode())
                    # S_ZEROFILL, S_GB_ZEROFILL and S_THREAD_LOCAL_ZEROFILL.
                    if section.flags & 0xff not in {1, 12, 18}:
                        stream.seek(header.offset + section.offset)
                        data = stream.read(section.size)
                        if len(data) != section.size:
                            raise RuntimeError("Truncated Mach-O section")
                        digest.update(data)
            if not count:
                raise RuntimeError("Native binary contains no attributable sections")
            slices[f"{header.header.cputype}:{header.header.cpusubtype}"] = digest.hexdigest()
    return slices


def native_files(app: Path) -> list[Path]:
    """Unique regular Mach-O files; bundle-internal aliases remain aliases."""
    paths = []
    for path in sorted(app.rglob("*")):
        if path.is_symlink():
            if not path.exists() or not path.resolve().is_relative_to(app.resolve()):
                raise RuntimeError("Native inventory refuses external/broken bundle symlink")
            continue
        if is_native(path):
            paths.append(path)
    return paths


def component_notices(component: dict, components: dict) -> list[str]:
    """Retain header/static-library notices as part of the consumer's mapping."""
    required = [component, *(components[normalize(name)]
                             for name in component.get("notice_dependencies", []))]
    return list(dict.fromkeys(n["path"] for c in required for n in c["notices"]))


def native_source_index(root: Path, inventory: dict, *, suite_binaries=(), distributions=None,
                        runtime_prefix=None, exact_runtime=False) -> dict:
    """Local-only source evidence. Absolute source paths are never serialized.

    RECORD ownership gives extensions/libraries their exact reviewed distribution.
    Standalone runtime ownership is confined to its verified build prefix, excluding
    site-packages. No basename or directory-name alone grants native ownership.
    The reviewed runtime must match its archive-pinned identities; an unreviewed local
    runtime's own files belong to the published, explicitly uncovered build-runtime record.
    """
    components = {normalize(c["name"]): c for c in inventory["components"]}
    index = {}
    verified_distribution_natives = set()

    def register(path, component):
        path = Path(path)
        if not is_native(path):
            return
        for architecture, identity in native_identity(path).items():
            key = (architecture, identity)
            owner = {"component": component["name"], "version": component["version"],
                     "notices": component_notices(component, components)}
            # Identical helpers may appear in multiple distributions; all static
            # notices are retained, and this records the selected verified owner.
            index[key] = owner

    if distributions is None:
        distributions = importlib.metadata.distributions()
    for distribution in distributions:
        name = normalize(distribution.metadata["Name"])
        if name in BUILD_ONLY:
            continue
        component = components.get(name)
        if component is None or distribution.version != component["version"]:
            raise RuntimeError(f"Unmapped native distribution: {name}")
        for relative in distribution.files or []:
            source_path = Path(distribution.locate_file(relative))
            if is_native(source_path):
                expected = relative.hash
                if expected is None or expected.mode != "sha256" or (
                    base64.urlsafe_b64encode(bytes.fromhex(sha256(source_path))).decode().rstrip("=")
                    != expected.value
                ):
                    raise RuntimeError(f"Native source differs from distribution RECORD: {name}")
                verified_distribution_natives.add(source_path.resolve())
                register(source_path, component)
    prefix = Path(runtime_prefix or sys.base_prefix).resolve()
    coverage = runtime_coverage(inventory, exact_runtime=exact_runtime)
    if coverage["status"] == REVIEWED_RUNTIME:
        runtime_sources = inventory.get("runtime_native_sources")
        if not runtime_sources:
            raise RuntimeError("No reviewed standalone runtime native identities")
        for relative, identity in runtime_sources.items():
            path = prefix / safe_path(relative)
            if not path.is_file() or native_identity(path) != identity:
                raise RuntimeError(
                    f"Standalone native source differs from pinned archive: {relative}"
                )
            register(path, components["cpython"])
    else:
        local = build_runtime_component(coverage, prefix)
        components[normalize(BUILD_RUNTIME)] = local
        series = f"python{sys.version_info.major}.{sys.version_info.minor}"
        site_packages = prefix / "lib" / series / "site-packages"
        for path in sorted(prefix.rglob("*")):
            if not path.is_symlink() and not path.is_relative_to(site_packages):
                register(path, local)
    # Bootloader is part of the PyInstaller RECORD but has explicit exception terms.
    import PyInstaller
    bootloaders = Path(PyInstaller.__file__).parent / "bootloader"
    for path in bootloaders.rglob("*"):
        if is_native(path):
            if path.resolve() not in verified_distribution_natives:
                raise RuntimeError("Unreviewed PyInstaller bootloader source outside RECORD")
            register(path, components["pyinstaller-bootloader"])
    mmseqs = root / ".tools/mmseqs"
    records = list(mmseqs.glob("*.json"))
    binary = mmseqs / "bin/mmseqs"
    matched = False
    for receipt in records:
        data = json.loads(receipt.read_text())
        if data.get("version") == components["mmseqs2"]["version"] and (
            data.get("binary_sha256") == sha256(binary)
            and data.get("archive_sha256") ==
            "158102d4c3e3ce721424f455106fac348cd01d424ed5e66cb03eac3de5064eb6"
        ):
            matched = True
            break
    if not matched:
        raise RuntimeError("MMseqs native source lacks matching pinned provisioning receipt")
    register(binary, components["mmseqs2"])
    from kbiosuite.version import __version__
    suite = {"name": "kbiosuite", "version": __version__, "notices": [{"path": "LICENSE"}]}
    for path in suite_binaries:
        register(path, suite)
    return index


def audit_native_inventory(app: Path, index: dict) -> dict:
    """Map the complete frozen closure before signing; unknown code fails closed."""
    base = app / "Contents/Resources/Notices"
    records = []
    for path in native_files(app):
        identity = native_identity(path)
        owners = []
        for architecture, digest in identity.items():
            owner = index.get((architecture, digest))
            if owner is None:
                raise RuntimeError(f"Unmapped frozen native component: {path.relative_to(app)}")
            if owner not in owners:
                owners.append(owner)
        records.append({"path": path.relative_to(app).as_posix(),
                        "identity": identity, "owners": owners})
    if not records:
        raise RuntimeError("App has no verified native closure")
    receipt = {"schema": 1, "check": NATIVE_CHECK, "status": "passed",
               "notice_inventory_sha256": sha256(base / "notice-inventory.json"),
               "suite_license_sha256": sha256(base / "LICENSE"),
               "native_files": records}
    (base / "native-inventory.json").write_text(json.dumps(receipt, indent=2) + "\n")
    validate_native_inventory(app)
    return receipt


def validate_native_inventory(app: Path) -> dict:
    """Validate exported or re-signed app without accessing its build environment."""
    base = app / "Contents/Resources/Notices"
    if not (base / "native-inventory.json").is_file():
        raise RuntimeError("Missing required frozen-native coverage receipt")
    inventory = json.loads((base / "notice-inventory.json").read_text())
    receipt = json.loads((base / "native-inventory.json").read_text())
    if inventory.get("open_entries") or receipt.get("schema") != 1 or (
        inventory.get("required_checks") != [NATIVE_CHECK]
    ) or (
        receipt.get("check") != NATIVE_CHECK or receipt.get("status") != "passed"
        or receipt.get("notice_inventory_sha256") != sha256(base / "notice-inventory.json")
        or receipt.get("suite_license_sha256") != sha256(base / "LICENSE")
    ):
        raise RuntimeError("Native coverage receipt does not match notice inventory")
    expected = {r["path"]: r for r in receipt["native_files"]}
    actual = {p.relative_to(app).as_posix(): p for p in native_files(app)}
    if set(actual) != set(expected) or len(expected) != len(receipt["native_files"]):
        raise RuntimeError("Frozen native file set differs from reviewed inventory")
    components = {normalize(c["name"]): c for c in inventory["components"]}
    for relative, path in actual.items():
        record = expected[relative]
        safe_path(relative)
        if native_identity(path) != record["identity"] or not record["owners"]:
            raise RuntimeError(f"Frozen native identity changed: {relative}")
        for owner in record["owners"]:
            if owner["component"] == "kbiosuite":
                allowed = ["LICENSE"]
            else:
                component = components.get(normalize(owner["component"]))
                if component is None or component["version"] != owner["version"]:
                    raise RuntimeError("Native owner lacks reviewed version mapping")
                allowed = component_notices(component, components)
            if owner["notices"] != allowed:
                raise RuntimeError("Native owner notice mapping changed")
            for notice in allowed:
                target = base / safe_path(notice)
                if not target.is_file():
                    raise RuntimeError(f"Missing native notice: {notice}")
                if owner["component"] != "kbiosuite":
                    expected_notice = next(n for c in components.values()
                                           for n in c["notices"] if n["path"] == notice)
                    if sha256(target) != expected_notice["sha256"]:
                        raise RuntimeError(f"Changed native notice: {notice}")
    return receipt


if __name__ == "__main__":
    # Release tooling asks a freshly provisioned interpreter (run with -I) for its identity.
    if sys.argv[1:] != ["--runtime-identity"]:
        raise SystemExit("Usage: notice_inventory.py --runtime-identity")
    print(json.dumps(runtime_identity(), sort_keys=True))
