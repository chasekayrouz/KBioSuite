"""Build/sign five independent apps, then promote a verified suite installation."""

import argparse
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import platform
import plistlib
import re
import shutil
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

from kbiosuite.version import COPYRIGHT, LICENSE_NOTICE, MINIMUM_MACOS, __version__

try:
    from notice_inventory import (
        REVIEWED_RUNTIME,
        audit_native_inventory,
        collect_notices,
        describe_runtime,
        native_source_index,
        runtime_coverage,
        validate_native_inventory,
    )
except ModuleNotFoundError:  # imported by file path in offline build tests
    _notice_spec = importlib.util.spec_from_file_location(
        "notice_inventory", Path(__file__).with_name("notice_inventory.py")
    )
    _notice_module = importlib.util.module_from_spec(_notice_spec)
    _notice_spec.loader.exec_module(_notice_module)
    REVIEWED_RUNTIME = _notice_module.REVIEWED_RUNTIME
    collect_notices = _notice_module.collect_notices
    audit_native_inventory = _notice_module.audit_native_inventory
    describe_runtime = _notice_module.describe_runtime
    native_source_index = _notice_module.native_source_index
    runtime_coverage = _notice_module.runtime_coverage
    validate_native_inventory = _notice_module.validate_native_inventory

ROOT = Path(__file__).resolve().parents[1]
BUILD_ROOT = ROOT / "build.noindex"
APPS = {
    "KBioSuite": "org.kbiosuite.launcher",
    "Kblast": "org.kbiosuite.kblast",
    "Kannotate": "org.kbiosuite.kannotate",
    "Kstack": "org.kbiosuite.kstack",
    "Kmatch": "org.kbiosuite.kmatch",
}
# Kmatch resolves these through Contents/Frameworks/kmatch (PyInstaller's link to Resources).
KMATCH_RESOURCES = (
    "engine_runner.py",
    "pymol_runner.py",
    "static/index.html",
    "static/vendor/3Dmol-min.js",
    # Upstream 3Dmol.js notices travel beside the minified build they cover.
    "static/vendor/3Dmol-LICENSE",
    "static/vendor/3Dmol-min.js.LICENSE.txt",
    "resources/boltz-requirements.txt",
)
# Notice texts copied from the checkout into bundles; source_identity() fingerprints them so a
# --allow-uncommitted build can never ship notices that its manifest does not identify.
SHARED_NOTICES = ("docs/kstack/LICENSE.biopython.rst", "docs/kstack/LICENSE.mmseqs.md")
KMATCH_NOTICES = ("docs/kmatch/af3.md",)
SUITE_LICENSE = "LICENSE"  # the suite's own MIT license, copied to every app's Notices/LICENSE
NOTICE_SOURCES = ("licenses", SUITE_LICENSE, *SHARED_NOTICES, *KMATCH_NOTICES)
# Thin arm64 MMseqs2 plus the record scripts/kstack/setup-mmseqs.sh writes beside it (release,
# archive SHA-256, thinned-binary SHA-256); copied into Kstack's frozen runtime.
MMSEQS_RUNTIME = "Contents/Resources/runtime/mmseqs"
MMSEQS_PROVENANCE = "provenance.json"
MANIFEST = ".kbiosuite-install.json"
# Every shipped Mach-O must be thin arm64, run on MINIMUM_MACOS and never reach into Homebrew
# (a frozen Homebrew OpenSSL would look for its certificates under the build machine's prefix).
MACHO_MAGICS = (b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe", b"\xfe\xed\xfa\xcf", b"\xfe\xed\xfa\xce")
FAT_MAGICS = (
    b"\xca\xfe\xba\xbe", b"\xca\xfe\xba\xbf", b"\xbe\xba\xfe\xca", b"\xbf\xba\xfe\xca"
)
HOMEBREW_PREFIXES = ("/opt/homebrew", "/usr/local/Cellar")
PATH_COMMANDS = {
    "LC_ID_DYLIB", "LC_LOAD_DYLIB", "LC_LOAD_WEAK_DYLIB", "LC_REEXPORT_DYLIB",
    "LC_LAZY_LOAD_DYLIB", "LC_LOAD_UPWARD_DYLIB", "LC_RPATH", "LC_LOAD_DYLINKER",
    "LC_DYLD_ENVIRONMENT",
}
OPENSSL_DIRECTORY = re.compile(
    rb'(?:OPENSSLDIR|ENGINESDIR|MODULESDIR): "([^"\x00]*)"'
    rb'|((?:/opt/homebrew|/usr/local/Cellar)/[^"\x00]*openssl[^"\x00]*)'
)
LSREGISTER = (
    "/System/Library/Frameworks/CoreServices.framework/Frameworks/"
    "LaunchServices.framework/Support/lsregister"
)


def run(*args):
    subprocess.run([str(a) for a in args], cwd=ROOT, check=True)


def metadata(name):
    info = dict(
        CFBundleName=name,
        CFBundleDisplayName=name,
        CFBundleIdentifier=APPS[name],
        CFBundleExecutable=name,
        CFBundlePackageType="APPL",
        CFBundleVersion=__version__,
        CFBundleShortVersionString=__version__,
        LSMinimumSystemVersion=MINIMUM_MACOS,
        NSHumanReadableCopyright=f"{COPYRIGHT}. {LICENSE_NOTICE}.",
        NSHighResolutionCapable=True,
        NSPrincipalClass="NSApplication",
        CFBundleIconFile=name,
    )
    if name == "Kstack":
        info["CFBundleDocumentTypes"] = [
            {
                "CFBundleTypeName": "Kstack HTML",
                "CFBundleTypeRole": "Editor",
                "LSHandlerRank": "Alternate",
                "LSItemContentTypes": ["public.html"],
            }
        ]
    return info


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def uv_python_root():
    """Where uv keeps managed interpreters: `uv python dir`, else its documented default."""
    uv = os.environ.get("UV") or shutil.which("uv")
    if uv:
        try:
            found = subprocess.run(
                [uv, "python", "dir"], capture_output=True, text=True, check=True, timeout=60
            ).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            found = ""
        if found:
            return Path(found).expanduser().resolve()
    configured = os.environ.get("UV_PYTHON_INSTALL_DIR")
    return Path(configured or Path.home() / ".local/share/uv/python").expanduser().resolve()


def build_python():
    """Identity of the running interpreter, which must be the uv-managed CPython selected by
    `.python-version` (any patch release of that minor version)."""
    required = (ROOT / ".python-version").read_text().strip()
    version = platform.python_version()
    base = Path(sys.base_prefix).resolve()
    if required not in {version, ".".join(version.split(".")[:2])}:
        raise RuntimeError(
            f"Build with Python {required} from .python-version, not {version}: "
            "run ./scripts/build-install.sh or `uv run --locked python scripts/build_apps.py`"
        )
    root = uv_python_root()
    if sys.implementation.name != "cpython" or not base.is_relative_to(root):
        raise RuntimeError(
            f"Build with a uv-managed CPython {required}; the running {sys.implementation.name} "
            f"{version} at {base} is not under uv's managed Python directory {root}. "
            "Run `uv python install` and `uv sync --locked` (python-preference = only-managed)."
        )
    marker = base / "BUILD"
    if not marker.is_file() or not marker.read_text().strip():
        raise RuntimeError("uv-managed CPython is missing its standalone BUILD identity")
    return {"implementation": "CPython", "version": version, "distribution": base.name,
            "standalone_build": marker.read_text().strip()}


def version_key(text):
    parts = [int(part) for part in text.split(".")]
    return tuple(parts + [0] * (3 - len(parts)))


def is_macho(path):
    with path.open("rb") as handle:
        head = handle.read(8)
    if head[:4] in MACHO_MAGICS:
        return True
    # Java class files share the fat magic; their version field is never a small slice count.
    byteorder = "little" if head[:4] in FAT_MAGICS[2:] else "big"
    return head[:4] in FAT_MAGICS and 0 < int.from_bytes(head[4:8], byteorder) < 20


def parse_otool(output, paths):
    """Per file from `otool -arch all -h -l -v`: slices, deployment minimums and load paths."""
    known = set(paths)
    records, current, command = {}, None, None
    for line in output.splitlines():
        if line.endswith(":") and not line[:1].isspace():
            head, sliced = line[:-1], re.fullmatch(r"(.+) \(architecture (\S+)\)", line[:-1])
            if sliced and sliced[1] in known:
                head = sliced[1]
            if head in known:
                current = records.setdefault(
                    head, {"fat": False, "architectures": [], "minimums": [], "load_paths": []}
                )
                current["fat"] |= head != line[:-1]
                command = None
                continue
        fields = line.split()
        if current is None or not fields:
            continue
        if fields[0] in {"MH_MAGIC", "MH_MAGIC_64", "MH_CIGAM", "MH_CIGAM_64"} and len(fields) > 2:
            architecture, subtype = fields[1].lower(), fields[2].lower()
            if architecture == "arm64" and subtype not in {"all", "v8"}:
                architecture += subtype  # arm64e, arm64e.x1: not the plain arm64 that ships
            current["architectures"].append(architecture)
        elif fields[0] == "cmd" and len(fields) > 1:
            command = fields[1]
        elif len(fields) > 1 and (command, fields[0]) in {
            ("LC_BUILD_VERSION", "minos"), ("LC_VERSION_MIN_MACOSX", "version")
        }:
            current["minimums"].append((command, fields[1]))
        elif command in PATH_COMMANDS and fields[0] in {"name", "path"} and len(fields) > 1:
            current["load_paths"].append(line.split(None, 1)[1].rsplit(" (offset ", 1)[0])
    return records


def inspect_macho(paths):
    """Inspect every slice of each Mach-O with otool (batched; the arm64 slice is not assumed)."""
    paths, records = [str(path) for path in paths], {}
    for start in range(0, len(paths), 200):
        batch = paths[start:start + 200]
        result = subprocess.run(
            ["otool", "-arch", "all", "-h", "-l", "-v", *batch],
            capture_output=True, text=True, errors="replace",
        )
        if result.returncode:
            raise RuntimeError(f"otool could not inspect Mach-O files: {result.stderr[:2000]}")
        records.update(parse_otool(result.stdout, batch))
    return records


def openssl_homebrew_directories(data):
    """OpenSSL default directories (OPENSSLDIR, ENGINESDIR, MODULESDIR, cert paths) that point
    into a Homebrew prefix; such a library would trust the build machine's certificate store."""
    if b"OPENSSLDIR" not in data and not any(p.encode() in data for p in HOMEBREW_PREFIXES):
        return []
    found = []
    for match in OPENSSL_DIRECTORY.finditer(data):
        value = (match[1] if match[1] is not None else match[2]).decode(errors="replace")
        if value.startswith(HOMEBREW_PREFIXES):
            found.append(value)
    return sorted(set(found))


def macho_files(directory):
    return [
        p for p in sorted(Path(directory).rglob("*"))
        if not p.is_symlink() and p.is_file() and is_macho(p)
    ]


def macho_problems(directory):
    """Every Mach-O under directory: thin arm64, minimums <= MINIMUM_MACOS, no Homebrew paths."""
    directory = Path(directory)
    files = macho_files(directory)
    records = inspect_macho(files) if files else {}
    floor, problems = version_key(MINIMUM_MACOS), []
    for path in files:
        relative = path.relative_to(directory).as_posix()
        record = records.get(str(path))
        if record is None:
            problems.append(f"{relative}: otool could not read this Mach-O")
            continue
        if record["fat"] or record["architectures"] != ["arm64"]:
            kind = "universal" if record["fat"] else "thin"
            problems.append(
                f"{relative}: {kind} {'+'.join(record['architectures']) or 'unknown'}, "
                "not thin arm64"
            )
        for command, minimum in record["minimums"]:
            try:
                too_new = version_key(minimum) > floor
            except ValueError:
                too_new = True
            if too_new:
                problems.append(
                    f"{relative}: {command} minimum {minimum} exceeds macOS {MINIMUM_MACOS}"
                )
        for load in record["load_paths"]:
            if load.startswith(HOMEBREW_PREFIXES):
                problems.append(f"{relative}: Homebrew load path {load}")
        for value in openssl_homebrew_directories(path.read_bytes()):
            problems.append(f"{relative}: Homebrew OpenSSL directory {value}")
    return problems


def absolute_paths(value):
    """Strings (keys included) in a JSON-like record that are absolute or home-relative paths."""
    if isinstance(value, dict):
        for key, item in value.items():
            yield from absolute_paths(key)
            yield from absolute_paths(item)
    elif isinstance(value, list):
        for item in value:
            yield from absolute_paths(item)
    elif isinstance(value, str) and value.startswith(("/", "~")):
        yield value


def validate(directory):
    manifest = directory / MANIFEST
    record = json.loads(manifest.read_text()) if manifest.is_file() else {}
    if absolute := sorted(set(absolute_paths(record))):
        raise RuntimeError(f"Installation manifest records absolute paths: {absolute[:5]}")
    source = record.get("source")
    for name, identity in APPS.items():
        app = directory / f"{name}.app"
        info = plistlib.loads((app / "Contents/Info.plist").read_bytes())
        if (
            info["CFBundleIdentifier"] != identity
            or info["CFBundleExecutable"] != name
            or info["CFBundleShortVersionString"] != __version__
            or info["CFBundleVersion"] != __version__
            or info.get("LSMinimumSystemVersion") != MINIMUM_MACOS
        ):
            raise RuntimeError(f"Invalid identity/version/minimum macOS: {name}")
        if source and (
            info.get("KBioSuiteSourceCommit") != source["commit"]
            or info.get("KBioSuiteSourceDirty") != source["dirty"]
        ):
            raise RuntimeError(f"Build source does not match installation manifest: {name}")
        for p in app.rglob("*"):
            if p.is_symlink() and (not p.exists() or not p.resolve().is_relative_to(app.resolve())):
                raise RuntimeError(f"Broken or external bundle resource: {p}")
            if p.name == "__pycache__" or p.suffix == ".pyc":
                raise RuntimeError(f"Unwanted bytecode resource: {p}")
        if name == "Kmatch":
            package = app / "Contents/Frameworks/kmatch"
            for relative in KMATCH_RESOURCES:
                path = package / relative
                if not path.is_file() or not path.resolve().is_relative_to(app.resolve()):
                    raise RuntimeError(f"Missing Kmatch bundle resource: {relative}")
        if not (app / "Contents/Resources/Notices" / SUITE_LICENSE).is_file():
            raise RuntimeError(f"Missing suite license notice: {name}")
        validate_native_inventory(app)
        for tool, details in record.get("bundled_tools", {}).get(name, {}).items():
            path = app / details["path"]
            if not path.is_file() or sha256(path) != details["sha256"]:
                raise RuntimeError(f"Bundled {tool} does not match installation manifest: {name}")
            recorded = [("license", entry) for entry in details.get("licenses", [])]
            if "provenance" in details:
                recorded.append(("provenance", details["provenance"]))
            for kind, entry in recorded:
                path = app / entry["path"]
                if not path.is_file() or sha256(path) != entry["sha256"]:
                    raise RuntimeError(
                        f"Bundled {tool} {kind} does not match installation manifest: {name}"
                    )
        if problems := macho_problems(app):
            raise RuntimeError(
                f"{name} Mach-O audit failed ({len(problems)}):\n" + "\n".join(problems[:20])
            )
        run("codesign", "--verify", "--deep", "--strict", app)


def source_identity():
    """Record the exact local source without contacting a remote."""
    names = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "src", "apps",
         "packaging", "resources", "scripts/build_apps.py", "scripts/check_source_state.py",
         "scripts/notice_inventory.py",
         "pyproject.toml", "uv.lock",
         ".python-version", *NOTICE_SOURCES], cwd=ROOT, text=True,
    ).splitlines()
    fingerprints = {
        name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        for name in sorted(set(names)) if (ROOT / name).is_file()
    }
    return {
        "files_sha256": fingerprints,
        "commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "dirty": bool(subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=normal"], cwd=ROOT, text=True
        ).strip()),
    }


def bundled_tools(products):
    """Copied tools are outside source_identity(); identify them from the products themselves.

    Hashes are recomputed from the current bytes, so a manifest regenerated after re-signing
    still validates; release, archive and thinned-binary hashes come from the setup record.
    """
    app = products / "Kstack.app"
    runtime = app / MMSEQS_RUNTIME
    provenance = runtime / MMSEQS_PROVENANCE
    setup = json.loads(provenance.read_text())
    binary = runtime / "bin/mmseqs"

    def identity(path):
        return {"path": path.relative_to(app).as_posix(), "sha256": sha256(path)}

    return {"Kstack": {"mmseqs": {
        **identity(binary),
        "version": setup["version"],
        "architecture": setup["architecture"],
        "archive_sha256": setup["archive_sha256"],
        "thinned_sha256": setup["binary_sha256"],
        "provenance": identity(provenance),
        "licenses": [identity(p) for p in sorted(runtime.glob("LICENSE*"))],
    }}}


def write_install_manifest(products, *, built_utc=None, source=None, python=None):
    """Write products/.kbiosuite-install.json; release tooling calls it again after re-signing.

    Omitted build facts are carried over from the existing manifest; bundled-tool hashes are
    always recomputed. Only relative paths are recorded.
    """
    path = products / MANIFEST
    previous = json.loads(path.read_text()) if path.is_file() else {}
    try:
        record = {
            "version": __version__,
            "built_utc": built_utc or previous["built_utc"],
            "apps": APPS,
            "source": source or previous["source"],
            "python": python or previous["python"],
            "bundled_tools": bundled_tools(products),
        }
    except KeyError as missing:
        raise RuntimeError(f"Cannot regenerate {path.name} without its {missing} record") from None
    if absolute := sorted(set(absolute_paths(record))):
        raise RuntimeError(f"Refusing to record absolute paths in {path.name}: {absolute[:5]}")
    path.write_text(json.dumps(record, indent=2) + "\n")
    return record


def retain_latest_build(products):
    """Retain one completed build; never adopt or prune unmarked local material.

    Older marked builds are pruned when their recorded apps are a subset of the suite, so a
    build from before the app set grew is still retired; failed builds lack the manifest.
    """
    stage = products.parent
    root = stage.parent
    marker = {"owner": "kbiosuite-build-v1", "root": str(root.resolve())}
    (stage / ".kbiosuite-build.json").write_text(json.dumps(marker) + "\n")
    temporary = root / ".latest-new"
    temporary.unlink(missing_ok=True)
    temporary.symlink_to(products.relative_to(root))
    temporary.replace(root / "latest")
    for previous in root.iterdir():
        if previous == stage or previous.is_symlink() or not previous.is_dir():
            continue
        try:
            owned = json.loads((previous / ".kbiosuite-build.json").read_text()) == marker
            installed = json.loads((previous / "apps/.kbiosuite-install.json").read_text())
        except (OSError, ValueError):
            continue
        recorded = installed.get("apps")
        if owned and isinstance(recorded, dict) and recorded and recorded.items() <= APPS.items():
            shutil.rmtree(previous)


def build(*, exact_runtime=False):
    """Build all five apps. Release tooling passes exact_runtime=True, which refuses (before
    freezing) any runtime other than the one licenses/inventory.json reviewed; local builds
    accept any uv-managed patch and mark the result as not covered by that review."""
    if sys.platform != "darwin":
        raise RuntimeError("macOS is required")
    python = build_python()
    inventory = json.loads((ROOT / "licenses/inventory.json").read_text())
    coverage = runtime_coverage(inventory, exact_runtime=exact_runtime)
    python = {**python, "redistribution_inventory": coverage["status"]}
    if coverage["status"] != REVIEWED_RUNTIME:
        print(
            f"Local build with {describe_runtime(coverage['identity'])}; the reviewed release "
            f"runtime is {describe_runtime(coverage['reviewed_identity'])}. Its Notices carry "
            "this runtime's own license texts and mark it as not covered by the reviewed "
            "redistribution inventory: for your own use, not redistribution.",
            flush=True,
        )
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    build_source = source_identity()
    stage = BUILD_ROOT / stamp
    products = stage / "apps"
    products.mkdir(parents=True)
    for name, package in [
        ("Kblast", "kblast"),
        ("Kannotate", "kannotate"),
        ("Kmatch", "kmatch"),
        ("KBioRuntime", "kstack"),
    ]:
        print(f"Building {name}…", flush=True)
        args = [
            sys.executable,
            "-m",
            "PyInstaller",
            "--noconfirm",
            "--clean",
            "--onedir",
            "--name",
            name,
            "--distpath",
            str(stage / "frozen"),
            "--workpath",
            str(stage / "work" / name),
            "--specpath",
            str(stage / "specs"),
            "--paths",
            str(ROOT / "src"),
            "--additional-hooks-dir",
            str(ROOT / "packaging/hooks"),
            "--collect-data",
            package,
            "--collect-data",
            "kbiosuite",
            "--exclude-module",
            "tkinter",
            "--exclude-module",
            "pytest",
        ]
        if name != "KBioRuntime":
            args += ["--windowed", "--hidden-import", "webview.platforms.cocoa"]
        if name == "Kannotate" and importlib.util.find_spec("pyhmmer") is not None:
            args += [
                "--copy-metadata",
                "pyhmmer",
                "--exclude-module",
                "pyhmmer.tests",
                "--hidden-import",
                "pyhmmer.platform.bsd",
            ]
        if name == "Kmatch":
            # collect-data skips .py; external engine/PyMOL interpreters run these files.
            args += [
                "--add-data",
                f"{ROOT / 'src/kmatch/engine_runner.py'}:kmatch",
                "--add-data",
                f"{ROOT / 'src/kmatch/pymol_runner.py'}:kmatch",
                "--exclude-module",
                "tmtools",
            ]
        args.append(str(ROOT / "packaging" / f"{name}.py"))
        run(*args)
    for name in ["Kblast", "Kannotate", "Kmatch"]:
        run("/usr/bin/ditto", stage / "frozen" / f"{name}.app", products / f"{name}.app")
    for name in ["KBioSuite", "Kstack"]:
        (products / f"{name}.app/Contents/MacOS").mkdir(parents=True)
        (products / f"{name}.app/Contents/Resources").mkdir()
    run(
        "swiftc",
        "-parse-as-library",
        "-O",
        "-target",
        f"arm64-apple-macos{MINIMUM_MACOS}",
        *sorted((ROOT / "apps/kbiosuite-macos").glob("*.swift")),
        "-o",
        products / "KBioSuite.app/Contents/MacOS/KBioSuite",
    )
    # Package.swift declares the same deployment target (checked by check_versions.py).
    run("swift", "build", "--package-path", ROOT / "apps/kstack-macos", "-c", "release")
    shutil.copy2(
        ROOT / "apps/kstack-macos/.build/release/Kstack",
        products / "Kstack.app/Contents/MacOS/Kstack",
    )
    resources = products / "Kstack.app/Contents/Resources"
    shutil.copytree(
        ROOT / "src/kstack/resources",
        resources / "kstack",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    shutil.copytree(stage / "frozen/KBioRuntime", resources / "runtime")
    mmseqs = ROOT / ".tools/mmseqs"
    setup = mmseqs / MMSEQS_PROVENANCE
    if not (mmseqs / "bin/mmseqs").is_file() or not setup.is_file():
        raise RuntimeError(
            "Thin arm64 MMseqs2 and its setup record are missing; explicitly run "
            "scripts/kstack/setup-mmseqs.sh"
        )
    if sha256(mmseqs / "bin/mmseqs") != json.loads(setup.read_text())["binary_sha256"]:
        raise RuntimeError(
            "MMseqs2 differs from its setup record; rerun scripts/kstack/setup-mmseqs.sh"
        )
    runtime_mmseqs = products / "Kstack.app" / MMSEQS_RUNTIME
    (runtime_mmseqs / "bin").mkdir(parents=True)
    shutil.copy2(mmseqs / "bin/mmseqs", runtime_mmseqs / "bin/mmseqs")
    shutil.copy2(setup, runtime_mmseqs / MMSEQS_PROVENANCE)
    for p in mmseqs.glob("LICENSE*"):
        shutil.copy2(p, runtime_mmseqs / p.name)
    native_sources = native_source_index(
        ROOT, inventory,
        suite_binaries=[products / f"{name}.app/Contents/MacOS/{name}"
                        for name in ("KBioSuite", "Kstack")],
        exact_runtime=exact_runtime,
    )
    run("swift", ROOT / "scripts/create-icons.swift", stage / "icons")
    for name in APPS:
        app = products / f"{name}.app"
        res = app / "Contents/Resources"
        res.mkdir(exist_ok=True)
        notice = res / "Notices"
        shutil.copytree(ROOT / "licenses", notice)
        if (notice / SUITE_LICENSE).exists():
            raise RuntimeError(f"licenses/{SUITE_LICENSE} would shadow the suite license")
        shutil.copy2(ROOT / SUITE_LICENSE, notice / SUITE_LICENSE)
        for relative in SHARED_NOTICES:
            shutil.copy2(ROOT / relative, notice / Path(relative).name)
        if name == "Kmatch":
            # AlphaFold 3 terms are referenced from the GUI.
            for relative in KMATCH_NOTICES:
                shutil.copy2(ROOT / relative, notice / Path(relative).name)
        collect_notices(ROOT, notice, exact_runtime=exact_runtime)
        plist = app / "Contents/Info.plist"
        info = plistlib.loads(plist.read_bytes()) if plist.exists() else {}
        info.update(metadata(name))
        info["KBioSuiteSourceCommit"] = build_source["commit"]
        info["KBioSuiteSourceDirty"] = build_source["dirty"]
        plist.write_bytes(plistlib.dumps(info))
        icon = (
            ROOT / "resources/kstack/Kstack.icns"
            if name == "Kstack"
            else stage / "icons" / f"{name}.icns"
        )
        shutil.copy2(icon, res / f"{name}.icns")
        audit_native_inventory(app, native_sources)
        # Sign each embedded Mach-O before the outer bundle; do not mutate resources afterward.
        run("codesign", "--force", "--deep", "--sign", "-", "--identifier", APPS[name], app)
    if source_identity() != build_source:
        raise RuntimeError("Source changed during build; rebuild before installing")
    write_install_manifest(products, built_utc=stamp, source=build_source, python=python)
    validate(products)
    # Verified signed products suffice for reinstalling; intermediate apps are disposable.
    for name in ("frozen", "work", "specs", "icons"):
        shutil.rmtree(stage / name)
    retain_latest_build(products)
    return products


def install(products, destination):
    validate(products)
    manifest = json.loads((products / MANIFEST).read_text())
    if manifest["source"] != source_identity():
        raise RuntimeError("Build products do not match current source; rebuild before installing")
    destination = destination.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not (destination / ".kbiosuite-install.json").is_file():
        raise RuntimeError(f"Destination is not a suite-owned installation: {destination}")
    # Refuse to replace running applications, including apps launched by Finder.
    listing = subprocess.check_output(["/bin/ps", "-axo", "command="], text=True)
    if any(
        str(destination / f"{name}.app/Contents/") in line
        for name in APPS
        for line in listing.splitlines()
    ):
        raise RuntimeError("Quit the installed suite apps before updating")
    stage = destination.parent / f".KBioSuite-install-{os.getpid()}.noindex"
    if stage.exists():
        raise RuntimeError(f"Inspect existing staging directory: {stage}")
    run("/usr/bin/ditto", products, stage)
    validate(stage)
    rollback = destination.parent / ".KBioSuite-previous.noindex"
    if rollback.exists():
        if not (rollback / ".kbiosuite-install.json").is_file():
            raise RuntimeError("Unowned rollback destination")
        shutil.rmtree(rollback)
    if destination.exists():
        destination.rename(rollback)
    try:
        stage.rename(destination)
        validate(destination)
        for name in APPS:
            run(LSREGISTER, "-f", destination / f"{name}.app")
    except Exception:
        if destination.exists():
            destination.rename(stage)
        if rollback.exists():
            rollback.rename(destination)
        raise
    # The previous installation is needed only until promotion/verification succeeds.
    if rollback.exists():
        shutil.rmtree(rollback)
    print(f"Installed all five apps at {destination}")


def main():
    from check_source_state import report

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--install", action="store_true")
    parser.add_argument("--install-only", type=Path)
    parser.add_argument("--destination", type=Path, default=Path.home() / "Applications/KBioSuite")
    parser.add_argument(
        "--reconcile-installed", help="Approved reason for replacing unproven installed source"
    )
    args = parser.parse_args()
    report(
        ROOT, args.destination, install=args.install, reconcile_installed=args.reconcile_installed
    )
    products = args.install_only or build()
    if args.install:
        report(ROOT, args.destination, install=True, reconcile_installed=args.reconcile_installed)
        install(products, args.destination)
    print(products)


if __name__ == "__main__":
    main()
