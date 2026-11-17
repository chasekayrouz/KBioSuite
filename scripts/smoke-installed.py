"""Exercise installed apps via LaunchServices, real GUI jobs and close/reopen."""

import argparse
import contextlib
import csv
import json
import os
import shutil
import subprocess
import time
from datetime import UTC, datetime
from pathlib import Path

import AppKit
import Foundation
import gemmi
from Bio import SeqIO
from build_apps import APPS, source_identity, validate
from check_versions import check

from kbiosuite.version import __version__

ROOT = Path(__file__).resolve().parents[1]
# Installed apps run without KMATCH_HOME, so they must report the suite default home.
KMATCH_HOME = Path.home() / "Library/Application Support/KBioSuite/Kmatch"
KMATCH_METHODS = {
    "choose_msa", "export_structure", "reveal_run", "export_alignment", "open_pdb"
}
KMATCH_DEPENDENCIES = {"boltz2", "foldseek", "pdb100"}
STRUCTURE_SUFFIXES = {".cif", ".mmcif", ".pdb", ".ent"}
# The installed Kmatch ends its own --smoke-dir run (and writes report.json) within its
# SMOKE_DEADLINE_SECONDS (2400 s, src/kmatch/gui.py); the harness waits longer so that the
# app's own deadline report, not a harness timeout, explains a slow job.
KMATCH_REPORT_TIMEOUT_SECONDS = 2700
NORMAL_QUIT_SECONDS = 20
FORCED_QUIT_SECONDS = 10
CHILD_EXIT_SECONDS = 10


def running_applications():
    # NSWorkspace observations update through the Cocoa run loop, even in a CLI harness.
    Foundation.NSRunLoop.currentRunLoop().runUntilDate_(
        Foundation.NSDate.dateWithTimeIntervalSinceNow_(0.05)
    )
    return list(AppKit.NSWorkspace.sharedWorkspace().runningApplications())


def bundle_path(application):
    url = application.bundleURL()
    path = url.path() if url is not None else None
    return Path(str(path)).resolve() if path else None


def installed_here(application, directory):
    """Only bundles inside the smoked installation count; same-ID apps elsewhere never do."""
    path = bundle_path(application)
    return path is not None and path.is_relative_to(directory)


def select(applications, identity, directory):
    return [
        a
        for a in applications
        if a.bundleIdentifier() == identity and installed_here(a, directory)
    ]


def foreign(applications, directory):
    """Running apps that share a suite bundle ID but live outside the smoked installation."""
    identities = set(APPS.values())
    return [
        a
        for a in applications
        if a.bundleIdentifier() in identities and not installed_here(a, directory)
    ]


def describe(application):
    return (
        f"{application.bundleIdentifier()} pid {application.processIdentifier()} at "
        f"{bundle_path(application) or 'an unknown location'}"
    )


def refuse_foreign(directory):
    """Never control a standalone app (e.g. ~/Applications/Kmatch.app) that shares an ID."""
    others = foreign(running_applications(), directory)
    if others:
        raise RuntimeError(
            "A non-suite app with a suite bundle identifier is running: "
            + "; ".join(map(describe, others))
            + f". The acceptance run only controls apps inside {directory} and never quits or "
            "cancels others. Let its jobs finish, quit it yourself, then rerun."
        )


def apps(identity, directory):
    return select(running_applications(), identity, directory)


def wait(check, timeout=150, what="the installed app"):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        value = check()
        if value:
            return value
        time.sleep(0.2)
    raise RuntimeError(f"Installed-app acceptance timed out after {timeout:g} s waiting for {what}")


def close(identity, directory):
    for app in apps(identity, directory):
        if not app.terminate():
            raise RuntimeError(f"Application refused normal quit: {identity}")
    try:
        wait(lambda: not apps(identity, directory), NORMAL_QUIT_SECONDS)
    except RuntimeError as exc:
        raise RuntimeError(f"Normal quit timed out: {identity}") from exc


def ready(identity, directory):
    current = apps(identity, directory)
    return current if current and all(app.isFinishedLaunching() for app in current) else None


def descendants(pid):
    import psutil

    try:
        return psutil.Process(pid).children(recursive=True)
    except psutil.Error:
        return []


def stop_children(children):
    """Terminate, then kill, child processes (e.g. a Boltz engine) that outlived their app."""
    import psutil

    _, alive = psutil.wait_procs(children, timeout=CHILD_EXIT_SECONDS)
    for child in alive:
        with contextlib.suppress(psutil.Error):
            child.terminate()
    _, stubborn = psutil.wait_procs(alive, timeout=5)
    for child in stubborn:
        with contextlib.suppress(psutil.Error):
            child.kill()
    return [child.pid for child in alive]


def stop_launched(directory, names):
    """Quit every app the smoke launched inside directory, then any engine left behind.

    Runs on every exit path so a failed acceptance never leaves a Kmatch/Boltz job running
    or blocks the next run. Apps outside directory are never touched. Never raises.
    """
    outcome = {"closed": [], "forced": [], "children_stopped": [], "errors": []}
    for name in dict.fromkeys(names):
        identity = APPS[name]
        try:
            current = apps(identity, directory)
            if not current:
                continue
            children = [child for app in current for child in descendants(app.processIdentifier())]
            for app in current:
                app.terminate()
            try:
                wait(lambda identity=identity: not apps(identity, directory), NORMAL_QUIT_SECONDS)
                outcome["closed"].append(name)
            except RuntimeError:
                for app in apps(identity, directory):
                    app.forceTerminate()
                try:
                    wait(
                        lambda identity=identity: not apps(identity, directory),
                        FORCED_QUIT_SECONDS,
                    )
                    outcome["forced"].append(name)
                except RuntimeError:
                    outcome["errors"].append(f"{name} is still running after a forced quit")
            outcome["children_stopped"] += [
                f"{name} child pid {pid}" for pid in stop_children(children)
            ]
        except Exception as exc:  # cleanup must never mask the original failure
            outcome["errors"].append(f"{name}: {exc}")
    return outcome


def launch(directory, name, *args):
    refuse_foreign(directory)
    # No activated shell, checkout cwd or developer environment variables.
    environment = {
        "HOME": str(Path.home()),
        "PATH": "/usr/bin:/bin",
        "TMPDIR": os.environ.get("TMPDIR", "/tmp"),
        "LANG": "en_US.UTF-8",
    }
    subprocess.run(
        ["/usr/bin/open", "-n", str(directory / f"{name}.app"), "--args", *map(str, args)],
        cwd="/tmp",
        env=environment,
        check=True,
    )
    current = wait(lambda: ready(APPS[name], directory), 30, f"{name} to finish launching")
    if not all(bundle_path(a) == (directory / f"{name}.app").resolve() for a in current):
        raise RuntimeError(f"Unexpected installed identity/path: {name}")


def read_report(path):
    if not path.is_file():
        return None
    return json.loads(path.read_text())


def running_report(path, identity, directory):
    """A long job's report, failing fast if the app quit without writing one."""
    result = read_report(path)
    if result is None and not apps(identity, directory):
        raise RuntimeError(f"{identity} quit before writing {path.name}")
    return result


def check_kmatch_bridge(bridge):
    """Installed Kmatch initialized pywebview, its native methods and local dependencies."""
    assert bridge["ok"] and not bridge.get("error"), bridge
    assert bridge["version"] == __version__, bridge
    assert Path(bridge["home"]).resolve() == KMATCH_HOME.resolve(), bridge
    assert bridge["bridge"]["event_seen"] and bridge["bridge"]["ready"], bridge
    assert bridge["bridge"]["desktop"] and not bridge["bridge"].get("error"), bridge
    methods = bridge["bridge"]["methods"]
    assert set(methods) >= KMATCH_METHODS and all(methods.values()), bridge
    dependencies = bridge["dependencies_ready"]
    assert set(dependencies) >= KMATCH_DEPENDENCIES and all(dependencies.values()), bridge


def check_kmatch_job(result, evidence):
    """One real MSA-free installed job: isolated run directory and parseable exports."""
    assert result["ok"] and not result.get("error"), result
    assert result["version"] == __version__, result
    assert result["checks"] and all(result["checks"].values()), result
    assert result["run_id"] and result["prediction"], result
    assert 1 <= result["hits"] <= 50 and result["top_hit"], result
    run_directory = Path(result["run_directory"]).resolve()
    assert run_directory.is_dir() and run_directory.is_relative_to(evidence / "runs"), result
    assert result["exports"], result
    structures = 0
    for name, path in result["exports"].items():
        path = Path(path).resolve()
        assert path.is_relative_to(evidence) and path.is_file() and path.stat().st_size, name
        if path.suffix.lower() in STRUCTURE_SUFFIXES:
            structure = gemmi.read_structure(str(path))
            assert len(structure) and structure[0].count_atom_sites(), name
            structures += 1
    assert structures, result["exports"]


def smoke(
    directory,
    output,
    *,
    bridge_launches=10,
    launch_cycles=3,
    check_eutils=False,
    kblast_rid=None,
    allow_uncommitted=False,
):
    directory = directory.resolve()
    output.mkdir(parents=True)
    validate(directory)
    check(directory)
    source = json.loads((directory / ".kbiosuite-install.json").read_text())["source"]
    assert source == source_identity(), "Installed build does not match current source"
    assert not source["dirty"] or (allow_uncommitted and source.get("files_sha256")), (
        "Build clean source, or use --allow-uncommitted with a fully fingerprinted build"
    )
    refuse_foreign(directory)
    if any(apps(identity, directory) for identity in APPS.values()):
        raise RuntimeError(f"Quit the suite apps in {directory} before the acceptance run")
    report = {"ok": False, "installed": str(directory), "source": source, "checks": []}
    launched = []  # every app this run started; all are quit on any exit path

    def start(name, *args):
        launched.append(name)
        launch(directory, name, *args)

    try:
        start(
            "Kblast",
            "--data-dir",
            output / "kblast/runs",
            "--smoke-dir",
            output / "kblast",
            *(["--check-eutils"] if check_eutils else []),
            *(["--qualify-rid", kblast_rid] if kblast_rid else []),
        )
        result = wait(
            lambda: read_report(output / "kblast/report.json"),
            1250 if kblast_rid else 150,
            "the Kblast report",
        )
        assert result["ok"], result
        report["kblast"] = result
        report["checks"].append("Kblast native GUI loaded from installed frozen resources")
        close(APPS["Kblast"], directory)
        report["kannotate_bridge_launches"] = []
        for index in range(bridge_launches):
            evidence = output / "native-bridge" / f"{index + 1:02d}"
            start("Kannotate", "--bridge-check-dir", evidence)
            bridge = wait(
                lambda evidence=evidence: read_report(evidence / "bridge.json"),
                30,
                "the Kannotate bridge report",
            )
            assert bridge["ok"] and bridge["bridge"]["event_seen"], bridge
            assert all(bridge["bridge"]["methods"].values()), bridge
            assert all(bridge["databases_ready"].values()), bridge
            report["kannotate_bridge_launches"].append(bridge)
            close(APPS["Kannotate"], directory)
        report["checks"].append(
            f"{bridge_launches} installed Kannotate native bridge initializations"
        )
        quick_input = output / "quick-search-five.faa"
        records = list(SeqIO.parse(ROOT / "tests/kannotate/fixtures/quick-search.faa", "fasta"))
        SeqIO.write(records[:5], quick_input, "fasta")
        start(
            "Kannotate",
            "--smoke-dir",
            output / "kannotate",
            "--quick-smoke-input",
            quick_input,
        )
        result = wait(
            lambda: read_report(output / "kannotate/report.json"), what="the Kannotate report"
        )
        assert result["ok"], result
        assert result["quick_search"]["ok"], result
        assert result["stop"]["ok"], result
        assert all(result["format_checks"].values()), result
        assert set(result["genbank_modes"]) == {"genbank", "csv"}, result
        assert len(result["state"]["results"]) == 2
        with (output / "kannotate/annotation outputs/public EF Tu.kannotate.csv").open(
            encoding="utf-8-sig", newline=""
        ) as handle:
            rows = list(csv.reader(handle))
        assert len(rows) == 2 and len(rows[0]) == 28
        assert "elongation factor" in rows[1][2].lower()
        genbank = SeqIO.read(
            output / "kannotate/annotation outputs/public control.kannotate.gbk", "genbank"
        )
        assert "elongation factor" in genbank.features[0].qualifiers["product"][0].lower()
        assert genbank.features[0].qualifiers["locus_tag"] == ["control"]
        assert genbank.features[0].location.start == 0
        with (output / "kannotate/annotation outputs/public control.kannotate.csv").open(
            encoding="utf-8-sig", newline=""
        ) as handle:
            rows = list(csv.reader(handle))
        assert len(rows) == 2 and all(len(row) == 31 for row in rows)
        assert rows[0][3:7] == ["From", "To", "Strand", "Length"]
        assert not {
            "Location note", "Kannotate version", "Scoring version", "Pfam version", "CDD version"
        }.intersection(rows[0])
        assert rows[1][3:6] == ["1", str(len(genbank)), "+"]
        translation = genbank.features[0].qualifiers["translation"][0]
        assert rows[1][6] == str(len(translation.replace("*", "")))
        assert rows[1][7] == "old test product"
        assert rows[1][9] == genbank.features[0].qualifiers["product"][0]
        with (output / "kannotate/csv only/public control.kannotate.csv").open(
            encoding="utf-8-sig", newline=""
        ) as handle:
            assert list(csv.reader(handle)) == rows
        assert [len(r["outputs"]) for r in result["state"]["results"]] == [1, 2]
        assert not any(r["warnings"] for r in result["state"]["results"])
        report["kannotate"] = result
        report["checks"].append(
            "Kannotate real GUI: all three GenBank output modes, "
            "31-column CSV with protein length, paired outputs, first-input folder and "
            "session override; both search engines"
        )
        report["checks"].append(
            "Installed Quick Search: GenBank file converted to five queries, domains, "
            "CDD filtering/alignment, "
            "bidirectional selection, navigation and zero-hit query"
        )
        report["checks"].append(
            "Installed Stop: decline/confirmation, both engine processes exit, "
            "temporary data removed, inputs preserved and a successful rerun"
        )
        close(APPS["Kannotate"], directory)
        fixtures = output / "kstack-fixture-root"
        shutil.copytree(ROOT / "tests/kstack/fixtures", fixtures / "tests/kstack/fixtures")
        (fixtures / "tests/kstack/renderer").mkdir(parents=True)
        shutil.copy2(
            ROOT / "tests/kstack/renderer/visual-harness.js", fixtures / "tests/kstack/renderer"
        )
        start("Kstack", "--smoke", fixtures)
        result = wait(
            lambda: read_report(fixtures / ".artifacts/native-smoke/report.json"),
            what="the Kstack report",
        )
        assert result["ok"], result
        report["kstack"] = result
        report["checks"].append(
            "Kstack native/WebKit acceptance passed with fixture-only checkout substitute"
        )
        close(APPS["Kstack"], directory)
        report["kmatch_bridge_launches"] = []
        for index in range(bridge_launches):
            evidence = output / "kmatch-bridge" / f"{index + 1:02d}"
            start("Kmatch", "--bridge-check-dir", evidence)
            bridge = wait(
                lambda evidence=evidence: read_report(evidence / "bridge.json"),
                60,
                "the Kmatch bridge report",
            )
            check_kmatch_bridge(bridge)
            report["kmatch_bridge_launches"].append(bridge)
            close(APPS["Kmatch"], directory)
        report["checks"].append(
            f"{bridge_launches} installed Kmatch native bridge initializations with the suite "
            "default home and ready Boltz-2, Foldseek and PDB100"
        )
        evidence = output / "kmatch"
        kmatch_input = output / "kmatch-1csp-A.fasta"
        shutil.copy2(ROOT / "tests/kmatch/fixtures/1csp-A.fasta", kmatch_input)
        start("Kmatch", "--smoke-dir", evidence, "--smoke-input", kmatch_input)
        check_kmatch_bridge(
            wait(lambda: read_report(evidence / "bridge.json"), 60, "the Kmatch bridge report")
        )
        # Real local MSA-free prediction plus PDB100 search; the app enforces its own deadline.
        result = wait(
            lambda: running_report(evidence / "report.json", APPS["Kmatch"], directory),
            KMATCH_REPORT_TIMEOUT_SECONDS,
            "the Kmatch job report",
        )
        check_kmatch_job(result, evidence)
        report["kmatch"] = result
        report["checks"].append(
            f"Kmatch real GUI job: MSA-free prediction of public 1CSP chain A, {result['hits']} "
            "structural hits, isolated run directory and gemmi-parsed coordinate exports"
        )
        close(APPS["Kmatch"], directory)
        report["launcher"] = []
        for cycle in range(launch_cycles):
            evidence = output / f"launcher-{cycle + 1}.json"
            start("KBioSuite", "--smoke", evidence)
            children = ["Kblast", "Kannotate", "Kstack", "Kmatch"]
            launched.extend(children)  # opened by the launcher, quit by cleanup on failure
            result = wait(lambda evidence=evidence: read_report(evidence), what="the launcher")
            assert result["ok"], result
            for name in children:
                wait(lambda name=name: ready(APPS[name], directory), 30, f"launched {name}")
            report["launcher"].append(result)
            for name in APPS:
                close(APPS[name], directory)
        report["checks"].append(
            "Native suite launcher opened all four independently installed child apps"
        )
        for _ in range(launch_cycles):
            for name in APPS:
                start(name)
                time.sleep(1.5)
                assert apps(APPS[name], directory), name
                close(APPS[name], directory)
        report["checks"].append("All five apps closed normally and reopened independently")
        validate(directory)
        report["checks"].append(
            "All strict signatures remain valid; no bytecode/resource mutations"
        )
        report["ok"] = True
    finally:
        # Also after an assertion, timeout or exception: no app or engine is left running.
        report["cleanup"] = stop_launched(directory, launched)
        if report["cleanup"]["errors"]:
            report["ok"] = False
        (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    if report["cleanup"]["errors"]:
        raise RuntimeError("Acceptance cleanup failed: " + "; ".join(report["cleanup"]["errors"]))
    print(
        json.dumps(
            {"ok": report["ok"], "checks": report["checks"], "output": str(output)}, indent=2
        )
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apps", type=Path, default=Path.home() / "Applications/KBioSuite")
    parser.add_argument("--bridge-launches", type=int, default=10)
    parser.add_argument("--launch-cycles", type=int, default=3)
    parser.add_argument("--check-eutils", action="store_true", help="Opt-in harmless live EInfo")
    parser.add_argument(
        "--kblast-rid", help="Opt-in browser RID for the fixed public CAC44203.1 installed gate"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / ".local/validation/installed" / datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ"),
    )
    parser.add_argument(
        "--allow-uncommitted",
        action="store_true",
        help="Qualify a fingerprinted build without claiming Git integration",
    )
    args = parser.parse_args()
    if args.bridge_launches < 1 or args.launch_cycles < 1:
        parser.error("Launch counts must be positive")
    smoke(
        args.apps,
        args.output.resolve(),
        bridge_launches=args.bridge_launches,
        launch_cycles=args.launch_cycles,
        check_eutils=args.check_eutils,
        kblast_rid=args.kblast_rid,
        allow_uncommitted=args.allow_uncommitted,
    )


if __name__ == "__main__":
    main()
