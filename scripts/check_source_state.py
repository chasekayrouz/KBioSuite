"""Report local source/default/installed mismatches and guard release installation.

Uses only local Git and signed-build records.
"""

import argparse
import ast
import json
import os
import subprocess
from pathlib import Path

# Private checkouts commit a push guard; its installed copy lives in the Git common directory.
GUARD_HOOK = "scripts/hooks/pre-push"
GUARD_FILES = {"pre-push": GUARD_HOOK, "export_public.py": "scripts/export_public.py"}


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


def version(text):
    for node in ast.parse(text).body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "__version__" for target in node.targets
        ):
            return ast.literal_eval(node.value)
    raise ValueError("Missing application version")


def ancestor(root, older, newer):
    result = subprocess.run(
        ["git", "-C", str(root), "merge-base", "--is-ancestor", older, newer],
        capture_output=True,
    )
    return result.returncode == 0


def committed(root, revision, path):
    result = subprocess.run(
        ["git", "-C", str(root), "cat-file", "blob", f"{revision}:{path}"], capture_output=True
    )
    return None if result.returncode else result.stdout


def hooks_path(worktree):
    result = subprocess.run(
        ["git", "-C", str(worktree), "config", "--get", "core.hooksPath"],
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() or None


def push_guard(root, worktrees):
    """Installed guard versus this revision's committed copy (None without a guard)."""
    if committed(root, "HEAD", GUARD_HOOK) is None:
        return None
    target = Path(git(root, "rev-parse", "--path-format=absolute", "--git-common-dir"))
    target = target / "kbio-hooks"
    files = {}
    for name, source in GUARD_FILES.items():
        expected = committed(root, "HEAD", source)
        installed = target / name
        files[name] = (
            expected is not None and installed.is_file() and installed.read_bytes() == expected
        )
    configured = {tree: hooks_path(tree) for tree in worktrees}
    return {
        "directory": str(target),
        "files": files,
        "worktrees": configured,
        "current": all(files.values())
        and os.access(target / "pre-push", os.X_OK)
        and bool(configured)
        and all(value == str(target) for value in configured.values()),
        "private_terms_present": (target / "private-terms.txt").is_file(),
    }


def inspect(root, destination=None):
    path = "src/kbiosuite/version.py"
    destination = (destination or Path.home() / "Applications/KBioSuite").expanduser()
    record = destination / ".kbiosuite-install.json"
    head, main = git(root, "rev-parse", "HEAD"), git(root, "rev-parse", "refs/heads/main")
    state = {
        "branch": git(root, "branch", "--show-current") or "(detached)",
        "revision": head,
        "version": version((root / path).read_text()),
        "dirty": bool(git(root, "status", "--porcelain", "--untracked-files=all")),
        "main_revision": main,
        "main_version": version(git(root, "show", f"refs/heads/main:{path}")),
        "contains_main": ancestor(root, main, head),
        "destination": str(destination),
        "installed": None,
        "worktrees": [],
    }
    for line in git(root, "worktree", "list", "--porcelain").splitlines():
        if line.startswith("worktree ") and Path(line[9:]).is_dir():
            worktree = Path(line[9:])
            state["worktrees"].append(
                {
                    "path": str(worktree),
                    "branch": git(worktree, "branch", "--show-current") or "(detached)",
                    "dirty": bool(git(worktree, "status", "--porcelain", "--untracked-files=all")),
                }
            )
    guard = push_guard(root, [tree["path"] for tree in state["worktrees"]])
    if guard is not None:
        state["push_guard"] = guard
    if destination.exists():
        try:
            data = json.loads(record.read_text())
            source = data["source"]
            revision = source.get("commit")
            state["installed"] = {
                "version": data["version"],
                "revision": revision,
                "dirty": source["dirty"],
                "ancestor_of_source": bool(revision and ancestor(root, revision, head)),
            }
        except (OSError, ValueError, KeyError, TypeError):
            state["installed"] = {"error": "Missing or invalid installed source provenance"}
    return state


def issues(state):
    result = []
    if state["dirty"]:
        result.append("Current source has uncommitted changes")
    if not state["contains_main"]:
        result.append("Current source does not contain local main")
    if state["version"] != state["main_version"]:
        result.append(
            f"Source version {state['version']} differs from main {state['main_version']}"
        )
    for tree in state["worktrees"]:
        if tree["dirty"]:
            result.append(f"Uncommitted work in {tree['path']} ({tree['branch']})")
    if "push_guard" in state and not state["push_guard"]["current"]:
        result.append("Push guard is missing or stale; run scripts/hooks/install_hooks.py")
    installed = state["installed"]
    if installed:
        if "error" in installed:
            result.append(installed["error"])
        else:
            if installed["dirty"]:
                result.append("Installed app was built from uncommitted source")
            if installed["version"] != state["version"]:
                result.append(
                    f"Installed version {installed['version']} "
                    f"differs from source {state['version']}"
                )
            if installed["revision"] != state["revision"]:
                result.append("Installed revision differs from current source")
            if not installed["ancestor_of_source"]:
                result.append("Installed source is not a verified ancestor of current source")
    return result


def require_installable(state, reconcile_installed=None):
    if state["dirty"] or state["branch"] != "main" or state["revision"] != state["main_revision"]:
        raise RuntimeError("Release installation requires clean, integrated local main")
    if any(tree["dirty"] for tree in state["worktrees"]):
        raise RuntimeError(
            "Resolve uncommitted work in other worktrees before release installation"
        )
    installed = state["installed"]
    if not installed:
        return
    if "error" not in installed:
        # Applications use numeric release versions; fail closed on unknown version formats.
        if tuple(map(int, installed["version"].split("."))) > tuple(
            map(int, state["version"].split("."))
        ):
            raise RuntimeError("Refusing to downgrade the installed application version")
        if installed["dirty"] is False and installed["ancestor_of_source"]:
            return
    if not reconcile_installed or not reconcile_installed.strip():
        raise RuntimeError(
            "Installed provenance needs reconciliation. Review it first, then record the "
            "approved reason with --reconcile-installed REASON; source checks remain mandatory."
        )
    print(f"Installed provenance reconciliation: {reconcile_installed}", flush=True)


def report(root, destination=None, *, install=False, reconcile_installed=None):
    state = inspect(root, destination)
    print(json.dumps(state, indent=2), flush=True)
    for issue in issues(state):
        print(f"MISMATCH: {issue}", flush=True)
    if install:
        require_installable(state, reconcile_installed)
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path)
    parser.add_argument("--require-installable", action="store_true")
    parser.add_argument("--reconcile-installed")
    args = parser.parse_args()
    report(
        Path(__file__).resolve().parents[1],
        args.destination,
        install=args.require_installable,
        reconcile_installed=args.reconcile_installed,
    )


if __name__ == "__main__":
    main()
