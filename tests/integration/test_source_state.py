"""Release guards use temporary repositories and never touch installed applications."""

import importlib.util
import json
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "source_state", ROOT / "scripts/check_source_state.py"
)
state = importlib.util.module_from_spec(spec)
spec.loader.exec_module(state)


def git(root, *args):
    return subprocess.check_output(["git", "-C", str(root), *args], text=True).strip()


@pytest.fixture
def repository(tmp_path):
    root = tmp_path / "repo"
    root.mkdir()
    git(root, "init", "-b", "main")
    git(root, "config", "user.name", "Release test")
    git(root, "config", "user.email", "test@example.invalid")
    version = root / "src/kbiosuite/version.py"
    version.parent.mkdir(parents=True)
    version.write_text('__version__ = "0.3.3"\n')
    git(root, "add", "src")
    git(root, "commit", "-m", "fixture")
    return root, tmp_path / "installed"


def install_record(root, destination, **overrides):
    record = {
        "version": "0.3.3",
        "source": {"commit": git(root, "rev-parse", "HEAD"), "dirty": False},
    }
    record.update(overrides)
    destination.mkdir(exist_ok=True)
    (destination / ".kbiosuite-install.json").write_text(json.dumps(record))


def test_clean_main_and_installed_ancestor(repository):
    root, destination = repository
    state.require_installable(state.inspect(root, destination))
    install_record(root, destination)
    (root / "next").write_text("next")
    git(root, "add", "next")
    git(root, "commit", "-m", "descendant")
    result = state.inspect(root, destination)
    state.require_installable(result)
    assert "Installed revision differs from current source" in state.issues(result)


@pytest.mark.parametrize("change", ["dirty", "feature", "detached", "other_dirty"])
def test_source_cannot_be_bypassed(repository, change):
    root, destination = repository
    if change == "dirty":
        (root / "uncommitted").write_text("pending")
    elif change == "feature":
        git(root, "switch", "-c", "feature")
    elif change == "detached":
        git(root, "switch", "--detach")
    else:
        other = root.parent / "other"
        git(root, "worktree", "add", "-b", "feature", str(other))
        (other / "pending").write_text("pending")
    with pytest.raises(RuntimeError):
        state.require_installable(state.inspect(root, destination), "authorized reconciliation")


@pytest.mark.parametrize("case", ["dirty", "unknown", "missing", "divergent", "standalone"])
def test_installed_provenance_requires_review(repository, case):
    root, destination = repository
    revision = git(root, "rev-parse", "HEAD")
    if case == "divergent":
        git(root, "switch", "-c", "divergent")
        (root / "branch").write_text("other")
        git(root, "add", "branch")
        git(root, "commit", "-m", "different release")
        revision = git(root, "rev-parse", "HEAD")
        git(root, "switch", "main")
    if case == "missing":
        destination.mkdir()
    elif case == "standalone":
        # The retired standalone Kmatch provenance format is not suite build provenance.
        install_record(root, destination, source={"revision": revision, "dirty": False})
    else:
        install_record(
            root,
            destination,
            source={
                "commit": "0" * 40 if case == "unknown" else revision,
                "dirty": case == "dirty",
            },
        )
    result = state.inspect(root, destination)
    with pytest.raises(RuntimeError, match="reconciliation"):
        state.require_installable(result)
    state.require_installable(result, "reviewed previous local source")


def test_newer_installed_version_never_bypassed(repository):
    root, destination = repository
    install_record(root, destination, version="0.3.4")
    with pytest.raises(RuntimeError, match="downgrade"):
        state.require_installable(state.inspect(root, destination), "reviewed")


def test_old_feature_reports_main_mismatch(repository):
    root, destination = repository
    git(root, "branch", "old-feature")
    (root / "src/kbiosuite/version.py").write_text('__version__ = "0.3.4"\n')
    git(root, "commit", "-am", "new release")
    git(root, "switch", "old-feature")
    result = state.inspect(root, destination)
    assert not result["contains_main"]
    assert "Source version 0.3.3 differs from main 0.3.4" in state.issues(result)


def test_default_destination_is_the_suite_folder(repository, monkeypatch):
    root, destination = repository
    monkeypatch.setenv("HOME", str(destination.parent / "home"))
    result = state.inspect(root)
    assert result["destination"] == str(destination.parent / "home/Applications/KBioSuite")
    assert result["installed"] is None
