"""Explicit Boltz/AF3 setup prerequisites, offline: the pinned uv download and the
Command Line Tools pre-check. No network, no real setup, no developer tools are run."""

import hashlib
import io
import json
import os
import shutil
import subprocess
import tarfile
from types import SimpleNamespace

import pytest

from kmatch import cli
from kmatch import setup as kmatch_setup
from kmatch.config import Config
from kmatch.process import Cancellation, Cancelled
from kmatch.setup import UV_MEMBER, UV_RECORD, UV_URL, UV_VERSION, install_uv, managed_uv

PAYLOAD = b"#!/bin/sh\necho 'uv fixture'\n"


def fake_archive(path, payload=PAYLOAD, member=UV_MEMBER, symlink=False):
    """A stand-in for uv-aarch64-apple-darwin.tar.gz (uv plus uvx); returns its SHA-256."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as tar:
        directory = tarfile.TarInfo("uv-aarch64-apple-darwin")
        directory.type, directory.mode = tarfile.DIRTYPE, 0o755
        tar.addfile(directory)
        for name, data in ((member, payload), ("uv-aarch64-apple-darwin/uvx", b"uvx fixture")):
            info = tarfile.TarInfo(name)
            if symlink and name == member:
                info.type, info.linkname = tarfile.SYMTYPE, "/bin/sh"
                tar.addfile(info)
                continue
            info.size, info.mode = len(data), 0o755
            tar.addfile(info, io.BytesIO(data))
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def release(tmp_path, monkeypatch):
    """The pinned release served from a local fixture archive; downloads are recorded."""
    source = tmp_path / "release" / "uv-aarch64-apple-darwin.tar.gz"
    monkeypatch.setattr(kmatch_setup, "UV_SHA256", fake_archive(source))
    calls = []

    def download(url, destination, cancel, progress, log, **kwargs):
        calls.append((url, destination))
        if not destination.exists():
            shutil.copyfile(source, destination)

    monkeypatch.setattr(kmatch_setup, "download", download)
    config = Config(tmp_path / "home", tmp_path / "cache")
    return SimpleNamespace(source=source, calls=calls, config=config, tools=config.home / "tools")


def leftovers(tools):
    return sorted(p.name for p in tools.iterdir() if p.name.startswith(".uv"))


def test_setup_downloads_verifies_and_installs_only_the_pinned_uv(release):
    events = []
    uv = install_uv(release.config, Cancellation(), events.append)
    assert uv == release.tools / "uv" and uv.read_bytes() == PAYLOAD
    assert uv.stat().st_mode & 0o777 == 0o755
    assert [url for url, _ in release.calls] == [UV_URL]
    assert release.calls[0][1].parent == release.tools / ".uv-download"
    assert not (release.tools / "uvx").exists() and leftovers(release.tools) == []
    record = json.loads((release.tools / UV_RECORD).read_text())
    assert record["version"] == UV_VERSION == "0.13.9" and record["url"] == UV_URL
    assert record["archive_sha256"] == kmatch_setup.UV_SHA256
    assert record["sha256"] == hashlib.sha256(PAYLOAD).hexdigest()
    assert managed_uv(release.config.home) == record
    assert any("uv 0.13.9 verified" in event.message for event in events)
    # Installed and still matching its record: later setups reuse it without a download.
    assert install_uv(release.config, Cancellation(), events.append) == uv
    assert len(release.calls) == 1


def test_pinned_release_is_the_published_arm64_archive():
    assert kmatch_setup.UV_URL == (
        "https://github.com/astral-sh/uv/releases/download/0.13.9/uv-aarch64-apple-darwin.tar.gz"
    )
    assert kmatch_setup.UV_SHA256 == (
        "5aeca67f2283f75235de40d0bd9c0a5bbc5fce23b6d48b3bf9fa1b2db08382ca"
    )


def test_hash_mismatch_discards_the_download_and_keeps_the_previous_uv(release, monkeypatch):
    release.tools.mkdir(parents=True)
    (release.tools / "uv").write_bytes(b"previous uv")
    monkeypatch.setattr(kmatch_setup, "UV_SHA256", "0" * 64)
    with pytest.raises(RuntimeError, match="does not match its published SHA-256"):
        install_uv(release.config, Cancellation(), lambda event: None)
    assert (release.tools / "uv").read_bytes() == b"previous uv"
    assert not (release.tools / UV_RECORD).exists()
    assert list((release.tools / ".uv-download").iterdir()) == []  # Retry downloads afresh.


@pytest.mark.parametrize("problem", ["missing", "symlink"])
def test_archive_without_a_regular_uv_executable_installs_nothing(release, monkeypatch, problem):
    fake = release.source.with_name("other.tar.gz")
    if problem == "missing":
        digest = fake_archive(fake, member="uv-aarch64-apple-darwin/not-uv")
    else:
        digest = fake_archive(fake, symlink=True)
    shutil.copyfile(fake, release.source)
    monkeypatch.setattr(kmatch_setup, "UV_SHA256", digest)
    with pytest.raises(RuntimeError, match="has no uv executable"):
        install_uv(release.config, Cancellation(), lambda event: None)
    assert not (release.tools / "uv").exists() and not (release.tools / UV_RECORD).exists()
    assert leftovers(release.tools) == [".uv-download"]


def test_cancelled_transfer_installs_nothing_and_resumes_on_retry(tmp_path, monkeypatch):
    """The real download path: curl's partial survives cancellation and is resumed."""
    source = tmp_path / "release.tar.gz"
    monkeypatch.setattr(kmatch_setup, "UV_SHA256", fake_archive(source))
    data = source.read_bytes()
    headers = (
        f'HTTP/2 302\ncontent-length: 0\n\nHTTP/2 200\netag: "a"\ncontent-length: {len(data)}\n'
    )
    monkeypatch.setattr(kmatch_setup.subprocess, "check_output", lambda *a, **k: headers)
    transfers = []

    def curl(command, log, cancel, progress, phase, **kwargs):
        assert command[0] == "/usr/bin/curl" and command[-1] == UV_URL
        partial = command[command.index("--output") + 1]
        transfers.append(partial.stat().st_size if partial.exists() else 0)
        if len(transfers) == 1:
            partial.write_bytes(data[:100])
            cancel.cancel()
            cancel.check()
        with partial.open("ab") as stream:
            stream.write(data[partial.stat().st_size :])

    monkeypatch.setattr(kmatch_setup, "run_process", curl)
    config = Config(tmp_path / "home", tmp_path / "cache")
    tools = config.home / "tools"
    with pytest.raises(Cancelled):
        install_uv(config, Cancellation(), lambda event: None)
    assert not (tools / "uv").exists() and not (tools / UV_RECORD).exists()
    assert install_uv(config, Cancellation(), lambda event: None).read_bytes() == PAYLOAD
    assert transfers == [0, 100] and leftovers(tools) == []


def test_partial_of_a_changed_remote_file_is_discarded_and_restarted_once(tmp_path, monkeypatch):
    """A stale partial (ETag/length changed) restarts from zero; the pinned hash still decides."""
    source = tmp_path / "release.tar.gz"
    monkeypatch.setattr(kmatch_setup, "UV_SHA256", fake_archive(source))
    data = source.read_bytes()
    etag = {"value": '"old"'}

    def head(*args, **kwargs):
        return f"HTTP/2 200\netag: {etag['value']}\ncontent-length: {len(data)}\n"

    monkeypatch.setattr(kmatch_setup.subprocess, "check_output", head)
    transfers = []

    def curl(command, log, cancel, progress, phase, **kwargs):
        partial = command[command.index("--output") + 1]
        transfers.append(partial.stat().st_size if partial.exists() else 0)
        if len(transfers) == 1:
            partial.write_bytes(b"x" * 100)  # Bytes of the earlier remote object.
            cancel.cancel()
            cancel.check()
        with partial.open("ab") as stream:
            stream.write(data[partial.stat().st_size :])

    monkeypatch.setattr(kmatch_setup, "run_process", curl)
    config = Config(tmp_path / "home", tmp_path / "cache")
    tools = config.home / "tools"
    with pytest.raises(Cancelled):
        install_uv(config, Cancellation(), lambda event: None)
    staging = tools / ".uv-download"
    assert sorted(p.name.rsplit(".tar.gz", 1)[1] for p in staging.iterdir()) == [
        ".part",
        ".part.source.json",
    ]
    etag["value"] = '"new"'
    events = []
    assert install_uv(config, Cancellation(), events.append).read_bytes() == PAYLOAD
    assert transfers == [0, 0] and leftovers(tools) == []
    assert any("changed on the server; restarting" in event.message for event in events)


def test_cancellation_after_download_keeps_the_previous_uv(release, monkeypatch):
    release.tools.mkdir(parents=True)
    (release.tools / "uv").write_bytes(b"previous uv")
    cancel, download = Cancellation(), kmatch_setup.download

    def download_then_cancel(*args, **kwargs):
        download(*args, **kwargs)
        cancel.cancel()

    monkeypatch.setattr(kmatch_setup, "download", download_then_cancel)
    with pytest.raises(Cancelled):
        install_uv(release.config, cancel, lambda event: None)
    assert (release.tools / "uv").read_bytes() == b"previous uv"
    assert not (release.tools / UV_RECORD).exists()
    archive = release.calls[0][1]
    assert archive.read_bytes() == release.source.read_bytes()
    # The completed archive is kept; the retry verifies and installs it.
    monkeypatch.setattr(kmatch_setup, "download", download)
    assert install_uv(release.config, Cancellation(), lambda event: None).read_bytes() == PAYLOAD
    assert not archive.exists() and leftovers(release.tools) == []


def test_extraction_cancelled_midway_leaves_target_and_no_temporary(tmp_path):
    archive = tmp_path / "uv.tar.gz"
    fake_archive(archive, payload=b"x" * (3 << 20))
    target = tmp_path / "tools" / "uv"
    target.parent.mkdir()
    target.write_bytes(b"previous uv")
    cancel = Cancellation()

    class Midway:
        calls = 0

        def check(self):
            self.calls += 1
            if self.calls == 2:
                cancel.cancel()
            cancel.check()

    with pytest.raises(Cancelled):
        kmatch_setup.replace_uv(archive, target, Midway())
    assert target.read_bytes() == b"previous uv" and os.listdir(target.parent) == ["uv"]


def test_outdated_or_altered_uv_is_replaced_atomically(release):
    release.tools.mkdir(parents=True)
    old = release.tools / "uv"
    old.write_bytes(b"older pinned uv")
    old.chmod(0o755)
    record = {
        "version": "0.0.1",
        "archive_sha256": kmatch_setup.UV_SHA256,
        "sha256": hashlib.sha256(b"older pinned uv").hexdigest(),
    }
    (release.tools / UV_RECORD).write_text(json.dumps(record))
    assert managed_uv(release.config.home) is None
    with old.open("rb") as reader:
        inode = os.stat(old).st_ino
        uv = install_uv(release.config, Cancellation(), lambda event: None)
        # A rename, never an in-place rewrite: an open reader still sees the complete old file.
        assert reader.read() == b"older pinned uv"
    assert os.stat(uv).st_ino != inode and uv.read_bytes() == PAYLOAD
    assert leftovers(release.tools) == []
    # A later change to the executable invalidates it; setup installs the pinned one again.
    uv.write_bytes(b"altered")
    assert managed_uv(release.config.home) is None
    install_uv(release.config, Cancellation(), lambda event: None)
    assert uv.read_bytes() == PAYLOAD and len(release.calls) == 2


def test_setup_checks_tools_first_and_uses_only_the_managed_uv(tmp_path, monkeypatch):
    path_uv = tmp_path / "bin" / "uv"
    path_uv.parent.mkdir()
    path_uv.write_text("#!/bin/sh\nexit 1\n")
    path_uv.chmod(0o755)
    monkeypatch.setenv("PATH", f"{path_uv.parent}:{os.environ['PATH']}")
    config = Config(tmp_path / "home", tmp_path / "cache")
    order = []
    monkeypatch.setattr(cli, "check_command_line_tools", lambda engine: order.append(engine))

    def uv(config, cancel, progress):
        order.append("uv")
        return config.home / "tools/uv"

    monkeypatch.setattr(cli, "install_uv", uv)
    monkeypatch.setattr(cli, "install_boltz", lambda *args: order.append(("boltz", args[3])))
    monkeypatch.setattr(cli, "install_af3", lambda *args: order.append(("af3", args[3])))
    cli.setup_component("boltz", config, Cancellation(), lambda event: None)
    cli.setup_component("af3-code", config, Cancellation(), lambda event: None)
    managed = str(config.home / "tools/uv")
    assert order == ["Boltz", "uv", ("boltz", managed), "AlphaFold 3", "uv", ("af3", managed)]


def test_status_and_launch_never_install_uv(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("KMATCH_HOME", str(tmp_path / "home"))
    monkeypatch.setenv("KMATCH_CACHE", str(tmp_path / "cache"))

    def never(*args, **kwargs):
        pytest.fail("status must not check tools or install uv")

    monkeypatch.setattr(cli, "install_uv", never)
    monkeypatch.setattr(cli, "check_command_line_tools", never)
    cli.main(["status"])
    assert len(json.loads(capsys.readouterr().out)) == 4
    assert not (tmp_path / "home/tools").exists()


@pytest.mark.parametrize("backend", ["boltz", "af3"])
def test_backend_installation_records_the_uv_that_built_it(release, monkeypatch, backend):
    uv = install_uv(release.config, Cancellation(), lambda event: None)
    monkeypatch.setattr(kmatch_setup, "run_process", lambda *args, **kwargs: None)
    home = release.config.home
    if backend == "boltz":
        (home / "models/boltz2/mols").mkdir(parents=True)
        (home / "models/boltz2/mols/.kmatch-complete").touch()
        kmatch_setup.install_boltz(release.config, Cancellation(), lambda e: None, str(uv))
        root = home / "backends/boltz-community"
    else:
        root = home / "backends/alphafold3"
        (root / ".git").mkdir(parents=True)
        fake = SimpleNamespace(check_output=lambda *a, **k: kmatch_setup.AF3_REV + "\n")
        monkeypatch.setattr(kmatch_setup, "subprocess", fake)
        kmatch_setup.install_af3(release.config, Cancellation(), lambda e: None, str(uv))
    record = json.loads((root / "installation.json").read_text())
    assert record["status"] == "complete"
    assert record["uv"] == {
        "version": UV_VERSION,
        "archive_sha256": kmatch_setup.UV_SHA256,
        "sha256": hashlib.sha256(PAYLOAD).hexdigest(),
    }


# Command Line Tools pre-check: subprocess is mocked, so no developer tool or stub runs.

TOOLS_MESSAGE = (
    "Boltz setup needs Apple's Command Line Tools. Run `xcode-select --install` in Terminal, "
    "then try again. Nothing was downloaded."
)


def developer_tools(monkeypatch, answers):
    """Answer each command by its executable name; record what ran."""
    ran = []

    def run(command, **kwargs):
        ran.append(command)
        assert kwargs["stdin"] == subprocess.DEVNULL and kwargs["timeout"]
        answer = answers[os.path.basename(command[0])]
        if isinstance(answer, BaseException):
            raise answer
        code, output = answer
        return subprocess.CompletedProcess(command, code, output, "")

    fake = SimpleNamespace(
        run=run, DEVNULL=subprocess.DEVNULL, SubprocessError=subprocess.SubprocessError
    )
    monkeypatch.setattr(kmatch_setup, "subprocess", fake)
    return ran


@pytest.mark.parametrize(
    "answer",
    [
        (2, "xcode-select: error: unable to get active developer directory"),
        (0, "/Library/Developer/CommandLineTools-removed\n"),
        (0, ""),
        FileNotFoundError("xcode-select"),
        subprocess.TimeoutExpired("xcode-select", 60),
    ],
)
def test_missing_command_line_tools_stop_setup_before_any_stub_runs(monkeypatch, answer):
    ran = developer_tools(monkeypatch, {"xcode-select": answer})
    with pytest.raises(RuntimeError) as failure:
        kmatch_setup.check_command_line_tools("Boltz")
    assert str(failure.value) == TOOLS_MESSAGE
    assert ran == [["/usr/bin/xcode-select", "-p"]]


USAGE = "Usage: /Library/Developer/CommandLineTools/usr/bin/install_name_tool [-change old new]"


def test_working_command_line_tools_pass(tmp_path, monkeypatch):
    ran = developer_tools(
        monkeypatch,
        {
            "xcode-select": (0, f"{tmp_path}\n"),
            "git": (0, "git version 2.50.1 (Apple Git-155)"),
            "otool": (0, "otool(1): Apple Inc. version cctools-1030"),
            "install_name_tool": (1, USAGE),
        },
    )
    kmatch_setup.check_command_line_tools("Boltz")
    assert ran == [
        ["/usr/bin/xcode-select", "-p"],
        ["/usr/bin/git", "--version"],
        ["/usr/bin/otool", "--version"],
        ["/usr/bin/install_name_tool"],
    ]


@pytest.mark.parametrize(
    "broken, answer",
    [
        ("git", (69, "You have not agreed to the Xcode license agreements.")),
        ("otool", FileNotFoundError("otool")),
        ("install_name_tool", (1, "xcrun: error: invalid active developer path")),
    ],
)
def test_a_broken_developer_tool_is_named_in_the_message(tmp_path, monkeypatch, broken, answer):
    answers = {
        "xcode-select": (0, str(tmp_path)),
        "git": (0, "git version 2.50.1"),
        "otool": (0, "otool(1)"),
        "install_name_tool": (1, USAGE),
        broken: answer,
    }
    ran = developer_tools(monkeypatch, answers)
    with pytest.raises(RuntimeError) as failure:
        kmatch_setup.check_command_line_tools("AlphaFold 3")
    message = str(failure.value)
    assert message.startswith(
        f"AlphaFold 3 setup needs Apple's Command Line Tools, but `{broken}` did not run"
    )
    assert "xcode-select --install" in message and message.endswith("Nothing was downloaded.")
    assert os.path.basename(ran[-1][0]) == broken


def test_cli_setup_without_command_line_tools_downloads_nothing(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("KMATCH_HOME", str(tmp_path / "home"))
    monkeypatch.setenv("KMATCH_CACHE", str(tmp_path / "cache"))
    developer_tools(monkeypatch, {"xcode-select": (2, "")})

    def offline(*args, **kwargs):
        pytest.fail("nothing may be downloaded or built")

    monkeypatch.setattr(kmatch_setup, "download", offline)
    monkeypatch.setattr(kmatch_setup, "run_process", offline)
    with pytest.raises(SystemExit) as stopped:
        cli.main(["setup", "boltz"])
    assert stopped.value.code == 1
    assert capsys.readouterr().err.strip() == TOOLS_MESSAGE
    assert not (tmp_path / "home").exists()
