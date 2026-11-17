"""Explicit, resumable installers; verified PDB updates retire old managed generations."""

import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path

from .database import database_lease, database_status, inspect_database
from .models import Progress, write_json
from .process import environment, run_process

PDB_URL = "https://foldseek.steineggerlab.workers.dev/pdb100.tar.gz"
PDB_VERSION = "https://foldseek.steineggerlab.workers.dev/pdb100.version"
FOLDSEEK_URL = "https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz"
AF3_REV = "85c4d20505fd5cef05eac22b534d4e793971ae69"
BOLTZ_REV = "f2e51834a6e7ef6801e8c5cf8131cfbd8d2ad3ab"
# Kmatch does not bundle uv. Explicit Boltz/AF3 setup downloads this pinned release, checks it
# against the SHA-256 published with the release, and installs only its uv executable as
# <home>/tools/uv. Launch and status never download it.
UV_VERSION = "0.13.9"
UV_URL = (
    f"https://github.com/astral-sh/uv/releases/download/{UV_VERSION}/uv-aarch64-apple-darwin.tar.gz"
)
UV_SHA256 = "5aeca67f2283f75235de40d0bd9c0a5bbc5fce23b6d48b3bf9fa1b2db08382ca"
UV_MEMBER = "uv-aarch64-apple-darwin/uv"
UV_RECORD = "uv-installation.json"
COMMAND_LINE_TOOLS = (
    ["/usr/bin/git", "--version"],
    ["/usr/bin/otool", "--version"],
    # install_name_tool has no version option: run bare, it prints its usage and exits 1.
    ["/usr/bin/install_name_tool"],
)


class RemoteChangedError(RuntimeError):
    """The server's file changed (ETag or length) since a partial transfer began."""


@contextmanager
def installation_lock(root):
    root.mkdir(parents=True, exist_ok=True)
    with (root / ".setup.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("A Kmatch setup operation is already running.") from exc
        try:
            yield
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def download(url, destination, cancel, progress, log, *, expected_md5=None):
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if expected_md5:
            with destination.open("rb") as stream:
                valid = hashlib.file_digest(stream, "md5").hexdigest() == expected_md5
            if not valid:
                raise RuntimeError(
                    "Cached archive checksum mismatch. Preserve it and choose a "
                    "fresh setup generation."
                )
        return
    partial = destination.with_name(destination.name + ".part")
    headers = subprocess.check_output(
        [
            "/usr/bin/curl",
            "--fail",
            "--silent",
            "--show-error",
            "--head",
            "--location",
            "--max-time",
            "45",
            url,
        ],
        text=True,
        timeout=50,
    )
    final_headers = headers.strip().split("\n\n")[-1]
    values = {}
    for line in final_headers.splitlines():
        if ":" in line:
            key, value = line.split(":", 1)
            values[key.lower()] = value.strip()
    size, etag = values.get("content-length"), values.get("etag")
    metadata = partial.with_name(partial.name + ".source.json")
    source = {"url": url, "etag": etag, "content_length": size}
    if partial.exists() and metadata.exists() and json.loads(metadata.read_text()) != source:
        raise RemoteChangedError(
            "Remote download changed during a partial transfer. Preserve the "
            "partial and start a fresh generation."
        )
    write_json(metadata, source)
    if not (partial.exists() and size and partial.stat().st_size == int(size)):
        progress(
            Progress("setup", f"Downloading {destination.name}; partial transfers resume on retry.")
        )
        run_process(
            [
                "/usr/bin/curl",
                "--fail",
                "--location",
                "--retry",
                "5",
                "--connect-timeout",
                "30",
                "--continue-at",
                "-",
                "--output",
                partial,
                url,
            ],
            log,
            cancel,
            progress,
            "download",
            timeout=7200,
        )
    if size and partial.stat().st_size != int(size):
        raise RuntimeError(f"Incomplete download: {destination.name}. Retry Setup to resume.")
    if expected_md5:
        with partial.open("rb") as stream:
            actual = hashlib.file_digest(stream, "md5").hexdigest()
        if actual != expected_md5:
            raise RuntimeError(
                "PDB100 archive does not match the published checksum; it has not been activated."
            )
    partial.replace(destination)


def extract(archive, destination, cancel, progress):
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as tar:
        for index, member in enumerate(tar):
            cancel.check()
            if index % 2000 == 0:
                progress(Progress("setup", f"Extracting {archive.name} · {index} entries"))
            tar.extract(member, destination, filter="data")


def check_command_line_tools(engine):
    """Before any Boltz/AF3 setup download, require Apple's Command Line Tools (or Xcode).

    Boltz setup needs git (its pinned GitHub source), otool and install_name_tool (OpenMP
    relinking); AF3 setup clones its source with git. xcode-select runs first: without the
    tools, /usr/bin/git, otool and install_name_tool are stubs that open an installation dialog,
    so they run only once xcode-select reports an existing developer directory.
    """
    advice = "Run `xcode-select --install` in Terminal, then try again. Nothing was downloaded."

    def run(command):
        return subprocess.run(
            command,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            timeout=60,
            env=environment(),
        )

    try:
        selected = run(["/usr/bin/xcode-select", "-p"])
    except (OSError, subprocess.SubprocessError) as exc:
        raise RuntimeError(f"{engine} setup needs Apple's Command Line Tools. {advice}") from exc
    developer = selected.stdout.strip()
    if selected.returncode or not developer or not Path(developer).is_dir():
        raise RuntimeError(f"{engine} setup needs Apple's Command Line Tools. {advice}")
    for command in COMMAND_LINE_TOOLS:
        name = Path(command[0]).name
        try:
            result = run(command)
            output = (result.stdout or "") + (result.stderr or "")
            usage = len(command) == 1 and "Usage:" in output and name in output
            works = result.returncode == 0 or usage
        except (OSError, subprocess.SubprocessError) as exc:
            output, works = str(exc), False
        if not works:
            detail = next((line.strip() for line in output.splitlines() if line.strip()), "")
            raise RuntimeError(
                f"{engine} setup needs Apple's Command Line Tools, but `{name}` did not run"
                + (f" ({detail[:200]})" if detail else "")
                + ". Run `xcode-select --install` in Terminal (or open Xcode once to finish its "
                "setup), then try again. Nothing was downloaded."
            )


def file_sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def managed_uv(home):
    """The verified record of HOME's pinned uv (HOME/tools/uv), or None if setup must install it.

    The record must name the pinned version and archive, and the executable must still match
    the SHA-256 recorded when it was extracted. Nothing is executed or downloaded here.
    """
    tools = Path(home) / "tools"
    binary = tools / "uv"
    try:
        record = json.loads((tools / UV_RECORD).read_text())
        if (
            not isinstance(record, dict)
            or record.get("version") != UV_VERSION
            or record.get("archive_sha256") != UV_SHA256
            or binary.is_symlink()
            or not binary.is_file()
            or not os.access(binary, os.X_OK)
            or file_sha256(binary) != record.get("sha256")
        ):
            return None
    except (OSError, ValueError):
        return None
    return record


def install_uv(config, cancel, progress):
    """HOME/tools/uv at the pinned version, downloaded only by explicit Boltz/AF3 setup.

    PATH is never searched: an app opened from Finder does not see a shell's PATH, and only the
    pinned release is known to build the engines' locked environments. The archive must match
    UV_SHA256 before anything is extracted; only its uv executable is installed, replacing any
    previous copy atomically.
    """
    tools = config.home / "tools"
    if managed_uv(config.home):
        return tools / "uv"
    with installation_lock(tools):
        if managed_uv(config.home):
            return tools / "uv"
        staging = tools / ".uv-download"
        if staging.is_symlink():
            raise RuntimeError("The uv staging path is a symlink; it was not modified.")
        staging.mkdir(exist_ok=True)
        archive = staging / f"uv-{UV_VERSION}-aarch64-apple-darwin.tar.gz"
        # Partial transfers of an earlier pinned release cannot be resumed into this one.
        for stale in staging.iterdir():
            if not stale.name.startswith(archive.name) and not stale.is_dir():
                stale.unlink()
        progress(Progress("setup", f"Downloading uv {UV_VERSION} for engine setup…"))
        try:
            download(UV_URL, archive, cancel, progress, tools / "setup.log")
        except RemoteChangedError:
            # The pinned SHA-256 below verifies whatever arrives, so a partial transfer of a
            # changed remote file is discarded and the download restarts once from scratch.
            partial = archive.with_name(archive.name + ".part")
            partial.unlink(missing_ok=True)
            partial.with_name(partial.name + ".source.json").unlink(missing_ok=True)
            progress(Progress("setup", "The uv download changed on the server; restarting it."))
            download(UV_URL, archive, cancel, progress, tools / "setup.log")
        cancel.check()
        if file_sha256(archive) != UV_SHA256:
            archive.unlink()
            raise RuntimeError(
                f"The uv {UV_VERSION} download does not match its published SHA-256. It was "
                "discarded and nothing was installed; retry setup."
            )
        digest = replace_uv(archive, tools / "uv", cancel)
        write_json(
            tools / UV_RECORD,
            {
                "version": UV_VERSION,
                "url": UV_URL,
                "archive_sha256": UV_SHA256,
                "sha256": digest,
                "installed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            },
        )
        shutil.rmtree(staging)
        progress(Progress("setup", f"uv {UV_VERSION} verified and installed."))
    return tools / "uv"


def replace_uv(archive, target, cancel):
    """Extract only the archive's uv executable; atomically replace TARGET (mode 0755)."""
    with tarfile.open(archive, "r:gz") as tar:
        try:
            member = tar.getmember(UV_MEMBER)
        except KeyError:
            member = None
        if member is None or not member.isreg():
            raise RuntimeError(
                "The verified uv archive has no uv executable; nothing was installed."
            )
        source = tar.extractfile(member)
        handle, temporary = tempfile.mkstemp(prefix=".uv-", dir=target.parent)
        temporary = Path(temporary)
        digest = hashlib.sha256()
        try:
            with os.fdopen(handle, "wb") as stream:
                while chunk := source.read(1 << 20):
                    cancel.check()
                    digest.update(chunk)
                    stream.write(chunk)
                stream.flush()
                os.fsync(stream.fileno())
            temporary.chmod(0o755)
            cancel.check()
            os.replace(temporary, target)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    return digest.hexdigest()


def uv_provenance(config, uv):
    """The uv that built an engine environment, for that backend's installation.json."""
    record = managed_uv(config.home)
    if record and Path(uv) == config.home / "tools/uv":
        return {key: record[key] for key in ("version", "archive_sha256", "sha256")}
    return {"path": str(uv), "version": "unrecorded"}


def dispatch_uv_check(argv=None):
    """Explicit clean-Mac qualification of uv alone, with a new isolated home.

    No engine or model is installed, and no default application home is read. This
    probes the curl/download/extraction path on a Mac without Command Line Tools.
    """
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or arguments[0] != "--kmatch-uv-check":
        return None
    if len(arguments) != 3:
        raise SystemExit("Usage: --kmatch-uv-check NEW_HOME NEW_REPORT_PATH")
    from .config import Config
    from .process import Cancellation

    home, report = Path(arguments[1]), Path(arguments[2])
    if home.exists() or home.is_symlink():
        raise RuntimeError("The uv qualification home must be new")
    with report.open("x") as output:
        home.mkdir(parents=True, exist_ok=False)
        try:
            uv = install_uv(Config(home, home / "cache"), Cancellation(), lambda event: None)
            version = subprocess.check_output(
                [str(uv), "--version"], text=True, timeout=20, env=environment()
            ).strip()
            if not version.startswith(f"uv {UV_VERSION} ") and version != f"uv {UV_VERSION}":
                raise RuntimeError("Downloaded uv did not report the pinned version")
            result = {"passed": True, "version_output": version,
                      "uv": managed_uv(home), "engine_setup_attempted": False}
        except Exception as exc:
            result = {"passed": False, "error": str(exc), "error_type": type(exc).__name__}
        json.dump(result, output, indent=2)
        output.write("\n")
    return 0 if result["passed"] else 1


def install_foldseek(config, cancel, progress):
    if config.shared_resources:
        from .search import FoldseekBackend

        state = FoldseekBackend(config).executable_status()
        if not state.ready:
            raise RuntimeError(state.message)
        progress(
            Progress("setup", "Shared Foldseek is ready. Manage versions with KBioSuite resources.")
        )
        return
    if config.foldseek_override:
        raise RuntimeError(
            "Foldseek uses an external override. Manage that executable separately, or remove "
            "foldseek_override to install Kmatch's managed copy. No files were downloaded or "
            "changed."
        )
    root = config.home / "tools"
    with installation_lock(root):
        archive = root / "foldseek-osx-universal.tar.gz"
        download(FOLDSEEK_URL, archive, cancel, progress, root / "setup.log")
        extract(archive, root, cancel, progress)
        version = subprocess.check_output(
            [str(config.home / "tools/foldseek/bin/foldseek"), "version"], text=True, timeout=10
        ).strip()
        with archive.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        write_json(
            root / "foldseek-installation.json",
            {"version": version, "url": FOLDSEEK_URL, "archive_sha256": digest},
        )


def install_pdb(config, cancel, progress):
    if config.shared_resources:
        inspected = inspect_database(config.pdb_prefix, cancel)
        progress(
            Progress(
                "setup", "Shared PDB100 indices verified. Manage versions with KBioSuite resources."
            )
        )
        return inspected
    if config.pdb_prefix_override:
        raise RuntimeError(
            "PDB100 uses an external override. Manage that database separately, or remove "
            "pdb_prefix_override to update Kmatch's managed database. "
            "No external files were changed."
        )
    root = config.home / "databases/pdb100"
    with installation_lock(root):
        progress(Progress("setup", "Checking the official PDB100 snapshot…"))
        version_text = subprocess.check_output(
            [
                "/usr/bin/curl",
                "--fail",
                "--silent",
                "--show-error",
                "--location",
                "--max-time",
                "45",
                PDB_VERSION,
            ],
            text=True,
            timeout=50,
        )
        expected = version_text.splitlines()[0].split()[0]
        if len(expected) != 32 or any(c not in "0123456789abcdef" for c in expected):
            raise RuntimeError("Unrecognized official PDB100 checksum metadata.")
        date = next(
            (line.split()[0] for line in version_text.splitlines() if "PDB_DATE" in line), "unknown"
        )
        if not re.fullmatch(r"\d{6}", date):
            raise RuntimeError("Unrecognized official PDB100 date metadata.")
        cancel.check()
        generation = root / f"{date}-{expected[:12]}"
        if generation.is_symlink():
            raise RuntimeError("The generation path is a symlink; it was not modified.")
        manifest = generation / "installation.json"
        if database_status(generation / "pdb").ready and manifest.exists():
            progress(Progress("setup", "Official PDB100 generation is installed; verifying it…"))
            verification = validate_pdb_generation(config, generation, cancel, progress, expected)
        else:
            if generation.exists():
                raise RuntimeError(
                    "A damaged generation exists; it was preserved. Inspect setup.log."
                )
            if shutil.disk_usage(root).free < 12 * 1024**3:
                raise RuntimeError("PDB100 setup needs 12 GB free for staged download/extraction.")
            stage = root / (".staging-" + generation.name)
            if stage.is_symlink():
                raise RuntimeError("The staging path is a symlink; it was not modified.")
            stage.mkdir(exist_ok=True)
            archive = stage / "pdb100.tar.gz"
            download(PDB_URL, archive, cancel, progress, root / "setup.log", expected_md5=expected)
            extract(archive, stage, cancel, progress)
            # Never write or remove the active generation during extraction/validation.
            verification = validate_pdb_generation(config, stage, cancel, progress, expected)
            write_json(
                stage / "installation.json",
                {
                    "status": "complete",
                    "pdb_date": date,
                    "source": PDB_URL,
                    "published_md5": expected,
                    "upstream_version": version_text,
                    "installed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                },
            )
            cancel.check()
            stage.rename(generation)
        verification.update(
            official_version_url=PDB_VERSION,
            official_version=version_text,
            checked_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        )
        write_json(generation / "verification.json", verification)
        removed = activate_pdb_generation(config, generation, cancel, progress)
        write_json(
            root / "last-update.json",
            {
                "active": generation.name,
                "removed_generations": removed,
                **verification,
            },
        )
        progress(
            Progress(
                "setup",
                f"PDB100 {date} verified and active; {len(removed)} old generations removed.",
            )
        )


def validate_pdb_generation(config, generation, cancel, progress, expected_md5):
    """Validate all index ranges and an actual extract/search before promotion."""
    prefix = generation / "pdb"
    progress(Progress("setup", "Validating PDB100 indices and a local extract/search…"))
    inspected = inspect_database(prefix, cancel)
    archive = generation / "pdb100.tar.gz"
    if archive.is_file():
        with archive.open("rb") as stream:
            checksum = hashlib.file_digest(stream, "md5").hexdigest()
        if checksum != expected_md5:
            raise RuntimeError("PDB100 cached archive differs from the official checksum.")
        inspected["verified_archive_md5"] = checksum
    cancel.check()
    with tempfile.TemporaryDirectory(prefix=".verify-", dir=generation) as temporary:
        root = Path(temporary)
        selection = root / "key"
        selection.write_text(str(inspected["sample_key"]) + "\n")
        query, subset, result = root / "query.cif", root / "subset", root / "self.tsv"
        for command in (
            [config.foldseek, "createsubdb", selection, prefix, subset, "--id-mode", "0"],
            [
                config.foldseek,
                "convert2pdb",
                subset,
                query,
                "--pdb-output-format",
                "1",
                "--threads",
                "1",
            ],
            [
                config.foldseek,
                "easy-search",
                query,
                prefix,
                result,
                root / "work",
                "--format-output",
                "target,qtmscore",
                "--alignment-type",
                "2",
                "--exact-tmscore",
                "1",
                "--threads",
                str(config.threads),
                "--max-seqs",
                "10",
            ],
        ):
            run_process(
                command + ["-v", "1"],
                generation / "verification.log",
                cancel,
                progress,
                "database validation",
                timeout=900,
                scratch_root=config.cache / "process",
            )
        scores = [float(line.split()[1]) for line in result.read_text().splitlines()]
        if not scores or max(scores) < 0.99:
            raise RuntimeError("PDB100 validation could not recover its sampled structure.")
        inspected["sample_search_max_qtm"] = max(scores)
    return inspected


def activate_pdb_generation(config, generation, cancel, progress):
    root = config.home / "databases/pdb100"
    if (
        generation.is_symlink()
        or generation.parent != root
        or not database_status(generation / "pdb").ready
    ):
        raise RuntimeError("Cannot activate an unverified or unmanaged PDB100 generation.")
    current = root / "current"
    if current.exists() and not current.is_symlink():
        raise RuntimeError("The current database path is unmanaged; it was not replaced.")
    progress(Progress("setup", "Activating PDB100; waiting for active searches/exports if needed…"))
    removed = []
    with database_lease(config, cancel, exclusive=True):
        temporary = root / ".current-next"
        temporary.unlink(missing_ok=True)
        temporary.symlink_to(generation.name, target_is_directory=True)
        temporary.replace(current)
        # After promotion, finish retirement even if cancellation arrives. Never remove external,
        # unmarked, incomplete, or symlinked data. Active readers have released their lease.
        for old in sorted(root.iterdir()):
            if old == generation or old.is_symlink() or not old.is_dir():
                continue
            try:
                record = json.loads((old / "installation.json").read_text())
                name = f"{record['pdb_date']}-{record['published_md5'][:12]}"
            except (OSError, ValueError, KeyError, TypeError):
                continue
            if (
                record.get("status") == "complete"
                and record.get("source") == PDB_URL
                and old.name == name
                and re.fullmatch(r"\d{6}-[0-9a-f]{12}", old.name)
            ):
                shutil.rmtree(old)
                removed.append(old.name)
    return removed


def install_boltz(config, cancel, progress, uv):
    root = config.home / "backends/boltz-community"
    with installation_lock(root):
        requirements = Path(__file__).parent / "resources/boltz-requirements.txt"
        run_process(
            [uv, "venv", "--python", "3.12.12", "--allow-existing", root / ".venv"],
            root / "setup.log",
            cancel,
            progress,
            "environment setup",
        )
        run_process(
            [uv, "pip", "install", "--python", root / ".venv/bin/python", "-r", requirements],
            root / "setup.log",
            cancel,
            progress,
            "environment setup",
            timeout=1800,
        )
        run_process(
            [root / ".venv/bin/boltz-fix-macos-libomp"],
            root / "setup.log",
            cancel,
            progress,
            "OpenMP setup",
        )
        models = config.home / "models/boltz2"
        models.mkdir(parents=True, exist_ok=True)
        with installation_lock(models):
            for filename, url in [
                (
                    "mols.tar",
                    "https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar",
                ),
                ("boltz2_conf.ckpt", "https://model-gateway.boltz.bio/boltz2_conf.ckpt"),
                ("boltz2_aff.ckpt", "https://model-gateway.boltz.bio/boltz2_aff.ckpt"),
            ]:
                download(url, models / filename, cancel, progress, root / "setup.log")
            if not (models / "mols/.kmatch-complete").exists():
                extract(models / "mols.tar", models, cancel, progress)
                (models / "mols/.kmatch-complete").touch()
        write_json(
            root / "installation.json",
            {
                "version": "2.11.0",
                "source_revision": BOLTZ_REV,
                "python": "3.12.12",
                "uv": uv_provenance(config, uv),
                "status": "complete",
            },
        )


def install_af3(config, cancel, progress, uv):
    root = config.home / "backends/alphafold3"
    if root.exists() and not (root / ".git").is_dir():
        raise RuntimeError("AF3 directory already exists and is not a managed source checkout.")
    if not root.exists():
        run_process(
            [
                "/usr/bin/git",
                "clone",
                "--branch",
                "v3.0.4",
                "--depth",
                "1",
                "https://github.com/google-deepmind/alphafold3.git",
                root,
            ],
            config.home / "af3-setup.log",
            cancel,
            progress,
            "AF3 source setup",
        )
    revision = subprocess.check_output(
        ["/usr/bin/git", "-C", str(root), "rev-parse", "HEAD"], text=True
    ).strip()
    if revision != AF3_REV:
        raise RuntimeError(
            "Existing AF3 source differs from the qualified revision; "
            "configure it as an external override."
        )
    with installation_lock(root):
        run_process(
            [uv, "sync", "--locked", "--python", "3.12.12"],
            root / "setup.log",
            cancel,
            progress,
            "AF3 environment setup",
            cwd=root,
            timeout=1800,
        )
        run_process(
            [uv, "run", "--locked", "build_data"],
            root / "setup.log",
            cancel,
            progress,
            "AF3 chemical data setup",
            cwd=root,
            timeout=600,
        )
        write_json(
            root / "installation.json",
            {
                "source_revision": AF3_REV,
                "python": "3.12.12",
                "uv": uv_provenance(config, uv),
                "status": "complete",
            },
        )
    progress(
        Progress(
            "setup",
            "AF3 code ready. Obtain authorized parameters separately; see AF3 setup documentation.",
        )
    )
