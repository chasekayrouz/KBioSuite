"""Bounded subprocess lifecycle shared by engines, search and setup."""

import os
import signal
import subprocess
import tempfile
import threading
import time
from pathlib import Path

from .models import Progress


class Cancelled(RuntimeError):
    pass


class ProcessError(RuntimeError):
    pass


class Cancellation:
    def __init__(self):
        self._event = threading.Event()

    def cancel(self):
        self._event.set()

    def check(self):
        if self._event.is_set():
            raise Cancelled("Run cancelled. Completed outputs and logs are retained.")


def environment(extra=None):
    keys = ("PATH", "HOME", "TMPDIR", "LANG", "LC_ALL", "SSL_CERT_FILE", "SSL_CERT_DIR")
    env = {key: os.environ[key] for key in keys if key in os.environ}
    env.update(
        PYTHONUNBUFFERED="1",
        PYTHONDONTWRITEBYTECODE="1",
        WANDB_MODE="disabled",
        DO_NOT_TRACK="1",
        PYTORCH_ENABLE_MPS_FALLBACK="1",
        OMP_NUM_THREADS="8",
    )
    env.update(extra or {})
    return env


def run_process(
    command,
    log: Path,
    cancel: Cancellation,
    progress,
    phase,
    *,
    cwd=None,
    timeout=3600,
    extra_env=None,
    scratch_root=None,
):
    cancel.check()
    log.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    # Library callers may supply their own cache; otherwise use ordinary OS temporary storage.
    # Never reload Kmatch configuration behind a library caller's back.
    if scratch_root is not None:
        scratch_root.mkdir(parents=True, exist_ok=True)
    with (
        log.open("a") as stream,
        tempfile.TemporaryDirectory(
            prefix="job-", dir=scratch_root, ignore_cleanup_errors=True
        ) as scratch,
    ):
        stream.write(f"\nStarting {Path(command[0]).name} · {phase}\n")
        stream.flush()
        try:
            process = subprocess.Popen(
                [str(item) for item in command],
                cwd=cwd,
                env=environment({"TMPDIR": scratch, **(extra_env or {})}),
                stdout=stream,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
        except OSError as exc:
            raise ProcessError(f"Cannot start {phase}: {exc}. See {log.name}.") from exc
        try:
            last_update = -5
            while process.poll() is None:
                cancel.check()
                elapsed = time.monotonic() - start
                if elapsed > timeout:
                    raise ProcessError(
                        f"{phase.capitalize()} exceeded {timeout:g} seconds. See {log.name}."
                    )
                if elapsed - last_update >= 2:
                    progress(
                        Progress(phase, f"{phase.capitalize()} running · {int(elapsed)} s", elapsed)
                    )
                    last_update = elapsed
                time.sleep(0.1)
            cancel.check()
            if process.returncode:
                lines = log.read_text(errors="replace").splitlines()
                detail = next(
                    (line.strip() for line in reversed(lines) if line.strip()),
                    "No engine diagnostics.",
                )
                raise ProcessError(
                    f"{phase.capitalize()} exited with code {process.returncode}: {detail[:400]} "
                    f"Open {log.name} in the run folder for details."
                )
        finally:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
            # The parent may exit while a child ignores SIGTERM. Its owned
            # process group must not leave an inference worker running.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    return time.monotonic() - start
