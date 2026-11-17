import os
import re
import signal
import subprocess
import time

from ..cancellation import check_cancelled


class SearchError(RuntimeError):
    pass


def local_environment():
    # This is process-local, never a persistent user setting.
    return {**os.environ, "BLAST_USAGE_REPORT": "false"}


def _terminate_search(process):
    """Stop the owned process group, including RPS worker subprocesses, then reap it."""

    def send(sig):
        try:
            if os.name == "posix":
                os.killpg(process.pid, sig)
            elif sig == signal.SIGTERM:
                process.terminate()
            else:
                process.kill()
        except ProcessLookupError:
            pass

    send(signal.SIGTERM)
    try:
        process.communicate(timeout=1)
    except subprocess.TimeoutExpired:
        pass
    finally:
        # The leader may have exited while a worker ignored TERM.
        send(signal.SIGKILL)
        process.communicate()


def run_search(command: list[str], *, timeout=None):
    check_cancelled()
    started = time.monotonic()
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=local_environment(),
            start_new_session=os.name == "posix",
        )
    except OSError as exc:
        raise SearchError(str(exc)) from exc
    try:
        while True:
            check_cancelled()
            if timeout is not None and time.monotonic() - started >= timeout:
                raise SearchError("Command timed out")
            try:
                stdout, stderr = process.communicate(timeout=0.1)
                break
            except subprocess.TimeoutExpired:
                # This is a cancellation polling interval, not a search deadline.
                continue
        check_cancelled()
    except BaseException:
        _terminate_search(process)
        raise
    if process.returncode or re.search(
        r"(?m)^(?:BLAST (?:engine|Database) error:|Error: \[rpsblast\])", stderr
    ):
        raise SearchError(stderr.strip()[-1000:] or f"Exit status {process.returncode}")
    return stdout


def write_queries(proteins, path):
    # Internal identifiers are independent of arbitrary input headers and duplicate IDs.
    with path.open("w") as handle:
        for protein in proteins:
            check_cancelled()
            handle.write(f">{protein.id}\n{protein.sequence}\n")


def blast_database_argument(path: str) -> str:
    """BLAST parses -db as a list, so preserve one pathname with literal quotes."""
    if '"' in path:
        raise ValueError("BLAST database paths cannot contain double quotes")
    return f'"{path}"' if any(c.isspace() for c in path) else path
