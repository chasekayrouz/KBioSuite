"""Lazy, process-isolated Pfam models for a qualified Apple Silicon runtime.

The cache sits outside this adapter. A cache hit never starts its child. Stop,
close, idle expiry and memory limits release the entire native model allocation.
No floating-point controls or reported scientific values are altered.
"""

import importlib.util
import json
import marshal
import os
import platform
import socket
import subprocess
import sys
import threading
import types
from importlib.metadata import version
from importlib.resources import files
from pathlib import Path

import psutil

from ..cancellation import RunCancelled, check_cancelled
from .common import SearchError, local_environment
from .hmmer import HMM_ARGUMENTS, parse_domtbl
from .identity import content_digest, file_identity, fingerprint, native_digest, reference_identity

PIN = "0.12.3"
IDLE_SECONDS = 300
MAX_RSS = 3 * 1024**3
EXPECTED_ARGUMENTS = ("--noali", "-E", "0.01", "--domE", "0.01")
PIPELINE_OPTIONS = {
    "E": 0.01,
    "domE": 0.01,
    "incE": 0.01,
    "incdomE": 0.01,
    "F1": 0.02,
    "F2": 0.001,
    "F3": 0.00001,
    "bias_filter": True,
    "null2": True,
    "seed": 42,
    "Z": None,
    "domZ": None,
    "T": None,
    "domT": None,
    "incT": None,
    "incdomT": None,
    "bit_cutoffs": None,
}


def package_identity(*, details=False):
    spec = importlib.util.find_spec("pyhmmer")
    if spec is None or version("pyhmmer") != PIN:
        raise ValueError("The qualified PyHMMER pin is unavailable")
    root = Path(spec.origin).parent
    assets = {}

    def normalize(code):
        return code.replace(
            co_filename="pyhmmer",
            co_consts=tuple(
                normalize(c) if isinstance(c, types.CodeType) else c for c in code.co_consts
            ),
        )

    # The pinned release's complete executable module set. Frozen and source
    # loaders return equivalent code; paths and appended signing blobs may differ.
    for name in (
        "pyhmmer",
        "pyhmmer._version",
        "pyhmmer.utils",
        "pyhmmer.hmmer",
        *(
            f"pyhmmer.hmmer.{n}"
            for n in (
                "_base",
                "_hmmpress",
                "_hmmscan",
                "_hmmsearch",
                "_jackhmmer",
                "_hmmalign",
                "_nhmmer",
                "_phmmer",
            )
        ),
    ):
        module = importlib.util.find_spec(name)
        assets[name] = fingerprint(marshal.dumps(normalize(module.loader.get_code(name))).hex())
    for name in (
        "plan7.abi3.so",
        "daemon.abi3.so",
        "easel.abi3.so",
        "errors.abi3.so",
        "platform/bsd.abi3.so",
    ):
        assets[name] = native_digest(root / name)
    return {
        "version": PIN,
        "content_sha256": fingerprint(assets),
        **({"assets": assets} if details else {}),
    }


def qualified_identity(adapter, reference):
    """Unknown builds/references use the CLI, without loading native models."""
    if sys.platform != "darwin" or platform.machine() != "arm64":
        return None
    try:
        receipt = json.loads(files("kannotate.resources").joinpath("resident.json").read_text())
        identity = {
            "backend": "pyhmmer-resident-mac-v1",
            "package": package_identity(),
            "cli_sha256": content_digest(tuple(file_identity(adapter.executable))),
            "reference_sha256": reference["installation_sha256"],
            "arguments": list(HMM_ARGUMENTS),
            "pipeline": PIPELINE_OPTIONS,
            "profiles": 30134,
        }
        if (
            receipt.get("status") == "passed"
            and receipt.get("qualified_identity") == identity
            and HMM_ARGUMENTS == EXPECTED_ARGUMENTS
        ):
            return identity
    except RunCancelled:
        raise
    except Exception:  # An unknown/partial package is ineligible, never a launch failure.
        pass
    return None


class ResidentHmmer:
    def __init__(self, fallback, identity, reference=None, *, minimum_queries=20):
        self.fallback = fallback
        self.executable, self.database = fallback.executable, fallback.database
        self.backend_identity = identity
        self.reference = reference
        self.runtime_options = {"minimum_resident_queries": minimum_queries}
        self.process = self.connection = None
        self.lock = threading.RLock()
        self.warnings = []
        self.last_call = {}
        self.last_backend = None
        self.disabled = False
        self.idle_timer = None
        self.lease = None

    def begin_run(self):
        self.warnings = []

    def close(self):
        with self.lock:
            self.lease = None
            if self.idle_timer is not None:
                self.idle_timer.cancel()
                self.idle_timer = None
            if self.process is not None:
                # This worker has native threads, never subprocesses. Popen's
                # signal methods first reap an already-exited child, avoiding a
                # signal to a recycled process-group ID after idle expiry.
                self.process.terminate()
                try:
                    if (
                        self.reference is not None
                        and reference_identity(self.reference["database"]) != self.reference
                    ):
                        raise SearchError("Resident Pfam reference changed")
                    self.process.wait(timeout=1)
                except subprocess.TimeoutExpired:
                    self.process.kill()
                    self.process.wait()
                self.process = None
            if self.connection is not None:
                self.connection.close()
                self.connection = None

    def _start(self):
        if self.idle_timer is not None:
            self.idle_timer.cancel()
        self.lease = object()
        if self.process is not None and self.process.poll() is None:
            return
        self.close()
        self.lease = object()
        parent, child = socket.socketpair()
        try:
            command = (
                [sys.executable, "--kannotate-resident-worker"]
                if getattr(sys, "frozen", False)
                else [sys.executable, "-m", "kannotate.search.resident"]
            )
            self.process = subprocess.Popen(
                [*command, str(child.fileno())],
                pass_fds=(child.fileno(),),
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                env=local_environment(),
                start_new_session=True,
            )
            self.connection = parent
            parent.settimeout(0.05)
        except BaseException:
            parent.close()
            raise
        finally:
            child.close()

    def _request(self, request):
        self._start()
        self.connection.sendall(json.dumps(request).encode() + b"\n")
        response = bytearray()
        while not response.endswith(b"\n"):
            check_cancelled()
            if psutil.Process(self.process.pid).memory_info().rss > MAX_RSS:
                raise SearchError("Resident Pfam exceeded its memory budget")
            try:
                part = self.connection.recv(4096)
            except TimeoutError:
                continue
            if not part or len(response) + len(part) > 65536:
                raise SearchError("Resident Pfam worker ended or returned invalid status")
            response.extend(part)
        result = json.loads(response)
        if not result.get("ok"):
            raise SearchError("Resident Pfam could not complete the search")
        self.last_call = result

    def _expire(self, lease):
        with self.lock:
            if self.lease is lease:
                self.close()

    def search(self, queries, workdir, threads):
        with self.lock:
            check_cancelled()
            from ..io.fasta import read_fasta

            if len(read_fasta(queries)) < self.runtime_options["minimum_resident_queries"]:
                # Mac whole-workflow measurements do not justify model loading for
                # small misses. Existing models keep their original idle expiry.
                self.last_backend = "canonical-cli"
                return self.fallback.search(queries, workdir, threads)
            if not self.disabled:
                try:
                    self._request(
                        {
                            "queries": str(queries),
                            "workdir": str(workdir),
                            "database": self.database,
                            "identity": self.backend_identity,
                            # PyHMMER >1 has an additional dispatcher thread.
                            "workers": 1 if threads < 3 else threads - 1,
                        }
                    )
                    check_cancelled()
                    hits = parse_domtbl(workdir / "resident.domtbl")
                    self.last_backend = "resident"
                    self.idle_timer = threading.Timer(IDLE_SECONDS, self._expire, (self.lease,))
                    self.idle_timer.daemon = True
                    self.idle_timer.start()
                    return hits
                except RunCancelled:
                    self.close()
                    raise
                except (OSError, ValueError, SearchError, psutil.Error):
                    self.close()
                    self.disabled = True
                    for name in ("resident.domtbl", "resident.json"):
                        (workdir / name).unlink(missing_ok=True)
                except BaseException:
                    self.close()
                    raise
            directory = workdir / "cli-fallback"
            directory.mkdir(exist_ok=True)
            self.warnings.append("Resident Pfam unavailable; used the canonical CLI.")
            self.last_backend = "canonical-cli"
            return self.fallback.search(queries, directory, threads)


def worker_main(fd):
    """Socket IPC also works in a windowed frozen app with no standard streams."""
    import time

    connection = socket.socket(fileno=int(fd))
    connection.settimeout(IDLE_SECONDS)
    stream = connection.makefile("rwb", buffering=0)
    profiles = alphabet = identity = None
    stage = "import"
    try:
        import pyhmmer

        from ..io.fasta import read_fasta

        while True:
            line = stream.readline(65536)
            if not line:
                return
            request = json.loads(line)
            stage = "runtime-identity"
            if request.get("inspect_identity"):
                stream.write(
                    json.dumps({"ok": True, "package": package_identity(details=True)}).encode()
                    + b"\n"
                )
                continue
            started = time.perf_counter()
            loaded = 0.0
            if identity is None:
                wanted = request["identity"]
                if package_identity() != wanted["package"] or pyhmmer.__version__ != PIN:
                    raise ValueError("Resident runtime changed")
                alphabet = pyhmmer.easel.Alphabet.amino()
                profiles = pyhmmer.plan7.OptimizedProfileBlock(alphabet)
                stage = "model-load"
                with pyhmmer.plan7.HMMFile(request["database"], alphabet=alphabet) as handle:
                    if not handle.is_pressed():
                        raise ValueError("Resident Pfam requires a pressed database")
                    for profile in handle.optimized_profiles():
                        profiles.append(profile)
                        if len(profiles) % 512 == 0:
                            if psutil.Process().memory_info().rss > MAX_RSS:
                                raise ValueError("Resident model memory budget exceeded")
                if len(profiles) != wanted["profiles"]:
                    raise ValueError("Resident profile count differs")
                identity = wanted
                loaded = time.perf_counter() - started
            if identity != request["identity"]:
                raise ValueError("Resident identity changed")
            proteins = read_fasta(Path(request["queries"]))
            stage = "digitalization"
            workers = min(request["workers"], len(proteins))
            digital = []
            for protein in proteins:
                query = pyhmmer.easel.TextSequence(name=protein.header, sequence=protein.sequence)
                item = query.digitize(alphabet)
                if item.textize().sequence != protein.sequence or len(item) != len(
                    protein.sequence
                ):
                    raise ValueError("Resident digitalization changed the sequence")
                digital.append(item)
            iterator = pyhmmer.hmmer.hmmscan(
                digital, profiles, cpus=workers, backend="threading", **PIPELINE_OPTIONS
            )
            count = 0
            stage = "search"
            try:
                with (Path(request["workdir"]) / "resident.domtbl").open("wb") as table:
                    for count, hits in enumerate(iterator, 1):
                        if hits.query.name != proteins[count - 1].header or hits.Z != len(profiles):
                            raise ValueError("Resident query order or search space changed")
                        hits.write(table, format="domains", header=count == 1)
                if count != len(proteins):
                    raise ValueError("Resident omitted a query")
            finally:
                iterator.close()
            report = {
                "ok": True,
                "pid": os.getpid(),
                "profiles": len(profiles),
                "load_seconds": loaded,
                "seconds": time.perf_counter() - started,
                "workers": workers,
                "rss_bytes": psutil.Process().memory_info().rss,
                "identity": identity,
            }
            (Path(request["workdir"]) / "resident.json").write_text(json.dumps(report, indent=2))
            stream.write(json.dumps(report).encode() + b"\n")
    except (Exception, KeyboardInterrupt) as error:
        # Parent falls back to CLI; never misrepresent a partial domain table as success.
        try:
            stream.write(
                json.dumps(
                    {
                        "ok": False,
                        "stage": stage,
                        "error": type(error).__name__,
                        "module": getattr(error, "name", None),
                    }
                ).encode()
                + b"\n"
            )
        except OSError:
            pass
        return
    finally:
        stream.close()
        connection.close()


if __name__ == "__main__":
    worker_main(sys.argv[1])
