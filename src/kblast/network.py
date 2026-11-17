"""Paced, cancel-aware HTTP. BLAST submissions are never automatically retried."""

import json
import re
import threading
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs

import requests

from . import __version__
from .config import Config
from .lifecycle import Clock
from .models import digest

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


class Cancelled(Exception):
    pass


class RemoteError(Exception):
    pass


class SubmissionUnknown(RemoteError):
    pass


def check_cancel(event: threading.Event) -> None:
    if event.is_set():
        raise Cancelled()


class RateLimiter:
    def __init__(self, interval: float, *, clock=None):
        self.interval = interval
        self.next_time = 0.0
        self.lock = threading.Lock()
        self.clock = clock or Clock()

    def wait(self, cancel: threading.Event) -> None:
        with self.lock:
            check_cancel(cancel)
            if self.clock.wait(cancel, max(0, self.next_time - self.clock.monotonic())):
                raise Cancelled()
            self.next_time = self.clock.monotonic() + self.interval


class Transport:
    def __init__(
        self,
        config: Config,
        cancel: threading.Event,
        log_path: Path,
        eutils_limiter: RateLimiter,
        blast_limiter: RateLimiter,
        *,
        clock=None,
    ):
        self.config, self.cancel, self.log_path = config, cancel, log_path
        self.eutils_limiter, self.blast_limiter = eutils_limiter, blast_limiter
        self._local = threading.local()
        self._sessions = []
        self._session_lock = threading.Lock()
        self._log_lock = threading.Lock()
        self.last_blast_poll = 0.0
        self.clock = clock or Clock()

    @property
    def session(self):
        # Requests sessions are owned by one thread; the service gates are shared.
        if not hasattr(self._local, "session"):
            session = requests.Session()
            session.headers["User-Agent"] = f"Kblast/{__version__} (local scientific client)"
            self._local.session = session
            with self._session_lock:
                self._sessions.append(session)
        return self._local.session

    def close(self) -> None:
        for session in self._sessions:
            session.close()

    def safe_message(self, message: str) -> str:
        if self.config.api_key:
            message = message.replace(self.config.api_key, "[redacted]")
        return re.sub(r"(?i)(api_key=)[^&\s]+", r"\1[redacted]", message)

    def log(self, entry: dict) -> None:
        with self._log_lock, self.log_path.open("a", encoding="utf-8") as handle:
            handle.write(self.safe_message(json.dumps(entry, ensure_ascii=False)) + "\n")

    def request(
        self, service: str, endpoint: str, params: dict[str, Any], *, submit=False, poll=False
    ) -> str:
        url = BLAST_URL if service == "blast" else EUTILS_URL + endpoint
        values = {**params, "tool": "Kblast"}
        if self.config.email:
            values["email"] = self.config.email
        if service == "eutils" and self.config.api_key:
            values["api_key"] = self.config.api_key
        safe_params = {k: v for k, v in values.items() if k not in {"api_key", "email", "QUERY"}}
        if "QUERY" in params:
            safe_params["query_sha256"] = digest(params["QUERY"])
        attempts = 1 if submit else self.config.attempts
        for attempt in range(1, attempts + 1):
            check_cancel(self.cancel)
            paced_at = self.clock.monotonic()
            if service == "blast":
                if poll:
                    if self.clock.wait(
                        self.cancel, max(0, self.last_blast_poll + 60 - self.clock.monotonic())
                    ):
                        raise Cancelled()
                self.blast_limiter.wait(self.cancel)
                if poll:
                    self.last_blast_poll = self.clock.monotonic()
            else:
                self.eutils_limiter.wait(self.cancel)
            entry = {
                "at": self.clock.now().isoformat(),
                "service": service,
                "endpoint": url,
                "method": "POST",
                "parameters": safe_params,
                "attempt": attempt,
                "poll": poll,
                "pace_seconds": round(self.clock.monotonic() - paced_at, 6),
            }
            self.log({**entry, "outcome": "started"})
            reason = ""
            delay = 2**attempt
            requested_at = self.clock.monotonic()
            try:
                # POST keeps API keys and unpublished queries out of URLs and access logs.
                response = self.session.post(url, data=values, timeout=self.config.timeout)
                entry["http_status"] = response.status_code
                entry["request_seconds"] = round(self.clock.monotonic() - requested_at, 6)
                prepared = getattr(response, "request", None)
                if prepared is not None and prepared.body:
                    wire = parse_qs(
                        prepared.body.decode()
                        if isinstance(prepared.body, bytes)
                        else prepared.body
                    )
                    # Evidence of the actual encoded POST, without credential/query values.
                    entry["sent_parameters"] = {
                        k: v[0] for k, v in wire.items() if k not in {"api_key", "email", "QUERY"}
                    }
                    entry["api_key_used"] = bool(wire.get("api_key"))
                    if "QUERY" in wire:
                        entry["query_sha256"] = digest(wire["QUERY"][0])
                        entry["query_length"] = len(wire["QUERY"][0])
                if response.ok:
                    result = response.text
                    # Only scientific payloads receive a response hash. An API error may echo a key.
                    is_scientific = (
                        "BlastOutput2" in result
                        or "BlastXML2" in result
                        or "<IPGReportSet" in result
                        or result.startswith(("LOCUS", ">", "Id\t"))
                        or (
                            endpoint == "esummary.fcgi"
                            and '"result"' in result
                            and '"error"' not in result.lower()
                        )
                    )
                    if self.config.api_key and self.config.api_key in result:
                        is_scientific = False
                    if re.search(r"<(?:ERROR|Error)\b", result):
                        is_scientific = False
                    if is_scientific:
                        entry["response_sha256"] = digest(result)
                    if service == "blast":
                        match = re.search(r"\bStatus\s*=\s*([A-Z]+)", result)
                        entry["ncbi_status"] = (
                            "READY"
                            if "BlastXML2" in result or "BlastOutput2" in result
                            else match.group(1)
                            if match
                            else "UNRECOGNIZED"
                        )
                    entry["received_at"] = self.clock.now().isoformat()
                    entry["outcome"] = "received"
                    self.log(entry)
                    return result
                retryable = response.status_code in {408, 429, 500, 502, 503, 504}
                reason = f"NCBI returned HTTP {response.status_code}."
                retry_after = response.headers.get("Retry-After", "")
                if retry_after:
                    try:
                        delay = max(delay, float(retry_after))
                    except ValueError:
                        try:
                            until = parsedate_to_datetime(retry_after)
                            delay = max(delay, (until - self.clock.now()).total_seconds())
                        except (ValueError, TypeError):
                            pass
            except requests.RequestException as exc:
                entry["request_seconds"] = round(self.clock.monotonic() - requested_at, 6)
                retryable = isinstance(
                    exc,
                    (
                        requests.Timeout,
                        requests.ConnectionError,
                        requests.exceptions.ChunkedEncodingError,
                        requests.exceptions.ContentDecodingError,
                    ),
                )
                reason = f"NCBI connection failed ({type(exc).__name__})."
            entry.update(outcome="failure", error=reason)
            self.log(entry)
            if submit:
                raise SubmissionUnknown(
                    reason + " Submission may have reached NCBI; it was not repeated. "
                    "Recover the RID if available, or explicitly start a new run."
                )
            if not retryable or attempt == attempts or delay > 120:
                raise RemoteError(reason + f" Stopped after {attempt} attempt(s).")
            if self.clock.wait(self.cancel, delay):
                raise Cancelled()
        raise RemoteError("NCBI request failed.")
