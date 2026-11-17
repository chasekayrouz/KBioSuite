"""Official BLAST URL API and E-utilities/IPG, with scientific response validation."""

import io
import json
import re
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path

from Bio import SeqIO

from .lifecycle import Clock, record
from .models import Hit, Occurrence, Run, utc_now
from .network import RemoteError, SubmissionUnknown, Transport, check_cancel
from .parsers import ResponseParseError, ipg_xml_reports
from .storage import atomic_write, safe_name


def accession_matches(requested: str, returned: str) -> bool:
    returned = returned.split(":", 1)[0]
    if "." in requested:
        return requested == returned
    return requested == returned.split(".", 1)[0]


def genbank_records(text: str) -> list[tuple[str, object]]:
    """Split multi-record EFetch without rewriting source bytes or feature locations."""
    records, lines = [], []
    for line in text.splitlines(keepends=True):
        if not lines and not line.strip():
            continue
        if not lines and not line.startswith("LOCUS "):
            raise ResponseParseError("NCBI did not return complete GenBank records.")
        lines.append(line)
        if line.strip() == "//":
            raw = "".join(lines)
            try:
                record = SeqIO.read(io.StringIO(raw), "genbank")
            except ValueError as exc:
                raise ResponseParseError("Cannot parse an NCBI GenBank record.") from exc
            records.append((raw, record))
            lines = []
    if lines or not records:
        raise ResponseParseError("Incomplete NCBI GenBank response.")
    return records


class MappingUnresolved(ValueError):
    """NCBI has nucleotide links, but their annotations cannot establish occurrences."""


class NcbiClient:
    def __init__(self, transport: Transport, directory: Path, *, clock=None):
        self.http, self.directory = transport, directory
        self.local_seconds = 0.0
        self._timing_lock = threading.Lock()
        self._unindexed = set()
        self.clock = clock or Clock()

    @contextmanager
    def local_work(self):
        started = time.monotonic()
        try:
            yield
        finally:
            with self._timing_lock:
                self.local_seconds += time.monotonic() - started

    def batch_failure(self, operation, exc):
        self.http.log(
            {
                "at": utc_now(),
                "operation": operation,
                "outcome": "individual_fallback",
                "error": self.http.safe_message(str(exc)),
            }
        )

    def prefetch_hits(self, hits: list[Hit]) -> None:
        """Batch scientific downloads; failures are isolated by the individual path."""
        pending = [
            h
            for h in hits
            if h.sequence_status != "success"
            and not (self.directory / "raw" / f"hit_{h.rank:04d}.fasta").is_file()
        ]
        for offset in range(0, len(pending), 100):
            batch = pending[offset : offset + 100]
            try:
                text = self.http.request(
                    "eutils",
                    "efetch.fcgi",
                    {
                        "db": "protein",
                        "id": ",".join(h.accession for h in batch),
                        "rettype": "fasta",
                        "retmode": "text",
                    },
                )
                with self.local_work():
                    records = list(SeqIO.parse(io.StringIO(text), "fasta"))
                    for hit in batch:
                        matches = [
                            r
                            for r in records
                            if any(accession_matches(hit.accession, i) for i in r.id.split("|"))
                            and len(r) == hit.length
                            and r.seq.defined
                        ]
                        if len(matches) == 1:
                            # The source record is reconstructed only for FASTA, never GenBank.
                            atomic_write(
                                self.directory / "raw" / f"hit_{hit.rank:04d}.fasta",
                                matches[0].format("fasta"),
                            )
            except (RemoteError, ValueError, OSError) as exc:
                self.batch_failure("protein_batch", exc)
        pending = [
            h
            for h in hits
            if h.mapping_status not in {"success", "no_mapping"}
            and not (self.directory / "raw" / f"hit_{h.rank:04d}.ipg.tsv").is_file()
        ]
        for offset in range(0, len(pending), 100):
            batch = pending[offset : offset + 100]
            path = (
                self.directory / "raw" / f"batch_{batch[0].rank:04d}_{batch[-1].rank:04d}.ipg.xml"
            )
            try:
                text = self.http.request(
                    "eutils",
                    "efetch.fcgi",
                    {
                        "db": "protein",
                        "id": ",".join(h.accession for h in batch),
                        "rettype": "ipg",
                        "retmode": "xml",
                    },
                )
                atomic_write(path, self.http.safe_message(text))
                with self.local_work():
                    reports = ipg_xml_reports(text)
                    for hit in batch:
                        # A report is emitted for each submitted accession, even within one IPG.
                        matches = [r for r in reports if accession_matches(hit.accession, r[0])]
                        if len(matches) == 1:
                            atomic_write(
                                self.directory / "raw" / f"hit_{hit.rank:04d}.ipg.tsv",
                                matches[0][2],
                            )
            except (RemoteError, ValueError, OSError) as exc:
                if path.is_file():
                    path.replace(path.with_suffix(".ipg-error.txt"))
                self.batch_failure("ipg_batch", exc)

    def submit_blast(self, run: Run) -> tuple[str, int]:
        text = self.http.request(
            "blast",
            "",
            {
                "CMD": "Put",
                "QUERY": run.query.sequence,
                **(run.submitted_blast_parameters or run.settings.blast_parameters()),
            },
            submit=True,
        )
        rid = re.search(r"\bRID\s*=\s*([A-Z0-9-]+)", text)
        rtoe = re.search(r"\bRTOE\s*=\s*(\d+)", text)
        if not rid:
            raise SubmissionUnknown(
                "NCBI did not return a RID. Submission was not repeated "
                "because it may have succeeded."
            )
        result = rid.group(1), int(rtoe.group(1)) if rtoe else 60
        atomic_write(
            self.directory / "raw" / "blast-submission.txt", f"RID={result[0]}\nRTOE={result[1]}\n"
        )
        atomic_write(
            self.directory / "raw" / "blast-submission-response.html",
            self.http.safe_message(text),
        )
        return result

    def wait_for_blast(self, run: Run, checkpoint) -> str:
        path = self.directory / "raw" / "blast.xml"
        if path.is_file():
            return path.read_text()
        clock = self.clock
        started = clock.monotonic()

        def event(name, *, at=None, **details):
            record(run, self.directory, name, at=at or clock.now(), **details)
            checkpoint()

        previous_get = run.last_blast_get_at
        log_path = self.directory / "requests.jsonl"
        if log_path.is_file():
            for line in log_path.read_text().splitlines():
                try:
                    entry = json.loads(line)
                    if (
                        entry.get("service") == "blast"
                        and entry["parameters"].get("CMD") == "Get"
                        and entry["parameters"].get("RID") == run.rid
                    ):
                        previous_get = max(previous_get or "", entry["at"])
                except (ValueError, KeyError):
                    continue
        last = previous_get or run.rid_adopted_at or run.submitted_at
        elapsed = (clock.now() - datetime.fromisoformat(last)).total_seconds() if last else 60
        # RTOE is advisory. Keep the 10-second submission gate and minute RID gate.
        delay = max(0, (60 if previous_get else 10) - elapsed)
        if clock.wait(self.http.cancel, delay):
            from .network import Cancelled

            raise Cancelled()
        invalid = 0
        while clock.monotonic() - started < self.http.config.blast_timeout_seconds:
            # Persist intent before I/O; transport logs each retry's actual wire time.
            run.last_blast_get_at = clock.now().isoformat()
            poll_started = clock.monotonic()
            run.message = f"Waiting at NCBI · RID {run.rid} · polling at most once/minute"
            event("status_poll_started", phase="waiting")
            text = self.http.request(
                "blast", "", {"CMD": "Get", "RID": run.rid, "FORMAT_TYPE": "XML2_S"}, poll=True
            )
            is_xml = "BlastOutput2" in text or "BlastXML2" in text
            match = re.search(r"\bStatus\s*=\s*([A-Z]+)", text)
            status = "READY" if is_xml else match.group(1) if match else "UNRECOGNIZED"
            if not is_xml:
                stamp = clock.now().strftime("%Y%m%dT%H%M%S%fZ")
                safe = self.http.safe_message(text) if hasattr(self.http, "safe_message") else text
                atomic_write(self.directory / "raw" / f"blast-status-{stamp}.txt", safe)
            event("status_received", status=status)
            if status != run.ncbi_status:
                run.ncbi_status = status
                event("status_transition", status=status)
            if status == "READY":
                event("ready_detected", phase="ready")
            if is_xml:
                event(
                    "result_retrieval_started",
                    at=datetime.fromisoformat(run.last_blast_get_at),
                    phase="retrieving",
                    delivery="status_request",
                )
                event("result_received", phase="retrieving", delivery="status_request")
                run.timings["blast_result_retrieval_parse"] = (
                    run.timings.get("blast_result_retrieval_parse", 0.0)
                    + clock.monotonic()
                    - poll_started
                )
                return text  # Pipeline validates the report before using any hits.
            if status in {"FAILED", "UNKNOWN"}:
                raise RemoteError(f"BLAST RID {run.rid} is {status.lower()} at NCBI.")
            if status not in {"WAITING", "READY"}:
                invalid += 1
                event("protocol_warning", status=status, consecutive_invalid=invalid)
                if invalid >= getattr(self.http.config, "attempts", 3):
                    raise RemoteError(
                        "NCBI repeatedly returned neither a BLAST status nor XML2 results."
                    )
            else:
                invalid = 0
            if status == "READY":
                # Result retrieval is not a status poll; the transport retains its 10s gate.
                run.message = f"Results ready · RID {run.rid} · retrieving XML2"
                event("result_retrieval_started", phase="retrieving")
                retrieving = clock.monotonic()
                text = self.http.request(
                    "blast", "", {"CMD": "Get", "RID": run.rid, "FORMAT_TYPE": "XML2_S"}
                )
                if "BlastOutput2" in text or "BlastXML2" in text:
                    event("result_received", delivery="result_request")
                    run.timings["blast_result_retrieval_parse"] = (
                        run.timings.get("blast_result_retrieval_parse", 0.0)
                        + clock.monotonic()
                        - retrieving
                    )
                    return text
                raise RemoteError("NCBI marked this RID READY but did not return XML2 results.")
            run.message = f"Waiting at NCBI · latest status {status} · RID {run.rid}"
            checkpoint()
            last_poll = getattr(self.http, "last_blast_poll", poll_started)
            if clock.wait(self.http.cancel, max(0, 60 - (clock.monotonic() - last_poll))):
                from .network import Cancelled

                raise Cancelled()
        raise RemoteError("BLAST polling timed out. Resume this run later using the saved RID.")

    def fetch_ipg(self, accession: str, rank: int) -> str:
        path = self.directory / "raw" / f"hit_{rank:04d}.ipg.tsv"
        if path.is_file():
            return path.read_text()
        return self.http.request(
            "eutils",
            "efetch.fcgi",
            {
                "db": "protein",
                "id": accession,
                "rettype": "ipg",
                "retmode": "text",
            },
        )

    def linked_occurrences(
        self, hit: Hit, flank: int, members: set[str], *, accessions=None
    ) -> list[Occurrence]:
        """Resolve all history-linked records and every exact annotated protein CDS."""
        prefix = f"hit_{hit.rank:04d}"
        text = self.http.request(
            "eutils",
            "elink.fcgi",
            {
                "dbfrom": "protein",
                "db": "nuccore",
                "id": ",".join(accessions or [hit.accession]),
                "cmd": "neighbor_history",
                "retmode": "xml",
            },
        )
        atomic_write(self.directory / "raw" / f"{prefix}.links.xml", self.http.safe_message(text))
        try:
            root = ET.fromstring(text)
            linksets = root.findall("LinkSet")
            if root.tag != "eLinkResult" or len(linksets) != 1 or root.find(".//ERROR") is not None:
                raise ResponseParseError("Invalid NCBI protein-to-nucleotide link response.")
            linkset = linksets[0]
            if linkset.findtext("DbFrom") != "protein" or not linkset.findall("IdList/Id"):
                raise ResponseParseError("NCBI did not resolve the requested protein for ELink.")
            histories = [
                e for e in linkset.findall("LinkSetDbHistory") if e.findtext("DbTo") == "nuccore"
            ]
            if not histories:
                return []
            keys = [e.findtext("QueryKey") for e in histories]
            query_key = keys[0]
            webenv = linkset.findtext("WebEnv")
            if not all(keys) or not webenv:
                raise ResponseParseError("Incomplete NCBI nucleotide link history.")
            count_text = self.http.request(
                "eutils",
                "esearch.fcgi",
                {
                    "db": "nuccore",
                    "term": " OR ".join(f"#{key}" for key in keys),
                    "WebEnv": webenv,
                    "usehistory": "y",
                    "retmax": 0,
                    "retmode": "xml",
                },
            )
            atomic_write(
                self.directory / "raw" / f"{prefix}.linked-search.xml",
                self.http.safe_message(count_text),
            )
            search = ET.fromstring(count_text)
            if search.tag != "eSearchResult" or search.find(".//ERROR") is not None:
                raise ResponseParseError("Invalid linked-nucleotide count response.")
            count = int(search.findtext("Count", "-1"))
            if count < 1:
                raise ResponseParseError("NCBI returned an empty or invalid link history.")
            if len(keys) > 1:
                query_key = search.findtext("QueryKey")
                webenv = search.findtext("WebEnv") or webenv
                if not query_key:
                    raise ResponseParseError("NCBI did not return the combined nucleotide history.")
        except (ET.ParseError, TypeError, ValueError) as exc:
            raise ResponseParseError(str(exc)) from exc
        occurrences = []
        for offset in range(0, count, 10):
            check_cancel(self.http.cancel)
            raw = self.http.request(
                "eutils",
                "efetch.fcgi",
                {
                    "db": "nuccore",
                    "query_key": query_key,
                    "WebEnv": webenv,
                    "retstart": offset,
                    "retmax": min(10, count - offset),
                    "rettype": "gbwithparts",
                    "retmode": "text",
                    "strand": 1,
                },
            )
            evidence = f"{prefix}_{offset:06d}.linked.gb"
            atomic_write(self.directory / "raw" / evidence, self.http.safe_message(raw))
            with self.local_work():
                records = genbank_records(raw)
                if len(records) != min(10, count - offset):
                    raise ResponseParseError(
                        "NCBI omitted records from the linked-nucleotide page."
                    )
                for record_raw, record in records:
                    self.cache_record_info(record.id, record, record_raw, "ELink / EFetch GenBank")
                    for index, feature in enumerate(record.features):
                        proteins = feature.qualifiers.get("protein_id", [])
                        if feature.type != "CDS" or not any(
                            accession_matches(m, p) for m in members for p in proteins
                        ):
                            continue
                        protein = next(
                            p for p in proteins if any(accession_matches(m, p) for m in members)
                        )
                        loc = feature.location
                        if (
                            loc is None
                            or any(p.ref for p in loc.parts)
                            or loc.strand not in {1, -1}
                        ):
                            raise MappingUnresolved(
                                "Linked target CDS has an unsupported location."
                            )
                        low, high = int(loc.start) + 1, int(loc.end)
                        if not 1 <= low <= high <= len(record):
                            raise MappingUnresolved(
                                "Linked target CDS has invalid native coordinates."
                            )
                        n = len(occurrences) + 1
                        occurrences.append(
                            Occurrence(
                                id=f"h{hit.rank:04d}_l{n:06d}",
                                hit_rank=hit.rank,
                                hit_accession=hit.accession,
                                ipg_row=0,
                                source="GenBank CDS",
                                protein_accession=protein,
                                protein_name=feature.qualifiers.get("product", [""])[0],
                                nucleotide_accession=record.id,
                                cds_start=low,
                                cds_end=high,
                                requested_start=low - flank,
                                requested_end=high + flank,
                                strand="+" if loc.strand == 1 else "-",
                                organism=record.annotations.get("organism", ""),
                                assembly=";".join(
                                    x.split(":", 1)[1]
                                    for x in record.dbxrefs
                                    if x.startswith("Assembly:")
                                ),
                                mapping_method="ELink + GenBank CDS",
                                raw_mapping={
                                    "mapping_method": "ELink + GenBank CDS",
                                    "evidence_file": evidence,
                                    "feature_index": str(index),
                                    "protein_id": protein,
                                    "location": str(loc),
                                },
                            )
                        )
        if not occurrences:
            raise MappingUnresolved(
                "NCBI returned nucleotide links but no matching annotated protein CDS."
            )
        return occurrences

    def cache_record_info(self, requested, record, raw, method):
        if not accession_matches(requested, record.id) or len(record) < 1:
            raise ResponseParseError("Cannot resolve this nucleotide accession to a record length.")
        info = {
            "accession": record.id,
            "length": len(record),
            "method": method,
            "retrieved_at": utc_now(),
        }
        atomic_write(self.directory / "raw" / f"{safe_name(requested)}.length.gb", raw)
        atomic_write(
            self.directory / "cache" / f"{safe_name(requested)}.record.json",
            json.dumps(info, indent=2),
        )
        return info

    def prefetch_nucleotide_info(self, accessions: list[str]) -> None:
        pending = [
            a
            for a in dict.fromkeys(accessions)
            if not (self.directory / "cache" / f"{safe_name(a)}.record.json").is_file()
        ]
        for offset in range(0, len(pending), 200):
            batch = pending[offset : offset + 200]
            try:
                text = self.http.request(
                    "eutils",
                    "esummary.fcgi",
                    {
                        "db": "nuccore",
                        "id": ",".join(batch),
                        "retmode": "json",
                    },
                )
                atomic_write(
                    self.directory / "raw" / f"nucleotide_{offset:06d}.summary.json",
                    self.http.safe_message(text),
                )
                with self.local_work():
                    result = json.loads(text)["result"]
                    for accession in batch:
                        matches = [
                            result[u]
                            for u in result["uids"]
                            if accession_matches(accession, result[u].get("accessionversion", ""))
                            and int(result[u].get("slen", 0)) > 0
                        ]
                        if len(matches) == 1:
                            info = {
                                "accession": matches[0]["accessionversion"],
                                "length": int(matches[0]["slen"]),
                                "method": "E-utilities ESummary (batch)",
                                "retrieved_at": utc_now(),
                            }
                            atomic_write(
                                self.directory / "cache" / f"{safe_name(accession)}.record.json",
                                json.dumps(info, indent=2),
                            )
            except (RemoteError, ValueError, KeyError, TypeError, OSError) as exc:
                self.batch_failure("nucleotide_summary_batch", exc)
        self._unindexed.update(
            a
            for a in pending
            if not (self.directory / "cache" / f"{safe_name(a)}.record.json").is_file()
        )

    def fetch_protein(self, accession: str, rank: int) -> tuple[str, str]:
        path = self.directory / "raw" / f"hit_{rank:04d}.fasta"
        text = (
            path.read_text()
            if path.is_file()
            else self.http.request(
                "eutils",
                "efetch.fcgi",
                {
                    "db": "protein",
                    "id": accession,
                    "rettype": "fasta",
                    "retmode": "text",
                },
            )
        )
        record = SeqIO.read(io.StringIO(text), "fasta")
        # NCBI FASTA can use accession, sp|accession|name, or pdb|entry|chain identifiers.
        identifiers = record.id.split("|")
        if not any(accession_matches(accession, i) for i in identifiers):
            raise ValueError("The fetched protein accession does not match this BLAST subject.")
        if not record.seq or not record.seq.defined:
            raise ValueError("NCBI returned no protein sequence.")
        atomic_write(path, text)
        return str(record.seq), record.id

    def nucleotide_info(self, accession: str, *, allow_unindexed=False) -> dict:
        path = self.directory / "cache" / f"{safe_name(accession)}.record.json"
        if path.is_file():
            return json.loads(path.read_text())
        if allow_unindexed and accession in self._unindexed:
            # No invented contig length: EFetch's native REGION/full-record header
            # will validate the requested slice and establish any right-edge clip.
            return {"accession": accession, "length": None}
        info = None
        try:
            text = self.http.request(
                "eutils",
                "esummary.fcgi",
                {
                    "db": "nuccore",
                    "id": accession,
                    "retmode": "json",
                },
            )
            data = json.loads(text)["result"]
            if len(data["uids"]) == 1:
                summary = data[data["uids"][0]]
                version = summary.get("accessionversion", "")
                if accession_matches(accession, version) and int(summary.get("slen", 0)) > 0:
                    info = {
                        "accession": version,
                        "length": int(summary["slen"]),
                        "method": "E-utilities ESummary",
                        "retrieved_at": utc_now(),
                    }
                    atomic_write(
                        self.directory / "raw" / f"{safe_name(accession)}.summary.json", text
                    )
        except (RemoteError, ValueError, KeyError, TypeError):
            pass
        if info is None:
            # ESummary cannot resolve every unindexed sequence accession. EFetch can.
            text = self.http.request(
                "eutils",
                "efetch.fcgi",
                {
                    "db": "nuccore",
                    "id": accession,
                    "rettype": "gb",
                    "retmode": "text",
                    "complexity": 1,
                },
            )
            record = SeqIO.read(io.StringIO(text), "genbank")
            if not accession_matches(accession, record.id) or len(record) < 1:
                raise ValueError("Cannot resolve this nucleotide accession to a record length.")
            info = {
                "accession": record.id,
                "length": len(record),
                "method": "E-utilities EFetch full GenBank header",
                "retrieved_at": utc_now(),
            }
            atomic_write(self.directory / "raw" / f"{safe_name(accession)}.length.gb", text)
        atomic_write(path, json.dumps(info, indent=2))
        return info

    def fetch_context(
        self, accession: str, start: int, end: int, *, allow_clipping=False
    ) -> tuple[str, object, str]:
        check_cancel(self.http.cancel)
        path = self.directory / "cache" / f"{safe_name(accession)}_{start}_{end}.gb"
        fetched_at = utc_now()
        stamp = path.with_suffix(".json")
        if path.is_file() and stamp.is_file():
            text = path.read_text()
            fetched_at = json.loads(stamp.read_text())["retrieved_at"]
        else:
            full_path = self.directory / "raw" / f"{safe_name(accession)}.length.gb"
            info_path = self.directory / "cache" / f"{safe_name(accession)}.record.json"
            full = full_path.is_file() and info_path.is_file() and start == 1
            info = json.loads(info_path.read_text()) if full else {}
            if full and info.get("length") == end:
                candidate = full_path.read_text()
                with self.local_work():
                    full_record = SeqIO.read(io.StringIO(candidate), "genbank")
                if full_record.seq.defined:
                    self.save_context(
                        accession, start, end, candidate, full_record, info["retrieved_at"]
                    )
                    return candidate, full_record, info["retrieved_at"]
            text = self.http.request(
                "eutils",
                "efetch.fcgi",
                {
                    "db": "nuccore",
                    "id": accession,
                    "rettype": "gbwithparts",
                    "retmode": "text",
                    "seq_start": start,
                    "seq_stop": end,
                    "strand": 1,
                    "complexity": 1,
                },
            )
        with self.local_work():
            record = SeqIO.read(io.StringIO(text), "genbank")
            self.save_context(
                accession, start, end, text, record, fetched_at, allow_clipping=allow_clipping
            )
        return text, record, fetched_at

    def save_context(
        self, accession, start, end, text, record, fetched_at, *, allow_clipping=False
    ):
        if not accession_matches(accession, record.id):
            raise ValueError("The GenBank response accession does not match the requested record.")
        actual_end = start + len(record) - 1
        if (
            not record.seq.defined
            or len(record) < 1
            or (actual_end != end and not (allow_clipping and actual_end < end))
        ):
            raise ValueError("The GenBank region is incomplete or has no retrievable sequence.")
        if not record.annotations.get("molecule_type"):
            raise ValueError("The returned region lacks GenBank molecule metadata.")
        # Validate the native REGION header as well as length: an offset error can
        # otherwise produce a same-length slice of the wrong part of a contig.
        region = re.search(r"^ACCESSION[^\n]*\bREGION:\s*(\d+)\.\.(\d+)\s*$", text, re.MULTILINE)
        if region and (int(region[1]), int(region[2])) != (start, actual_end):
            raise ValueError("The GenBank REGION coordinates do not match the requested range.")
        if allow_clipping and not region and start != 1:
            raise ValueError("Unindexed context lacks a verifiable native REGION header.")
        if allow_clipping and (actual_end < end or not region):
            info = {
                "accession": record.id.split(":", 1)[0],
                "length": actual_end,
                "method": "E-utilities EFetch native record boundary",
                "retrieved_at": fetched_at,
            }
            atomic_write(
                self.directory / "cache" / f"{safe_name(accession)}.record.json",
                json.dumps(info, indent=2),
            )
        path = self.directory / "cache" / f"{safe_name(accession)}_{start}_{end}.gb"
        atomic_write(path, text)
        atomic_write(path.with_suffix(".json"), json.dumps({"retrieved_at": fetched_at}))

    def context_results(self, requests, whole_records, unindexed=()):
        """At most three in-flight slice requests, sharing the application's rate gate."""
        pending = list(dict.fromkeys(requests))
        # Full contigs share one EFetch parameter set and can be fetched together.
        full = [
            key
            for key in pending
            if key in whole_records
            and not (
                self.directory / "cache" / f"{safe_name(key[0])}_{key[1]}_{key[2]}.gb"
            ).is_file()
        ]
        batches, batch, bases = [], [], 0
        for key in full:
            if batch and (len(batch) == 20 or bases + key[2] > 2_000_000):
                batches.append(batch)
                batch, bases = [], 0
            batch.append(key)
            bases += key[2]
        if batch:
            batches.append(batch)
        for batch in batches:
            if len(batch) < 2:
                continue
            try:
                text = self.http.request(
                    "eutils",
                    "efetch.fcgi",
                    {
                        "db": "nuccore",
                        "id": ",".join(k[0] for k in batch),
                        "rettype": "gbwithparts",
                        "retmode": "text",
                        "strand": 1,
                        "complexity": 1,
                    },
                )
                fetched_at = utc_now()
                completed = []
                with self.local_work():
                    records = genbank_records(text)
                    for key in batch:
                        matches = [
                            (raw, r) for raw, r in records if accession_matches(key[0], r.id)
                        ]
                        if len(matches) == 1:
                            raw, record = matches[0]
                            try:
                                self.save_context(*key, raw, record, fetched_at)
                            except (ValueError, OSError) as exc:
                                self.batch_failure("whole_context_record", exc)
                                continue
                            completed.append((key, (raw, record, fetched_at)))
                for key, result in completed:
                    pending.remove(key)
                    yield key, result
            except (RemoteError, ValueError, OSError) as exc:
                self.batch_failure("whole_context_batch", exc)
        with ThreadPoolExecutor(max_workers=3, thread_name_prefix="ncbi-context") as pool:
            futures = {
                pool.submit(self.fetch_context, *key, allow_clipping=key in unindexed): key
                for key in pending
            }
            try:
                for future in as_completed(futures):
                    try:
                        result = future.result()
                    except Exception as exc:
                        result = exc
                    yield futures[future], result
            finally:
                for future in futures:
                    future.cancel()
