"""Literal subject-to-IPG expansion. Each independent item has its own failure state."""

import csv
import io
import re
from contextlib import contextmanager
from datetime import UTC, datetime

from .lifecycle import Clock, record, remember_rid, update_rid_result
from .models import NO_MAPPING_MESSAGE, TARGET_CDS_ERROR, Run, digest
from .ncbi import MappingUnresolved, NcbiClient
from .network import Cancelled, RemoteError, check_cancel
from .parsers import ResponseParseError, parse_blast_xml, parse_ipg
from .storage import RunStore, atomic_write

ITEM_ERRORS = (RemoteError, ValueError, KeyError, TypeError, OSError)


def target_cds_present(record, occurrence) -> bool:
    expected_strand = {"+": 1, "-": -1}.get(occurrence.strand)
    if expected_strand is None:
        return False
    target_start = occurrence.cds_start - occurrence.context_start
    target_end = occurrence.cds_end - occurrence.context_start + 1
    for feature in record.features:
        if feature.type == "CDS" and feature.location is not None:
            if any(part.ref for part in feature.location.parts):
                continue
            try:
                low, high = int(feature.location.start), int(feature.location.end)
            except (TypeError, ValueError):
                continue
            if low == target_start and high == target_end:
                if expected_strand == feature.location.strand:
                    return True
    return False


class Pipeline:
    def __init__(self, store: RunStore, run: Run, client: NcbiClient):
        self.store, self.run, self.client = store, run, client
        self.directory = store.directory(run.id)
        self.checkpoint_seconds = 0.0
        self.clock = getattr(client, "clock", Clock())

    @contextmanager
    def timed(self, stage):
        started = self.clock.monotonic()
        try:
            yield
        finally:
            self.run.timings[stage] = (
                self.run.timings.get(stage, 0.0) + self.clock.monotonic() - started
            )

    def checkpoint(self, message: str | None = None) -> None:
        if message:
            self.run.message = message
        started = self.clock.monotonic()
        self.store.save(self.run)
        self.checkpoint_seconds += self.clock.monotonic() - started

    def blast_event(self, name, **details):
        record(self.run, self.directory, name, at=self.clock.now(), **details)

    def execute(self) -> None:
        if self.run.rid:
            self.blast_event("resumed", phase="waiting")
        run = self.run
        prior_client_local = getattr(self.client, "local_seconds", 0.0)
        prior_checkpoints = self.checkpoint_seconds
        for stage in (
            "blast_submission_wait",
            "blast_result_retrieval_parse",
            "protein_mapping",
            "nucleotide_contexts",
            "local_processing",
        ):
            run.timings.setdefault(stage, 0.0)
        run.status, run.error = "running", ""
        self.checkpoint("Starting gene-neighborhood retrieval")
        try:
            self.blast()
            if not run.hits:
                run.finish()
                self.checkpoint()
                return
            with self.timed("protein_mapping"):
                self.map_hits()
            with self.timed("nucleotide_contexts"):
                self.retrieve_contexts()
            run.finish()
            self.checkpoint()
        except Cancelled:
            self.blast_event("cancelled", phase="cancelled")
            run.status = "cancelled"
            for h in run.hits:
                if h.mapping_status in {"pending", "running"}:
                    h.mapping_status = "cancelled"
                if h.sequence_status in {"pending", "running"}:
                    h.sequence_status = "cancelled"
            for o in run.occurrences:
                if o.status in {"pending", "running"}:
                    o.status = "cancelled"
            self.checkpoint(
                "Cancelled locally. Successful results and the BLAST RID are preserved."
            )
        except Exception as exc:
            run.status = "failure"
            run.error = self.client.http.safe_message(f"{type(exc).__name__}: {exc}")
            self.blast_event("failed", phase="failed", error=run.error)
            self.checkpoint("Run stopped. Saved results remain available for export or resume.")
        finally:
            self.client.http.close()
            self.run.timings["local_processing"] = (
                self.run.timings.get("local_processing", 0.0)
                + getattr(self.client, "local_seconds", 0.0)
                - prior_client_local
                + self.checkpoint_seconds
                - prior_checkpoints
            )
            self.checkpoint()

    def blast(self) -> None:
        run = self.run
        path = self.directory / "raw" / "blast.xml"
        if path.is_file():
            self.blast_event(
                "local_result_reused", phase="complete" if run.blast_metadata else "parsing"
            )
            if not run.blast_metadata:
                self.blast_event("parser_started", phase="parsing")
                with self.timed("blast_result_retrieval_parse"), self.timed("local_processing"):
                    self.parse_blast(path.read_text())
                run.blast_completed_at = self.clock.now().isoformat()
                self.blast_event("parser_completed", phase="complete", hits=len(run.hits))
                self.checkpoint()
            return
        run.stage = "blast"
        started = self.clock.monotonic()
        prior_retrieval = run.timings.get("blast_result_retrieval_parse", 0.0)
        try:
            self.ensure_blast_submission()
            text = self.client.wait_for_blast(run, self.checkpoint)
        finally:
            run.timings["blast_submission_wait"] = run.timings.get(
                "blast_submission_wait", 0.0
            ) + max(
                0.0,
                self.clock.monotonic()
                - started
                - (run.timings.get("blast_result_retrieval_parse", 0.0) - prior_retrieval),
            )
        # Keep the actual server report even when its program/database is rejected.
        self.blast_event("parser_started", phase="parsing")
        self.checkpoint("Parsing NCBI BLAST XML2 results")
        with self.timed("blast_result_retrieval_parse"), self.timed("local_processing"):
            atomic_write(path, text)
            self.parse_blast(text)
        run.blast_completed_at = self.clock.now().isoformat()
        self.blast_event("parser_completed", phase="complete", hits=len(run.hits))
        self.checkpoint(f"Received {len(run.hits)} BLAST subject hits in NCBI rank order")

    def parse_blast(self, text):
        run = self.run
        try:
            hits, metadata = parse_blast_xml(
                text,
                run.settings.top_n,
                len(run.query.sequence),
                "nr" if run.rid_source in {"external", "browser"} else run.settings.database,
                expected_query_sequence=run.query.sequence,
            )
        except ValueError as exc:
            run.rid_validation = {"status": "rejected", "error": str(exc)}
            update_rid_result(run)
            self.blast_event("rid_result_rejected", error=str(exc))
            raise
        run.hits, run.blast_metadata = hits, metadata
        run.rid_validation = {"status": "validated", **metadata["query_validation"]}
        update_rid_result(run)
        if metadata["returned_hit_count"] < run.settings.top_n:
            self.blast_event(
                "fewer_hits_available",
                requested=run.settings.top_n,
                returned=metadata["returned_hit_count"],
            )

    def ensure_blast_submission(self) -> None:
        run = self.run
        if not run.rid:
            receipt = self.directory / "raw" / "blast-submission.txt"
            if receipt.is_file():
                saved = receipt.read_text()
                rid = re.search(r"^RID=([A-Z0-9-]+)$", saved, re.MULTILINE)
                rtoe = re.search(r"^RTOE=(\d+)$", saved, re.MULTILINE)
                if rid:
                    run.rid = rid.group(1)
                    run.rtoe = int(rtoe.group(1)) if rtoe else 60
                    run.submitted_at = datetime.fromtimestamp(
                        receipt.stat().st_mtime, UTC
                    ).isoformat()
                    run.submission_attempted = True
                    run.rid_source = "kblast_url_api"
                    run.original_rid = run.original_rid or run.rid
                    remember_rid(run)
                    self.blast_event("rid_recovered", phase="waiting", rtoe_seconds=run.rtoe)
                    self.checkpoint("Recovered the saved NCBI submission receipt")
        if not run.rid:
            if run.submission_attempted:
                raise RemoteError(
                    "A prior submission has no saved RID; automatic resubmission is disabled. "
                    "Recover its RID or explicitly start a new run."
                )
            run.submission_attempted = True
            run.submitted_blast_parameters = run.settings.blast_parameters()
            run.submitted_at = self.clock.now().isoformat()
            self.blast_event(
                "submission_started", phase="submitting", parameters=run.submitted_blast_parameters
            )
            self.checkpoint(f"Submitting BLASTP to NCBI {run.settings.database}")
            run.rid, run.rtoe = self.client.submit_blast(run)
            run.rid_source = "kblast_url_api"
            run.original_rid = run.original_rid or run.rid
            remember_rid(run)
            self.blast_event("rid_received", phase="waiting", rtoe_seconds=run.rtoe)
            self.checkpoint(f"Submitted · RID {run.rid} · NCBI estimated {run.rtoe} seconds")

    def map_hits(self) -> None:
        run = self.run
        run.stage = "mapping"
        self.checkpoint("Batching protein and IPG retrieval")
        self.client.prefetch_hits(run.hits)
        for hit in run.hits:
            check_cancel(self.client.http.cancel)
            if hit.sequence_status != "success":
                hit.sequence_status, hit.sequence_error = "running", ""
                self.checkpoint(
                    f"Fetching full protein · {hit.rank}/{len(run.hits)} · {hit.accession}"
                )
                try:
                    hit.sequence, hit.sequence_accession = self.client.fetch_protein(
                        hit.accession, hit.rank
                    )
                    if len(hit.sequence) != hit.length:
                        raise ValueError("Protein length differs from the BLAST subject length.")
                    hit.sequence_status = "success"
                except ITEM_ERRORS as exc:
                    hit.sequence, hit.sequence_status = "", "retrieval_failure"
                    hit.sequence_error = self.client.http.safe_message(str(exc))
                self.checkpoint()
            check_cancel(self.client.http.cancel)
            if hit.mapping_status in {"success", "no_mapping"}:
                continue
            hit.mapping_status, hit.mapping_error = "running", ""
            self.checkpoint(f"Expanding IPG · {hit.rank}/{len(run.hits)} · {hit.accession}")
            try:
                text = self.client.fetch_ipg(hit.accession, hit.rank)
                # Preserve the complete response even if parsing fails.
                path = self.directory / "raw" / f"hit_{hit.rank:04d}.ipg.tsv"
                atomic_write(path, self.client.http.safe_message(text))
                try:
                    with self.timed("local_processing"):
                        rows = parse_ipg(text, hit, run.settings.flank)
                except ValueError:
                    # An error response is evidence, not a usable cache on resume.
                    path.replace(path.with_name(f"hit_{hit.rank:04d}.ipg-error.txt"))
                    raise
                # Publish authoritative IPG rows before a supplementary request so
                # a link failure does not discard these independent resolved facts.
                run.occurrences = [o for o in run.occurrences if o.hit_rank != hit.rank] + rows
                raw_rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
                missing = list(
                    dict.fromkeys(
                        r["Protein"]
                        for r in raw_rows
                        if r.get("Protein")
                        and r.get("Source") in {"RefSeq", "INSDC"}
                        and r.get("Nucleotide Accession", "").strip() in {"", "N/A", "-"}
                    )
                )
                if not rows or missing:
                    members = {hit.accession} | {r["Protein"] for r in raw_rows if r.get("Protein")}
                    self.checkpoint(f"Checking NCBI nucleotide links · {hit.accession}")
                    additional = self.client.linked_occurrences(
                        hit,
                        run.settings.flank,
                        members,
                        accessions=list(
                            dict.fromkeys(([hit.accession] if not rows else []) + missing)
                        ),
                    )

                    # Reconcile duplicate evidence paths, preserving every original
                    # IPG row and every distinct CDS feature returned by the fallback.
                    def location_key(o):
                        return (
                            o.protein_accession,
                            o.nucleotide_accession,
                            o.cds_start,
                            o.cds_end,
                            o.strand,
                        )

                    represented = {location_key(o) for o in rows}
                    rows += [o for o in additional if location_key(o) not in represented]
                if not rows:
                    hit.mapping_status = "no_mapping"
                    hit.mapping_error = NO_MAPPING_MESSAGE
                else:
                    hit.mapping_status = "success"
                hit.mapping_resolution_version = 3
                # Replace an incomplete prior expansion without biological dereplication.
                run.occurrences = [o for o in run.occurrences if o.hit_rank != hit.rank] + rows
            except MappingUnresolved as exc:
                hit.mapping_status = "unresolved_mapping"
                hit.mapping_error = self.client.http.safe_message(str(exc))
            except ResponseParseError as exc:
                hit.mapping_status = "parse_failure"
                hit.mapping_error = self.client.http.safe_message(str(exc))
            except ITEM_ERRORS as exc:
                hit.mapping_status = "retrieval_failure"
                hit.mapping_error = self.client.http.safe_message(str(exc))
            self.checkpoint()
        run.occurrences.sort(key=lambda o: (o.hit_rank, o.ipg_row))
        self.checkpoint(
            f"{len(run.hits)} BLAST hits expanded to {len(run.occurrences)} genomic occurrences"
        )

    def retrieve_contexts(self) -> None:
        run = self.run
        run.stage = "contexts"
        self.checkpoint("Batching nucleotide record metadata")
        self.client.prefetch_nucleotide_info(
            [o.nucleotide_accession for o in run.occurrences if o.status != "unresolved_mapping"]
        )
        pending, whole, unindexed = {}, set(), set()
        for index, o in enumerate(run.occurrences, 1):
            check_cancel(self.client.http.cancel)
            if o.status == "unresolved_mapping":
                continue
            if o.status == "success":
                path = self.directory / o.context_file
                if (
                    o.target_cds_found is True
                    and path.is_file()
                    and digest(path.read_bytes()) == o.context_sha256
                ):
                    continue
            o.status, o.error = "running", ""
            self.checkpoint(
                f"Retrieving context · {index}/{len(run.occurrences)} · {o.nucleotide_accession}"
            )
            try:
                info = self.client.nucleotide_info(o.nucleotide_accession, allow_unindexed=True)
                o.record_accession, o.record_length = info["accession"], info["length"]
                if o.cds_start < 1 or (o.record_length is not None and o.cds_end > o.record_length):
                    o.status = "unresolved_mapping"
                    o.error = "Mapped CDS coordinates lie outside the current nucleotide record."
                    self.checkpoint()
                    continue
                o.context_start = max(1, o.requested_start)
                o.context_end = (
                    min(o.record_length, o.requested_end) if o.record_length else o.requested_end
                )
                o.clipped_left = o.context_start != o.requested_start
                o.clipped_right = o.context_end != o.requested_end
                key = o.record_accession, o.context_start, o.context_end
                pending.setdefault(key, []).append(o)
                if o.context_start == 1 and o.context_end == o.record_length:
                    whole.add(key)
                if o.record_length is None:
                    unindexed.add(key)
            except ITEM_ERRORS as exc:
                o.status = "retrieval_failure"
                o.error = self.client.http.safe_message(str(exc))
            self.checkpoint()
        self.checkpoint(f"Retrieving {len(pending)} distinct context downloads")
        for key, result in self.client.context_results(list(pending), whole, unindexed):
            for o in pending[key]:
                # Apply the first completed response even if cancellation arrived during I/O.
                # The following occurrence observes cancellation, preserving completed work.
                try:
                    if isinstance(result, Exception):
                        raise result
                    with self.timed("local_processing"):
                        self.save_context(o, *result)
                except ITEM_ERRORS as exc:
                    o.status = "retrieval_failure"
                    o.error = self.client.http.safe_message(str(exc))
                self.checkpoint()
                check_cancel(self.client.http.cancel)

    def save_context(self, o, text, record, retrieved_at):
        actual_end = o.context_start + len(record) - 1
        if o.record_length is None:
            o.record_accession = record.id.split(":", 1)[0]
            info = self.client.nucleotide_info(o.nucleotide_accession, allow_unindexed=True)
            o.record_length = info["length"]
            o.context_end = actual_end
            o.clipped_right = o.context_end != o.requested_end
        if o.cds_end > actual_end:
            o.status = "unresolved_mapping"
            o.error = "Mapped CDS coordinates lie outside the retrieved nucleotide record."
            return
        o.feature_count = len(record.features)
        o.target_cds_found = target_cds_present(record, o)
        o.context_file = f"contexts/{o.id}.gb"
        atomic_write(self.directory / o.context_file, text)
        o.context_sha256, o.retrieved_at = digest(text), retrieved_at
        if o.target_cds_found:
            o.status = "success"
        else:
            o.status, o.error = "integrity_failure", TARGET_CDS_ERROR
