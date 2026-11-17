"""Explicit live NCBI validation; inputs and raw stores remain local."""

import argparse
import csv
import io
import json
import threading
import uuid
import zipfile
from pathlib import Path

from Bio import Blast, SeqIO

from kblast.app import InstanceLock
from kblast.config import Config
from kblast.export import context_filename, export_bundle
from kblast.models import Run, Settings
from kblast.ncbi import NcbiClient
from kblast.network import RateLimiter, Transport
from kblast.pipeline import Pipeline, target_cds_present
from kblast.sequence import parse_query
from kblast.storage import RunStore, atomic_write


def audit_outputs(store, run):
    directory = store.directory(run.id)
    bundle = export_bundle(directory, run)
    atomic_write(directory / "validation_export.zip", bundle)
    successful = [o for o in run.occurrences if o.status == "success"]
    proteins = [h for h in run.hits if h.sequence_status == "success"]
    raw_ipg_rows = 0
    target_translation_checks = 0
    target_translation_checks_skipped = 0
    with zipfile.ZipFile(io.BytesIO(bundle)) as archive:
        combined = list(SeqIO.parse(io.StringIO(archive.read("contexts.gb").decode()), "genbank"))
        assert len(combined) == len(successful)
        for record, o in zip(combined, successful, strict=True):
            raw = (directory / o.context_file).read_bytes()
            assert archive.read(context_filename(o)) == raw
            source = SeqIO.read(io.StringIO(raw.decode()), "genbank")
            assert record.seq == source.seq
            assert [str(f.location) for f in record.features] == [
                str(f.location) for f in source.features
            ]
            assert o.target_cds_found is True and target_cds_present(record, o)
            hit = next(h for h in run.hits if h.rank == o.hit_rank)
            targets = [
                f
                for f in record.features
                if f.type == "CDS"
                and f.location is not None
                and int(f.location.start) + o.context_start == o.cds_start
                and int(f.location.end) + o.context_start - 1 == o.cds_end
                and f.location.strand == {"+": 1, "-": -1}[o.strand]
            ]
            if hit.sequence_status == "success":
                assert any(f.qualifiers.get("translation") == [hit.sequence] for f in targets)
                target_translation_checks += 1
            else:
                target_translation_checks_skipped += 1
        fasta = list(SeqIO.parse(io.StringIO(archive.read("blast_hits.fasta").decode()), "fasta"))
        assert [str(r.seq) for r in fasta] == [h.sequence for h in proteins]
        rows = list(
            csv.DictReader(io.StringIO(archive.read("metadata.tsv").decode()), delimiter="\t")
        )
        assert len(rows) == len(run.hits) + len(run.occurrences)
        assert sum(r["context_exported"] == "True" for r in rows) == len(successful)
        for hit in run.hits:
            path = directory / "raw" / f"hit_{hit.rank:04d}.ipg.tsv"
            if not path.is_file():
                assert hit.mapping_status in {"retrieval_failure", "parse_failure", "cancelled"}
                continue
            raw = path.read_bytes()
            assert archive.read(f"provenance/{path.name}") == raw
            raw_rows = list(csv.DictReader(io.StringIO(raw.decode()), delimiter="\t"))
            raw_ipg_rows += len(raw_rows)
            for o in (o for o in run.occurrences if o.hit_rank == hit.rank):
                if o.mapping_method == "IPG":
                    assert o.raw_mapping == raw_rows[o.ipg_row - 1]
                else:
                    assert o.ipg_row == 0
                    assert archive.read(f"provenance/{o.raw_mapping['evidence_file']}")
    reference = Blast.read(io.BytesIO((directory / "raw" / "blast.xml").read_bytes()))
    assert run.blast_metadata["program"] == "blastp"
    assert run.blast_metadata["database"] == run.settings.database
    query_checks = 0
    assert len(run.hits) == min(run.settings.top_n, len(reference))
    for hit, ref in zip(run.hits, reference[: run.settings.top_n], strict=True):
        assert hit.subject_id == ref.target.id
        assert hit.bit_score == max(h.annotations["bit score"] for h in ref)
        assert hit.evalue == min(h.annotations["evalue"] for h in ref)
        assert hit.percent_identity == max(
            100 * h.annotations["identity"] / h.shape[1] for h in ref
        )
        positions = set()
        for hsp in ref:
            low, high = sorted(map(int, (hsp.coordinates[1, 0], hsp.coordinates[1, -1])))
            positions.update(range(low, high))
        assert hit.query_coverage == 100 * len(positions) / len(run.query.sequence)
        assert len(hit.hsps) == len(ref)
        for hsp in hit.hsps:
            assert (
                hsp["qseq"].replace("-", "").upper()
                == run.query.sequence[hsp["query-from"] - 1 : hsp["query-to"]]
            )
            query_checks += 1
    from datetime import datetime

    events = [json.loads(line) for line in (directory / "requests.jsonl").read_text().splitlines()]
    submissions = [
        e
        for e in events
        if e.get("parameters", {}).get("CMD") == "Put" and e["outcome"] == "started"
    ]
    assert len(submissions) == 1
    sent = [
        e
        for e in events
        if e.get("parameters", {}).get("CMD") == "Put" and e["outcome"] == "received"
    ]
    assert sent[0]["sent_parameters"]["PROGRAM"] == "blastp"
    assert sent[0]["sent_parameters"]["DATABASE"] == run.settings.database
    assert int(sent[0]["sent_parameters"]["HITLIST_SIZE"]) == run.settings.top_n
    assert sent[0]["query_sha256"] == run.to_dict()["query_sha256"]
    assert sent[0]["api_key_used"] is False
    eutils = [e for e in events if e.get("service") == "eutils" and e["outcome"] == "received"]
    assert run.status in {"success", "partial_success", "no_hits"}
    audit = {
        "passed": True,
        "run_id": run.id,
        "rid": run.rid,
        "program": run.blast_metadata["program"],
        "database": run.blast_metadata["database"],
        "blast_version": run.blast_metadata["version"],
        "submitted_parameters": run.submitted_blast_parameters,
        "returned_parameters": run.blast_metadata["parameters"],
        "blast_elapsed_seconds": (
            datetime.fromisoformat(run.blast_completed_at)
            - datetime.fromisoformat(run.submitted_at)
        ).total_seconds(),
        "submission_count": len(submissions),
        "query_alignment_checks": query_checks,
        "raw_xml_subjects": len(reference),
        "hits": len(run.hits),
        "raw_ipg_rows": raw_ipg_rows,
        "genomic_occurrences": len(run.occurrences),
        "non_genomic_ipg_rows_excluded": raw_ipg_rows
        - sum(o.mapping_method == "IPG" for o in run.occurrences),
        "linked_cds_occurrences_added": sum(
            o.mapping_method != "IPG" for o in run.occurrences
        ),
        "contexts": len(successful),
        "protein_sequences": len(proteins),
        "unresolved_rows": sum(o.status == "unresolved_mapping" for o in run.occurrences),
        "retrieval_failures": sum(o.status == "retrieval_failure" for o in run.occurrences),
        "integrity_failures": sum(o.status == "integrity_failure" for o in run.occurrences),
        "target_cds_translation_checks": target_translation_checks,
        "target_cds_translation_checks_skipped": target_translation_checks_skipped,
        "run_status": run.status,
        "top_n": run.settings.top_n,
        "flank": run.settings.flank,
        "stage_seconds": run.timings,
        "per_hit": [
            {
                "rank": h.rank,
                "accession": h.accession,
                "mapping_status": h.mapping_status,
                "mapping_error": h.mapping_error,
                "occurrences": sum(o.hit_rank == h.rank for o in run.occurrences),
                "contexts": sum(
                    o.hit_rank == h.rank and o.status == "success" for o in run.occurrences
                ),
                "context_failures": [
                    {"id": o.id, "status": o.status, "error": o.error}
                    for o in run.occurrences
                    if o.hit_rank == h.rank and o.status != "success"
                ],
            }
            for h in run.hits
        ],
        "eutils_attempts": sum(
            e["service"] == "eutils" and e["outcome"] == "started" for e in events if "service" in e
        ),
        "eutils_key_used_on_all_received_requests": bool(eutils)
        and all(e.get("api_key_used") for e in eutils),
        "checks": "Source-byte exports, GenBank orientation/features, FASTA, all audit rows, "
        "and BLAST subject order/statistics against Biopython",
    }
    atomic_write(directory / "audit.json", json.dumps(audit, indent=2))
    return audit


def record_completion(run, audit):
    path = Path("docs/kblast/validation.md")
    if not path.is_file():
        return
    start, end = "<!-- LIVE_BLAST_BEGIN -->", "<!-- LIVE_BLAST_END -->"
    text = path.read_text()
    if start not in text or end not in text:
        return
    result = (
        f"{start}\nThe live {run.settings.database} run `{run.id}` "
        f"finished with status **{run.status}**, "
        f"RID `{run.rid}`.\n\n"
    )
    if audit and audit.get("passed"):
        result += (
            f"The completed export audit passed: {audit['hits']} subject hits, "
            f"{audit['genomic_occurrences']} genomic occurrences, "
            f"{audit['contexts']} GenBank contexts, "
            f"and {audit['protein_sequences']} full protein sequences. "
            f"{audit['non_genomic_ipg_rows_excluded']} non-genomic rows are preserved only "
            "in raw IPG provenance. "
            f"BLAST completed in {audit['blast_elapsed_seconds']:.1f} seconds. "
            f"Returned metadata: {audit['program']}, {audit['database']}, "
            f"{audit['blast_version']}. Target CDS/translation checks: "
            f"{audit['target_cds_translation_checks']}. "
            "See the run’s `audit.json` and `validation_export.zip`.\n"
        )
    else:
        result += "The live end-to-end check did not complete; inspect the saved error/status.\n"
    atomic_write(path, text.split(start)[0] + result + end + text.split(end, 1)[1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume")
    parser.add_argument("--data-dir", type=Path, default=Path("runs/nr_validation"))
    parser.add_argument("--query", type=Path, default=Path("examples/kblast/actinorhodin_cyclase.fasta"))
    parser.add_argument("--database", choices=["nr", "nr_cluster_seq"], default="nr")
    parser.add_argument("--top-n", type=int, default=2)
    parser.add_argument("--flank", type=int, default=20_000)
    parser.add_argument("--update-docs", action="store_true")
    args = parser.parse_args()
    config = Config.discover(args.data_dir)
    lock = InstanceLock(config.data_dir)
    store = RunStore(config.data_dir)
    run = (
        store.load(args.resume)
        if args.resume
        else Run(
            uuid.uuid4().hex,
            parse_query(args.query.read_text()),
            Settings(top_n=args.top_n, flank=args.flank, database=args.database),
        )
    )
    store.save(run)
    print(f"Live validation run: {run.id}", flush=True)
    transport = Transport(
        config,
        threading.Event(),
        store.directory(run.id) / "requests.jsonl",
        RateLimiter(0.11 if config.api_key else 0.36),
        RateLimiter(10),
    )

    class ReportingPipeline(Pipeline):
        def checkpoint(self, message=None):
            super().checkpoint(message)
            if message and not message.startswith("Retrieving context ·"):
                print(message, flush=True)

    ReportingPipeline(store, run, NcbiClient(transport, store.directory(run.id))).execute()
    print(f"Final status: {run.status}", flush=True)
    if run.error:
        print(run.error, flush=True)
    try:
        audit = audit_outputs(store, run) if run.blast_metadata else None
    except Exception as exc:
        audit = {"passed": False, "error": transport.safe_message(f"{type(exc).__name__}: {exc}")}
        atomic_write(store.directory(run.id) / "audit.json", json.dumps(audit, indent=2))
        print("The live export audit failed; inspect audit.json.", flush=True)
    if args.update_docs:
        record_completion(run, audit)
    lock.close()
    if not audit or not audit.get("passed"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
