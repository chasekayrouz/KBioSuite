"""Live mapping/context validation independent of the potentially long BLAST queue."""

import json
import argparse
import uuid
import threading
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace

from kblast.config import Config
from kblast.models import TARGET_CDS_ERROR, digest, utc_now
from kblast.ncbi import NcbiClient
from kblast.network import RateLimiter, Transport
from kblast.parsers import parse_ipg
from kblast.pipeline import target_cds_present
from kblast.storage import atomic_write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=Path("runs/kblast") / ("mapping-" + uuid.uuid4().hex))
    parser.add_argument("--limit", type=int, help="Explicitly test only the first N occurrences")
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    directory = args.output_dir
    directory.mkdir(parents=True, exist_ok=True)
    config = Config.discover()
    transport = Transport(
        config,
        threading.Event(),
        directory / "requests.jsonl",
        RateLimiter(0.11 if config.api_key else 0.36),
        RateLimiter(10),
    )
    client = NcbiClient(transport, directory)
    accession = "CAC44203.1"
    text = client.fetch_ipg(accession, 1)
    rows = parse_ipg(text, SimpleNamespace(rank=1, accession=accession), 20_000)
    if args.limit:
        rows = rows[:args.limit]
    atomic_write(directory / "raw" / "hit_0001.ipg.tsv", text)
    print(f"Live mapping check: {accession}, {len(rows)} genomic occurrences", flush=True)
    for index, o in enumerate(rows, 1):
        if o.status == "unresolved_mapping":
            continue
        try:
            info = client.nucleotide_info(o.nucleotide_accession)
            o.record_accession, o.record_length = info["accession"], info["length"]
            if o.cds_end > o.record_length:
                o.status, o.error = (
                    "unresolved_mapping",
                    "CDS coordinates exceed the nucleotide record",
                )
                continue
            o.context_start, o.context_end = (
                max(1, o.requested_start),
                min(o.record_length, o.requested_end),
            )
            o.clipped_left = o.context_start != o.requested_start
            o.clipped_right = o.context_end != o.requested_end
            text, record, stamp = client.fetch_context(
                o.record_accession, o.context_start, o.context_end
            )
            o.context_file = f"contexts/{o.id}.gb"
            atomic_write(directory / o.context_file, text)
            o.context_sha256, o.retrieved_at = digest(text), stamp
            o.feature_count = len(record.features)
            o.target_cds_found = target_cds_present(record, o)
            if o.target_cds_found:
                o.status = "success"
            else:
                o.status, o.error = "integrity_failure", TARGET_CDS_ERROR
        except Exception as exc:
            o.status, o.error = "retrieval_failure", transport.safe_message(str(exc))
        print(f"{index}/{len(rows)} {o.nucleotide_accession}: {o.status}", flush=True)
        atomic_write(
            directory / "validation.json",
            json.dumps(
                {
                    "scope": "Public-protein IPG and context retrieval, independently of BLAST",
                    "protein_accession": accession,
                    "flank": 20_000,
                    "at": utc_now(),
                    "occurrences": [asdict(row) for row in rows],
                },
                indent=2,
            ),
        )
    transport.close()
    print(
        {
            "rows": len(rows),
            "success": sum(o.status == "success" for o in rows),
            "unresolved": sum(o.status == "unresolved_mapping" for o in rows),
            "failed": sum(o.status == "retrieval_failure" for o in rows),
            "integrity_failures": sum(o.status == "integrity_failure" for o in rows),
        },
        flush=True,
    )


if __name__ == "__main__":
    main()
