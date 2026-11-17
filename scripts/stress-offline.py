"""Repeat deterministic workload/error lifecycles; no network or app launches."""

import argparse
import json
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iterations", type=int, default=5)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.iterations < 1:
        parser.error("iterations must be positive")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    suites = {
        "kblast": [
            "tests/kblast/test_blast_lifecycle.py",
            "tests/kblast/test_rid_import.py",
            "tests/kblast/test_retrieval.py",
            "tests/kblast/test_pipeline.py",
        ],
        "kannotate": [
            "tests/kannotate/test_gui.py",
            "tests/kannotate/test_desktop_api.py",
            "tests/kannotate/test_io.py",
        ],
        # Engines are faked: export/cancellation and database activation lifecycles only.
        "kmatch": [
            "tests/kmatch/test_kmatch_core.py",
            "tests/kmatch/test_alignment.py",
            "tests/kmatch/test_stabilization.py",
        ],
    }
    report = {
        "status": "passed",
        "layer": "offline stress",
        "started_at": datetime.now(UTC).isoformat(),
        "results": [],
    }
    for name, paths in suites.items():
        for iteration in range(args.iterations):
            stamp = f"{name}-{iteration + 1:02d}"
            started = time.perf_counter()
            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "pytest",
                    "-q",
                    *paths,
                    f"--junitxml={output / (stamp + '.xml')}",
                ],
                cwd=ROOT,
                capture_output=True,
                text=True,
            )
            (output / (stamp + ".log")).write_text(result.stdout + result.stderr)
            entry = {
                "component": name,
                "iteration": iteration + 1,
                "status": "passed" if result.returncode == 0 else "failed",
                "seconds": round(time.perf_counter() - started, 3),
                "exit_code": result.returncode,
            }
            if (output / (stamp + ".xml")).exists():
                suite = ET.parse(output / (stamp + ".xml")).getroot().find("testsuite")
                entry.update(
                    {
                        key: int(suite.get(key, "0"))
                        for key in ["tests", "failures", "errors", "skipped"]
                    }
                )
            report["results"].append(entry)
            if result.returncode:
                report["status"] = "failed"
            (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if report["status"] == "failed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
