"""One JSON request on stdin; one JSON response on stdout. Diagnostics on stderr."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from . import __version__
from .model import create_project, recompute, validate
from .parser import inspect_inputs, parse_genbank
from .persistence import load, to_html


def handle(request: dict) -> dict:
    action = request["action"]
    if action == "inspect":
        return {"inputs": inspect_inputs(request["paths"])}
    if action == "generate":
        start = time.perf_counter()
        if "clusters" in request:
            clusters = request["clusters"]
            if any(c.get("topology") == "circular" and not c.get("linearized_at_origin") for c in clusters):
                raise ValueError("Circular input requires explicit Linearize at Origin or Omit.")
        else:
            inputs = request["inputs"]
            if not 1 <= len(inputs) <= 50:
                raise ValueError("Select 1–50 inputs.")
            clusters = []
            for i, item in enumerate(inputs):
                cluster = parse_genbank(item["path"], item.get("name"), slot=i)
                if cluster["topology"] == "circular":
                    if item.get("circular") == "omit":
                        continue
                    if item.get("circular") != "linearize":
                        raise ValueError("Circular input requires Linearize at Origin or Omit.")
                    cluster["linearized_at_origin"] = True
                clusters.append(cluster)
        project = create_project(clusters, request.get("settings"), request.get("preserve", False), request.get("connectors", True))
        if "title" in request:
            if not isinstance(request["title"], str):
                raise ValueError("Project title must be a string.")
            project["metadata"]["title"] = request["title"]
        html = to_html(project)
        return {"project": project, "html": html, "elapsed_seconds": time.perf_counter() - start}
    if action == "recompute":
        project = recompute(request["project"], request["settings"], request.get("relayout", False))
        return {"project": project}
    if action == "open":
        return {"project": load(request["path"])}
    if action == "serialize":
        return {"html": to_html(request["project"])}
    if action == "validate":
        return {"project": validate(request["project"])}
    raise ValueError(f"Unknown backend action: {action}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Kstack JSON biological backend")
    parser.add_argument("--request", type=Path, help="Read request JSON from file instead of stdin")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()
    try:
        request = json.loads(args.request.read_text() if args.request else sys.stdin.read())
        response = {"ok": True, "result": handle(request)}
    except Exception as exc:
        response = {"ok": False, "error": {"type": type(exc).__name__, "message": str(exc)}}
    print(json.dumps(response, ensure_ascii=False, allow_nan=False))
    return 0 if response["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
