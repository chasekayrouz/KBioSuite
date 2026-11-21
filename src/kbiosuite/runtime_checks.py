"""Explicit frozen-runtime qualification; ordinary application launch does not call it.

No network stack lives here. Every frozen Python app imports this module and PyInstaller
bundles whatever it can see imported, function-level imports included, so each entry point
passes the probes for the client it ships: Kmatch and Kannotate use `kbiosuite.tls` (stdlib
urllib with macOS trust) against their own hosts, and Kblast's `requests` probe is
`kbiosuite.requests_probe`, imported only by packaging/Kblast.py.
"""

import json
import os
import sys
from pathlib import Path

NETWORK_FLAG = "--kbiosuite-network-check"
PUBLIC_URL = "https://www.ncbi.nlm.nih.gov/"
CERTIFICATE_OVERRIDES = ("SSL_CERT_FILE", "SSL_CERT_DIR", "REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE")


def network_check(stack, probes):
    """One verified public request per endpoint, each recorded whether or not it passed.

    probes maps a URL to probe(url) -> HTTP status; a probe verifies certificates, ignores
    proxy settings and raises when the endpoint does not answer as the app needs.
    """
    if any(name in os.environ for name in CERTIFICATE_OVERRIDES):
        raise RuntimeError("Remove certificate override variables before qualification")
    endpoints = []
    for url, probe in probes.items():
        try:
            status = probe(url)
            if status != 200:
                raise RuntimeError(f"Qualification endpoint returned HTTP {status}")
            endpoints.append({"url": url, "passed": True, "http_status": status})
        except Exception as exc:
            endpoints.append({"url": url, "passed": False,
                              "error": str(exc), "error_type": type(exc).__name__})
    passed = bool(endpoints) and all(e["passed"] for e in endpoints)
    return {"stack": stack, "endpoints": endpoints, "tls_verified": passed}


def dispatch_network_check(tool, stack, probes, argv=None):
    """Return None for a normal launch, or a process status for the explicit probe.

    probes() -> the tool's {url: probe}; it is called only for the explicit probe.
    """
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or arguments[0] != NETWORK_FLAG:
        return None
    if len(arguments) != 2:
        raise SystemExit(f"Usage: {NETWORK_FLAG} NEW_REPORT_PATH")
    report = Path(arguments[1])
    # Reserve the report before any network access; never overwrite earlier evidence.
    with report.open("x") as output:
        try:
            checked = network_check(stack, probes())
            result = {"tool": tool, "passed": checked["tls_verified"], **checked}
        except Exception as exc:
            result = {"tool": tool, "stack": stack, "passed": False,
                      "error": str(exc), "error_type": type(exc).__name__}
        json.dump(result, output, indent=2)
        output.write("\n")
    return 0 if result["passed"] else 1
