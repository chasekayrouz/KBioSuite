"""Explicit frozen-runtime qualification; ordinary application launch does not call it.

Only the stdlib `urllib` stack lives here. Every frozen Python app imports this module and
PyInstaller bundles whatever it can see imported, function-level imports included, so a
third-party stack is passed in by the one entry point that ships it: Kblast's `requests`
probe is `kbiosuite.requests_probe`, imported only by packaging/Kblast.py.
"""

import json
import os
import ssl
import sys
import urllib.request
from pathlib import Path

NETWORK_FLAG = "--kbiosuite-network-check"
PUBLIC_URL = "https://www.ncbi.nlm.nih.gov/"
CERTIFICATE_OVERRIDES = ("SSL_CERT_FILE", "SSL_CERT_DIR", "REQUESTS_CA_BUNDLE", "CURL_CA_BUNDLE")


def urllib_head(url):
    """HEAD through stdlib urllib (Kannotate downloads, Kmatch RCSB access); no proxies."""
    context = ssl.create_default_context()
    if context.verify_mode != ssl.CERT_REQUIRED or not context.check_hostname:
        raise RuntimeError("TLS certificate and hostname verification must be enabled")
    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({}), urllib.request.HTTPSHandler(context=context)
    )
    with opener.open(urllib.request.Request(url, method="HEAD"), timeout=20) as response:
        return response.status


def network_check(stack, head=None):
    """One public HEAD request, certificate verification enabled and proxies disabled.

    head(url) -> HTTP status performs the request for a stack other than stdlib urllib.
    """
    if any(name in os.environ for name in CERTIFICATE_OVERRIDES):
        raise RuntimeError("Remove certificate override variables before qualification")
    if head is None:
        if stack != "urllib":
            raise ValueError("Unsupported qualification network stack")
        head = urllib_head
    status = head(PUBLIC_URL)
    if status != 200:
        raise RuntimeError(f"Qualification endpoint returned HTTP {status}")
    return {"stack": stack, "url": PUBLIC_URL, "http_status": status, "tls_verified": True}


def dispatch_network_check(tool, stack, argv=None, head=None):
    """Return None for a normal launch, or a process status for the explicit probe."""
    arguments = list(sys.argv[1:] if argv is None else argv)
    if not arguments or arguments[0] != NETWORK_FLAG:
        return None
    if len(arguments) != 2:
        raise SystemExit(f"Usage: {NETWORK_FLAG} NEW_REPORT_PATH")
    report = Path(arguments[1])
    # Reserve the report before any network access; never overwrite earlier evidence.
    with report.open("x") as output:
        try:
            result = {"tool": tool, "passed": True, **network_check(stack, head)}
        except Exception as exc:
            result = {"tool": tool, "stack": stack, "passed": False,
                      "error": str(exc), "error_type": type(exc).__name__}
        json.dump(result, output, indent=2)
        output.write("\n")
    return 0 if result["passed"] else 1
