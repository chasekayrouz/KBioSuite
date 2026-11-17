"""Kblast's `requests` qualification stack; imported only by packaging/Kblast.py.

Kept out of kbiosuite.runtime_checks so PyInstaller bundles requests, urllib3 and certifi
only into Kblast, whose NCBI client uses them; Kannotate and Kmatch use stdlib urllib.
"""

import requests


def head(url):
    """HEAD with certificate verification; no proxy, netrc or CA-bundle environment settings."""
    with requests.Session() as session:
        session.trust_env = False
        with session.head(url, timeout=20, verify=True, allow_redirects=False) as response:
            response.raise_for_status()
            return response.status_code
