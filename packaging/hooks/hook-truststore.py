"""truststore's optional `inject_into_ssl` helpers import requests and urllib3 inside functions.

The apps never call them: `kbiosuite.tls` passes a truststore context straight to urllib.
Without this, PyInstaller would follow those imports and bundle requests, urllib3 and certifi
into Kannotate and Kmatch. Kblast imports requests itself, so it keeps its own copy.
"""

excludedimports = ["requests", "urllib3"]
