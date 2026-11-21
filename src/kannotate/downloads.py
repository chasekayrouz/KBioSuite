"""Checksummed downloads with bounded parallel HTTP ranges for large assets."""

import hashlib
import os
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.error import HTTPError

from kbiosuite import tls


def digest(path: Path, algorithm="sha256") -> str:
    value = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def fetch(url: str, path: Path, *, sha256=None, md5=None, workers=8):
    def valid():
        return path.is_file() and (
            sha256 and digest(path) == sha256 or not sha256 and md5 and digest(path, "md5") == md5
        )

    if valid():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    request = urllib.request.Request(url, method="HEAD")
    try:
        with tls.urlopen(request, timeout=60) as response:
            size = int(response.headers.get("Content-Length", "0"))
            ranges = response.headers.get("Accept-Ranges") == "bytes"
    except HTTPError as exc:
        if exc.code not in (405, 501):
            raise
        size, ranges = 0, False

    def part(start, end):
        for retry in range(3):
            try:
                req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
                with (
                    tls.urlopen(req, timeout=120) as response,
                    temporary.open("r+b") as out,
                ):
                    if response.status != 206:
                        raise RuntimeError("Server did not honor HTTP Range request")
                    out.seek(start)
                    count = 0
                    while block := response.read(1024 * 1024):
                        count += len(block)
                        if count > end - start + 1:
                            raise RuntimeError("Oversized range response")
                        out.write(block)
                    if count != end - start + 1:
                        raise RuntimeError("Incomplete range response")
                return
            except (OSError, RuntimeError):
                if retry == 2:
                    raise
                time.sleep(1 + retry)

    if size > 32 * 1024 * 1024 and ranges and workers > 1:
        with temporary.open("wb") as handle:
            handle.truncate(size)
        chunk = (size + workers - 1) // workers
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(part, start, min(start + chunk - 1, size - 1))
                for start in range(0, size, chunk)
            ]
            for future in futures:
                future.result()
    else:
        with tls.urlopen(url, timeout=120) as response, temporary.open("wb") as handle:
            while block := response.read(1024 * 1024):
                handle.write(block)
    if sha256 and digest(temporary) != sha256 or md5 and digest(temporary, "md5") != md5:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"Checksum mismatch for {url}; refusing changed database content")
    os.replace(temporary, path)
    return path
