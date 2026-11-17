"""Local runtime identities; stat keys detect replacement between operations."""

import hashlib
import json
import struct
from functools import lru_cache
from pathlib import Path

from ..dependencies import executable_version


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def file_identity(path):
    path = Path(path).resolve()
    stat = path.stat()
    return [str(path), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns]


def digest(path):
    from ..cancellation import check_cancelled

    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            check_cancelled()
            value.update(block)
    return value.hexdigest()


@lru_cache(maxsize=32)
def tool_version(identity):
    return executable_version(identity[0])


@lru_cache(maxsize=64)
def content_digest(identity):
    return digest(identity[0])


def native_digest(path):
    """Hash all arm64 Mach-O sections, independent of appended code signatures.

    Signing changes the header's LC_CODE_SIGNATURE and __LINKEDIT. HMM computation
    lives in the sections; bind their names, flags and every byte, without copying
    or modifying the library. Unknown binary formats use their complete digest.
    """
    data = Path(path).read_bytes()
    if data[:4] != b"\xcf\xfa\xed\xfe":
        return hashlib.sha256(data).hexdigest()
    _, cpu, _, _, commands, _, _, _ = struct.unpack_from("<8I", data)
    if cpu != 0x100000C:
        raise ValueError("Resident qualification requires an arm64 Mach-O")
    value = hashlib.sha256()
    offset, sections = 32, 0
    for _ in range(commands):
        command, size = struct.unpack_from("<II", data, offset)
        if size < 8 or offset + size > len(data):
            raise ValueError("Invalid Mach-O command")
        if command == 0x19:  # LC_SEGMENT_64
            count = struct.unpack_from("<I", data, offset + 64)[0]
            for i in range(count):
                section = offset + 72 + i * 80
                name, segment, _, length, position, _, _, _, flags, _, _, _ = struct.unpack_from(
                    "<16s16sQQ8I", data, section
                )
                value.update(name + segment + struct.pack("<QI", length, flags))
                # Zerofill sections have no file contents.
                if flags & 0xFF not in {1, 0xC, 0x12}:
                    if position + length > len(data):
                        raise ValueError("Invalid Mach-O section")
                    value.update(data[position : position + length])
                sections += 1
        offset += size
    if not sections:
        raise ValueError("Mach-O has no qualified sections")
    return value.hexdigest()


def reference_identity(database):
    parent = Path(database["path"]).parent
    manifest = parent / "installation.json"
    metadata = json.loads(manifest.read_text())
    assets = [file_identity(parent / name) for name in sorted(metadata["files"])]
    if database.get("descriptions"):
        description = file_identity(database["descriptions"])
        if description not in assets:
            assets.append(description)
    return {
        "database": database,
        "installation_sha256": digest(manifest),
        "installation_identity": file_identity(manifest),
        "assets": assets,
    }
