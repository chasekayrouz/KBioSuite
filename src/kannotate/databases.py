import gzip
import json
import shutil
import subprocess
import sys
import tarfile
import tempfile
from dataclasses import replace
from datetime import UTC, datetime
from importlib.resources import files
from pathlib import Path

from . import __version__
from .config import default_root, load_config, save_config
from .dependencies import REQUIRED, executable_version, resolve_dependencies
from .downloads import digest, fetch
from .io.fasta import read_fasta
from .search.common import blast_database_argument, local_environment, write_queries
from .search.hmmer import HmmerAdapter
from .search.rpsblast import RpsBlastAdapter


def specification():
    return json.loads(files("kannotate.resources").joinpath("databases.json").read_text())


def disk_usage(directory: Path) -> int:
    return sum(p.stat().st_size for p in directory.rglob("*") if p.is_file())


def safe_extract(archive: Path, destination: Path):
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            relative = Path(member.name)
            if (
                relative.is_absolute()
                or ".." in relative.parts
                or not (member.isfile() or member.isdir())
            ):
                raise RuntimeError("Unsafe member in database archive")
        tar.extractall(destination, filter="data")


def verify_databases(config: dict):
    tools = config["executables"]
    db = config["databases"]
    info = subprocess.run(
        [tools["blastdbcmd"], "-db", blast_database_argument(db["cdd"]["path"]), "-info"],
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
        env=local_environment(),
    )
    if info.returncode:
        raise RuntimeError("CDD validation failed: " + info.stderr.strip())
    expected = specification()["cdd"]["expected_profiles"]
    if f"{expected:,} sequences" not in info.stdout or "CDD.v3.21" not in info.stdout:
        raise RuntimeError("CDD profile count differs from pinned comprehensive release")
    with tempfile.TemporaryDirectory(prefix="kannotate-smoke-") as temp:
        directory = Path(temp)
        query = directory / "smoke.faa"
        query.write_text(files("kannotate.resources").joinpath("smoke.faa").read_text())
        proteins = read_fasta(query)
        write_queries([replace(p, id="smoke") for p in proteins], query)
        hmm = HmmerAdapter(tools["hmmscan"], db["pfam"]["path"]).search(query, directory, 2)
        rps = RpsBlastAdapter(
            tools["rpsblast"], db["cdd"]["path"], db["cdd"]["descriptions"]
        ).search(query, directory, 2)
    if not any(h.evalue < 1e-10 for h in hmm) or not any(h.evalue < 1e-10 for h in rps):
        raise RuntimeError("Positive-control protein did not match both databases")
    return {
        "blastdb_info": info.stdout.strip(),
        "hmmer_hits": len(hmm),
        "rps_hits": len(rps),
        "control": "UniProt P0CE47 (E. coli elongation factor Tu)",
        "passed": True,
    }


def setup(root: Path | None = None, *, report=print):
    root = (root or default_root()).resolve()
    root.mkdir(parents=True, exist_ok=True)
    spec = specification()
    try:
        previous = load_config(root)
    except (RuntimeError, OSError, ValueError):
        config_path = root / ".kannotate/config.json"
        if config_path.is_file():
            # Never rebuild legacy copies when an explicitly shared pin is broken.
            raise
        previous = {}
    if previous.get("shared_resources"):
        from kbiosuite.resource_client import consumer
        from kbiosuite.resources import verify

        pins = consumer(previous["shared_resources"])["bindings"]
        for identity in sorted({pin["resource"] for pin in pins.values()}):
            verify(identity)
        validation = verify_databases(previous)
        report("Shared resources verified; both local search controls passed.")
        return dict(previous, validation=validation)
    tools = resolve_dependencies(previous.get("executables"))
    missing = [name for name in REQUIRED if not tools[name]]
    if missing:
        raise RuntimeError(
            "Missing local tools: "
            + ", ".join(missing)
            + ". Install HMMER/BLAST+ or place existing executables on PATH."
        )
    if sys.byteorder != "little":
        raise RuntimeError("The pinned CDD distribution requires a little-endian host")
    storage = root / "data" / "databases"
    storage.mkdir(parents=True, exist_ok=True)
    archives = storage / ".downloads"
    archives.mkdir(exist_ok=True)
    databases = {}
    for key in ("pfam", "cdd"):
        entry = spec[key]
        target = storage / f"{key}-{entry['version']}"
        marker = target / "installation.json"
        if marker.is_file():
            metadata = json.loads(marker.read_text())
            if metadata.get("specification") != entry:
                raise RuntimeError(f"Installed {key} metadata differs from pinned specification")
            for name, checksum in metadata["files"].items():
                if not (target / name).is_file() or digest(target / name) != checksum:
                    raise RuntimeError(f"Installed {key} file is missing or changed: {name}")
            report(f"{key}: verified existing {entry['version']} installation")
        else:
            if shutil.disk_usage(storage).free < 12 * 1024**3:
                raise RuntimeError(
                    "Setup needs at least 12 GiB free for temporary database preparation"
                )
            if target.exists():
                raise RuntimeError(f"Incomplete installation at {target}; inspect before removing")
            report(f"Installing {key} {entry['version']}…")
            assets = []
            for asset in entry["assets"]:
                archive = archives / asset["filename"]
                fetch(asset["url"], archive, sha256=asset["sha256"], md5=asset.get("md5"))
                assets.append(archive)
            with tempfile.TemporaryDirectory(prefix=f".{key}-", dir=storage) as temporary:
                staging = Path(temporary)
                if key == "pfam":
                    with (
                        gzip.open(assets[0], "rb") as source,
                        (staging / "Pfam-A.hmm").open("wb") as out,
                    ):
                        shutil.copyfileobj(source, out)
                    count = 0
                    with (staging / "Pfam-A.hmm").open() as handle:
                        count = sum(line.startswith("HMMER3/") for line in handle)
                    if count != entry["expected_profiles"]:
                        raise RuntimeError("Pfam profile count differs from pinned release")
                    run = subprocess.run(
                        [tools["hmmpress"], str(staging / "Pfam-A.hmm")],
                        capture_output=True,
                        text=True,
                        check=False,
                    )
                    if run.returncode:
                        raise RuntimeError("hmmpress failed: " + run.stderr.strip())
                else:
                    safe_extract(assets[0], staging)
                    with (
                        gzip.open(assets[1], "rb") as source,
                        (staging / "cddid.tbl").open("wb") as out,
                    ):
                        shutil.copyfileobj(source, out)
                file_hashes = {
                    str(p.relative_to(staging)): digest(p)
                    for p in sorted(staging.rglob("*"))
                    if p.is_file()
                }
                metadata = {
                    "schema_version": 1,
                    "kannotate_version": __version__,
                    "installed_at": datetime.now(UTC).isoformat(),
                    "specification": entry,
                    "files": file_hashes,
                    "disk_bytes": disk_usage(staging),
                    "executables": {
                        n: {"path": p, "version": executable_version(p)}
                        for n, p in tools.items()
                        if p
                    },
                }
                (staging / "installation.json").write_text(json.dumps(metadata, indent=2) + "\n")
                staging.rename(target)
        databases[key] = {
            "version": entry["version"],
            "path": str(target / entry["database_prefix"]),
        }
        if key == "cdd":
            databases[key]["descriptions"] = str(target / "cddid.tbl")
    config = {
        "schema_version": 1,
        "kannotate_version": __version__,
        "executables": tools,
        "databases": databases,
    }
    report("Validating both databases with a public positive-control protein…")
    validation = verify_databases(config)
    config["validation"] = validation
    save_config(config, root)
    # Remove only the exact downloaded database assets after both controls pass.
    for key in ("pfam", "cdd"):
        for asset in spec[key]["assets"]:
            (archives / asset["filename"]).unlink(missing_ok=True)
    used = disk_usage(storage)
    report(
        f"Ready: Pfam {spec['pfam']['version']}; CDD {spec['cdd']['version']}; "
        f"{used / 1024**3:.2f} GiB"
    )
    return config


def status(root: Path | None = None):
    root = root or default_root()
    try:
        config = load_config(root)
    except (RuntimeError, OSError, ValueError):
        config = {}
    tools = resolve_dependencies(config.get("executables"))
    db_status = {}
    for key, entry in config.get("databases", {}).items():
        prefix = Path(entry["path"])
        marker = prefix.parent / "installation.json"
        try:
            metadata = json.loads(marker.read_text())
            required = [prefix.parent / name for name in metadata["files"]]
            matches = metadata["specification"] == specification()[key]
        except (OSError, ValueError, KeyError):
            required, matches = [], False
        db_status[key] = {
            "version": entry["version"],
            "ready": matches and bool(required) and all(p.is_file() for p in required),
        }
    return {
        "executables": tools,
        "databases": db_status,
        "disk_bytes": sum(
            disk_usage(directory) for directory in
            {Path(entry["path"]).parent for entry in config.get("databases", {}).values()}
        ),
    }
