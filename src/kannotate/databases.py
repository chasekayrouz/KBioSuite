import gzip
import json
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request
from dataclasses import replace
from datetime import UTC, datetime
from importlib.resources import files
from pathlib import Path

from . import __version__
from .config import cdd_release, default_root, load_config, save_config
from .dependencies import REQUIRED, executable_version, resolve_dependencies
from .downloads import digest, fetch
from .io.fasta import read_fasta
from .search.common import blast_database_argument, local_environment, write_queries
from .search.hmmer import HmmerAdapter
from .search.rpsblast import RpsBlastAdapter, load_descriptions


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


def read_release_info(url: str) -> dict:
    """Read NCBI's small cdd.info note naming the release it currently serves."""
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            text = response.read(64 * 1024).decode("utf-8", "replace")
    except OSError as exc:
        raise RuntimeError(f"Could not reach NCBI for CDD release information: {exc}") from exc
    # The release names the install directory, so accept only one plain path component.
    version = re.search(r"(?m)^cdd version (\w[\w.-]*)\s*$", text)
    profiles = re.search(r"(?m)^(\d+) models in the default search set", text)
    if not version or not profiles:
        raise RuntimeError("Could not read NCBI's CDD release information at " + url)
    return {"version": version[1], "profiles": int(profiles[1])}


def blastdbcmd(executable, prefix, *arguments):
    run = subprocess.run(
        [executable, "-db", blast_database_argument(prefix), *arguments],
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
        env=local_environment(),
    )
    if run.returncode:
        raise RuntimeError("CDD validation failed: " + run.stderr.strip())
    return run.stdout


def verify_cdd(executable, prefix, descriptions, release):
    """Check a CDD database against the release recorded for it."""
    info = blastdbcmd(executable, prefix, "-info")
    count = re.search(r"([\d,]+) sequences", info)
    if not count or int(count[1].replace(",", "")) != release["profiles"]:
        raise RuntimeError(
            f"CDD profile count differs from release {release['version']} "
            f"({release['profiles']:,} expected)"
        )
    title = re.search(r"CDD\.v(\S+)", info)
    if title and title[1] != release["version"]:
        raise RuntimeError(f"CDD database {title[0]} is not release {release['version']}")
    # NCBI's own tables omit an odd model (3.21 lacks one), so allow 1% before
    # treating the description table as belonging to another release.
    described = load_descriptions(descriptions)
    models = [
        i.rsplit("|", 1)[-1]
        for i in blastdbcmd(executable, prefix, "-entry", "all", "-outfmt", "%i").split()
    ]
    if (
        len(models) != release["profiles"]
        or sum(m not in described for m in models) > len(models) // 100
    ):
        raise RuntimeError(f"CDD description table does not match release {release['version']}")
    return info.strip()


def verify_databases(config: dict):
    tools = config["executables"]
    db = config["databases"]
    receipt = Path(db["cdd"]["path"]).parent / "installation.json"
    if not receipt.is_file():
        raise RuntimeError(f"CDD installation record is missing: {receipt}")
    release = cdd_release(json.loads(receipt.read_text()))
    info = verify_cdd(tools["blastdbcmd"], db["cdd"]["path"], db["cdd"]["descriptions"], release)
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
        "blastdb_info": info,
        "hmmer_hits": len(hmm),
        "rps_hits": len(rps),
        "control": "UniProt P0CE47 (E. coli elongation factor Tu)",
        "passed": True,
    }


def write_receipt(staging: Path, entry: dict, tools: dict, **extra):
    metadata = {
        "schema_version": 1,
        "kannotate_version": __version__,
        "installed_at": datetime.now(UTC).isoformat(),
        "specification": entry,
        **extra,
        "files": {
            str(p.relative_to(staging)): digest(p)
            for p in sorted(staging.rglob("*"))
            if p.is_file()
        },
        "disk_bytes": disk_usage(staging),
        "executables": {
            n: {"path": p, "version": executable_version(p)} for n, p in tools.items() if p
        },
    }
    (staging / "installation.json").write_text(json.dumps(metadata, indent=2) + "\n")


def cdd_target(previous: dict, storage: Path, entry: dict):
    """Keep a configured CDD installation; otherwise target NCBI's current release."""
    configured = previous.get("databases", {}).get("cdd", {}).get("path")
    if configured and Path(configured).parent.is_dir():
        if not (Path(configured).parent / "installation.json").is_file():
            raise RuntimeError(
                f"Incomplete CDD installation at {Path(configured).parent}; inspect before removing"
            )
        return Path(configured).parent, None
    release = read_release_info(entry["release_info"])
    return storage / f"cdd-{release['version']}", release


def install_cdd(entry, release, target, tools, databases):
    """Download, check and only then publish NCBI's current CDD release."""
    with tempfile.TemporaryDirectory(prefix=".cdd-", dir=target.parent) as temporary:
        work = Path(temporary)
        downloads = [fetch(asset["url"], work / asset["filename"]) for asset in entry["assets"]]
        if read_release_info(entry["release_info"]) != release:
            raise RuntimeError(
                "NCBI published a new CDD release during the download; run setup again."
            )
        observed = {path.name: digest(path) for path in downloads}
        staging = work / "cdd"
        staging.mkdir()
        safe_extract(downloads[0], staging)
        with gzip.open(downloads[1], "rb") as source, (staging / "cddid.tbl").open("wb") as out:
            shutil.copyfileobj(source, out)
        for path in downloads:
            path.unlink()
        write_receipt(staging, entry, tools, release={**release, "assets": observed})
        staged = {
            "path": str(staging / entry["database_prefix"]),
            "descriptions": str(staging / "cddid.tbl"),
        }
        verify_databases({"executables": tools, "databases": {**databases, "cdd": staged}})
        staging.rename(target)


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
        raise RuntimeError("The CDD distribution requires a little-endian host")
    storage = root / "data" / "databases"
    storage.mkdir(parents=True, exist_ok=True)
    archives = storage / ".downloads"
    archives.mkdir(exist_ok=True)
    cdd, release = cdd_target(previous, storage, spec["cdd"])
    databases = {}
    for key in ("pfam", "cdd"):
        entry = spec[key]
        target = cdd if key == "cdd" else storage / f"pfam-{entry['version']}"
        marker = target / "installation.json"
        if marker.is_file():
            metadata = json.loads(marker.read_text())
            if key == "pfam" and metadata.get("specification") != entry:
                raise RuntimeError(f"Installed {key} metadata differs from pinned specification")
            for name, checksum in metadata["files"].items():
                if not (target / name).is_file() or digest(target / name) != checksum:
                    raise RuntimeError(f"Installed {key} file is missing or changed: {name}")
            report(f"{key}: verified existing installation at {target.name}")
        else:
            if shutil.disk_usage(storage).free < 12 * 1024**3:
                raise RuntimeError(
                    "Setup needs at least 12 GiB free for temporary database preparation"
                )
            if target.exists():
                raise RuntimeError(f"Incomplete installation at {target}; inspect before removing")
            if key == "cdd":
                report(f"Installing NCBI's current CDD release {release['version']}…")
                install_cdd(entry, release, target, tools, databases)
            else:
                report(f"Installing {key} {entry['version']}…")
                asset = entry["assets"][0]
                archive = archives / asset["filename"]
                fetch(asset["url"], archive, sha256=asset["sha256"], md5=asset.get("md5"))
                with tempfile.TemporaryDirectory(prefix=f".{key}-", dir=storage) as temporary:
                    staging = Path(temporary)
                    with (
                        gzip.open(archive, "rb") as source,
                        (staging / "Pfam-A.hmm").open("wb") as out,
                    ):
                        shutil.copyfileobj(source, out)
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
                    write_receipt(staging, entry, tools)
                    staging.rename(target)
        version = (
            entry["version"]
            if key == "pfam"
            else cdd_release(json.loads(marker.read_text()))["version"]
        )
        databases[key] = {"version": version, "path": str(target / entry["database_prefix"])}
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
        f"Ready: Pfam {databases['pfam']['version']}; CDD {databases['cdd']['version']}; "
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
            # CDD follows NCBI's releases; any complete recorded installation is usable.
            matches = key == "cdd" or metadata["specification"] == specification()[key]
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
