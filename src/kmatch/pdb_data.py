"""Cached public PDB annotations/assemblies. Only public entry IDs leave this Mac."""

import hashlib
import json
import os
import re
import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

from . import __version__
from .models import Progress
from .process import Cancelled

TARGET = re.compile(r"^([0-9][A-Za-z0-9]{3})(?:-assembly([0-9]+))?_([A-Za-z0-9_.-]+)$")
GRAPHQL = """query { entries(entry_ids: IDS) {
 rcsb_id struct { title }
 polymer_entities {
  rcsb_polymer_entity { pdbx_description }
  rcsb_polymer_entity_container_identifiers { auth_asym_ids asym_ids entity_id }
 }
 nonpolymer_entities {
  pdbx_entity_nonpoly { comp_id name }
  rcsb_nonpolymer_entity_container_identifiers { auth_asym_ids asym_ids entity_id }
  nonpolymer_comp { chem_comp { id name formula } }
 }
} }"""


def pdb_identity(target):
    match = TARGET.fullmatch(target)
    if not match:
        raise ValueError("This hit has no supported PDB entry/chain identifier.")
    return match[1].upper(), match[2], match[3]


def pdb_url(target):
    return "https://www.rcsb.org/structure/" + pdb_identity(target)[0]


def fetch(url, cancel, limit):
    """Bounded reads; cancelled/partial responses never become cache entries."""
    cancel.check()
    request = urllib.request.Request(
        url, headers={"User-Agent": f"Kmatch/{__version__} public-PDB-reader"}
    )
    deadline = time.monotonic() + 60
    chunks, size = [], 0
    try:
        with urllib.request.urlopen(request, timeout=8) as response:
            while True:
                cancel.check()
                chunk = response.read(256 * 1024)
                if not chunk:
                    break
                size += len(chunk)
                if size > limit or time.monotonic() > deadline:
                    raise RuntimeError("Public PDB response exceeded the size/time limit.")
                chunks.append(chunk)
    except OSError:
        cancel.check()
        raise
    cancel.check()
    return b"".join(chunks)


def cache_write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".pdb-", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def normalize_entry(entry):
    polymers, ligands = [], []
    for entity in entry.get("polymer_entities") or []:
        identifiers = entity["rcsb_polymer_entity_container_identifiers"]
        polymers.append(
            {
                "name": (entity.get("rcsb_polymer_entity") or {}).get("pdbx_description"),
                "author_chains": identifiers.get("auth_asym_ids", []),
                "label_chains": identifiers.get("asym_ids", []),
                "entity_id": identifiers["entity_id"],
            }
        )
    for entity in entry.get("nonpolymer_entities") or []:
        component = (entity.get("nonpolymer_comp") or {}).get("chem_comp") or {}
        fallback = entity.get("pdbx_entity_nonpoly") or {}
        code = component.get("id") or fallback.get("comp_id")
        if not code or code in {"HOH", "DOD", "WAT"}:
            continue
        identifiers = entity["rcsb_nonpolymer_entity_container_identifiers"]
        ligands.append(
            {
                "id": code,
                "name": component.get("name") or fallback.get("name") or code,
                "formula": component.get("formula"),
                "author_chains": identifiers.get("auth_asym_ids", []),
                "copies": len(identifiers.get("asym_ids", [])),
            }
        )
    return {
        "pdb_id": entry["rcsb_id"],
        "title": (entry.get("struct") or {}).get("title"),
        "polymers": polymers,
        "ligands": ligands,
    }


def entry_for_hit(entry, target):
    pdb_id, _, chain = pdb_identity(target)
    # Biological-assembly copy suffixes are not new author chain identifiers.
    chains = [chain, re.sub(r"-\d+$", "", chain)]
    matches = [p for p in entry["polymers"] if any(c in p["author_chains"] for c in chains)]
    if not matches:
        matches = [p for p in entry["polymers"] if chain in p["label_chains"]]
    names = sorted({p["name"] for p in matches if p["name"]})
    return {
        **entry,
        "protein_name": "; ".join(names) or None,
        "url": "https://www.rcsb.org/structure/" + pdb_id,
        "status": "ready",
    }


def load_details(config, targets, allow_network, cancel):
    entries, missing, failures = {}, set(), {}
    for target in targets:
        try:
            identifier = pdb_identity(target)[0]
        except ValueError:
            continue
        path = config.pdb_cache / "metadata" / (identifier + ".json")
        try:
            entry = json.loads(path.read_text())
            if entry["pdb_id"] != identifier or not isinstance(entry["polymers"], list):
                raise ValueError("Invalid cached PDB metadata")
            entries[identifier] = entry
        except (OSError, ValueError, KeyError, TypeError):
            missing.add(identifier)
    if missing and allow_network:
        try:
            query = GRAPHQL.replace("IDS", json.dumps(sorted(missing)))
            url = "https://data.rcsb.org/graphql?" + urllib.parse.urlencode({"query": query})
            response = json.loads(fetch(url, cancel, 8 * 1024 * 1024))
            if response.get("errors"):
                raise RuntimeError("RCSB could not return the requested annotations.")
            for raw in response["data"]["entries"] or []:
                if not raw:
                    continue
                entry = normalize_entry(raw)
                identifier = entry["pdb_id"]
                if identifier not in missing:
                    continue
                entry.update(
                    source="https://data.rcsb.org/graphql",
                    retrieved_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                )
                cancel.check()
                cache_write(
                    config.pdb_cache / "metadata" / (identifier + ".json"),
                    json.dumps(entry, indent=2).encode(),
                )
                entries[identifier] = entry
        except Cancelled:
            raise
        except (OSError, ValueError, KeyError, TypeError, RuntimeError) as exc:
            failures = dict.fromkeys(missing, "RCSB details unavailable: " + str(exc))
    result = {}
    for target in targets:
        cancel.check()
        try:
            identifier = pdb_identity(target)[0]
            result[target] = (
                entry_for_hit(entries[identifier], target)
                if identifier in entries
                else {
                    "status": "unavailable",
                    "url": pdb_url(target),
                    "message": failures.get(
                        identifier,
                        "PDB details are not cached."
                        if not allow_network
                        else "This entry was not returned by RCSB.",
                    ),
                }
            )
        except ValueError as exc:
            result[target] = {"status": "unavailable", "message": str(exc)}
    return result


def assembly_file(config, target, allow_network, cancel, progress):
    identifier, assembly, _ = pdb_identity(target)
    filename = identifier + (f"-assembly{assembly}" if assembly else "") + ".cif"
    path = config.pdb_cache / "coordinates" / filename
    url = "https://files.rcsb.org/download/" + filename
    if not path.is_file():
        if not allow_network:
            raise RuntimeError(
                "Full-atom PDB coordinates are not cached. Enable public PDB downloads "
                "or turn off full-atom hit export to use the local Cα trace."
            )
        progress(
            Progress("export", f"Downloading public {identifier} coordinates from RCSB (cached)…")
        )
        contents = fetch(url, cancel, 128 * 1024 * 1024)
        import gemmi

        block = gemmi.cif.read_string(contents.decode()).sole_block()
        entry_id = (block.find_value("_entry.id") or "").strip("'\"").upper()
        # RCSB-generated biological assemblies can carry the generated-file placeholder XXXX.
        # Their identity/frame is further checked against the local hit trace before export.
        if entry_id != identifier and not (assembly and entry_id == "XXXX"):
            raise RuntimeError(
                "RCSB coordinate response does not identify the requested PDB entry."
            )
        if not len(gemmi.make_structure_from_block(block)):
            raise RuntimeError("RCSB returned no coordinate models.")
        cancel.check()
        cache_write(path, contents)
    return path, {
        "url": url,
        "path": str(path),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
