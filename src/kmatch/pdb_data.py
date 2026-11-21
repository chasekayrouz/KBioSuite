"""Cached public PDB annotations/assemblies. Only public entry IDs leave this Mac."""

import hashlib
import http.client
import json
import os
import re
import ssl
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from kbiosuite import tls

from . import __version__
from .models import Progress
from .process import Cancelled

AGENT = {"User-Agent": f"Kmatch/{__version__} public-PDB-reader"}

TARGET = re.compile(r"^([0-9][A-Za-z0-9]{3})(?:-assembly([0-9]+))?_([A-Za-z0-9_.-]+)$")
GRAPHQL = """query { entries(entry_ids: IDS) {
 rcsb_id struct { title }
 polymer_entities {
  rcsb_polymer_entity { pdbx_description }
  rcsb_polymer_entity_container_identifiers { auth_asym_ids asym_ids entity_id }
 }
 nonpolymer_entities {
  pdbx_entity_nonpoly { comp_id name }
  rcsb_nonpolymer_entity_container_identifiers { auth_asym_ids asym_ids }
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


def failure_message(exc):
    """Why a public RCSB request failed, for the user; the caller says how to retry."""
    if isinstance(exc, urllib.error.HTTPError):
        return f"RCSB returned HTTP {exc.code}."
    if isinstance(exc, http.client.HTTPException):
        return INTERRUPTED
    reason = exc.reason if isinstance(exc, urllib.error.URLError) else exc
    if isinstance(reason, ssl.SSLCertVerificationError):
        return (
            "macOS could not verify RCSB's security certificate. Check this Mac's date and "
            "time and any network filter."
        )
    return "RCSB could not be reached. Check the network connection."


INTERRUPTED = "RCSB's response was interrupted."
UNEXPECTED = "RCSB returned an unexpected response."
RELOAD = " Choose Reload details to try again."


class RequestFailed(RuntimeError):
    """A public RCSB request failed in transit. Nothing was cached, so a retry can succeed."""


def declared_length(headers):
    """The body length a response declares, or None when its framing declares none.

    An intermediary may repeat or join identical values ("5, 5"). Anything else that is not
    one byte count is refused as a failed transfer, never read as "no length declared".
    """
    if "chunked" in (headers.get("Transfer-Encoding") or "").lower():
        return None  # chunk decoding, and any error in it, is the HTTP reader's
    declared = headers.get_all("Content-Length") or []
    values = {part.strip() for value in declared for part in value.split(",")}
    if not values:
        return None
    if len(values) != 1 or not all(v.isascii() and v.isdigit() for v in values):
        raise http.client.HTTPException("Unusable Content-Length: " + ", ".join(sorted(values)))
    return int(values.pop())


def fetch(url, cancel, limit):
    """Bounded reads of one public RCSB response.

    Cancellation, errors reported by the HTTP reader, and detected Content-Length mismatches
    abort the transfer. Chunk decoding relies on the HTTP reader. For a close-delimited
    response, EOF alone cannot establish completeness.
    """
    cancel.check()
    request = urllib.request.Request(url, headers=AGENT)
    deadline = time.monotonic() + 60
    chunks, size = [], 0
    try:
        with tls.urlopen(request, timeout=8) as response:
            while True:
                cancel.check()
                chunk = response.read(256 * 1024)
                if not chunk:
                    break
                size += len(chunk)
                if size > limit or time.monotonic() > deadline:
                    raise RuntimeError("Public PDB response exceeded the size/time limit.")
                chunks.append(chunk)
            # A fixed-length body that ends early reads as an ordinary end of file.
            declared = declared_length(response.headers)
            if declared is not None and size != declared:
                raise http.client.IncompleteRead(b"", max(declared - size, 0))
    except (OSError, http.client.HTTPException) as exc:
        cancel.check()
        raise RequestFailed(failure_message(exc)) from exc
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


def normalize_entry(raw):
    """One answered entry as the record Kmatch caches and displays.

    Every field GRAPHQL asks for must be present with a type RCSB documents: a collection or
    text field may be null, but a field left out, or one of another type, raises ValueError. An
    incomplete answer must never be saved as "no ligands" or a wrong instance count, and an
    instance list RCSB leaves null is kept as an unknown count.
    """

    def field(mapping, key, kinds):
        if not isinstance(mapping, dict) or key not in mapping:
            raise ValueError(f"RCSB's answer has no {key}")
        if not isinstance(mapping[key], kinds):
            raise ValueError(f"RCSB's answer has an unexpected {key}")
        return mapping[key]

    def optional(mapping, key, inner, kinds):
        """mapping[key][inner], where mapping[key] may be null."""
        outer = field(mapping, key, dict | None)
        return None if outer is None else field(outer, inner, kinds)

    def chains(identifiers, key):
        """A list of chain IDs, or None when RCSB reports none for this entity."""
        value = field(identifiers, key, list | None)
        if value is not None and not all(isinstance(chain, str) for chain in value):
            raise ValueError(f"RCSB's answer has an unexpected {key}")
        return value

    text = str | None
    polymers, ligands = [], []
    for entity in field(raw, "polymer_entities", list | None) or []:
        identifiers = field(entity, "rcsb_polymer_entity_container_identifiers", dict)
        polymers.append(
            {
                "name": optional(entity, "rcsb_polymer_entity", "pdbx_description", text),
                # No chains reported: nothing for a hit's chain to match.
                "author_chains": chains(identifiers, "auth_asym_ids") or [],
                "label_chains": chains(identifiers, "asym_ids") or [],
                "entity_id": field(identifiers, "entity_id", str),
            }
        )
    for entity in field(raw, "nonpolymer_entities", list | None) or []:
        # Either description may be null; one that is present must be complete.
        component = optional(entity, "nonpolymer_comp", "chem_comp", dict | None)
        fallback = field(entity, "pdbx_entity_nonpoly", dict | None)
        described = {key: None for key in ("id", "name", "formula", "comp_id", "fallback_name")}
        if component is not None:
            for key in ("id", "name", "formula"):
                described[key] = field(component, key, text)
        if fallback is not None:
            described["comp_id"] = field(fallback, "comp_id", text)
            described["fallback_name"] = field(fallback, "name", text)
        code = described["id"] or described["comp_id"]
        if not code:
            raise ValueError("RCSB's answer has a small molecule without a component ID")
        identifiers = field(entity, "rcsb_nonpolymer_entity_container_identifiers", dict)
        instances = chains(identifiers, "asym_ids")
        ligand = {
            "id": code,
            "name": described["name"] or described["fallback_name"] or code,
            "formula": described["formula"],
            "author_chains": chains(identifiers, "auth_asym_ids") or [],
            # The number of instances RCSB lists; None when it lists none, which is not zero.
            "copies": None if instances is None else len(instances),
        }
        if code not in {"HOH", "DOD", "WAT"}:
            ligands.append(ligand)
    return {
        "pdb_id": field(raw, "rcsb_id", str),
        "title": optional(raw, "struct", "title", text),
        "polymers": polymers,
        "ligands": ligands,
    }


def usable(entry, identifier):
    """Whether a saved record for IDENTIFIER has every field the hit table and ligand panel read.

    Checked before a record is cached and again when it is read back, so a record that could
    not be displayed is requested again instead of failing on every reload.
    """

    def text(value, optional=False):
        return isinstance(value, str) or (optional and value is None)

    def chains(value):
        return isinstance(value, list) and all(text(chain) for chain in value)

    try:
        return (
            entry["pdb_id"] == identifier
            and text(entry.get("title"), optional=True)
            and text(entry.get("retrieved_utc"), optional=True)
            and isinstance(entry["polymers"], list)
            and isinstance(entry["ligands"], list)
            and all(
                text(polymer["name"], optional=True)
                and chains(polymer["author_chains"])
                and chains(polymer["label_chains"])
                for polymer in entry["polymers"]
            )
            and all(
                text(ligand["id"])
                and text(ligand["name"])
                and text(ligand["formula"], optional=True)
                and (
                    ligand["copies"] is None
                    or (type(ligand["copies"]) is int and ligand["copies"] >= 0)
                )
                for ligand in entry["ligands"]
            )
        )
    except (KeyError, TypeError, AttributeError):
        return False


def cached_entry(path, identifier):
    """A saved record ready for display, or None. The saved file is never rewritten here."""
    try:
        entry = json.loads(path.read_text())
        # 0.6.1 and earlier saved RCSB's null chain lists as they arrived, and displayed them.
        for item in [*entry["polymers"], *entry["ligands"]]:
            for key in ("author_chains", "label_chains"):
                if isinstance(item, dict) and key in item and item[key] is None:
                    item[key] = []
    except (OSError, ValueError, KeyError, TypeError):
        return None
    return entry if usable(entry, identifier) else None


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


def metadata_url(identifiers):
    query = GRAPHQL.replace("IDS", json.dumps(sorted(identifiers)))
    return "https://data.rcsb.org/graphql?" + urllib.parse.urlencode({"query": query})


def load_details(config, targets, allow_network, cancel):
    # failures: entries worth asking for again, as identifier -> (message, technical detail).
    entries, missing, failures = {}, set(), {}
    for target in targets:
        try:
            identifier = pdb_identity(target)[0]
        except ValueError:
            continue
        entry = cached_entry(config.pdb_cache / "metadata" / (identifier + ".json"), identifier)
        if entry is None:
            missing.add(identifier)
        else:
            entries[identifier] = entry
    if missing and allow_network:
        failed = None  # the whole request, or a record that names no requested entry
        try:
            response = json.loads(fetch(metadata_url(missing), cancel, 8 * 1024 * 1024))
            if not isinstance(response, dict):
                raise ValueError("RCSB's answer is not a JSON object")
            if response.get("errors"):
                raise RuntimeError("RCSB could not return the requested annotations.")
            returned = response["data"]["entries"]
            if not isinstance(returned, list | None):
                raise ValueError("RCSB's entries are not a list")
            # Each record stands alone: one that cannot be used or saved never costs another.
            for raw in returned or []:
                if raw is None:
                    continue  # RCSB's answer for an ID it does not know
                identifier = raw.get("rcsb_id") if isinstance(raw, dict) else None
                if not isinstance(identifier, str):
                    failed = (UNEXPECTED + RELOAD, "")
                    continue
                if identifier not in missing or identifier in entries:
                    continue
                try:
                    entry = normalize_entry(raw)
                except ValueError as exc:
                    failures[identifier] = (UNEXPECTED + RELOAD, str(exc))
                    continue
                entry.update(
                    source="https://data.rcsb.org/graphql",
                    retrieved_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                )
                cancel.check()
                try:
                    cache_write(
                        config.pdb_cache / "metadata" / (identifier + ".json"),
                        json.dumps(entry, indent=2).encode(),
                    )
                except OSError as exc:
                    failures[identifier] = (
                        f"PDB details for {identifier} could not be saved to Kmatch's cache.",
                        str(exc),
                    )
                    continue
                failures.pop(identifier, None)
                entries[identifier] = entry
        except Cancelled:
            raise
        except RuntimeError as exc:
            # The technical cause stays out of the sentence and is kept for the tooltip.
            failed = (str(exc) + RELOAD, str(exc.__cause__ or ""))
        except (ValueError, KeyError, TypeError):
            failed = (UNEXPECTED + RELOAD, "")
        if failed:
            for identifier in missing - entries.keys():
                failures.setdefault(identifier, failed)
    result = {}
    for target in targets:
        cancel.check()
        try:
            identifier = pdb_identity(target)[0]
            if identifier in entries:
                result[target] = entry_for_hit(entries[identifier], target)
            elif identifier in failures:
                message, detail = failures[identifier]
                result[target] = {
                    "status": "unavailable",
                    "url": pdb_url(target),
                    "message": message,
                    "detail": detail,
                    "retryable": True,
                }
            else:
                result[target] = {
                    "status": "unavailable",
                    "url": pdb_url(target),
                    "message": "PDB details are not cached."
                    if not allow_network
                    else "This entry was not returned by RCSB.",
                }
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
        progress(Progress("export", f"Downloading public {identifier} coordinates from RCSB…"))
        try:
            contents = fetch(url, cancel, 128 * 1024 * 1024)
        except RequestFailed as exc:
            raise RuntimeError(
                f"{exc} Export the 3D alignment again to retry, or turn off full-atom hit "
                "export to use the local Cα trace."
            ) from exc
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


def network_probes():
    """Explicit frozen qualification only: Kmatch's two public requests, for 1CSP chain A."""

    def metadata(url):
        request = urllib.request.Request(url, headers=AGENT)
        with tls.urlopen(request, timeout=20, direct=True) as response:
            status, body = response.status, response.read(1024 * 1024)
        answer = json.loads(body)
        if answer.get("errors"):
            raise RuntimeError("RCSB reported an error for public 1CSP")
        entry = normalize_entry(answer["data"]["entries"][0])
        if entry["pdb_id"] != "1CSP" or not entry_for_hit(entry, "1csp_A")["protein_name"]:
            raise RuntimeError("RCSB did not return public 1CSP with a chain A protein name")
        return status

    return {
        metadata_url(["1CSP"]): metadata,
        "https://files.rcsb.org/download/1CSP.cif": tls.head,
    }
