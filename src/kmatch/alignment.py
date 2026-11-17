"""File-based Foldseek superposition, with optional cached public full-atom targets."""

import hashlib
import json
import os
import re
import shutil
import tempfile
import time
import uuid
from dataclasses import asdict, dataclass
from pathlib import Path

import gemmi
import numpy as np

from .database import database_lease, database_status
from .models import PredictionResult, Progress, SearchResult, write_json
from .process import Cancellation, Cancelled, run_process
from .target_traces import load_target_trace

FORMATS = ("pse", "cif", "pdb")
PYMOL_RUNNER = Path(__file__).with_name("pymol_runner.py")
AF3_TERMS = "https://github.com/google-deepmind/alphafold3/blob/main/OUTPUT_TERMS_OF_USE.md"


def pymol_executable(config):
    """Detect only; probing a licensed program never blocks GUI polling."""
    candidates = (
        [Path(config.pymol_override).expanduser()]
        if config.pymol_override
        else [
            Path.home() / "Applications/PyMOL.app/Contents/bin/pymol",
            Path("/Applications/PyMOL.app/Contents/bin/pymol"),
            Path(shutil.which("pymol") or "/nonexistent/kmatch-pymol"),
        ]
    )
    for path in candidates:
        if path.suffix == ".app":
            path /= "Contents/bin/pymol"
        if path.is_file() and os.access(path, os.X_OK):
            return path.absolute()
    return None


def export_support(config):
    executable = pymol_executable(config)
    return {
        "pse": bool(executable),
        "message": (
            "Local PyMOL detected; session generation is verified on export."
            if executable
            else "PyMOL is unavailable. Use mmCIF/PDB, or install PyMOL locally and set "
            "pymol_override in config.json to its executable or .app."
        ),
        "target_representation": (
            "Full-atom export uses cached public PDB coordinates; "
            "local Cα export remains available."
        ),
    }


def run_structure(run_directory, recorded):
    """The predicted structure a saved run should read.

    A recorded path outside RUN_DIRECTORY (for example after the run was moved or copied)
    is replaced by the run's own prediction/<same name> copy when that file exists. Otherwise
    the recorded path is returned unchanged, so a run without its own copy still reads its
    recorded (possibly old-root) location; a missing file then fails when it is read.
    """
    recorded, run_directory = Path(recorded), Path(run_directory)
    own = run_directory / "prediction" / recorded.name
    if not recorded.resolve().is_relative_to(run_directory.resolve()) and own.is_file():
        return own
    return recorded


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_]", "_", str(value))[:64] or "structure"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ca_chain(structure, label):
    if len(structure) != 1 or len(structure[0]) != 1:
        raise ValueError(f"{label} must contain exactly one model and one protein chain.")
    chain = structure[0][0]
    atoms, sequence = [], []
    for residue in chain:
        ca = [atom for atom in residue if atom.name == "CA"]
        if not ca:
            continue
        if len(ca) != 1:
            raise ValueError(f"{label} has ambiguous alternate Cα coordinates.")
        atoms.append(ca[0])
        # Gemmi uses a blank code for unknown residues; Foldseek represents these as X.
        sequence.append(
            gemmi.find_tabulated_residue(residue.name).one_letter_code.upper().strip() or "X"
        )
    xyz = np.array([[atom.pos.x, atom.pos.y, atom.pos.z] for atom in atoms], dtype=float)
    if not len(atoms) or not np.isfinite(xyz).all():
        raise ValueError(f"{label} has missing or nonfinite Cα coordinates.")
    return chain, "".join(sequence), xyz


def aligned_pairs(hit, query_sequence, target_sequence):
    if len(query_sequence) != hit.query_length or len(target_sequence) != hit.target_length:
        raise ValueError("Coordinate lengths do not match this saved Foldseek hit.")
    if len(hit.query_alignment) != len(hit.target_alignment):
        raise ValueError("The saved Foldseek alignment is incomplete.")
    qi, ti = hit.query_start - 1, hit.target_start - 1
    pairs = []
    for query, target in zip(
        hit.query_alignment.upper(), hit.target_alignment.upper(), strict=True
    ):
        for letter, index, sequence in ((query, qi, query_sequence), (target, ti, target_sequence)):
            if letter != "-" and (not 0 <= index < len(sequence) or letter != sequence[index]):
                raise ValueError(
                    "Coordinate sequence does not match the saved alignment. "
                    "The query or database entry may have changed."
                )
        if query != "-" and target != "-":
            pairs.append((qi, ti))
        qi += query != "-"
        ti += target != "-"
    if not pairs or qi != hit.query_end or ti != hit.target_end:
        raise ValueError("Saved Foldseek alignment endpoints do not match its aligned sequences.")
    return np.array(pairs, dtype=int)


def fixed_tm(distances, length):
    d0 = max(0.5, 1.24 * np.cbrt(length - 15) - 1.8) if length > 21 else 0.5
    return float(np.sum(1 / (1 + (distances / d0) ** 2)) / length)


def superpose(query, target, hit):
    """Apply the saved target→query transform once. Never optimize/refit coordinates."""
    qchain, qseq, qxyz = ca_chain(query, "Query")
    tchain, tseq, txyz = ca_chain(target, "Database hit")
    pairs = aligned_pairs(hit, qseq, tseq)
    rotation, translation = np.asarray(hit.rotation), np.asarray(hit.translation)
    if (
        rotation.shape != (3, 3)
        or translation.shape != (3,)
        or not np.isfinite(rotation).all()
        or not np.isfinite(translation).all()
        or not np.allclose(rotation @ rotation.T, np.eye(3), atol=0.002)
        or not np.isclose(np.linalg.det(rotation), 1.0, atol=0.002)
    ):
        raise ValueError("The saved Foldseek transformation is not a finite proper rigid rotation.")
    transformed = txyz @ rotation.T + translation
    qi, ti = pairs.T
    before = np.linalg.norm(qxyz[qi] - txyz[ti], axis=1)
    after = np.linalg.norm(qxyz[qi] - transformed[ti], axis=1)
    metrics = {
        "aligned_ca_pairs": len(pairs),
        "aligned_ca_rmsd_before_angstrom": float(np.sqrt(np.mean(before**2))),
        "aligned_ca_rmsd_after_angstrom": float(np.sqrt(np.mean(after**2))),
        "fixed_transform_query_tm": fixed_tm(after, hit.query_length),
        "fixed_transform_target_tm": fixed_tm(after, hit.target_length),
        "note": "Measured over saved aligned residue pairs without refitting. Foldseek's "
        "individually optimized TM normalizations can differ from these fixed-transform scores.",
    }
    combined = gemmi.Structure()
    combined.name = "Kmatch_superposition"
    model = gemmi.Model("1")
    query_copy, target_copy = qchain.clone(), tchain.clone()
    query_copy.name, target_copy.name = "Q", "T"
    transform = gemmi.Transform()
    transform.mat.fromlist(rotation.tolist())
    transform.vec.fromlist(translation.tolist())
    for residue in target_copy:
        for atom in residue:
            atom.pos = list(transform.apply(atom.pos))
    model.add_chain(query_copy)
    model.add_chain(target_copy)
    combined.add_model(model)
    combined.setup_entities()
    combined.assign_label_seq_id()
    return combined, metrics, {"Q": qchain.name, "T": tchain.name}


def write_coordinates(structure, path, provenance):
    encoded = json.dumps(provenance, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
    if path.suffix == ".cif":
        document = structure.make_mmcif_document()
        document.sole_block().set_pair(
            "_kmatch_alignment.provenance_json", gemmi.cif.quote(encoded)
        )
        document.write_file(str(path))
    else:
        atoms = [atom for chain in structure[0] for residue in chain for atom in residue]
        if (
            len(atoms) > 99999
            or any(
                not -999 <= residue.seqid.num <= 9999 for chain in structure[0] for residue in chain
            )
            or any(not -999.999 <= value <= 9999.999 for atom in atoms for value in atom.pos)
        ):
            raise ValueError("This superposition exceeds PDB field limits. Export mmCIF instead.")
        remarks = "REMARK 950 KMATCH SUPERPOSITION: QUERY Q; HIT T; OPTIONAL LIGANDS L\n"
        remarks += "REMARK 950 PROVENANCE JSON: JOIN REMARK 951 PAYLOADS IN ORDER\n"
        remarks += "".join(
            f"REMARK 951 {encoded[i : i + 68]}\n" for i in range(0, len(encoded), 68)
        )
        path.write_text(remarks + structure.make_pdb_string())
    # Verify a format round-trip before making the artifact available.
    reread = gemmi.read_structure(str(path))
    original = np.array([list(a.pos) for c in structure[0] for r in c for a in r])
    restored = np.array([list(a.pos) for c in reread[0] for r in c for a in r])
    if (
        len(reread) != 1
        or [c.name for c in reread[0]] != [c.name for c in structure[0]]
        or (
            original.shape != restored.shape
            or not np.allclose(original, restored, atol=0.001, rtol=0)
        )
    ):
        raise RuntimeError(
            "Exported coordinates did not round-trip correctly; no file was published."
        )


def materialize_target(config, metadata, hit, work, log, cancel, progress):
    database = Path(metadata["database"])
    state = database_status(database)
    if not state.ready or not Path(str(database) + ".lookup").is_file():
        raise RuntimeError(
            "This legacy hit has no retained target trace, and its original PDB100 generation "
            f"is unavailable: {database}. Run a new search with the saved prediction to export "
            "this hit. An RCSB download alone cannot safely recover the original Foldseek frame."
        )
    if metadata.get("database_version") and state.version != metadata["database_version"]:
        raise RuntimeError(
            "The database generation differs from the saved search; restore it first."
        )
    if not hit.target or any(character.isspace() or character == "\0" for character in hit.target):
        raise ValueError("The saved target identifier cannot be materialized safely.")
    selected, subset = work / "selected.txt", work / "target"
    selected.write_text(hit.target + "\n")
    progress(Progress("export", "Extracting the selected hit's local PDB100 Cα coordinates…"))
    for command in (
        [config.foldseek, "createsubdb", selected, database, subset, "--id-mode", "1", "-v", "1"],
        [
            config.foldseek,
            "convert2pdb",
            subset,
            work / "target.cif",
            "--pdb-output-format",
            "1",
            "--threads",
            "1",
            "-v",
            "1",
        ],
    ):
        run_process(
            command,
            log,
            cancel,
            progress,
            "target extraction",
            timeout=120,
            scratch_root=config.cache / "process",
        )
    entries = subset.with_suffix(".index").read_text().splitlines()
    if len(entries) != 1:
        raise RuntimeError("The selected identifier did not resolve to exactly one database entry.")
    target = gemmi.read_structure(str(work / "target.cif"))
    if any(
        atom.name != "CA" for model in target for chain in model for res in chain for atom in res
    ):
        raise RuntimeError("Unexpected coordinate representation in the local PDB100 hit.")
    return target, entries[0].split()[0]


@dataclass(frozen=True)
class AlignmentExport:
    export_id: str
    directory: Path
    artifact: Path
    format: str
    target: str
    query: str
    metrics: dict
    saved_to: str | None = None


def save_export(result, destination, config, cancel):
    destination = Path(destination).expanduser().absolute()
    if destination.suffix.lower() != "." + result.format:
        raise ValueError(f"Use a .{result.format} filename for this export format.")
    provenance = json.loads((result.directory / "provenance.json").read_text())
    protected = [config.home, config.cache, Path(provenance["run_directory"])]
    if (
        destination.is_symlink()
        or destination.resolve() == Path(provenance["query_structure"]).resolve()
        or any(destination.resolve().is_relative_to(root.resolve()) for root in protected)
    ):
        raise ValueError("Choose a destination outside Kmatch's managed assets and run data.")
    # Same-directory staging means cancellation/failure never truncates an existing export.
    descriptor, temporary = tempfile.mkstemp(prefix=".kmatch-export-", dir=destination.parent)
    os.close(descriptor)
    try:
        shutil.copyfile(result.artifact, temporary)
        cancel.check()
        os.replace(temporary, destination)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return str(destination)


def export_alignment(
    config,
    run_directory,
    hit_index,
    format,
    *,
    cancel=None,
    progress=lambda event: None,
    destination=None,
    output_root=None,
    expected_target=None,
    full_atom=False,
    allow_network=False,
):
    if type(full_atom) is not bool or type(allow_network) is not bool:
        raise ValueError("Full-atom and public PDB download options must be true or false.")
    cancel = cancel or Cancellation()
    cancel.check()
    if format not in FORMATS:
        raise ValueError("Choose PyMOL session, mmCIF or PDB for the 3D alignment.")
    pymol = pymol_executable(config) if format == "pse" else None
    if format == "pse" and not pymol:
        raise RuntimeError(export_support(config)["message"])
    run_directory = Path(run_directory)
    metadata = json.loads((run_directory / "search/search.json").read_text())
    if type(hit_index) is not int or not 0 <= hit_index < len(metadata["hits"]):
        raise ValueError("Select an existing structural-search hit.")
    hit = SearchResult(**metadata["hits"][hit_index])
    if expected_target is not None and hit.target != expected_target:
        raise ValueError("The selected search hit changed. Select it again before exporting.")
    prediction_data = json.loads((run_directory / "prediction/prediction.json").read_text())
    # The saved search's SHA-256 still gates whichever copy is read.
    prediction_data["structure_path"] = run_structure(
        run_directory, prediction_data["structure_path"]
    )
    prediction = PredictionResult(**prediction_data)
    query_hash = digest(prediction.structure_path)
    if metadata.get("query_structure_sha256") and metadata["query_structure_sha256"] != query_hash:
        raise ValueError("The query coordinates changed since this search. Run a new search first.")
    query = gemmi.read_structure(str(prediction.structure_path))
    run = json.loads((run_directory / "run.json").read_text())
    query_name = run.get("sequence_identifier", hit.query)
    export_id = time.strftime("%Y%m%d-%H%M%S") + "-" + uuid.uuid4().hex[:8]
    directory = Path(output_root or run_directory / "exports") / export_id
    directory.mkdir(parents=True, mode=0o700)
    work = directory / ".work"
    work.mkdir(mode=0o700)
    log = directory / "export.log"
    status = {"export_id": export_id, "status": "running", "format": format, "target": hit.target}
    write_json(directory / "export.json", status)
    try:
        retained = load_target_trace(run_directory / "search", metadata, hit)
        if retained is not None:
            target, key = retained
            progress(Progress("export", "Using this run's retained target Cα coordinates…"))
        else:
            with database_lease(config, cancel):
                target, key = materialize_target(config, metadata, hit, work, log, cancel, progress)
        cancel.check()
        progress(Progress("export", "Applying the saved Foldseek transformation; no refitting…"))
        combined, metrics, chain_mapping = superpose(query, target, hit)
        full_atom_details = None
        if full_atom:
            from .hit_structure import full_atom_superposition

            full_atom_details = full_atom_superposition(
                config, combined, target, hit, allow_network, cancel, progress
            )
        provenance = {
            "schema_version": 1,
            "export_id": export_id,
            "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "run_id": run_directory.name,
            "run_directory": str(run_directory.resolve()),
            "query_identifier": query_name,
            "query_structure": str(prediction.structure_path),
            "query_sha256": query_hash,
            "query_hash_recorded_at_search": bool(metadata.get("query_structure_sha256")),
            "predictor": prediction.predictor,
            "predictor_version": prediction.predictor_version,
            "database": metadata["database"],
            "database_version": metadata.get("database_version"),
            "database_key": key,
            "target_coordinate_source": "run-retained trace" if retained else "original database",
            "target_trace": metadata.get("target_traces", {}).get(hit.target),
            "foldseek_version": metadata.get("foldseek_version"),
            "hit": asdict(hit),
            "transform": "x_query = u @ x_target + t; raw saved Foldseek u/t",
            "refitted": False,
            "target_representation": "Cα only from local Foldseek database",
            "target_numbering": "1-based database sequence positions, not original PDB author IDs",
            "coordinate_chains": {"Q": "query", "T": "superposed target"},
            "source_chain_names": chain_mapping,
            "metrics": metrics,
            "pymol_objects": {
                "query": "query_" + safe_name(query_name),
                "hit": "hit_" + safe_name(hit.target),
            },
        }
        if full_atom_details is not None:
            provenance["full_atom"] = full_atom_details
            provenance["target_representation"] = full_atom_details["representation"]
            provenance["target_numbering"] = "original PDB author numbers and insertion codes"
            provenance["source_chain_names"]["T"] = full_atom_details["source_chain"]
            if full_atom_details["ligands"]:
                provenance["coordinate_chains"]["L"] = "contacting small molecules / ions"
        if prediction.predictor == "af3":
            provenance["query_output_terms"] = AF3_TERMS
            provenance["query_output_note"] = (
                "AF3 query remains in its original coordinate frame; serialization may round "
                "coordinates and the exported chain label is Q. The target was superposed "
                "separately. AF3 output terms apply."
            )
        write_json(directory / "provenance.json", provenance)
        coordinate_path = work / ("alignment.pdb" if format == "pdb" else "alignment.cif")
        write_coordinates(combined, coordinate_path, provenance)
        artifact = work / ("alignment." + format)
        if format == "pse":
            request = {
                "coordinates": str(coordinate_path),
                "output": str(artifact),
                "provenance": provenance,
                "verification": str(work / "session-check.json"),
            }
            write_json(work / "session-request.json", request)
            progress(Progress("export", "Creating and reopening the local PyMOL session…"))
            try:
                run_process(
                    [pymol, "-cqk", "-r", PYMOL_RUNNER],
                    log,
                    cancel,
                    progress,
                    "PyMOL session export",
                    timeout=120,
                    extra_env={"KMATCH_ALIGNMENT_REQUEST": str(work / "session-request.json")},
                    scratch_root=config.cache / "process",
                )
                verification = json.loads((work / "session-check.json").read_text())
            except Cancelled:
                raise
            except (OSError, ValueError, RuntimeError) as exc:
                raise RuntimeError(
                    f"Local PyMOL could not create a verified session: {exc} "
                    "Use mmCIF/PDB export, or repair the configured PyMOL installation. "
                    f"Diagnostics: {log}"
                ) from exc
            if not verification.get("roundtrip_verified") or not artifact.is_file():
                raise RuntimeError("PyMOL did not produce a verified session. Try mmCIF export.")
            write_json(directory / "session-check.json", verification)
        cancel.check()
        artifact = artifact.replace(directory / artifact.name)
        result = AlignmentExport(
            export_id, directory, artifact, format, hit.target, query_name, metrics
        )
        if destination:
            from dataclasses import replace

            result = replace(result, saved_to=save_export(result, destination, config, cancel))
        status.update(
            status="complete", artifact=str(artifact), saved_to=result.saved_to, metrics=metrics
        )
        write_json(directory / "export.json", status)
        progress(Progress("complete", f"3D alignment exported · {hit.target} · {format.upper()}"))
        return result
    except BaseException as exc:
        status.update(
            status="cancelled" if isinstance(exc, (Cancelled, KeyboardInterrupt)) else "failed",
            error=str(exc),
        )
        write_json(directory / "export.json", status)
        raise
    finally:
        shutil.rmtree(work, ignore_errors=True)
