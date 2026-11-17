"""Asymmetric fixtures expose inverse/transposed transforms and gap/numbering errors."""

import json
import re
import threading
import time
from dataclasses import asdict, replace
from pathlib import Path

import gemmi
import numpy as np
import pytest

from kmatch import alignment
from kmatch.config import Config
from kmatch.gui import DesktopApi, Job, create_app
from kmatch.models import DependencyStatus, PredictionResult, SearchResult, write_json
from kmatch.process import Cancellation, Cancelled

QXYZ = np.array([[0, 0, 1], [3, 1, 2], [4, 5, 0], [6, 7, 4], [8, 5, 6], [11, 4, 8]])
ROTATION = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
TRANSLATION = np.array([12, -4, 7])
NAMES = dict(A="ALA", C="CYS", D="ASP", E="GLU", F="PHE", G="GLY", K="LYS")


def make_structure(sequence, xyz, full=False):
    structure = gemmi.Structure()
    model, chain = gemmi.Model("1"), gemmi.Chain("A")
    for i, (letter, position) in enumerate(zip(sequence, xyz, strict=True)):
        residue = gemmi.Residue()
        residue.name = NAMES[letter]
        residue.seqid = gemmi.SeqId(i + 10 if full else i + 1, "A" if full and i == 2 else " ")
        residue.het_flag = "A"
        for name, delta in (("N", -1), ("CA", 0), ("C", 1)) if full else (("CA", 0),):
            atom = gemmi.Atom()
            atom.name, atom.element = name, gemmi.Element("N" if name == "N" else "C")
            atom.pos = list(np.asarray(position, dtype=float) + [delta, 0, 0])
            atom.b_iso = 88 if full else 0
            residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    return structure


@pytest.fixture
def case(tmp_path, monkeypatch):
    cfg = Config(tmp_path / "home", tmp_path / "cache", pymol_override=str(tmp_path / "absent"))
    run = cfg.runs / "fixture-run"
    (run / "prediction").mkdir(parents=True)
    (run / "search").mkdir()
    query = make_structure("ACDEFG", QXYZ, full=True)
    target_xyz = np.insert((QXYZ - TRANSLATION) @ ROTATION, 3, [20, 30, 40], axis=0)
    target = make_structure("ACDKEFG", target_xyz)
    source = run / "prediction/predicted.pdb"
    query.write_pdb(str(source))
    prediction = PredictionResult(source, "fixture", "test", 6, {"ptm": 0.9}, [88] * 6, 1)
    write_json(run / "prediction/prediction.json", prediction)
    hit = SearchResult(
        query="query",
        target="1abc-assembly1_B",
        evalue=1e-10,
        probability=1,
        bits=100,
        sequence_identity=6 / 7,
        alignment_length=7,
        query_length=6,
        target_length=7,
        query_coverage=1,
        target_coverage=1,
        query_tm=1,
        target_tm=6 / 7,
        alignment_tm=1,
        lddt=1,
        query_start=1,
        query_end=6,
        target_start=1,
        target_end=7,
        query_alignment="ACD-EFG",
        target_alignment="ACDKEFG",
        rotation=ROTATION.tolist(),
        translation=TRANSLATION.tolist(),
        raw={"u": "0,-1,0,1,0,0,0,0,1", "t": "12,-4,7"},
    )
    metadata = dict(
        database=str(tmp_path / "old-generation/pdb"),
        database_version="test-date",
        foldseek_version="test",
        hits=[asdict(hit)],
        query_structure_sha256=alignment.digest(source),
    )
    write_json(run / "search/search.json", metadata)
    write_json(run / "run.json", {"sequence_identifier": "fixture_query", "status": "complete"})
    monkeypatch.setattr(alignment, "materialize_target", lambda *a: (target.clone(), "42"))
    return cfg, run, query, target, hit, metadata


def test_transform_direction_gaps_insertion_codes_and_originals(case):
    _, _, query, target, hit, _ = case
    before = query.make_pdb_string(), target.make_pdb_string()
    combined, metrics, _ = alignment.superpose(query, target, hit)
    assert metrics["aligned_ca_pairs"] == 6
    assert metrics["aligned_ca_rmsd_before_angstrom"] > 10
    assert metrics["aligned_ca_rmsd_after_angstrom"] < 1e-10
    assert metrics["fixed_transform_query_tm"] == pytest.approx(1)
    assert [c.name for c in combined[0]] == ["Q", "T"]
    assert combined[0]["Q"][2].seqid.icode == "A"
    assert combined[0]["Q"][0][0].b_iso == 88
    assert before == (query.make_pdb_string(), target.make_pdb_string())


@pytest.mark.parametrize(
    "changes",
    [
        {"rotation": [[-1, 0, 0], [0, 1, 0], [0, 0, 1]]},
        {"translation": [float("nan"), 0, 0]},
        {"query_alignment": "ACC-EFG"},
        {"target_end": 6},
        {"target_length": 8},
    ],
)
def test_rejects_invalid_transform_or_coordinate_correspondence(case, changes):
    _, _, query, target, hit, _ = case
    with pytest.raises(ValueError):
        alignment.superpose(query, target, replace(hit, **changes))


@pytest.mark.parametrize("fmt", ["cif", "pdb"])
def test_coordinate_export_roundtrip_provenance_and_repeat(case, fmt):
    cfg, run, _, _, hit, _ = case
    originals = {p: p.read_bytes() for p in run.rglob("*") if p.is_file()}
    first = alignment.export_alignment(cfg, run, 0, fmt)
    contents = first.artifact.read_bytes()
    second = alignment.export_alignment(cfg, run, 0, fmt)
    assert first.directory != second.directory and first.artifact.read_bytes() == contents
    assert all(path.read_bytes() == data for path, data in originals.items())
    assert not first.directory.joinpath(".work").exists()
    structure = gemmi.read_structure(str(first.artifact))
    assert len(structure) == 1 and [c.name for c in structure[0]] == ["Q", "T"]
    pairs = alignment.aligned_pairs(hit, "ACDEFG", "ACDKEFG")
    q = np.array([list(r["CA"][0].pos) for r in structure[0]["Q"]])
    t = np.array([list(r["CA"][0].pos) for r in structure[0]["T"]])
    assert np.allclose(q[pairs[:, 0]], t[pairs[:, 1]], atol=0.001)
    if fmt == "cif":
        block = gemmi.cif.read(str(first.artifact)).sole_block()
        provenance = json.loads(
            gemmi.cif.as_string(block.find_value("_kmatch_alignment.provenance_json"))
        )
    else:
        provenance = json.loads(
            "".join(
                line[11:]
                for line in first.artifact.read_text().splitlines()
                if line.startswith("REMARK 951 ")
            )
        )
    assert provenance["hit"]["target"] == hit.target and provenance["hit"]["raw"] == hit.raw
    assert provenance["refitted"] is False and provenance["query_hash_recorded_at_search"]


def test_legacy_runs_and_af3_terms_remain_exportable(case):
    cfg, run, _, _, _, metadata = case
    del metadata["query_structure_sha256"]
    write_json(run / "search/search.json", metadata)
    path = run / "prediction/prediction.json"
    prediction = json.loads(path.read_text())
    prediction["predictor"] = "af3"
    write_json(path, prediction)
    result = alignment.export_alignment(cfg, run, 0, "cif")
    text = result.artifact.read_text()
    assert "OUTPUT_TERMS_OF_USE" in text and '"query_hash_recorded_at_search":false' in text


def test_moved_run_reads_its_own_prediction_not_the_recorded_old_root(case, tmp_path):
    cfg, run, _, _, _, _ = case
    path = run / "prediction/prediction.json"
    prediction = json.loads(path.read_text())
    # An adopted run still records its original absolute location in the legacy home.
    legacy = tmp_path / "legacy-home/runs/fixture-run/prediction/predicted.pdb"
    legacy.parent.mkdir(parents=True)
    legacy.write_text("REMARK legacy coordinates must never be read\n")
    prediction["structure_path"] = str(legacy)
    write_json(path, prediction)
    result = alignment.export_alignment(cfg, run, 0, "cif")
    provenance = json.loads((result.directory / "provenance.json").read_text())
    assert provenance["query_structure"] == str(run / "prediction/predicted.pdb")
    assert result.metrics["aligned_ca_pairs"] == 6
    legacy.unlink()
    assert alignment.export_alignment(cfg, run, 0, "pdb").artifact.is_file()
    # The search-time SHA-256 still protects the run's own copy.
    source = run / "prediction/predicted.pdb"
    source.write_text(source.read_text() + "REMARK changed\n")
    with pytest.raises(ValueError, match="coordinates changed"):
        alignment.export_alignment(cfg, run, 0, "cif")


def test_recorded_structure_path_is_kept_inside_run_or_without_own_copy(tmp_path):
    run = tmp_path / "run"
    inside = run / "prediction/predicted.pdb"
    assert alignment.run_structure(run, inside) == inside
    outside = tmp_path / "elsewhere/predicted.pdb"
    # Without a same-named run copy the recorded path is used (and fails loudly if missing).
    assert alignment.run_structure(run, outside) == outside
    inside.parent.mkdir(parents=True)
    inside.write_text("MODEL\n")
    assert alignment.run_structure(run, str(outside)) == inside


def test_changed_query_or_stale_selection_is_rejected(case):
    cfg, run, _, _, _, _ = case
    with pytest.raises(ValueError, match="changed"):
        alignment.export_alignment(cfg, run, 0, "cif", expected_target="different")
    source = run / "prediction/predicted.pdb"
    source.write_text(source.read_text() + "REMARK changed\n")
    with pytest.raises(ValueError, match="coordinates changed"):
        alignment.export_alignment(cfg, run, 0, "cif")
    assert not (run / "exports").exists()


def test_cancellation_cleans_staging_and_preserves_existing_destination(case, tmp_path):
    cfg, run, _, _, _, _ = case
    target = tmp_path / "previous.cif"
    target.write_bytes(b"previous export")
    cancel = Cancellation()

    def progress(event):
        if event.message.startswith("Applying"):
            cancel.cancel()

    with pytest.raises(Cancelled):
        alignment.export_alignment(
            cfg, run, 0, "cif", cancel=cancel, progress=progress, destination=target
        )
    assert target.read_bytes() == b"previous export"
    folder = next((run / "exports").iterdir())
    assert not (folder / ".work").exists() and not (folder / "alignment.cif").exists()
    assert json.loads((folder / "export.json").read_text())["status"] == "cancelled"
    assert alignment.export_alignment(cfg, run, 0, "cif").artifact.exists()


def test_atomic_save_and_managed_original_protection(case, tmp_path):
    cfg, run, _, _, _, _ = case
    result = alignment.export_alignment(cfg, run, 0, "pdb")
    source = run / "prediction/predicted.pdb"
    original = source.read_bytes()
    with pytest.raises(ValueError, match="managed"):
        alignment.save_export(result, source, cfg, Cancellation())
    assert source.read_bytes() == original
    destination = tmp_path / "existing.pdb"
    destination.write_bytes(b"previous export")
    cancel = Cancellation()
    cancel.cancel()
    with pytest.raises(Cancelled):
        alignment.save_export(result, destination, cfg, cancel)
    assert destination.read_bytes() == b"previous export"
    alignment.save_export(result, destination, cfg, Cancellation())
    assert destination.read_bytes() == result.artifact.read_bytes()
    assert not list(tmp_path.glob(".kmatch-export-*"))


def test_missing_or_failing_pymol_has_portable_fallback(case, tmp_path):
    cfg, run, _, _, _, _ = case
    assert not alignment.export_support(cfg)["pse"]
    with pytest.raises(RuntimeError, match="PyMOL is unavailable"):
        alignment.export_alignment(cfg, run, 0, "pse")
    program = tmp_path / "pymol"
    program.write_text("#!/bin/sh\nexit 7\n")
    program.chmod(0o755)
    with pytest.raises(RuntimeError, match="Use mmCIF/PDB"):
        alignment.export_alignment(replace(cfg, pymol_override=str(program)), run, 0, "pse")
    folder = next((run / "exports").iterdir())
    assert (folder / "export.log").is_file() and not (folder / ".work").exists()
    assert not (folder / "alignment.pse").exists()
    assert alignment.export_alignment(cfg, run, 0, "cif").artifact.exists()


def test_materialization_pins_original_generation(tmp_path, monkeypatch):
    # Exercise the real selection/command boundary with disposable database fixtures.
    cfg = Config(tmp_path, tmp_path / "cache", pdb_prefix_override=str(tmp_path / "new/pdb"))
    prefix = tmp_path / "old/pdb"
    prefix.parent.mkdir()
    Path(str(prefix) + ".lookup").touch()
    monkeypatch.setattr(
        alignment, "database_status", lambda p: DependencyStatus("PDB100", True, "ready", "old")
    )
    target = make_structure("ACDEFG", QXYZ)
    commands = []

    def process(command, *args, **kwargs):
        commands.append(command)
        if command[1] == "createsubdb":
            assert command[3] == prefix and command[6] == "1"
            Path(str(command[4]) + ".index").write_text("42\t0\t1\n")
        else:
            target.make_mmcif_document().write_file(str(command[3]))

    monkeypatch.setattr(alignment, "run_process", process)
    from types import SimpleNamespace

    work = tmp_path / "work"
    work.mkdir()
    result, key = alignment.materialize_target(
        cfg,
        dict(database=str(prefix), database_version="old"),
        SimpleNamespace(target="1abc_A"),
        work,
        tmp_path / "log",
        Cancellation(),
        lambda e: None,
    )
    assert key == "42" and len(result[0][0]) == 6 and len(commands) == 2
    Path(str(prefix) + ".lookup").unlink()
    with pytest.raises(RuntimeError, match="legacy hit has no retained target trace"):
        alignment.materialize_target(
            cfg,
            dict(database=str(prefix)),
            SimpleNamespace(target="1abc_A"),
            work,
            tmp_path / "log",
            Cancellation(),
            lambda e: None,
        )


def test_gui_exports_off_thread_recovers_and_rejects_stale_requests(case, monkeypatch):
    cfg, run, _, _, hit, _ = case
    job = Job(cfg)
    job.update(
        directory=str(run),
        run_id=run.name,
        hits=[asdict(hit)],
        prediction=json.loads((run / "prediction/prediction.json").read_text()),
    )
    client = create_app(job).test_client()
    token = re.search('name="kmatch-token" content="([^"]+)"', client.get("/").text)[1]
    headers = {"X-Kmatch-Token": token}
    assert (
        client.post(
            "/api/alignment-export",
            headers=headers,
            json={"run_id": "previous", "hit_index": 0, "format": "cif"},
        ).status_code
        == 400
    )
    entered = threading.Event()

    def waiting(*a, **kwargs):
        entered.set()
        while True:
            kwargs["cancel"].check()
            time.sleep(0.01)

    monkeypatch.setattr("kmatch.gui.export_alignment", waiting)
    body = {"run_id": run.name, "hit_index": 0, "format": "cif"}
    assert client.post("/api/alignment-export", headers=headers, json=body).status_code == 200
    assert entered.wait(1)
    start = time.monotonic()
    assert client.get("/api/status", headers=headers).json["job"]["busy"]
    assert time.monotonic() - start < 0.5
    client.post("/api/cancel", headers=headers, json={})
    job.thread.join(2)
    assert job.snapshot()["phase"] == "cancelled" and job.snapshot()["hits"] == [asdict(hit)]
    monkeypatch.setattr("kmatch.gui.export_alignment", alignment.export_alignment)
    client.post("/api/alignment-export", headers=headers, json=body)
    job.thread.join(2)
    state = job.snapshot()
    assert not state["busy"] and state["phase"] == "complete"
    identifier = state["alignment_export"]["export_id"]
    assert client.get("/api/alignment-export/not-current", headers=headers).status_code == 400
    response = client.get("/api/alignment-export/" + identifier, headers=headers)
    assert response.status_code == 200 and b"_atom_site." in response.data


def test_native_save_dialog_cancel_starts_no_worker(case):
    cfg, run, _, _, hit, _ = case
    job = Job(cfg)
    job.update(directory=str(run), run_id=run.name, hits=[asdict(hit)], prediction={"ready": True})
    api = DesktopApi(job)

    class Window:
        def create_file_dialog(self, *args, **kwargs):
            return None

    api._window = Window()
    assert api.export_alignment(run.name, 0, "cif") is False
    assert job.thread is None and not (run / "exports").exists()
