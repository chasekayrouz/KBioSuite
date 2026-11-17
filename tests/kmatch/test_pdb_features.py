"""Public metadata, frame reconciliation, ligand filtering and residue correspondence."""

import json
from dataclasses import asdict
from pathlib import Path

import gemmi
import numpy as np
import pytest
from test_alignment import case as case
from test_alignment import make_structure

from kmatch import alignment, hit_structure, pdb_data
from kmatch.config import Config
from kmatch.gui import Job, create_app
from kmatch.process import Cancellation, Cancelled


def entry():
    return {
        "rcsb_id": "1ABC",
        "struct": {"title": "A public structure"},
        "polymer_entities": [
            {
                "rcsb_polymer_entity": {"pdbx_description": "Enzyme alpha"},
                "rcsb_polymer_entity_container_identifiers": {
                    "auth_asym_ids": ["B"],
                    "asym_ids": ["X"],
                    "entity_id": "1",
                },
            }
        ],
        "nonpolymer_entities": [
            {
                "pdbx_entity_nonpoly": {"comp_id": code, "name": name},
                "rcsb_nonpolymer_entity_container_identifiers": {
                    "auth_asym_ids": ["B"],
                    "asym_ids": ["Y"],
                    "entity_id": "2",
                },
            }
            for code, name in [("ZN", "Zinc ion"), ("HOH", "Water"), ("ATP", "ATP")]
        ],
    }


def test_metadata_name_chain_ligands_cache_and_public_only_requests(tmp_path, monkeypatch):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    calls = []

    def fetch(url, *args):
        calls.append(url)
        return json.dumps({"data": {"entries": [entry()]}}).encode()

    monkeypatch.setattr(pdb_data, "fetch", fetch)
    target = "1abc-assembly1_B"
    result = pdb_data.load_details(cfg, [target], True, Cancellation())[target]
    assert result["protein_name"] == "Enzyme alpha"
    assert [ligand["id"] for ligand in result["ligands"]] == ["ZN", "ATP"]
    assert result["url"] == "https://www.rcsb.org/structure/1ABC"
    assert pdb_data.load_details(cfg, [target], False, Cancellation())[target] == result
    assert len(calls) == 1 and "1ABC" in calls[0] and "sequence" not in calls[0]
    assert (
        pdb_data.entry_for_hit(pdb_data.normalize_entry(entry()), "1abc-assembly1_B-2")[
            "protein_name"
        ]
        == "Enzyme alpha"
    )
    assert (
        pdb_data.entry_for_hit(pdb_data.normalize_entry(entry()), "1abc-assembly1_C")[
            "protein_name"
        ]
        is None
    )


@pytest.mark.parametrize(
    "target", ["../secret", "1abc_A?x=secret", "AF-Q1234", "https://evil/1abc_A"]
)
def test_identifiers_cannot_supply_arbitrary_urls(target):
    with pytest.raises(ValueError):
        pdb_data.pdb_url(target)


def test_network_failure_is_unavailable_not_no_ligands(tmp_path, monkeypatch):
    cfg = Config(tmp_path, tmp_path / "cache")

    def fail(*args):
        raise OSError("offline")

    monkeypatch.setattr(pdb_data, "fetch", fail)
    result = pdb_data.load_details(cfg, ["1abc_A"], True, Cancellation())["1abc_A"]
    assert result["status"] == "unavailable" and "offline" in result["message"]
    assert "ligands" not in result
    cancel = Cancellation()
    cancel.cancel()
    with pytest.raises(Cancelled):
        pdb_data.load_details(cfg, ["1abc_A"], False, cancel)
    with pytest.raises(RuntimeError, match="not cached"):
        pdb_data.assembly_file(cfg, "1abc-assembly1_A", False, Cancellation(), lambda e: None)


def test_invalid_coordinate_response_is_not_cached(tmp_path, monkeypatch):
    cfg = Config(tmp_path, tmp_path / "cache")
    monkeypatch.setattr(pdb_data, "fetch", lambda *a: b"data_wrong\n_entry.id 9XYZ\n")
    with pytest.raises(RuntimeError, match="does not identify"):
        pdb_data.assembly_file(cfg, "1abc-assembly1_A", True, Cancellation(), lambda e: None)
    assert not list(cfg.pdb_cache.rglob("*.cif"))


def add_ligand(model, name, element, position, entity_type, number):
    chain = gemmi.Chain("B")
    residue = gemmi.Residue()
    residue.name = name
    residue.seqid = gemmi.SeqId(number, " ")
    residue.het_flag = "H"
    residue.entity_type = entity_type
    residue.subchain = name
    atom = gemmi.Atom()
    atom.name = element
    atom.element = gemmi.Element(element)
    atom.pos = list(position)
    atom.occ = 1
    residue.add_atom(atom)
    chain.add_residue(residue)
    model.add_chain(chain)


@pytest.fixture
def full_case(case, monkeypatch, tmp_path):
    cfg, run, query, trace, hit, _ = case
    residues = hit_structure.ca_residues(trace[0][0])
    # Independent source→database transform in addition to nonidentity Foldseek u/t.
    source_u = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    source_t = np.array([40, 25, -15])
    source_xyz = (hit_structure.xyz(residues) - source_t) @ source_u
    full = make_structure("ACDKEFG", source_xyz, full=True)
    full[0][0].name = "B"
    for residue in full[0][0]:
        residue.entity_type = gemmi.EntityType.Polymer
        residue.subchain = "B"
    # Close to the first target residue only; the alignment contains an insertion later.
    contact = source_xyz[0] + [0, 0, 2]
    add_ligand(full[0], "ZN", "Zn", contact, gemmi.EntityType.NonPolymer, 200)
    add_ligand(full[0], "HOH", "O", contact, gemmi.EntityType.Water, 201)
    add_ligand(full[0], "NAG", "C", contact, gemmi.EntityType.Branched, 202)
    add_ligand(full[0], "SO4", "S", [999, 999, 999], gemmi.EntityType.NonPolymer, 203)
    # Supply parsed source to isolate scientific logic; fixture never changes normal cached data.
    path = tmp_path / "assembly.cif"
    path.touch()
    monkeypatch.setattr(hit_structure, "assembly_file", lambda *a: (path, {"sha256": "fixture"}))
    original_reader = gemmi.read_structure
    monkeypatch.setattr(
        gemmi, "read_structure", lambda p: full.clone() if Path(p) == path else original_reader(p)
    )
    return cfg, run, query, trace, hit, full


def test_full_atom_frame_and_ligand_contacts_follow_foldseek_alignment(full_case):
    cfg, _, query, trace, hit, full = full_case
    original = full.make_pdb_string()
    combined, _, _ = alignment.superpose(query, trace, hit)
    info = hit_structure.full_atom_superposition(
        cfg, combined, trace, hit, False, Cancellation(), lambda e: None
    )
    assert info["registration"]["ca_rmsd_angstrom"] < 1e-9
    assert [x["component_id"] for x in info["ligands"]] == ["ZN"]
    assert info["aligned_query_contacts"][0]["query_position"] == 1
    assert info["aligned_query_contacts"][0]["database_target_position"] == 1
    assert info["aligned_query_contacts"][0]["hit_residue"]["number"] == 10
    assert [c.name for c in combined[0]] == ["Q", "T", "L"]
    q = hit_structure.xyz(hit_structure.ca_residues(combined[0]["Q"]))
    t = hit_structure.xyz(hit_structure.ca_residues(combined[0]["T"]))
    pairs = alignment.aligned_pairs(hit, "ACDEFG", "ACDKEFG")
    assert np.allclose(q[pairs[:, 0]], t[pairs[:, 1]])
    assert full.make_pdb_string() == original


@pytest.mark.parametrize("fmt", ["pdb", "cif"])
def test_full_atom_coordinate_export_preserves_ligands_and_provenance(full_case, fmt):
    cfg, run, _, _, _, _ = full_case
    result = alignment.export_alignment(cfg, run, 0, fmt, full_atom=True)
    structure = gemmi.read_structure(str(result.artifact))
    assert [c.name for c in structure[0]] == ["Q", "T", "L"]
    assert structure[0]["L"][0].name == "ZN"
    info = json.loads((result.directory / "provenance.json").read_text())
    assert info["full_atom"]["contact_cutoff_angstrom"] == 4
    assert info["source_chain_names"]["T"] == "B"
    assert info["refitted"] is False


def test_wrong_conformation_rejected_instead_of_moving_ligands(full_case):
    cfg, _, query, trace, hit, full = full_case
    full[0][0][2]["CA"][0].pos.x += 9
    combined, _, _ = alignment.superpose(query, trace, hit)
    with pytest.raises(RuntimeError, match="does not reproduce"):
        hit_structure.full_atom_superposition(
            cfg, combined, trace, hit, False, Cancellation(), lambda e: None
        )


def test_mapping_extra_residues_missing_residues_and_alternate_conformer():
    assert hit_structure.trace_mapping("ACDE", "XACDE") == [1, 2, 3, 4]
    with pytest.raises(ValueError):
        hit_structure.trace_mapping("ACDE", "ACE")
    with pytest.raises(ValueError):
        hit_structure.trace_mapping("ACDE", "ACKE")
    r = gemmi.Residue()
    for label, occ in [("A", 0.2), ("B", 0.8)]:
        a = gemmi.Atom()
        a.name = "CA"
        a.element = gemmi.Element("C")
        a.altloc = label
        a.occ = occ
        r.add_atom(a)
    best = hit_structure.major_conformer(r)
    assert len(best) == 1 and best[0].occ == pytest.approx(0.8) and not best[0].has_altloc()


def test_metadata_job_is_independent_and_error_does_not_destroy_hits(case, monkeypatch):
    cfg, run, _, _, hit, _ = case
    job = Job(cfg)
    job.update(run_id=run.name, directory=str(run), prediction={"ready": True}, hits=[asdict(hit)])

    def fail(*args):
        raise OSError("fixture offline")

    monkeypatch.setattr("kmatch.gui.load_details", fail)
    job.load_pdb_details(run.name, False)
    job.details_thread.join(2)
    assert not job.snapshot()["busy"] and len(job.snapshot()["hits"]) == 1
    assert "fixture offline" in job.snapshot()["pdb_details_status"]
    with pytest.raises(ValueError):
        job.load_pdb_details("stale", True)
    assert (
        create_app(job)
        .test_client()
        .post("/api/pdb-details", json={"run_id": run.name})
        .status_code
        == 403
    )


def test_cancellation_during_network_error_remains_cancellation(monkeypatch):
    cancel = Cancellation()

    def timeout(*args, **kwargs):
        cancel.cancel()
        raise TimeoutError("fixture request timed out")

    monkeypatch.setattr(pdb_data.urllib.request, "urlopen", timeout)
    with pytest.raises(Cancelled):
        pdb_data.fetch("https://data.rcsb.org/graphql", cancel, 1024)


def test_query_contact_mapping_accounts_for_target_insertion(full_case):
    cfg, _, query, trace, hit, full = full_case
    position = list(full[0][0][4]["CA"][0].pos)
    position[2] += 1
    add_ligand(full[0], "LIG", "C", position, gemmi.EntityType.NonPolymer, 204)
    combined, _, _ = alignment.superpose(query, trace, hit)
    info = hit_structure.full_atom_superposition(
        cfg, combined, trace, hit, False, Cancellation(), lambda e: None
    )
    pair = next(p for p in info["aligned_query_contacts"] if p["database_target_position"] == 5)
    assert pair["query_position"] == 4
    assert pair["query_residue"]["number"] == 13
    assert pair["hit_residue"]["number"] == 14


def test_native_export_forwards_full_atom_and_network_options(case, tmp_path, monkeypatch):
    from kmatch.gui import DesktopApi

    cfg, run, _, _, hit, _ = case
    job = Job(cfg)
    job.update(run_id=run.name, directory=str(run), prediction={"ready": True}, hits=[asdict(hit)])
    api = DesktopApi(job)

    class Window:
        def create_file_dialog(self, *args, **kwargs):
            return str(tmp_path / "alignment.pse")

    api._window = Window()
    calls = []
    monkeypatch.setattr(job, "export_alignment", lambda *a, **kw: calls.append((a, kw)))
    assert api.export_alignment(run.name, 0, "pse", True, False)
    assert calls[0][1] == {"full_atom": True, "allow_network": False}


def test_backbone_alternate_follows_database_not_occupancy(full_case):
    cfg, _, query, trace, hit, full = full_case
    residue = full[0][0][2]
    original = residue["CA"][0]
    original.altloc = "A"
    original.occ = 0.2
    alternate = original.clone()
    alternate.altloc = "B"
    alternate.occ = 0.8
    alternate.pos.x += 2
    residue.add_atom(alternate)
    combined, _, _ = alignment.superpose(query, trace, hit)
    info = hit_structure.full_atom_superposition(
        cfg, combined, trace, hit, False, Cancellation(), lambda e: None
    )
    assert info["registration"]["ca_rmsd_angstrom"] < 1e-8
    assert len([a for a in combined[0]["T"][2] if a.name == "CA"]) == 1
