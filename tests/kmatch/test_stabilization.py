import json
import shutil
import struct
import threading
import time
from dataclasses import replace
from pathlib import Path

import pytest
from test_alignment import case as case
from test_pdb_features import full_case as full_case

from kmatch import alignment, setup, target_traces
from kmatch.confidence import boltz_warnings
from kmatch.config import Config
from kmatch.database import REQUIRED, database_lease, inspect_database
from kmatch.process import Cancellation, Cancelled
from kmatch.search import hit_rankings, retained_hit_indices


def ranked_fixture(hit):
    # Opposing qTM and significance order; 20 overlapping candidates and 30 new ones.
    return [
        replace(hit, target=f"1abc_{i}", query_tm=1 - i / 100, evalue=1 / (i + 1))
        for i in range(80)
    ]


def test_rankings_independently_select_top_50_and_union(full_case):
    hits = ranked_fixture(full_case[4])
    rankings = hit_rankings(hits)
    assert rankings["qtm"] == list(range(50))
    assert rankings["significance"] == list(range(79, 29, -1))
    assert retained_hit_indices(rankings) == list(range(80))
    assert hit_rankings([]) == {"qtm": [], "significance": []}
    assert hit_rankings(hits[:3])["significance"] == [2, 1, 0]
    disjoint = ranked_fixture(full_case[4]) + [
        replace(full_case[4], target=f"1abc_{i}", query_tm=0, evalue=1e-30) for i in range(80, 110)
    ]
    assert len(retained_hit_indices(hit_rankings(disjoint))) == 100
    # Probability precedes E-value; bit score breaks equal probability/E-value ties.
    tied = [
        replace(hits[0], probability=0.9, evalue=1e-30),
        replace(hits[0], probability=1, evalue=1e-10, bits=10),
        replace(hits[0], probability=1, evalue=1e-10, bits=20),
        replace(hits[0], probability=1, evalue=1e-10, bits=20),
    ]
    assert hit_rankings(tied)["significance"] == [2, 3, 1, 0]


def test_search_retains_union_and_keeps_full_saved_indices(full_case, monkeypatch):
    from kmatch.models import DependencyStatus, PredictionResult
    from kmatch.search import FoldseekBackend

    cfg, run, _, _, hit, _ = full_case
    hits = ranked_fixture(hit)
    backend = FoldseekBackend(cfg)
    monkeypatch.setattr(backend, "status", lambda: DependencyStatus("test", True, "ready"))
    monkeypatch.setattr("kmatch.search.run_process", lambda *a, **k: 0.1)
    monkeypatch.setattr("kmatch.search.parse_foldseek", lambda *a: hits)
    retained = []

    def retain(config, database, candidates, *args):
        retained.extend(candidates)
        return {h.target: {} for h in candidates}

    monkeypatch.setattr("kmatch.search.retain_target_traces", retain)
    p = json.loads((run / "prediction/prediction.json").read_text())
    p["structure_path"] = Path(p["structure_path"])
    output = run / "new-search"
    returned = backend.search(PredictionResult(**p), output, Cancellation(), lambda e: None)
    saved = json.loads((output / "search.json").read_text())
    assert returned == hits and retained == hits
    assert saved["retained_target_count"] == 80
    assert saved["hit_rankings"]["significance"][0] == 79
    assert saved["hits"][79]["target"] == returned[79].target


def test_gui_annotations_use_candidate_union_and_export_keeps_saved_index(full_case, monkeypatch):
    from dataclasses import asdict
    from types import SimpleNamespace

    from kmatch.gui import Job

    cfg, run, _, _, hit, _ = full_case
    hits = ranked_fixture(hit) + [replace(hit, target="1abc_outside", query_tm=0, evalue=1)]
    job = Job(cfg)
    job.update(
        directory=str(run),
        run_id=run.name,
        prediction={"ready": True},
        hits=[asdict(h) for h in hits],
        hit_rankings=hit_rankings(hits),
    )
    requested = []
    monkeypatch.setattr(
        "kmatch.gui.load_details",
        lambda config, targets, *args: requested.extend(targets) or {},
    )
    job.load_pdb_details(run.name, False)
    job.details_thread.join(3)
    assert requested == [h.target for h in hits[:80]]
    exported = []

    def export(config, directory, index, format, **kwargs):
        exported.append((index, kwargs["expected_target"]))
        raise RuntimeError("stop after verifying export arguments")

    monkeypatch.setattr("kmatch.gui.export_alignment", export)
    job.export_alignment(run.name, 79, "cif")
    job.thread.join(3)
    assert exported == [(79, hits[79].target)]
    assert not job.snapshot()["busy"]
    # Production Job publishes rankings for all results, rather than truncating its hit list.
    monkeypatch.setattr(
        "kmatch.gui.Pipeline",
        lambda config: SimpleNamespace(
            run=lambda *a, **k: SimpleNamespace(prediction=hit, hits=hits)
        ),
    )
    job.run("ACDEFGHIKLMNPQRSTVWY", "boltz2", False)
    job.thread.join(3)
    assert len(job.snapshot()["hits"]) == 81
    assert job.snapshot()["hit_rankings"] == hit_rankings(hits)


@pytest.mark.parametrize(
    ("plddt", "ptm", "warn"),
    [
        (58.4, 0.439, True),
        (69.9, 0.8, True),
        (92, 0.49, True),
        (70, 0.5, False),
        (91.3, 0.705, False),
        (95.3, 0.877, False),
    ],
)
def test_confidence_warning_boundary_and_msa_scope(plddt, ptm, warn):
    confidence = dict(mean_plddt=plddt, ptm=ptm)
    assert bool(boltz_warnings("boltz2", confidence, {"msa": "single-sequence"})) is warn
    assert not boltz_warnings("boltz2", confidence, {"msa": "custom-a3m"})
    assert not boltz_warnings("af3", confidence, {"msa": "single-sequence"})
    assert not boltz_warnings("boltz2", confidence, {})


def retain_fixture(full_case, monkeypatch):
    cfg, run, query, trace, hit, full = full_case
    database = cfg.home / "old/pdb"
    database.parent.mkdir()
    Path(str(database) + ".lookup").write_text("42\t" + hit.target + "\t0\n")

    def extract(command, *args, **kwargs):
        if command[1] == "createsubdb":
            Path(str(command[4]) + ".index").write_text("42\t0\t1\n")
            # Model the external links which must NEVER be retained with a saved run.
            Path(str(command[4]) + ".lookup").symlink_to(Path(str(database) + ".lookup"))
        else:
            destination = Path(command[3])
            destination.mkdir()
            trace.make_mmcif_document().write_file(str(destination / (hit.target + ".cif")))

    monkeypatch.setattr(target_traces, "run_process", extract)
    records = target_traces.retain_target_traces(
        cfg, database, [hit], run / "search", Cancellation(), lambda e: None
    )
    metadata_path = run / "search/search.json"
    metadata = json.loads(metadata_path.read_text())
    metadata.update(target_traces=records, database=str(database), schema_version=2)
    metadata_path.write_text(json.dumps(metadata))
    shutil.rmtree(database.parent)
    return records, metadata


def test_saved_trace_exports_after_database_and_foldseek_removed(full_case, monkeypatch):
    cfg, run, _, _, hit, _ = full_case
    records, metadata = retain_fixture(full_case, monkeypatch)
    assert len(records) == 1 and not list((run / "search").glob(".traces-*"))
    assert not any(p.is_symlink() for p in (run / "search/target-traces").rglob("*"))

    def unavailable(*args):
        raise AssertionError("An export with retained coordinates must not need the database")

    monkeypatch.setattr(alignment, "materialize_target", unavailable)
    for full_atom in (False, True):
        result = alignment.export_alignment(cfg, run, 0, "cif", full_atom=full_atom)
        assert result.artifact.is_file()
        info = json.loads((result.directory / "provenance.json").read_text())
        assert info["target_coordinate_source"] == "run-retained trace"
        assert info["refitted"] is False
        if full_atom:
            assert info["full_atom"]["ligands"][0]["component_id"] == "ZN"
    saved, key = target_traces.load_target_trace(run / "search", metadata, hit)
    assert key == "42" and len(saved) == 1
    path = run / "search" / records[hit.target]["file"]
    path.write_text(path.read_text() + "\n# changed\n")
    with pytest.raises(RuntimeError, match="missing or changed"):
        alignment.export_alignment(cfg, run, 0, "cif")


def test_trace_paths_cannot_escape_run(full_case):
    _, run, _, _, hit, _ = full_case
    metadata = {"target_traces": {hit.target: {"file": "../../elsewhere.cif"}}}
    with pytest.raises(ValueError, match="outside this run"):
        target_traces.load_target_trace(run / "search", metadata, hit)


def test_unknown_database_residue_is_retained_as_x(full_case, monkeypatch):
    from kmatch.hit_structure import sequence

    cfg, run, query, trace, hit, full = full_case
    trace[0][0][1].name = "UNK"
    changed = replace(hit, target_alignment=hit.target_alignment.replace("C", "X"))
    assert alignment.ca_chain(trace, "test")[1][1] == "X"
    assert sequence(list(trace[0][0]))[1] == "X"
    records, _ = retain_fixture((cfg, run, query, trace, changed, full), monkeypatch)
    assert changed.target in records


def test_cancelled_retention_does_not_publish_partial_snapshot(full_case, monkeypatch):
    cfg, run, _, _, hit, _ = full_case
    monkeypatch.setattr(
        target_traces, "run_process", lambda *a, **k: (_ for _ in ()).throw(Cancelled())
    )
    with pytest.raises(Cancelled):
        target_traces.retain_target_traces(
            cfg, Path("not-used"), [hit], run / "search", Cancellation(), lambda e: None
        )
    assert not (run / "search/target-traces").exists()
    assert not list((run / "search").glob(".traces-*"))


def generation(config, date, checksum):
    path = config.home / "databases/pdb100" / f"{date}-{checksum[:12]}"
    path.mkdir(parents=True)
    for suffix in REQUIRED:
        dest = path / ("pdb" + suffix)
        if suffix.endswith(".dbtype"):
            dest.write_bytes(
                struct.pack(
                    "<I",
                    {".dbtype": 0, "_ss.dbtype": 0, "_ca.dbtype": 101, "_h.dbtype": 12}[suffix],
                )
            )
        elif suffix.endswith(".index"):
            dest.write_text("42\t0\t4\n")
        elif suffix == ".lookup":
            dest.write_text("42\t1abc_A\t0\n")
        else:
            dest.write_bytes(b"ABC\0")
    (path / "installation.json").write_text(
        json.dumps(
            dict(status="complete", pdb_date=date, published_md5=checksum, source=setup.PDB_URL)
        )
    )
    return path


def test_index_validation_checks_bounds_and_matching_keys(tmp_path):
    cfg = Config(tmp_path, tmp_path / "cache")
    path = generation(cfg, "250101", "a" * 32)
    assert inspect_database(path / "pdb", Cancellation()) == {"entries": 1, "sample_key": 42}
    (path / "pdb_ca.index").write_text("42\t1\t4\n")
    with pytest.raises(RuntimeError, match="inconsistent"):
        inspect_database(path / "pdb", Cancellation())
    (path / "pdb_ca.index").write_text("43\t0\t4\n")
    with pytest.raises(RuntimeError, match="same hits"):
        inspect_database(path / "pdb", Cancellation())


def test_activation_waits_for_readers_then_retires_only_managed_generation(tmp_path):
    cfg = Config(tmp_path, tmp_path / "cache")
    old = generation(cfg, "250101", "a" * 32)
    new = generation(cfg, "260101", "b" * 32)
    root = old.parent
    (root / "current").symlink_to(old.name)
    external = tmp_path / "external"
    external.mkdir()
    (external / "original").write_text("preserve")
    (root / "linked").symlink_to(external)
    (root / "unmanaged").mkdir()
    cancel = Cancellation()
    results = []
    with database_lease(cfg, cancel):
        worker = threading.Thread(
            target=lambda: results.append(
                setup.activate_pdb_generation(cfg, new, cancel, lambda e: None)
            )
        )
        worker.start()
        time.sleep(0.2)
        assert worker.is_alive() and (root / "current").resolve() == old and old.exists()
    worker.join(3)
    assert not worker.is_alive() and results == [[old.name]]
    assert (root / "current").resolve() == new and not old.exists()
    assert (external / "original").read_text() == "preserve" and (root / "unmanaged").is_dir()


def test_cancelled_activation_and_failed_stage_preserve_active_database(tmp_path, monkeypatch):
    cfg = Config(tmp_path, tmp_path / "cache")
    old = generation(cfg, "250101", "a" * 32)
    root = old.parent
    (root / "current").symlink_to(old.name)
    new = generation(cfg, "260101", "b" * 32)
    cancel = Cancellation()
    cancel.cancel()
    with pytest.raises(Cancelled):
        setup.activate_pdb_generation(cfg, new, cancel, lambda e: None)
    assert (root / "current").resolve() == old and old.exists()
    shutil.rmtree(new)
    monkeypatch.setattr(
        setup.subprocess,
        "check_output",
        lambda *a, **k: "b" * 32 + "  pdb100.tar.gz\n260101\tPDB_DATE\n",
    )
    monkeypatch.setattr(setup, "download", lambda *a, **k: None)
    monkeypatch.setattr(setup, "extract", lambda *a, **k: None)

    def fail(*args):
        raise RuntimeError("fixture validation failure")

    monkeypatch.setattr(setup, "validate_pdb_generation", fail)
    with pytest.raises(RuntimeError, match="fixture validation failure"):
        setup.install_pdb(cfg, Cancellation(), lambda e: None)
    assert (root / "current").resolve() == old and old.exists() and not new.exists()
    with pytest.raises(RuntimeError, match="external override"):
        setup.install_pdb(
            replace(cfg, pdb_prefix_override="/external/pdb"), Cancellation(), lambda e: None
        )


def test_staged_update_validates_before_promotion_then_removes_old(tmp_path, monkeypatch):
    cfg = Config(tmp_path, tmp_path / "cache")
    old = generation(cfg, "250101", "a" * 32)
    root = old.parent
    (root / "current").symlink_to(old.name)
    payload = generation(Config(tmp_path / "fixture", tmp_path / "cache"), "260101", "b" * 32)
    monkeypatch.setattr(
        setup.subprocess,
        "check_output",
        lambda *a, **k: "b" * 32 + "  pdb100.tar.gz\n260101\tPDB_DATE\n",
    )
    monkeypatch.setattr(setup, "download", lambda *a, **k: None)

    def extract(archive, stage, *args):
        assert stage.name.startswith(".staging-") and (root / "current").resolve() == old
        for path in payload.glob("pdb*"):
            shutil.copy2(path, stage / path.name)

    def validate(config, stage, cancel, progress, checksum):
        assert (root / "current").resolve() == old and old.exists()
        assert checksum == "b" * 32 and stage.name.startswith(".staging-")
        return inspect_database(stage / "pdb", cancel)

    monkeypatch.setattr(setup, "extract", extract)
    monkeypatch.setattr(setup, "validate_pdb_generation", validate)
    setup.install_pdb(cfg, Cancellation(), lambda e: None)
    assert not old.exists()
    assert (root / "current").resolve().name == "260101-" + "b" * 12
    assert not list(root.glob(".staging-*"))
    assert json.loads((root / "last-update.json").read_text())["removed_generations"] == [old.name]
