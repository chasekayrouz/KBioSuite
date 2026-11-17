"""External jobs enter the same scientifically validated pipeline without a Put."""

import io
import json
import zipfile
from datetime import datetime, timedelta
from types import SimpleNamespace

import pytest
from conftest import PROTEIN, FixtureClient, blast_xml
from test_blast_lifecycle import FakeClock, rig

from kblast.app import create_app
from kblast.config import Config
from kblast.export import export_bundle
from kblast.pipeline import Pipeline
from kblast.service import RunManager


def import_rig(tmp_path, monkeypatch, replies, *, top_n=2, source="browser"):
    manager = RunManager(Config(data_dir=tmp_path / "runs"))
    monkeypatch.setattr(manager, "launch", lambda *_: None)
    run = manager.start(PROTEIN, top_n, 100, rid="BROWSER-RID", rid_source=source)
    clock = FakeClock()
    clock.base = datetime.fromisoformat(run.rid_adopted_at) - timedelta(seconds=clock.seconds)
    clock, store, run, client, sent = rig(
        tmp_path, monkeypatch, list(replies), clock=clock, run=run
    )
    return manager, clock, store, run, client, sent


@pytest.mark.parametrize("source", ["browser", "external"])
@pytest.mark.parametrize("waiting", [False, True])
def test_fresh_external_ready_or_waiting_never_submits(tmp_path, monkeypatch, source, waiting):
    replies = (["Status=WAITING"] if waiting else []) + [blast_xml()]
    _, _, store, run, client, sent = import_rig(tmp_path, monkeypatch, replies, source=source)
    Pipeline(store, run, client).blast()
    assert len(run.hits) == 2 and run.rid_validation["status"] == "validated"
    assert run.rid_source == source and run.rid_adopted_at
    assert not run.submission_attempted and not run.submitted_blast_parameters
    assert all(p["CMD"] == "Get" and p["RID"] == "BROWSER-RID" for _, p in sent)
    saved = store.load(run.id)
    assert saved.rid_history[-1]["returned_metadata"]["database"] == "nr"
    assert saved.to_dict()["blast_parameters"] == saved.blast_metadata["parameters"]


def test_switch_waiting_local_rid_retains_original_and_uses_adopted_job(tmp_path, monkeypatch):
    clock, store, run, client, original_sent = rig(tmp_path, monkeypatch, ["Status=WAITING"])
    run.rid_source = "kblast_url_api"
    run.original_rid = run.rid
    run.submission_attempted = True
    run.submitted_blast_parameters = run.settings.blast_parameters()
    store.save(run)
    clock.cancel_at = clock.seconds + 20
    Pipeline(store, run, client).execute()
    assert run.status == "cancelled" and run.ncbi_status == "WAITING"
    original_at = run.submitted_at
    manager = RunManager(Config(data_dir=store.root))
    monkeypatch.setattr(manager, "launch", lambda *_: None)
    monkeypatch.setattr("kblast.service.utc_now", lambda: clock.now().isoformat())
    adopted = manager.resume(run.id, rid="BROWSER-RID", rid_source="browser")
    assert adopted.original_rid == "TEST-RID" and adopted.rid == "BROWSER-RID"
    assert adopted.submitted_at == original_at
    assert adopted.submitted_blast_parameters == run.submitted_blast_parameters
    assert adopted.rid_history[0]["rid"] == "TEST-RID"
    assert adopted.rid_history[0]["superseded_by"] == "BROWSER-RID"
    clock.cancel_at = None
    clock, store, adopted, client, sent = rig(
        tmp_path, monkeypatch, [blast_xml()], clock=clock, run=adopted
    )
    Pipeline(store, adopted, client).blast()
    assert adopted.rid_source == "browser" and len(adopted.hits) == 2
    assert all(p["CMD"] == "Get" for _, p in original_sent + sent)
    assert sent[0][1]["RID"] == "BROWSER-RID"
    assert store.load(adopted.id).original_rid == "TEST-RID"


@pytest.mark.parametrize(
    "mutate,match",
    [
        (lambda x: x.replace("<db>nr</db>", "<db>nr_cluster_seq</db>"), "database mismatch"),
        (
            lambda x: x.replace("<program>blastp</program>", "<program>blastn</program>"),
            "program/database mismatch",
        ),
        (lambda x: x.replace("<Search>", "<Search/><Search>", 1), "exactly one"),
        (
            lambda x: x.replace("<query-len>20</query-len>", "<query-len>21</query-len>"),
            "query length",
        ),
        (lambda x: x[:-30], "Malformed NCBI BLAST"),
    ],
)
def test_incompatible_or_malformed_result_blocks_downstream(tmp_path, monkeypatch, mutate, match):
    _, _, store, run, client, sent = import_rig(tmp_path, monkeypatch, [mutate(blast_xml())])
    Pipeline(store, run, client).execute()
    assert run.status == "failure" and match in run.error
    assert not run.hits and not run.occurrences
    assert run.rid_validation["status"] == "rejected"
    assert all(p["CMD"] == "Get" for _, p in sent)
    assert (store.directory(run.id) / "raw/blast.xml").exists()


@pytest.mark.parametrize(
    "rid", ["bad", "lowercase", "ABC 123", "../ABC123", "A" * 65, "ABC123&CMD=Put", 0, None, False]
)
def test_malformed_rid_rejected_before_any_remote_work(tmp_path, monkeypatch, rid):
    manager = RunManager(Config(data_dir=tmp_path / "runs"))
    launched = []
    monkeypatch.setattr(manager, "launch", lambda *args: launched.append(args))
    with pytest.raises(ValueError, match="valid NCBI RID"):
        manager.start(PROTEIN, 2, 100, rid=rid)
    assert not launched and manager.store.list_runs() == []


@pytest.mark.parametrize("status", ["UNKNOWN", "FAILED"])
def test_expired_or_failed_rid_never_submits_or_maps(tmp_path, monkeypatch, status):
    _, _, store, run, client, sent = import_rig(tmp_path, monkeypatch, [f"Status={status}"])
    Pipeline(store, run, client).execute()
    assert run.status == "failure" and status.lower() in run.error
    assert run.ncbi_status == status and not run.hits and not run.occurrences
    assert len(sent) == 1 and sent[0][1]["CMD"] == "Get"


def test_reactivate_original_job_preserves_origin_and_every_adoption(tmp_path, monkeypatch):
    _, _, store, run, _, _ = import_rig(tmp_path, monkeypatch, [])
    run.rid = "ORIGINAL-RID"
    run.rid_source = "kblast_url_api"
    run.rid_adopted_at = ""
    run.original_rid = run.rid
    run.submitted_blast_parameters = run.settings.blast_parameters()
    run.rid_history = []
    store.save(run)
    manager = RunManager(Config(data_dir=store.root))
    monkeypatch.setattr(manager, "launch", lambda *_: None)
    switched = manager.resume(run.id, rid="BROWSER-RID", rid_source="browser")
    restored = manager.resume(run.id, rid="ORIGINAL-RID", rid_source="external")
    assert restored.original_rid == restored.rid == "ORIGINAL-RID"
    assert restored.rid_source == "kblast_url_api"
    assert [e["rid"] for e in restored.rid_history] == [
        "ORIGINAL-RID",
        "BROWSER-RID",
        "ORIGINAL-RID",
    ]
    assert restored.rid_history[0]["submitted_parameters"] == run.settings.blast_parameters()
    assert restored.rid_history[1]["superseded_by"] == "ORIGINAL-RID"
    assert switched.original_rid == "ORIGINAL-RID"


def test_rejected_report_is_archived_before_another_rid_is_adopted(tmp_path, monkeypatch):
    _, _, store, run, client, _ = import_rig(tmp_path, monkeypatch, [blast_xml("nr_cluster_seq")])
    Pipeline(store, run, client).execute()
    rejected = (store.directory(run.id) / "raw/blast.xml").read_bytes()
    manager = RunManager(Config(data_dir=store.root))
    monkeypatch.setattr(manager, "launch", lambda *_: None)
    adopted = manager.resume(run.id, rid="SECOND-RID")
    assert not (store.directory(run.id) / "raw/blast.xml").exists()
    archived = list((store.directory(run.id) / "raw").glob("blast-BROWSER-RID-*.xml"))
    assert len(archived) == 1 and archived[0].read_bytes() == rejected
    assert adopted.rid_history[0]["validation"]["status"] == "rejected"
    with zipfile.ZipFile(io.BytesIO(export_bundle(store.directory(run.id), adopted))) as bundle:
        assert bundle.read("provenance/" + archived[0].name) == rejected


@pytest.mark.parametrize("top_n,expected", [(1, 1), (3, 3), (5, 3)])
def test_available_hits_consumed_in_original_rank_order(tmp_path, monkeypatch, top_n, expected):
    _, _, store, run, client, _ = import_rig(tmp_path, monkeypatch, [blast_xml()], top_n=top_n)
    Pipeline(store, run, client).blast()
    assert len(run.hits) == expected
    assert [h.accession for h in run.hits] == ["PROTA.1", "PROTB.1", "PROTC.1"][:expected]
    assert run.blast_metadata["returned_hit_count"] == 3
    assert run.blast_metadata["consumed_hit_count"] == expected
    assert run.blast_metadata["requested_top_n"] == top_n


def test_zero_hits_is_normal_completion_without_mapping(tmp_path, monkeypatch):
    xml = blast_xml()
    xml = xml.replace(xml.split("<hits>")[1].split("</hits>")[0], "")
    _, _, store, run, client, sent = import_rig(tmp_path, monkeypatch, [xml])
    Pipeline(store, run, client).execute()
    assert run.status == "no_hits" and not run.occurrences
    assert run.rid_validation["status"] == "validated" and len(sent) == 1
    assert run.blast_metadata["query_validation"]["method"] == "length_only"
    assert run.blast_metadata["query_validation"]["residual_risk"]


def test_different_valid_parameters_are_preserved_not_claimed_as_defaults(tmp_path, monkeypatch):
    xml = blast_xml().replace("<expect>0.05</expect>", "<expect>10</expect>")
    xml = xml.replace("<matrix>BLOSUM62</matrix>", "<matrix>PAM30</matrix>")
    _, _, store, run, client, _ = import_rig(tmp_path, monkeypatch, [xml])
    Pipeline(store, run, client).blast()
    assert run.to_dict()["blast_parameters"] == {"matrix": "PAM30", "expect": "10"}
    assert run.to_dict()["planned_blast_parameters"]["EXPECT"] == 0.05
    assert run.rid_history[0]["submitted_parameters"] == {}


def test_cancel_and_restart_resume_adopted_rid_without_put(tmp_path, monkeypatch):
    _, clock, store, run, client, sent = import_rig(tmp_path, monkeypatch, ["Status=WAITING"])
    clock.cancel_at = clock.seconds + 20
    Pipeline(store, run, client).execute()
    assert run.status == "cancelled" and run.rid == "BROWSER-RID"
    run.status = "running"
    store.save(run)
    manager = RunManager(Config(data_dir=store.root))
    restored = store.load(run.id)
    assert restored.status == "interrupted" and restored.rid_source == "browser"
    monkeypatch.setattr(manager, "launch", lambda *_: None)
    restored = manager.resume(run.id)
    clock.cancel_at = None
    clock, store, restored, client, later = rig(
        tmp_path, monkeypatch, [blast_xml()], clock=clock, run=restored
    )
    Pipeline(store, restored, client).blast()
    assert all(p["CMD"] == "Get" and p["RID"] == "BROWSER-RID" for _, p in sent + later)
    assert later[0][0] - sent[0][0] >= 60
    assert restored.rid_history[0]["adopted_at"] == run.rid_adopted_at


def test_partial_query_residues_checked_even_beyond_consumed_top_n(tmp_path, monkeypatch):
    xml = blast_xml().replace(
        "<positive>20</positive>", f"<qseq>{PROTEIN}</qseq><positive>20</positive>", 1
    )
    _, _, store, run, client, _ = import_rig(tmp_path, monkeypatch, [xml], top_n=1)
    Pipeline(store, run, client).blast()
    assert run.rid_validation["checked_positions"] == 20
    assert not run.rid_validation["residual_risk"]
    # A wrong same-length query in a later hit must invalidate the entire report.
    bad = blast_xml().replace(
        "<bit-score>97</bit-score>",
        "<qseq>AAAAAAAAAAAAAAAAAAAA</qseq><bit-score>97</bit-score>",
    )
    _, _, store, run, client, _ = import_rig(tmp_path / "bad", monkeypatch, [bad], top_n=1)
    Pipeline(store, run, client).execute()
    assert run.status == "failure" and not run.hits


def test_imported_and_local_results_have_identical_downstream_and_exports(
    run_fixture,
    tmp_path,
    monkeypatch,
):
    local_store, local, fixture, local_pipeline = run_fixture
    fixture.fail_context = False
    local_pipeline.execute()
    _, _, store, run, client, sent = import_rig(tmp_path / "external", monkeypatch, [blast_xml()])
    downstream = FixtureClient()
    downstream.fail_context = False
    for name in (
        "prefetch_hits",
        "fetch_protein",
        "fetch_ipg",
        "linked_occurrences",
        "prefetch_nucleotide_info",
        "nucleotide_info",
        "context_results",
    ):
        setattr(client, name, getattr(downstream, name))
    Pipeline(store, run, client).execute()
    assert run.status == local.status and run.hits == local.hits
    assert [
        (o.id, o.status, o.context_start, o.context_end, o.strand) for o in run.occurrences
    ] == [(o.id, o.status, o.context_start, o.context_end, o.strand) for o in local.occurrences]
    assert len(sent) == 1 and sent[0][1]["CMD"] == "Get"
    with zipfile.ZipFile(io.BytesIO(export_bundle(store.directory(run.id), run))) as adopted:
        with zipfile.ZipFile(
            io.BytesIO(export_bundle(local_store.directory(local.id), local))
        ) as own:
            for name in ("contexts.gb", "blast_hits.fasta"):
                assert adopted.read(name) == own.read(name)
        provenance = json.loads(adopted.read("provenance/run.json"))
        assert provenance["rid_source"] == "browser"
        assert provenance["rid_history"][0]["validation"]["status"] == "validated"


def test_active_worker_must_finish_before_adoption_and_completed_data_is_preserved(
    run_fixture,
    monkeypatch,
):
    store, run, _, pipeline = run_fixture
    manager = RunManager(Config(data_dir=store.root))
    manager.worker = SimpleNamespace(is_alive=lambda: True)
    with pytest.raises(ValueError, match="already active"):
        manager.resume(run.id, rid="BROWSER-RID")
    manager.worker = None
    pipeline.execute()
    original = store.load(run.id).to_dict()
    with pytest.raises(ValueError, match="already has parsed"):
        manager.resume(run.id, rid="BROWSER-RID")
    assert store.load(run.id).to_dict() == original


def test_flask_fresh_import_switch_and_official_page_without_query(tmp_path, monkeypatch):
    manager = RunManager(Config(data_dir=tmp_path / "runs"))
    monkeypatch.setattr(manager, "launch", lambda *_: None)
    app = create_app(manager)
    client = app.test_client()
    headers = {"X-Kblast-Token": app.config["KBLAST_TOKEN"]}
    response = client.post(
        "/api/runs",
        headers=headers,
        json={
            "sequence": PROTEIN,
            "top_n": 2,
            "flank": 100,
            "rid": "BROWSER-RID",
            "rid_source": "browser",
        },
    )
    assert response.status_code == 201
    run_id = response.json["id"]
    response = client.post(
        f"/api/runs/{run_id}/resume",
        headers=headers,
        json={
            "rid": "SECOND-RID",
            "rid_source": "external",
        },
    )
    assert response.status_code == 200
    data = client.get(f"/api/runs/{run_id}", headers=headers).json
    assert data["rid"] == "SECOND-RID" and data["rid_history"][0]["rid"] == "BROWSER-RID"
    opened = []
    monkeypatch.setattr("kblast.app.webbrowser.open", lambda url: opened.append(url))
    assert client.post("/api/open-protein-blast", headers=headers, json={}).status_code == 200
    assert opened == ["https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins"]
    assert PROTEIN not in opened[0]
