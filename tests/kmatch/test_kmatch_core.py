import json
import struct
import sys
import threading
import time
from pathlib import Path

import gemmi
import pytest

from kmatch.config import Config
from kmatch.database import REQUIRED, database_status
from kmatch.gui import Job, create_app
from kmatch.models import DependencyStatus, PredictionResult
from kmatch.pipeline import Pipeline
from kmatch.predictors import normalize_af3, normalize_boltz, validate_a3m
from kmatch.process import Cancellation, Cancelled, ProcessError, run_process
from kmatch.search import FIELDS, parse_foldseek
from kmatch.sequence import parse_sequence

SEQ = "ACDEFGHIKLMNPQRSTVWY"


@pytest.mark.parametrize(
    "text", [SEQ, SEQ.lower(), f">query name\n{SEQ[:10]}\n {SEQ[10:]}\n", f">\n{SEQ}"]
)
def test_sequence_normalization(text):
    assert parse_sequence(text).sequence == SEQ


@pytest.mark.parametrize(
    "text",
    [
        "",
        ">empty",
        "ACDEFGHIKX",
        "ACDEFGHIK*",
        "ACDEFGHIK1",
        "ACDEFGHIK-",
        "ACDEFGHIKU",
        "ACDEFGHIKB",
        f">one\n{SEQ}\n>two\n{SEQ}",
        "ACD",
    ],
)
def test_invalid_sequences_are_never_silently_rewritten(text):
    with pytest.raises(ValueError):
        parse_sequence(text)


def test_a3m_sequence_validation(tmp_path):
    msa = tmp_path / "query.a3m"
    msa.write_text(f">q\n{SEQ}\n>hit\nACDEaaFGHIKLMNPQRSTVWY\n")
    assert validate_a3m(msa, SEQ).startswith(">q")
    msa.write_text(f">q\n{SEQ[:-1]}\n")
    with pytest.raises(ValueError, match="match"):
        validate_a3m(msa, SEQ)


def structure(path, length=20, confidence=85):
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for i in range(length):
        residue = gemmi.Residue()
        residue.name = "ALA"
        residue.seqid = gemmi.SeqId(i + 1, " ")
        residue.het_flag = "A"
        atom = gemmi.Atom()
        atom.name = "CA"
        atom.element = gemmi.Element("C")
        atom.pos = gemmi.Position(i * 3.8, 0, 0)
        atom.b_iso = confidence
        residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    obj = gemmi.Structure()
    obj.add_model(model)
    if path.suffix == ".cif":
        obj.make_mmcif_document().write_file(str(path))
    else:
        obj.write_pdb(str(path))


def test_boltz_normalization_ranks_by_model_confidence(tmp_path):
    raw = tmp_path / "raw/boltz_results_query/predictions/query"
    raw.mkdir(parents=True)
    for index, score in enumerate([0.8, 0.9]):
        structure(raw / f"query_model_{index}.pdb", confidence=80 + index * 10)
        (raw / f"confidence_query_model_{index}.json").write_text(
            json.dumps({"confidence_score": score, "ptm": 0.8})
        )
    result = normalize_boltz(
        tmp_path / "raw", tmp_path, "boltz2", "2.11.0", 20, 42, {"msa": "single-sequence"}
    )
    assert result.structure_path.exists() and result.confidence["mean_plddt"] == 90
    assert result.per_residue_confidence == [90] * 20
    assert result.provenance["selected_model"] == "query_model_1.pdb"


def test_af3_normalization_preserves_confidence_and_restricted_provenance(tmp_path):
    raw = tmp_path / "raw/query"
    raw.mkdir(parents=True)
    structure(raw / "query_model.cif")
    (raw / "query_summary_confidences.json").write_text(
        json.dumps({"ptm": 0.9, "ranking_score": 0.85})
    )
    result = normalize_af3(tmp_path / "raw", tmp_path, "3.0.4", 20, 100, {})
    assert result.confidence["ptm"] == 0.9 and result.per_residue_confidence == [85] * 20
    assert (
        result.provenance["experimental_mps"]
        and "OUTPUT_TERMS" in result.provenance["output_terms"]
    )


def foldseek_row(**overrides):
    data = dict(
        zip(
            FIELDS,
            [
                "q",
                "1csp_A",
                "1e-12",
                "1",
                "400",
                "0.5",
                "20",
                "20",
                "40",
                "1",
                "0.5",
                "0.95",
                "0.5",
                "0.96",
                "0.9",
                "1",
                "20",
                "2",
                "21",
                SEQ,
                SEQ,
                "1,0,0,0,1,0,0,0,1",
                "1,2,3",
            ],
            strict=True,
        )
    )
    data.update(overrides)
    return "\t".join(str(data[key]) for key in FIELDS) + "\n"


def test_foldseek_retains_every_field_and_transform(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text(foldseek_row(lddt="nan"))
    hit = parse_foldseek(path)[0]
    assert hit.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]] and hit.translation == [1, 2, 3]
    assert hit.lddt is None and hit.raw["lddt"] == "nan" and len(hit.raw) == 23
    assert hit.sequence_identity == 0.5 and hit.target_start == 2 and hit.query_alignment == SEQ


@pytest.mark.parametrize(
    "values",
    [{"u": "1,2"}, {"t": "1,2"}, {"prob": "nan"}, {"qcov": "2"}, {"qend": "25"}, {"qaln": "ABC"}],
)
def test_foldseek_rejects_incomplete_or_invalid_results(tmp_path, values):
    path = tmp_path / "hits.tsv"
    path.write_text(foldseek_row(**values))
    with pytest.raises(ValueError):
        parse_foldseek(path)


def test_database_missing_partial_and_complete(tmp_path):
    prefix = tmp_path / "pdb"
    assert not database_status(prefix).ready
    for suffix in REQUIRED:
        Path(str(prefix) + suffix).write_bytes(
            struct.pack(
                "<I", {".dbtype": 0, "_ss.dbtype": 0, "_ca.dbtype": 101, "_h.dbtype": 12}[suffix]
            )
            if suffix.endswith(".dbtype")
            else b"data\n"
        )
    assert database_status(prefix).ready
    (tmp_path / "installation.json").write_text('{"status":"installing"}')
    assert not database_status(prefix).ready
    (tmp_path / "installation.json").write_text('{"status":"complete","pdb_date":"250101"}')
    assert database_status(prefix).version == "250101"
    Path(str(prefix) + "_ca").unlink()
    assert not database_status(prefix).ready


def test_subprocess_errors_and_cancellation(tmp_path):
    cancel = Cancellation()
    with pytest.raises(ProcessError, match="code 7"):
        run_process(
            [sys.executable, "-c", 'print("diagnostic");raise SystemExit(7)'],
            tmp_path / "error.log",
            cancel,
            lambda e: None,
            "prediction",
        )
    cancel = Cancellation()
    timer = threading.Timer(0.25, cancel.cancel)
    timer.start()
    start = time.monotonic()
    with pytest.raises(Cancelled):
        run_process(
            [sys.executable, "-c", "import time;time.sleep(30)"],
            tmp_path / "cancel.log",
            cancel,
            lambda e: None,
            "prediction",
        )
    timer.join()
    assert time.monotonic() - start < 4


class FakePrediction:
    identity = "fixture"

    def __init__(self, fail=False):
        self.calls = 0
        self.fail = fail

    def status(self):
        return DependencyStatus("fixture", True, "ready", "test")

    def predict(self, sequence, output, cancel, progress, **kwargs):
        self.calls += 1
        if self.fail:
            raise RuntimeError("fixture predictor failure")
        output.mkdir(parents=True)
        path = output / "predicted.pdb"
        structure(path, length=sequence.length)
        return PredictionResult(
            path, "fixture", "test", sequence.length, {}, [85] * sequence.length, 1, {}
        )


class FakeSearch:
    def __init__(self, ready=True, fail=False):
        self.ready = ready
        self.fail = fail
        self.calls = 0

    def status(self):
        return DependencyStatus(
            "fixture", self.ready, "database missing" if not self.ready else "ready"
        )

    def search(self, prediction, output, cancel, progress):
        self.calls += 1
        assert prediction.structure_path.exists()
        if self.fail:
            raise RuntimeError("fixture search failure")
        output.mkdir()
        (output / "work").mkdir()
        (output / "work/scratch").write_text("temporary")
        return []


def test_pipeline_success_failure_cleanup_and_recovery(tmp_path):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    pred, search = FakePrediction(), FakeSearch()
    calls = []
    result = Pipeline(cfg, pred, search).run(SEQ, on_prediction=lambda r: calls.append(r))
    assert (
        pred.calls == search.calls == len(calls) == 1
        and not (result.directory / "search/work").exists()
    )
    assert json.loads((result.directory / "run.json").read_text())["status"] == "complete"
    with pytest.raises(RuntimeError, match="search failure"):
        Pipeline(cfg, FakePrediction(), FakeSearch(fail=True)).run(SEQ)
    failed = [
        p
        for p in cfg.runs.iterdir()
        if json.loads((p / "run.json").read_text())["status"] == "failed"
    ][0]
    assert (failed / "prediction/predicted.pdb").exists()
    assert Pipeline(cfg, FakePrediction(), FakeSearch()).run(SEQ).prediction.structure_path.exists()


def test_missing_database_prevents_expensive_folding(tmp_path):
    pred = FakePrediction()
    with pytest.raises(RuntimeError, match="database missing"):
        Pipeline(Config(tmp_path, tmp_path / "cache"), pred, FakeSearch(ready=False)).run(SEQ)
    assert pred.calls == 0


def test_gui_token_and_offline_startup(tmp_path):
    app = create_app(Job(Config(tmp_path, tmp_path / "cache")))
    client = app.test_client()
    page = client.get("/").text
    import re

    token = re.search(r'name="kmatch-token" content="([^"]+)"', page)[1]
    assert client.get("/api/status").status_code == 403
    headers = {"X-Kmatch-Token": token}
    assert client.get("/api/status", headers=headers).status_code == 200
    assert (
        client.get("/api/status", headers={**headers, "Origin": "https://example.com"}).status_code
        == 403
    )
    assert (
        client.post("/api/run", headers=headers, json={"sequence": "INVALID!"}).status_code == 400
    )
    assert client.get("/api/structure", headers=headers).status_code == 400


def test_network_connections_are_blocked_in_inference_subprocess(tmp_path):
    from kmatch.predictors import RUNNER

    command = [
        sys.executable,
        "-c",
        "import runpy,socket;ns=runpy.run_path(" + repr(str(RUNNER)) + ");"
        'ns["block_network"]();socket.create_connection(("127.0.0.1",9),timeout=.1)',
    ]
    with pytest.raises(ProcessError, match="cannot contact the network"):
        run_process(command, tmp_path / "blocked.log", Cancellation(), lambda e: None, "prediction")


def test_export_preserves_af3_output_notice(tmp_path):

    from kmatch.gui import DesktopApi

    source = tmp_path / "source.cif"
    source.write_text("data_query\n")
    target = tmp_path / "saved.cif"
    job = Job(Config(tmp_path, tmp_path / "cache"))
    job.update(prediction={"structure_path": str(source), "predictor": "af3"})
    api = DesktopApi(job)

    class Window:
        def create_file_dialog(self, *args, **kwargs):
            return (str(target),)

    api._window = Window()
    assert api.export_structure() == str(target)
    assert target.read_bytes() == source.read_bytes()
    assert "OUTPUT_TERMS_OF_USE" in target.with_name("saved.cif.terms.txt").read_text()


@pytest.mark.parametrize(
    "settings",
    [[], {"threads": "8"}, {"seed": -1}, {"pdb_prefix_override": 7}, {"runs_root": "/tmp"}],
)
def test_invalid_configuration_has_explicit_errors(tmp_path, monkeypatch, settings):
    monkeypatch.setenv("KMATCH_HOME", str(tmp_path))
    monkeypatch.setenv("KMATCH_CACHE", str(tmp_path / "cache"))
    (tmp_path / "config.json").write_text(json.dumps(settings))
    with pytest.raises(ValueError):
        Config.load()


def test_malformed_configuration_names_the_file(tmp_path, monkeypatch):
    monkeypatch.setenv("KMATCH_HOME", str(tmp_path))
    monkeypatch.setenv("KMATCH_CACHE", str(tmp_path / "cache"))
    (tmp_path / "config.json").write_text("{not json")
    with pytest.raises(ValueError, match="config.json is not valid JSON"):
        Config.load()


def test_foldseek_executable_readiness_is_independent_of_database(tmp_path, monkeypatch):
    from kmatch.search import FoldseekBackend

    executable = tmp_path / "foldseek"
    executable.touch()
    monkeypatch.setattr(
        "kmatch.search.subprocess.check_output", lambda *a, **k: "qualified-version"
    )
    backend = FoldseekBackend(
        Config(tmp_path, tmp_path / "cache", foldseek_override=str(executable))
    )
    assert backend.executable_status().ready
    assert not backend.status().ready


def test_job_cancellation_and_recovery(tmp_path):
    job = Job(Config(tmp_path, tmp_path / "cache"))
    entered = threading.Event()

    def operation():
        entered.set()
        while True:
            job.cancel_token.check()
            time.sleep(0.01)

    job.start(operation)
    assert entered.wait(timeout=1)
    with pytest.raises(ValueError, match="already running"):
        job.start(lambda: None)
    job.cancel()
    job.thread.join(timeout=2)
    assert job.snapshot()["phase"] == "cancelled" and not job.snapshot()["busy"]
    job.start(lambda: job.update(phase="complete"))
    job.thread.join(timeout=2)
    assert job.snapshot()["phase"] == "complete"


def test_completed_partial_download_is_verified_before_activation(tmp_path, monkeypatch):
    import hashlib

    from kmatch.setup import download

    target = tmp_path / "archive.tar"
    partial = tmp_path / "archive.tar.part"
    payload = b"public fixture archive bytes"
    partial.write_bytes(payload)
    monkeypatch.setattr(
        "kmatch.setup.subprocess.check_output",
        lambda *a, **k: f"HTTP/2 200\ncontent-length: {len(payload)}\n\n",
    )
    monkeypatch.setattr(
        "kmatch.setup.run_process",
        lambda *a, **k: pytest.fail("A complete partial must not be redownloaded"),
    )
    with pytest.raises(RuntimeError, match="checksum"):
        download(
            "https://example.com/fixture",
            target,
            Cancellation(),
            lambda e: None,
            tmp_path / "setup.log",
            expected_md5="0" * 32,
        )
    assert partial.read_bytes() == payload and not target.exists()
    download(
        "https://example.com/fixture",
        target,
        Cancellation(),
        lambda e: None,
        tmp_path / "setup.log",
        expected_md5=hashlib.md5(payload).hexdigest(),
    )
    assert target.read_bytes() == payload and not partial.exists()


def test_failed_database_update_preserves_active_generation(tmp_path, monkeypatch):
    from kmatch.setup import install_pdb

    cfg = Config(tmp_path, tmp_path / "cache")
    root = tmp_path / "databases/pdb100"
    (root / "previous").mkdir(parents=True)
    (root / "current").symlink_to("previous")
    monkeypatch.setattr(
        "kmatch.setup.subprocess.check_output",
        lambda *a, **k: "a" * 32 + "  pdb100.tar.gz\n260101 PDB_DATE\n",
    )

    def fail(*a, **k):
        raise RuntimeError("fixture interrupted download")

    monkeypatch.setattr("kmatch.setup.download", fail)
    with pytest.raises(RuntimeError, match="interrupted"):
        install_pdb(cfg, Cancellation(), lambda e: None)
    assert (root / "current").resolve() == root / "previous"
    assert not database_status(root / ("260101-" + "a" * 12) / "pdb").ready


def test_search_requests_exact_tm_scores(tmp_path, monkeypatch):
    from kmatch.search import FoldseekBackend

    generation = tmp_path / "generation"
    generation.mkdir()
    (generation / "pdb").touch()
    current = tmp_path / "current"
    current.symlink_to(generation, target_is_directory=True)
    cfg = Config(tmp_path, tmp_path / "cache", pdb_prefix_override=str(current / "pdb"))
    backend = FoldseekBackend(cfg)
    monkeypatch.setattr(
        backend, "status", lambda: DependencyStatus("fixture", True, "ready", "test")
    )

    def search_process(command, *args, **kwargs):
        assert command[command.index("--exact-tmscore") + 1] == "1"
        assert command[3] == generation / "pdb"
        Path(command[4]).write_text(foldseek_row())
        return 0.1

    monkeypatch.setattr("kmatch.search.run_process", search_process)
    monkeypatch.setattr("kmatch.search.retain_target_traces", lambda *args: {})
    structure(tmp_path / "query.pdb")
    pred = PredictionResult(tmp_path / "query.pdb", "fixture", "test", 20, {}, [], 1, {})
    hits = backend.search(pred, tmp_path / "search", Cancellation(), lambda e: None)
    assert len(hits) == 1 and hits[0].rotation[0] == [1, 0, 0]
    assert json.loads((tmp_path / "search/search.json").read_text())["query_structure_sha256"]


def test_cancellation_stops_child_that_ignores_term(tmp_path):
    import psutil

    pidfile = tmp_path / "child.pid"
    child = (
        "import signal,time,pathlib,os;signal.signal(signal.SIGTERM,signal.SIG_IGN);"
        f"pathlib.Path({str(pidfile)!r}).write_text(str(os.getpid()));time.sleep(30)"
    )
    parent = (
        "import subprocess,sys,time;"
        f'subprocess.Popen([sys.executable,"-c",{child!r}]);time.sleep(30)'
    )
    cancel = Cancellation()

    def stop_after_spawn():
        deadline = time.monotonic() + 3
        while not pidfile.exists() and time.monotonic() < deadline:
            time.sleep(0.01)
        cancel.cancel()

    thread = threading.Thread(target=stop_after_spawn)
    thread.start()
    with pytest.raises(Cancelled):
        run_process(
            [sys.executable, "-c", parent],
            tmp_path / "cancel.log",
            cancel,
            lambda e: None,
            "fixture",
            scratch_root=tmp_path / "scratch",
        )
    thread.join(timeout=4)
    assert pidfile.exists()
    try:
        child_process = psutil.Process(int(pidfile.read_text()))
        # SIGKILL delivery and orphan reaping are asynchronous on macOS.
        deadline = time.monotonic() + 1
        while child_process.is_running() and child_process.status() != psutil.STATUS_ZOMBIE:
            assert time.monotonic() < deadline
            time.sleep(0.01)
    except psutil.NoSuchProcess:
        pass
    assert not list((tmp_path / "scratch").iterdir())


def test_browser_af3_export_contains_unmodified_structure_and_notice(tmp_path):
    import io
    import re
    import zipfile

    source = tmp_path / "predicted.cif"
    source.write_text("data_query\n")
    job = Job(Config(tmp_path, tmp_path / "cache"))
    job.update(prediction={"structure_path": str(source), "predictor": "af3"})
    client = create_app(job).test_client()
    token = re.search(r'name="kmatch-token" content="([^"]+)"', client.get("/").text)[1]
    response = client.get("/api/export", headers={"X-Kmatch-Token": token})
    with zipfile.ZipFile(io.BytesIO(response.data)) as archive:
        assert archive.read("predicted.cif") == source.read_bytes()
        assert b"OUTPUT_TERMS_OF_USE" in archive.read("predicted.cif.terms.txt")


def test_registered_engine_runs_through_gui_without_frontend_changes(tmp_path, monkeypatch):
    import re

    from kmatch.predictors import BACKENDS, PredictorDefinition

    monkeypatch.setitem(
        BACKENDS,
        "fixture",
        PredictorDefinition(
            "Fixture predictor", "Fixture predictor", "Local fixture", lambda cfg: FakePrediction()
        ),
    )
    monkeypatch.setattr("kmatch.pipeline.FoldseekBackend", lambda cfg: FakeSearch())
    job = Job(Config(tmp_path, tmp_path / "cache", predictor="fixture"))
    client = create_app(job).test_client()
    token = re.search(r'name="kmatch-token" content="([^"]+)"', client.get("/").text)[1]
    headers = {"X-Kmatch-Token": token}
    data = client.get("/api/status", headers=headers).json
    assert any(item["id"] == "fixture" for item in data["backends"])
    response = client.post(
        "/api/run", headers=headers, json={"sequence": SEQ, "backend": "fixture"}
    )
    assert response.status_code == 200
    job.thread.join(timeout=2)
    result = client.get("/api/status", headers=headers).json["job"]
    assert result["phase"] == "complete" and result["prediction"]["predictor"] == "fixture"
