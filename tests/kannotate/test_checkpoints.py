"""Disposable sources and evidence exercise restart/publication integrity."""

import json
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from kannotate import AnnotationEngine, Hit
from kannotate.batch import annotate_batch
from kannotate.cancellation import RunCancelled, RunControl
from kannotate.checkpoints import CheckpointBatch
from kannotate.io.fasta import read_fasta
from kannotate.search.common import SearchError


class FixtureAdapter:
    def __init__(self, executable, database, engine):
        self.executable, self.database, self.engine = str(executable), str(database), engine
        self.calls = []
        self.fail = False

    def search(self, query, workdir, threads):
        proteins = read_fasta(query)
        self.calls.append([p.sequence for p in proteins])
        if self.fail:
            raise SearchError("database fixture failure")
        (workdir / "raw.tsv").write_text("synthetic raw evidence")
        return [
            Hit(
                self.engine,
                p.header,
                "fixture",
                "fixture",
                "Fixture family",
                1e-15,
                1,
                len(p.sequence),
                len(p.sequence),
                bitscore=60,
            )
            for p in proteins
        ]


def engine_at(tmp_path):
    runtime = tmp_path / "runtime"
    runtime.mkdir()
    executable = runtime / "hmmscan"
    executable.write_text('#!/bin/sh\necho "HMMER 3.4 fixture"\n')
    executable.chmod(0o755)
    databases = {}
    for source in ["pfam", "cdd"]:
        directory = runtime / source
        directory.mkdir()
        (directory / "db.bin").write_bytes(b"synthetic reference")
        (directory / "installation.json").write_text(json.dumps({"files": {"db.bin": {}}}))
        databases[source] = {"version": "fixture", "path": str(directory / "db.bin")}
    engine = AnnotationEngine(
        FixtureAdapter(executable, databases["pfam"]["path"], "hmm"),
        FixtureAdapter(executable, databases["cdd"]["path"], "rps"),
        database_versions={"pfam": "fixture", "cdd": "fixture"},
    )
    engine.config = {"databases": databases}
    return engine


def genbank_source(path):
    records = []
    for i, count in enumerate([260, 260, 0, 513, 2]):
        record = SeqRecord(Seq("ATG" * 8), id=f"record{i}", annotations={"molecule_type": "DNA"})
        record.features = [
            SeqFeature(
                SimpleLocation(0, 21, strand=-1 if j % 2 else 1),
                type="CDS",
                qualifiers={
                    "locus_tag": [f"000{j}"],
                    "product": ['old, "product" β'],
                    **({"translation": ["MACDEFG"]} if j % 31 else {"pseudo": [""]}),
                },
            )
            for j in range(count)
        ]
        records.append(record)
    SeqIO.write(records, path, "genbank")


def run(source, out, engine, **kwargs):
    return annotate_batch(
        [source], out, engine, threads=3, min_score=6, retain_annotations=False, **kwargs
    )[0]


@pytest.mark.parametrize("kind", ["fasta", "genbank"])
@pytest.mark.parametrize("mode", ["csv", "genbank", "both"])
def test_interrupt_resume_equals_unchunked_bytes_and_keeps_completed_groups(tmp_path, kind, mode):
    engine = engine_at(tmp_path)
    source = tmp_path / ("same.name.faa" if kind == "fasta" else "same.name.gbff")
    if kind == "fasta":
        source.write_text('>duplicate, "β"\nMACDEFG\n' * 1050)
    else:
        genbank_source(source)
    before = source.read_bytes()
    baseline = run(source, tmp_path / "baseline", engine, genbank_output=mode)
    assert not baseline.error
    engine.hmm_adapter.calls.clear()
    control = RunControl()
    checkpoint = tmp_path / "state"

    def stop(path, count):
        if count:
            control.request_stop()

    with pytest.raises(RunCancelled):
        run(
            source,
            tmp_path / "out",
            engine,
            genbank_output=mode,
            checkpoint_dir=checkpoint,
            control=control,
            on_progress=stop,
        )
    saved = {str(p): p.read_bytes() for p in checkpoint.rglob("*") if p.is_file()}
    assert len(list(checkpoint.glob("file-*/group-*"))) == 1
    calls = len(engine.hmm_adapter.calls)
    resumed = run(
        source,
        tmp_path / "out",
        engine,
        genbank_output=mode,
        checkpoint_dir=checkpoint,
        resume=True,
    )
    assert not resumed.error
    assert source.read_bytes() == before
    assert [p.name for p in resumed.outputs] == [p.name for p in baseline.outputs]
    assert [p.read_bytes() for p in resumed.outputs] == [p.read_bytes() for p in baseline.outputs]
    assert all(Path(p).read_bytes() == b for p, b in saved.items())
    assert len(engine.hmm_adapter.calls) > calls
    calls = len(engine.hmm_adapter.calls)
    again = run(
        source,
        tmp_path / "out",
        engine,
        genbank_output=mode,
        checkpoint_dir=checkpoint,
        resume=True,
    )
    assert not again.error and len(engine.hmm_adapter.calls) == calls


@pytest.mark.parametrize("mutation", ["source", "reference", "tool", "options", "format", "code"])
def test_resume_rejects_identity_changes_without_searching(tmp_path, monkeypatch, mutation):
    engine = engine_at(tmp_path)
    source = tmp_path / "input.faa"
    source.write_text(">p\nMACDEFG\n")
    state = tmp_path / "state"
    first = run(source, tmp_path / "out", engine, checkpoint_dir=state)
    assert not first.error
    calls = len(engine.hmm_adapter.calls)
    kwargs = {}
    if mutation == "source":
        source.write_text(">q\nMACDEFG\n")
    elif mutation == "reference":
        Path(engine.hmm_adapter.database).write_bytes(b"changed reference")
    elif mutation == "tool":
        Path(engine.hmm_adapter.executable).write_text('#!/bin/sh\necho "HMMER changed"\n')
    elif mutation == "options":
        kwargs["threads"] = 4
    elif mutation == "format":
        kwargs["genbank_output"] = "both"
    else:
        monkeypatch.setattr("kannotate.checkpoints.implementation_identity", lambda: {"changed": 1})
    with pytest.raises(ValueError, match="identity changed"):
        annotate_batch(
            [source],
            tmp_path / "out",
            engine,
            threads=kwargs.pop("threads", 3),
            min_score=6,
            retain_annotations=False,
            checkpoint_dir=state,
            resume=True,
            **kwargs,
        )
    assert len(engine.hmm_adapter.calls) == calls


def test_corrupt_checkpoint_and_no_clobber_are_explicit(tmp_path):
    engine = engine_at(tmp_path)
    source = tmp_path / "input.faa"
    source.write_text(">p\nMACDEFG\n")
    state = tmp_path / "state"
    first = run(source, tmp_path / "out", engine, checkpoint_dir=state)
    expected = first.output.read_bytes()
    evidence = next(state.glob("file-*/group-*/evidence/*/hmm/*/raw.tsv"))
    evidence.write_text("corrupt")
    result = run(source, tmp_path / "out", engine, checkpoint_dir=state, resume=True)
    assert "evidence" in result.error
    assert first.output.read_bytes() == expected
    first.output.write_bytes(b"user output")
    result = run(source, tmp_path / "out", engine, checkpoint_dir=state, resume=True)
    assert "Output exists" in result.error
    assert first.output.read_bytes() == b"user output"


def test_warned_group_cannot_commit_but_parser_warnings_are_retained(tmp_path):
    engine = engine_at(tmp_path)
    source = tmp_path / "input.faa"
    source.write_text(">p\nMACDEFG\n")
    engine.hmm_adapter.fail = True
    state = tmp_path / "state"
    first = run(source, tmp_path / "out", engine, checkpoint_dir=state)
    assert "warnings prevent" in first.error
    assert not list(state.glob("file-*/group-*"))
    engine.hmm_adapter.fail = False
    resumed = run(source, tmp_path / "out", engine, checkpoint_dir=state, resume=True)
    assert not resumed.error


def test_exclusive_lock_and_crash_between_pair_publication(tmp_path, monkeypatch):
    import kannotate.batch as batch

    engine = engine_at(tmp_path)
    source = tmp_path / "input.gbk"
    genbank_source(source)
    state = tmp_path / "state"
    parameters = {"threads": 3, "min_score": 6, "genbank_output": "both", "table_format": "csv"}
    with CheckpointBatch.open(
        state, [source], tmp_path / "out", engine, resume=False, parameters=parameters
    ):
        with pytest.raises(ValueError, match="already in use"):
            run(
                source,
                tmp_path / "out",
                engine,
                genbank_output="both",
                checkpoint_dir=state,
                resume=True,
            )
    original = batch._publish

    def crash(staged, **kwargs):
        original(staged[:1], **kwargs)
        raise KeyboardInterrupt("simulated process loss after first link")

    with monkeypatch.context() as context:
        context.setattr(batch, "_publish", crash)
        with pytest.raises(KeyboardInterrupt):
            run(
                source,
                tmp_path / "out",
                engine,
                genbank_output="both",
                checkpoint_dir=state,
                resume=True,
            )
    kept = next((tmp_path / "out").glob("*.gbk"))
    identity, contents = kept.stat().st_ino, kept.read_bytes()
    calls = len(engine.hmm_adapter.calls)
    result = run(
        source, tmp_path / "out", engine, genbank_output="both", checkpoint_dir=state, resume=True
    )
    assert not result.error and len(result.outputs) == 2
    assert len(engine.hmm_adapter.calls) == calls
    assert kept.stat().st_ino == identity and kept.read_bytes() == contents


def test_parser_warning_empty_rows_and_original_alias_name_survive(tmp_path):
    engine = engine_at(tmp_path)
    target = tmp_path / "original.gbk"
    record = SeqRecord(Seq("ATG" * 8), id="record", annotations={"molecule_type": "DNA"})
    record.features = [
        SeqFeature(
            SimpleLocation(0, 21),
            type="CDS",
            qualifiers={"translation": [""], "product": ["original"]},
        )
    ]
    SeqIO.write(record, target, "genbank")
    target.write_text(
        target.read_text().replace("CDS             1..21", "CDS             broken_location")
    )
    alias = tmp_path / "input.original.name.gbff"
    alias.symlink_to(target)
    before = target.read_bytes()
    state = tmp_path / "state"
    baseline = run(alias, tmp_path / "baseline", engine, genbank_output="both")
    first = run(alias, tmp_path / "out", engine, genbank_output="both", checkpoint_dir=state)
    assert not first.error and first.warnings
    assert [p.read_bytes() for p in first.outputs] == [p.read_bytes() for p in baseline.outputs]
    assert "input.original.name.gbff" in first.outputs[1].read_text(encoding="utf-8-sig")
    repeat = run(
        alias, tmp_path / "out", engine, genbank_output="both", checkpoint_dir=state, resume=True
    )
    assert not repeat.error and repeat.warnings == first.warnings
    assert target.read_bytes() == before


def test_source_mutation_during_group_never_commits_mixed_evidence(tmp_path):
    engine = engine_at(tmp_path)
    source = tmp_path / "input.faa"
    source.write_text(">p\nMACDEFG\n")
    search = engine.hmm_adapter.search

    def change(*args):
        result = search(*args)
        source.write_text(">changed\nMACDEFG\n")
        return result

    engine.hmm_adapter.search = change
    state = tmp_path / "state"
    result = run(source, tmp_path / "out", engine, checkpoint_dir=state)
    assert "Input changed" in result.error
    assert not list(state.glob("file-*/group-*"))
    assert not list((tmp_path / "out").iterdir())


def test_stop_during_initial_verification_can_restart_without_mixing_a_plan(tmp_path, monkeypatch):
    import kannotate.checkpoints as checkpoints

    engine = engine_at(tmp_path)
    source = tmp_path / "input.faa"
    source.write_text(">p\nMACDEFG\n")
    state = tmp_path / "state"
    original = checkpoints.engine_identity

    def stopped(engine):
        raise RunCancelled("stopped during reference hashing")

    monkeypatch.setattr(checkpoints, "engine_identity", stopped)
    with pytest.raises(RunCancelled):
        run(source, tmp_path / "out", engine, checkpoint_dir=state)
    assert (state / "initializing.json").is_file() and not (state / "plan.json").exists()
    monkeypatch.setattr(checkpoints, "engine_identity", original)
    result = run(source, tmp_path / "out", engine, checkpoint_dir=state, resume=True)
    assert not result.error and result.run.protein_count == 1
