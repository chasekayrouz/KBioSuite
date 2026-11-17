import csv
import io
import json
import zipfile
from pathlib import Path

import pytest
from Bio import Blast, SeqIO
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from conftest import HEADER, PROTEIN, FixtureClient, blast_xml, hsp_xml, ipg_text

from kblast.export import export_bundle
from kblast.models import Run, Settings, digest
from kblast.parsers import parse_blast_xml, parse_ipg
from kblast.pipeline import Pipeline
from kblast.preview import context_preview
from kblast.sequence import parse_query


def test_single_protein_normalization():
    query = parse_query(">enzyme description\n mAvkd lGa\nTNRELVKATGQV*\n")
    assert query.sequence == PROTEIN
    assert query.identifier == "enzyme"
    assert "terminal" in query.warnings[0]


@pytest.mark.parametrize("text", ["", ">a\nMAAA\n>b\nMAAA", "MA-A", "MAAA**", "M123", ">\nMAAA"])
def test_invalid_protein_is_rejected(text):
    with pytest.raises(ValueError):
        parse_query(text)


def test_ambiguous_amino_acids_and_nucleotide_like_warning():
    assert parse_query("MBXZJUO").sequence == "MBXZJUO"
    assert any("nucleotide-like" in w for w in parse_query("ACGT").warnings)


@pytest.mark.parametrize("n,flank", [(0, 20), (2, -1), (True, 2), (2, 2.5), (5001, 2)])
def test_settings_validation(n, flank):
    with pytest.raises(ValueError):
        Settings(n, flank)


def test_intentional_kblast_settings_and_legacy_submission_provenance():
    settings = Settings(2)
    params = settings.blast_parameters()
    assert params == {
        "PROGRAM": "blastp",
        "DATABASE": "nr",
        "HITLIST_SIZE": 2,
        "EXPECT": 0.05,
        "GAPCOSTS": "11 1",
        "WORD_SIZE": 5,
        "SHORT_QUERY_ADJUST": "true",
    }
    data = Run("a" * 32, parse_query(PROTEIN), settings).to_dict()
    data["schema_version"] = 3
    data["settings"].pop("database")
    data.pop("submitted_blast_parameters")
    data["rid"], data["submission_attempted"] = "LEGACY-RID", True
    original = {
        **params,
        "DATABASE": "nr",
        "EXPECT": 10,
        "WORD_SIZE": 3,
        "FILTER": "F",
        "SHORT_QUERY_ADJUST": "false",
    }
    data["blast_parameters"] = original
    restored = Run.from_dict(data)
    assert restored.settings.database == "nr"
    assert restored.to_dict()["blast_parameters"] == original
    assert restored.submitted_blast_parameters == original


def test_blast_subject_order_aliases_and_all_hsp_summary():
    hits, metadata = parse_blast_xml(blast_xml(), 2, 20)
    assert [h.accession for h in hits] == ["PROTA.1", "PROTB.1"]
    assert len(hits[0].descriptions) == 2
    assert len(hits[0].hsps) == 2
    assert hits[0].bit_score == 200
    assert hits[0].percent_identity == 90
    assert hits[0].query_coverage == 100
    assert hits[0].hsps[1]["bit-score"] == 200
    assert metadata["statistics"]["db-num"] == "123"


def test_subject_metrics_union_overlap_disjoint_nested_and_extrema():
    xml = blast_xml()
    first_hsps = xml.split("<hsps>", 1)[1].split("</hsps>", 1)[0]
    hsps = (
        hsp_xml(10, 2, 6, 3, "1e-5")
        + hsp_xml(40, 4, 9, 4, "1e-10")
        + hsp_xml(20, 5, 7, 3, "1e-20")
        + hsp_xml(30, 15, 18, 2, "1e-8")
        + hsp_xml(30, 15, 18, 2, "1e-8")
    )
    hit = parse_blast_xml(xml.replace(first_hsps, hsps, 1), 1, 20)[0][0]
    assert hit.query_coverage == 60  # Positions 2..9 and 15..18, counted once.
    assert hit.percent_identity == 100  # Third HSP, independent of max score.
    assert hit.bit_score == 40 and hit.evalue == 1e-20
    assert len(hit.hsps) == 5 and hit.hsps[3] == hit.hsps[4]
    assert hit.hsps[0]["bit-score"] == 10


def test_blast_checks_program_database_and_query():
    with pytest.raises(ValueError):
        parse_blast_xml(blast_xml("swissprot"), 2, 20)
    with pytest.raises(ValueError, match="program/database mismatch"):
        parse_blast_xml(
            blast_xml().replace("<program>blastp</program>", "<program>blastn</program>"), 2, 20
        )
    with pytest.raises(ValueError):
        parse_blast_xml(blast_xml(), 2, 21)


def test_actual_program_mismatch_preserves_report_and_stops_before_ipg(run_fixture):
    store, run, client, pipeline = run_fixture
    client.wait_for_blast = lambda *_: blast_xml().replace(
        "<program>blastp</program>", "<program>blastn</program>"
    )
    pipeline.execute()
    assert run.status == "failure" and "program/database mismatch" in run.error
    assert not run.hits and not run.occurrences
    assert not any(call[0] == "ipg" for call in client.calls)
    assert "<program>blastn</program>" in (store.directory(run.id) / "raw/blast.xml").read_text()


def test_blast_keeps_subject_accession_version_from_primary_identifier():
    xml = blast_xml().replace("<accession>PROTA.1</accession>", "<accession>PROTA</accession>")
    assert parse_blast_xml(xml, 1, 20)[0][0].accession == "PROTA.1"


def test_public_xml2_fixture_matches_biopython_reference_parser():
    raw = (Path(__file__).parent / "fixtures" / "biopython_blastp.xml").read_bytes()
    reference = Blast.read(io.BytesIO(raw))
    hits, metadata = parse_blast_xml(raw.decode(), 2, len(reference.query.seq), "nr")
    assert [h.accession for h in hits] == ["WP_003221446.1", "BAI85158.2"]
    assert metadata["version"] == "BLASTP 2.15.0+"
    for hit, ref_hit in zip(hits, reference[:2], strict=True):
        assert hit.bit_score == max(h.annotations["bit score"] for h in ref_hit)
        assert hit.evalue == min(h.annotations["evalue"] for h in ref_hit)
        assert hit.percent_identity == max(
            100 * h.annotations["identity"] / h.shape[1] for h in ref_hit
        )
        positions = set()
        for hsp in ref_hit:
            low, high = sorted(map(int, (hsp.coordinates[1, 0], hsp.coordinates[1, -1])))
            positions.update(range(low, high))
        assert hit.query_coverage == 100 * len(positions) / len(reference.query.seq)


def test_resume_recovers_receipt_saved_before_json_checkpoint(run_fixture):
    store, run, client, pipeline = run_fixture
    run.submission_attempted = True
    receipt = store.directory(run.id) / "raw" / "blast-submission.txt"
    receipt.parent.mkdir()
    receipt.write_text("RID=RECOVERED-RID\nRTOE=60\n")
    pipeline.execute()
    assert run.rid == "RECOVERED-RID"
    assert run.status == "partial_success"
    assert not any(call[0] == "submit" for call in client.calls)


def test_ipg_expands_rows_without_deduplication():
    hit = parse_blast_xml(blast_xml(), 1, 20)[0][0]
    rows = parse_ipg(ipg_text(), hit, 100)
    assert len(rows) == 5
    assert len({o.id for o in rows}) == 5
    assert [o.ipg_row for o in rows] == [1, 2, 3, 5, 6]
    assert rows[0].nucleotide_accession == rows[1].nucleotide_accession
    assert rows[0].requested_start == -90
    assert (rows[2].cds_start, rows[2].cds_end, rows[2].strand) == (900, 959, "-")
    assert all(o.status == "pending" for o in rows)
    identical = ipg_text().splitlines()[1]
    duplicates = parse_ipg(HEADER + identical + "\n" + identical + "\n", hit, 100)
    assert len(duplicates) == 2 and duplicates[0].id != duplicates[1].id
    assert duplicates[0].raw_mapping == duplicates[1].raw_mapping
    assert parse_ipg(HEADER, hit, 100) == []
    for malformed in ["", "Error: request failed", "Id\tStart\n1\t2"]:
        with pytest.raises(ValueError):
            parse_ipg(malformed, hit, 100)


@pytest.mark.parametrize(
    "accession,start,end",
    [
        ("", "", ""),
        ("-", "10", "69"),
        ("N/A", "", ""),
    ],
)
def test_non_genomic_ipg_rows_are_excluded(accession, start, end):
    hit = parse_blast_xml(blast_xml(), 1, 20)[0][0]
    row = f"1\tRefSeq\t{accession}\t{start}\t{end}\t+\tPROTA.1\tname\torganism\t\t\n"
    assert parse_ipg(HEADER + row, hit, 100) == []


@pytest.mark.parametrize("with_nucleotide_rows", [False, True])
def test_protein_only_rows_do_not_cause_partial_success(run_fixture, with_nucleotide_rows):
    store, run, client, pipeline = run_fixture
    run.settings.top_n = 1
    lines = ipg_text().splitlines()
    text = (
        "\n".join([lines[0], *lines[1:4], lines[4]]) + "\n"
        if with_nucleotide_rows
        else (lines[0] + "\n" + lines[4] + "\n")
    )
    client.fetch_ipg = lambda *_: text
    pipeline.execute()
    assert run.status == "success" and run.error == ""
    assert len(run.occurrences) == (3 if with_nucleotide_rows else 0)
    assert run.hits[0].mapping_status == ("success" if with_nucleotide_rows else "no_mapping")
    if not with_nucleotide_rows:
        assert "No nucleotide occurrence" in run.hits[0].mapping_error
    assert (store.directory(run.id) / "raw/hit_0001.ipg.tsv").read_text() == text
    with zipfile.ZipFile(io.BytesIO(export_bundle(store.directory(run.id), run))) as archive:
        assert archive.read("provenance/hit_0001.ipg.tsv").decode() == text
    assert run.to_dict()["counts"]["failed"] == 0


def test_invalid_ipg_response_is_preserved_as_retrieval_failure(run_fixture):
    store, run, client, pipeline = run_fixture
    run.settings.top_n = 1
    client.fetch_ipg = lambda *_: "Error: synthetic malformed IPG response"
    pipeline.execute()
    assert run.hits[0].mapping_status == "parse_failure"
    assert run.status == "partial_success" and not run.occurrences
    directory = store.directory(run.id)
    assert (directory / "raw/hit_0001.ipg-error.txt").read_text().startswith("Error:")
    assert not (directory / "raw/hit_0001.ipg.tsv").exists()
    with zipfile.ZipFile(io.BytesIO(export_bundle(directory, run))) as archive:
        assert archive.read("provenance/hit_0001.ipg-error.txt").startswith(b"Error:")
    client.fetch_ipg = lambda *_: "\n".join(ipg_text().splitlines()[:4]) + "\n"
    pipeline.execute()
    assert run.status == "success" and len(run.occurrences) == 3
    assert (directory / "raw/hit_0001.ipg.tsv").is_file()
    assert (directory / "raw/hit_0001.ipg-error.txt").is_file()


def test_partial_failures_clipping_and_native_orientation(run_fixture):
    store, run, client, pipeline = run_fixture
    pipeline.execute()
    assert run.status == "partial_success"
    assert [o.status for o in run.occurrences] == [
        "success",
        "success",
        "success",
        "retrieval_failure",
        "unresolved_mapping",
    ]
    assert run.hits[1].mapping_status == "no_mapping"
    assert run.hits[1].sequence_status == "retrieval_failure"
    left, duplicate, minus = run.occurrences[:3]
    assert (left.context_start, left.context_end, left.clipped_left, left.clipped_right) == (
        1,
        169,
        True,
        False,
    )
    assert (minus.context_start, minus.context_end, minus.clipped_left, minus.clipped_right) == (
        800,
        1000,
        False,
        True,
    )
    assert left.context_file != duplicate.context_file
    assert all(o.target_cds_found for o in [left, duplicate, minus])
    record = SeqIO.read(store.directory(run.id) / minus.context_file, "genbank")
    assert str(record.seq) == ("ACGTTAGC" * 125)[799:1000]
    assert record.features[1].location.strand == -1
    assert store.load(run.id).to_dict() == run.to_dict()


def test_export_all_and_selected_preserve_record_order_and_failures(run_fixture):
    store, run, _, pipeline = run_fixture
    pipeline.execute()
    directory = store.directory(run.id)
    # Secrets and unrelated files placed beside the run must never be swept into a bundle.
    (directory / ".env").write_text("NCBI_API_KEY=fixture_secret")
    bundle = export_bundle(directory, run)
    with zipfile.ZipFile(io.BytesIO(bundle)) as archive:
        assert not any(".env" in n for n in archive.namelist())
        records = list(SeqIO.parse(io.StringIO(archive.read("contexts.gb").decode()), "genbank"))
        assert len(records) == 3
        assert records[0].id == records[1].id == "NC_000001.1"
        assert records[2].id == "NC_000002.1"
        proteins = list(
            SeqIO.parse(io.StringIO(archive.read("blast_hits.fasta").decode()), "fasta")
        )
        assert len(proteins) == 1 and str(proteins[0].seq) == PROTEIN
        rows = list(
            csv.DictReader(io.StringIO(archive.read("metadata.tsv").decode()), delimiter="\t")
        )
        assert len(rows) == 7  # 2 hit rows + 5 qualifying genomic rows, including failures.
        assert archive.read("provenance/hit_0001.ipg.tsv").decode() == ipg_text()
        assert sum(r["context_exported"] == "True" for r in rows) == 3
        manifest = json.loads(archive.read("provenance/run.json"))
        assert manifest["export"]["protein_hits_missing"] == [2]
        for o in run.occurrences[:3]:
            name = next(n for n in archive.namelist() if n.startswith("contexts/" + o.id))
            assert archive.read(name) == (directory / o.context_file).read_bytes()
    selected = export_bundle(directory, run, [run.occurrences[2].id])
    with zipfile.ZipFile(io.BytesIO(selected)) as archive:
        assert (
            len(list(SeqIO.parse(io.StringIO(archive.read("contexts.gb").decode()), "genbank")))
            == 1
        )
        assert (
            len(
                list(
                    csv.DictReader(
                        io.StringIO(archive.read("metadata.tsv").decode()), delimiter="\t"
                    )
                )
            )
            == 7
        )
    with pytest.raises(ValueError):
        export_bundle(directory, run, [run.occurrences[3].id])


def test_resume_preserves_successes_and_never_resubmits(run_fixture):
    store, run, _, pipeline = run_fixture
    pipeline.execute()
    resumed = store.load(run.id)
    client = FixtureClient()
    client.fail_context = False
    Pipeline(store, resumed, client).execute()
    assert resumed.occurrences[3].status == "success"
    assert not any(call[0] in {"submit", "blast", "ipg"} for call in client.calls)
    contexts = [call for call in client.calls if call[0] == "context"]
    assert contexts == [("context", "NC_000003.1", 150, 409)]
    assert len(resumed.occurrences) == 5


def test_cancellation_preserves_completed_items_and_can_resume(run_fixture):
    store, run, client, pipeline = run_fixture
    client.cancel_after_context = True
    pipeline.execute()
    assert run.status == "cancelled"
    assert run.occurrences[0].status == "success"
    assert run.occurrences[1].status == "cancelled"
    resumed = store.load(run.id)
    client = FixtureClient()
    Pipeline(store, resumed, client).execute()
    assert resumed.occurrences[1].status == "success"
    assert not any(call[0] == "submit" for call in client.calls)


def test_unknown_submission_cannot_automatically_resubmit(run_fixture):
    _, run, client, pipeline = run_fixture
    run.submission_attempted = True
    pipeline.execute()
    assert run.status == "failure"
    assert "resubmission is disabled" in run.error
    assert client.calls == []


def test_preview_uses_native_coordinates_and_annotation(run_fixture):
    store, run, _, pipeline = run_fixture
    pipeline.execute()
    preview = context_preview(store.directory(run.id), run.occurrences[2])
    target = next(f for f in preview["features"] if f["target"])
    assert (target["start"], target["end"], target["strand"]) == (900, 959, -1)
    assert target["product"] == "test synthase"


def test_integrity_failure_blocks_preview_and_export(run_fixture):
    store, run, _, pipeline = run_fixture
    pipeline.execute()
    path = store.directory(run.id) / run.occurrences[0].context_file
    path.write_text("corrupted fixture")
    with pytest.raises(ValueError, match="integrity"):
        context_preview(store.directory(run.id), run.occurrences[0])
    with pytest.raises(ValueError, match="integrity"):
        export_bundle(store.directory(run.id), run)


def test_preview_preserves_joined_parts_and_omits_remote_coordinates(run_fixture):
    store, run, _, pipeline = run_fixture
    pipeline.execute()
    o = run.occurrences[2]
    path = store.directory(run.id) / o.context_file
    record = SeqIO.read(path, "genbank")
    record.features.extend(
        [
            SeqFeature(
                CompoundLocation(
                    [SimpleLocation(0, 20, strand=1), SimpleLocation(50, 80, strand=1)]
                ),
                type="CDS",
                qualifiers={"gene": ["joined"]},
            ),
            SeqFeature(
                SimpleLocation(0, 10, strand=1, ref="NC_REMOTE.1"),
                type="CDS",
                qualifiers={"gene": ["remote"]},
            ),
        ]
    )
    SeqIO.write(record, path, "genbank")
    o.context_sha256 = digest(path.read_bytes())
    preview = context_preview(store.directory(run.id), o)
    joined = next(f for f in preview["features"] if f["gene"] == "joined")
    assert [(p["start"], p["end"]) for p in joined["parts"]] == [(800, 819), (850, 879)]
    assert preview["omitted_features"] == 1
    assert "NC_REMOTE.1" in path.read_text()


def test_json_round_trip_does_not_store_credentials(run_fixture):
    _, run, _, _ = run_fixture
    serialized = run.to_dict()
    assert Run.from_dict(serialized).settings == run.settings
    assert "api_key" not in json.dumps(serialized)


@pytest.mark.parametrize("mismatch", ["missing", "coordinates", "strand", "remote", "unknown"])
def test_target_cds_integrity_failure_preserves_evidence_and_other_successes(run_fixture, mismatch):
    store, run, client, pipeline = run_fixture
    run.settings.top_n = 1
    text = "\n".join(ipg_text().splitlines()[:4]) + "\n"
    client.fetch_ipg = lambda *_: text
    fetch = client.fetch_context

    def altered_context(accession, start, end):
        raw, record, stamp = fetch(accession, start, end)
        if accession == "NC_000002.1":
            feature = record.features[1]
            low, high = int(feature.location.start), int(feature.location.end)
            if mismatch == "missing":
                record.features.remove(feature)
            elif mismatch == "coordinates":
                feature.location = SimpleLocation(low + 1, high, strand=-1)
            elif mismatch == "strand":
                feature.location = SimpleLocation(low, high, strand=1)
            elif mismatch == "remote":
                feature.location = SimpleLocation(low, high, strand=-1, ref="NC_REMOTE.1")
            else:
                run.occurrences[2].strand = ""
            output = io.StringIO()
            SeqIO.write(record, output, "genbank")
            raw = output.getvalue()
        return raw, record, stamp

    client.fetch_context = altered_context
    pipeline.execute()
    assert run.status == "partial_success"
    assert [o.status for o in run.occurrences] == ["success", "success", "integrity_failure"]
    bad = run.occurrences[2]
    assert bad.target_cds_found is False and "resolved coordinates and strand" in bad.error
    path = store.directory(run.id) / bad.context_file
    assert path.is_file() and digest(path.read_bytes()) == bad.context_sha256
    assert bad.retrieved_at and bad.feature_count > 0
    assert run.to_dict()["counts"]["integrity_failures"] == 1
    assert store.load(run.id).occurrences[2].status == "integrity_failure"
    with pytest.raises(ValueError, match="verified target CDS"):
        context_preview(store.directory(run.id), bad)
    with pytest.raises(ValueError, match="successful"):
        export_bundle(store.directory(run.id), run, [bad.id])
    with zipfile.ZipFile(io.BytesIO(export_bundle(store.directory(run.id), run))) as archive:
        assert (
            len(list(SeqIO.parse(io.StringIO(archive.read("contexts.gb").decode()), "genbank")))
            == 2
        )
        rows = list(
            csv.DictReader(io.StringIO(archive.read("metadata.tsv").decode()), delimiter="\t")
        )
        audit_row = next(row for row in rows if row["occurrence_id"] == bad.id)
        assert (
            audit_row["status"] == "integrity_failure" and audit_row["context_exported"] == "False"
        )
        assert audit_row["context_sha256"] == bad.context_sha256


def test_completed_zero_hit_report_is_normal_terminal_state(run_fixture):
    store, run, client, pipeline = run_fixture
    xml = blast_xml()
    subjects = xml.split("<hits>", 1)[1].split("</hits>", 1)[0]
    empty = xml.replace(subjects, "", 1)
    client.wait_for_blast = lambda *_: empty
    pipeline.execute()
    assert run.status == "no_hits" and run.stage == "complete" and run.error == ""
    assert "BLAST completed normally" in run.message and "no protein subject hits" in run.message
    assert not run.hits and not run.occurrences
    assert run.blast_metadata["program"] == "blastp" and run.blast_completed_at
    assert client.calls == [("submit",)]
    assert (store.directory(run.id) / "raw/blast.xml").read_text() == empty
    Pipeline(store, store.load(run.id), client).execute()
    assert client.calls == [("submit",)]  # A saved no-hit report needs no remote work.
    with zipfile.ZipFile(io.BytesIO(export_bundle(store.directory(run.id), run))) as archive:
        assert json.loads(archive.read("provenance/run.json"))["status"] == "no_hits"


def test_legacy_checkpoint_normalization_keeps_raw_evidence_and_successful_contexts(run_fixture):
    store, run, client, pipeline = run_fixture
    run.settings.top_n = 1
    client.fetch_ipg = lambda *_: "\n".join(ipg_text().splitlines()[:5]) + "\n"
    pipeline.execute()
    data = run.to_dict()
    data["schema_version"], data["status"] = 1, "partial_success"
    data["hits"][0]["bit_score"] = 99
    protein_only = {
        **data["occurrences"][0],
        "id": "h0001_o000004",
        "ipg_row": 4,
        "nucleotide_accession": "",
        "cds_start": None,
        "cds_end": None,
        "status": "unresolved_mapping",
        "context_file": "",
        "target_cds_found": None,
    }
    data["occurrences"].append(protein_only)
    evidence = {p: p.read_bytes() for p in store.directory(run.id).rglob("*") if p.is_file()}
    restored = Run.from_dict(data)
    assert restored.status == "success" and restored.schema_version == 4
    assert len(restored.occurrences) == 3 and restored.hits[0].bit_score == 200
    assert restored.to_dict()["blast_summary_definitions"]["query_coverage"].startswith("union")
    assert all(p.read_bytes() == raw for p, raw in evidence.items())
    data["occurrences"][0]["target_cds_found"] = False
    rejected = Run.from_dict(data)
    assert rejected.occurrences[0].status == "integrity_failure"
    assert rejected.status == "partial_success"
