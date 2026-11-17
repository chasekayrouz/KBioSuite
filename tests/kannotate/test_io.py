import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from openpyxl import load_workbook

from kannotate import AnnotationEngine, Protein, __version__
from kannotate.batch import annotate_batch
from kannotate.io.detection import detect_format
from kannotate.io.fasta import read_fasta
from kannotate.io.genbank import read_genbank, write_genbank
from kannotate.io.xlsx import GENBANK_HEADERS, HEADERS, write_xlsx
from kannotate.selection import select_annotation


class NoHits:
    def search(self, *args):
        return []


def sample_records():
    record = SeqRecord(Seq("ATGAAATTT" * 10), id="one", name="one", description="test")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["comment"] = "Original reference comment"
    record.features = [
        SeqFeature(
            SimpleLocation(0, 90), type="source", qualifiers={"organism": ["Test bacterium"]}
        ),
        SeqFeature(
            CompoundLocation([SimpleLocation(0, 9, -1), SimpleLocation(18, 27, -1)]),
            type="CDS",
            qualifiers={
                "locus_tag": ["A"],
                "translation": ["ACDEFGHIK*QQ"],
                "product": ["old product"],
                "gene": ["old"],
                "EC_number": ["1.2.3.4"],
                "note": ["old note"],
                "inference": ["old"],
                "db_xref": ["old"],
                "protein_id": ["old"],
                "custom_function": ["old"],
            },
        ),
        SeqFeature(
            SimpleLocation(0, 9, 1),
            type="CDS",
            qualifiers={"locus_tag": ["B"], "transl_table": ["11"]},
        ),
        SeqFeature(SimpleLocation(3, 5, 1), type="CDS", qualifiers={"locus_tag": ["C"]}),
    ]
    second = SeqRecord(Seq("ATGAAA"), id="two", name="two", description="second")
    second.annotations["molecule_type"] = "DNA"
    second.features = [
        SeqFeature(SimpleLocation(0, 6, 1), type="CDS", qualifiers={"translation": [""]})
    ]
    return [record, second]


def test_permissive_fasta_one_row_per_header(tmp_path):
    path = tmp_path / "proteins.gbk"
    path.write_text("\ufeff>same\n ac dEf!GG\n>same\n>\n M\n>last\nbjzoux\n")
    proteins = read_fasta(path)
    assert len(proteins) == 4
    assert [p.sequence for p in proteins] == ["ACDEF", "", "M", "BJZOUX"]
    assert detect_format(path) == "fasta"
    out = tmp_path / "output.xlsx"
    write_xlsx([select_annotation(p, [], []) for p in proteins], out)
    workbook = load_workbook(out)
    sheet = workbook.active
    assert sheet.max_row == 5
    assert sheet.max_column == 24
    assert list(next(sheet.values)) == HEADERS
    assert sheet.freeze_panes == "A2"
    assert sheet.auto_filter.ref == "A1:X5"
    assert all(row[2] == "hypothetical protein" for row in list(sheet.values)[1:])


def test_genbank_clean_rewrite_multirecord_and_translation_truth(tmp_path):
    path = tmp_path / "sample.faa"
    reference = sample_records()
    SeqIO.write(reference, path, "genbank")
    assert detect_format(path) == "genbank"
    data = read_genbank(path)
    assert [p.sequence for p in data.proteins] == ["ACDEFGHIK*QQ", "MKF", "", ""]
    engine = AnnotationEngine(NoHits(), NoHits())
    run = engine.annotate(data.proteins)
    output = tmp_path / "fresh.gbk"
    write_genbank(data, run.annotations, output)
    written = list(SeqIO.parse(output, "genbank"))
    assert len(written) == 2
    assert written[0].features[0] == data.records[0].features[0]
    assert written[0].features[1].location == reference[0].features[1].location
    assert written[0].features[1].qualifiers == {
        "locus_tag": ["A"],
        "translation": ["ACDEFGHIK*QQ"],
        "product": ["hypothetical protein"],
    }
    assert written[0].features[2].qualifiers["translation"] == ["MKF"]
    assert output.read_text().count(f"Annotated with Kannotate v{__version__}") == 1
    data2 = read_genbank(output)
    write_genbank(data2, engine.annotate(data2.proteins).annotations, output)
    assert output.read_text().count(f"Annotated with Kannotate v{__version__}") == 1


def test_batch_failures_no_clobber_and_one_to_one(tmp_path):
    a = tmp_path / "a.faa"
    a.write_text(">a\nMKT\n")
    bad = tmp_path / "bad.faa"
    bad.write_text("not fasta")
    gb = tmp_path / "b.gbk"
    SeqIO.write(sample_records(), gb, "genbank")
    results = tmp_path / "results"
    engine = AnnotationEngine(NoHits(), NoHits())
    batch = annotate_batch(
        [bad, a, gb],
        results,
        engine,
        table_format="xlsx",
    )
    assert [r.error is None for r in batch] == [False, True, True]
    assert sorted(p.name for p in results.iterdir()) == ["a.kannotate.xlsx", "b.kannotate.gbk"]
    original = (results / "a.kannotate.xlsx").read_bytes()
    assert annotate_batch(
        [a],
        results,
        engine,
        table_format="xlsx",
    )[0].error
    assert (results / "a.kannotate.xlsx").read_bytes() == original
    assert not annotate_batch(
        [a],
        results,
        engine,
        force=True,
        table_format="xlsx",
    )[0].error
    assert annotate_batch(
        [a, a],
        results,
        engine,
        force=True,
        table_format="xlsx",
    )[1].error


def test_literal_excel_header_is_not_formula(tmp_path):
    out = tmp_path / "test.xlsx"
    write_xlsx([select_annotation(Protein("x", "=1+1", ""), [], [])], out)
    assert load_workbook(out).active["A2"].data_type == "s"


def test_workbook_top_two_candidates_preserve_original_descriptions(tmp_path):
    from kannotate.models import Hit

    protein = Protein("q", "header", "M" * 100)
    hmm = [
        Hit(
            "hmm",
            "q",
            f"PF{i}",
            f"HMM{i}",
            f"Full HMM description {i}",
            10 ** (-10 - i),
            2 + i,
            90,
            100,
        )
        for i in range(3)
    ]
    rps = [
        Hit(
            "rps",
            "q",
            f"cd{i}",
            f"RPS{i}",
            f"Full RPS description {i}; Reviewed",
            10 ** (-20 - i),
            4 + i,
            95,
            100,
        )
        for i in range(3)
    ]
    annotation = select_annotation(protein, hmm, rps)
    path = tmp_path / "hits.xlsx"
    write_xlsx([annotation], path)
    sheet = load_workbook(path).active
    assert sheet["E2"].value == "HMM2 (PF2)"
    assert sheet["J2"].value == "HMM1 (PF1)"
    assert sheet["O2"].value == "RPS2 (cd2)"
    assert sheet["T2"].value == "RPS1 (cd1)"
    assert sheet["F2"].value == "Full HMM description 2"
    assert sheet["P2"].value == "Full RPS description 2; Reviewed"
    assert sheet["H2"].value == "4-90"
    assert sheet["R2"].value == "6-95"
    assert sheet["G2"].value == 1e-12
    assert sheet["Q2"].value == 1e-22
    from kannotate.io.genbank import GenBankInput

    record = SeqRecord(Seq("N" * 300), id="hit_record")
    record.features = [SeqFeature(SimpleLocation(0, 300, 1), type="CDS")]
    data = GenBankInput([record], [protein], [(0, 0)], [])
    gb_path = tmp_path / "genbank-hits.xlsx"
    write_xlsx([annotation], gb_path, genbank=data, input_name="hits.gb")
    gb_row = list(load_workbook(gb_path).active.values)[1]
    assert gb_row[7:10] == list(sheet.values)[1][1:4]
    assert gb_row[10:] == list(sheet.values)[1][4:]


def test_boundary_truncated_cds_derives_without_noisy_warnings(tmp_path):
    path = tmp_path / "boundary.gbk"
    records = sample_records()[:1]
    records[0].features.append(
        SeqFeature(SimpleLocation(0, 300), type="CDS", qualifiers={"locus_tag": ["boundary"]})
    )
    SeqIO.write(records, path, "genbank")
    data = read_genbank(path)
    assert len(data.proteins[-1].sequence) == 30
    assert data.warnings == []


def test_unreadable_location_is_retained_verbatim(tmp_path):
    path = tmp_path / "malformed.gbk"
    SeqIO.write(sample_records()[:1], path, "genbank")
    text = path.read_text()
    import re

    text = re.sub(r"(     CDS             )[^\n]+", r"\g<1>unreadable_location", text, count=1)
    path.write_text(text)
    data = read_genbank(path)
    assert data.records[0].features[1].location is None
    out = tmp_path / "output.gbk"
    write_genbank(data, [select_annotation(p, [], []) for p in data.proteins], out)
    assert "     CDS             unreadable_location\n" in out.read_text()
    assert out.read_text().count("     CDS ") == 3


def test_genbank_workbook_retains_source_context_and_all_cds_rows(tmp_path):
    source = tmp_path / "=sample.gbk"
    records = sample_records()
    records[0].features[1].qualifiers["product"] = ["=original product"]
    SeqIO.write(records, source, "genbank")
    original = source.read_bytes()
    data = read_genbank(source)
    run = AnnotationEngine(NoHits(), NoHits()).annotate(data.proteins)
    output = tmp_path / "annotations.xlsx"
    write_xlsx(run.annotations, output, genbank=data, input_name=source.name)
    sheet = load_workbook(output).active
    assert list(next(sheet.values)) == GENBANK_HEADERS
    assert GENBANK_HEADERS[:10] == [
        "Input filename",
        "Record ID",
        "Locus tag",
        "From",
        "To",
        "Strand",
        "Original product",
        "Normalized/search sequence",
        "Kannotate annotation",
        "Kannotate score",
    ]
    assert GENBANK_HEADERS[10:] == HEADERS[4:]
    assert (sheet.max_row, sheet.max_column) == (5, 30)
    rows = list(sheet.values)[1:]
    assert rows[0][:10] == (
        "=sample.gbk",
        "one",
        "A",
        1,
        27,
        "-",
        "=original product",
        "ACDEFGHIK*QQ",
        "hypothetical protein",
        0,
    )
    assert rows[1][3:6] == (1, 9, "+")
    assert rows[2][7] is None and rows[3][7] is None
    assert rows[3][1:7] == ("two", None, 1, 6, "+", None)
    assert sheet["D2"].comment and "join" in sheet["D2"].comment.text
    assert sheet["A2"].data_type == sheet["G2"].data_type == "s"
    assert sheet["J2"].number_format == "0.00"
    assert sheet["M2"].number_format == "0.00E+00"
    assert sheet.freeze_panes == "A2" and sheet.auto_filter.ref == "A1:AD5"
    assert source.read_bytes() == original


def test_genbank_workbook_uncertain_remote_and_unreadable_locations(tmp_path):
    from Bio.SeqFeature import AfterPosition, BeforePosition

    from kannotate.io.genbank import GenBankInput

    locations = [
        SimpleLocation(BeforePosition(3), AfterPosition(21), -1),
        CompoundLocation([SimpleLocation(80, 90, 1), SimpleLocation(0, 9, 1)]),
        CompoundLocation([SimpleLocation(0, 9, 1), SimpleLocation(18, 27, -1)]),
        SimpleLocation(0, 9, 1, ref="another_record"),
        None,
        SimpleLocation(9, 9, 1),
    ]
    record = SeqRecord(Seq("N" * 90), id="locations")
    record.features = [SeqFeature(location, type="CDS") for location in locations]
    proteins = [Protein(str(i), str(i), "") for i in range(len(locations))]
    data = GenBankInput([record], proteins, [(0, i) for i in range(len(locations))], [])
    out = tmp_path / "locations.xlsx"
    write_xlsx(
        [select_annotation(p, [], []) for p in proteins],
        out,
        genbank=data,
        input_name="locations.gb",
    )
    sheet = load_workbook(out).active
    rows = list(sheet.values)[1:]
    assert [r[3:6] for r in rows] == [
        (4, 21, "-"),
        (1, 90, "+"),
        (1, 27, None),
        (None, None, "+"),
        (None, None, None),
        (None, None, "+"),
    ]
    assert all(sheet.cell(i, 4).comment for i in range(2, 8))


@pytest.mark.parametrize(
    "mode,extensions",
    [
        ("genbank", [".gbk"]),
        ("xlsx", [".xlsx"]),
        ("both", [".gbk", ".xlsx"]),
    ],
)
def test_genbank_output_modes_search_once_and_report_every_file(tmp_path, mode, extensions):
    source = tmp_path / "input.gbk"
    SeqIO.write(sample_records(), source, "genbank")
    original = source.read_bytes()
    engine = AnnotationEngine(NoHits(), NoHits())
    calls = []

    class CountingEngine:
        def annotate(self, proteins, **kwargs):
            calls.append(proteins)
            return engine.annotate(proteins, **kwargs)

    reports = []
    result = annotate_batch(
        [source],
        None,
        CountingEngine(),
        genbank_output=mode,
        on_result=reports.append,
        table_format="xlsx",
    )[0]
    assert not result.error and len(calls) == 1 and reports == [result]
    assert [p.suffix for p in result.outputs] == extensions
    assert result.output == result.outputs[0]
    assert all(p.is_file() and p.parent == tmp_path for p in result.outputs)
    assert source.read_bytes() == original
    assert not list(tmp_path.glob(".kannotate-*"))
    if mode == "both":
        sheet = load_workbook(result.outputs[1]).active
        assert sheet["G2"].value == "old product"
        records = list(SeqIO.parse(result.outputs[0], "genbank"))
        assert records[0].features[1].qualifiers["product"] == ["hypothetical protein"]


def test_genbank_without_cds_still_writes_header_only_workbook(tmp_path):
    source = tmp_path / "empty.gbk"
    record = sample_records()[0]
    record.features = []
    SeqIO.write(record, source, "genbank")
    result = annotate_batch(
        [source],
        tmp_path / "out",
        AnnotationEngine(NoHits(), NoHits()),
        genbank_output="both",
        table_format="xlsx",
    )[0]
    assert not result.error
    sheet = load_workbook(result.outputs[1]).active
    assert (sheet.max_row, sheet.max_column) == (1, 30)


def test_both_preflights_existing_second_output_and_protects_every_input(tmp_path):
    source = tmp_path / "input.gbk"
    SeqIO.write(sample_records(), source, "genbank")
    output = tmp_path / "out"
    output.mkdir()
    existing = output / "input.kannotate.xlsx"
    existing.write_bytes(b"preexisting result")

    class MustNotSearch:
        def annotate(self, *args, **kwargs):
            raise AssertionError("Preflight must precede searching")

    result = annotate_batch(
        [source],
        output,
        MustNotSearch(),
        genbank_output="both",
        table_format="xlsx",
    )[0]
    assert "Output exists" in result.error
    assert list(output.iterdir()) == [existing]
    result = annotate_batch(
        [source, existing],
        output,
        MustNotSearch(),
        genbank_output="both",
        force=True,
        table_format="xlsx",
    )[0]
    assert "cannot replace an input" in result.error
    assert existing.read_bytes() == b"preexisting result"


@pytest.mark.parametrize("force", [False, True])
def test_both_rolls_back_partial_publication(tmp_path, monkeypatch, force):
    import os

    source = tmp_path / "input.gbk"
    SeqIO.write(sample_records(), source, "genbank")
    output = tmp_path / "out"
    output.mkdir()
    gb, xlsx = [output / f"input.kannotate{ext}" for ext in (".gbk", ".xlsx")]
    if force:
        gb.write_bytes(b"old genbank")
        xlsx.write_bytes(b"old workbook")
        replace = os.replace

        def fail_second(src, dst):
            if dst == xlsx and src.suffix == ".xlsx":
                raise OSError("injected publication failure")
            return replace(src, dst)

        monkeypatch.setattr("kannotate.batch.os.replace", fail_second)
    else:
        link = os.link

        def racing_second(src, dst):
            if dst == xlsx:
                xlsx.write_bytes(b"created by another writer during search")
            return link(src, dst)

        monkeypatch.setattr("kannotate.batch.os.link", racing_second)
    result = annotate_batch(
        [source],
        output,
        AnnotationEngine(NoHits(), NoHits()),
        genbank_output="both",
        force=force,
        table_format="xlsx",
    )[0]
    assert result.error and result.outputs == ()
    assert not list(output.glob(".kannotate-*"))
    if force:
        assert gb.read_bytes() == b"old genbank" and xlsx.read_bytes() == b"old workbook"
    else:
        assert not gb.exists()
        assert xlsx.read_bytes() == b"created by another writer during search"


@pytest.mark.parametrize("force", [False, True])
def test_stop_before_both_outputs_publish_preserves_originals(tmp_path, monkeypatch, force):
    from kannotate.cancellation import RunCancelled, RunControl

    source = tmp_path / "input.gbk"
    SeqIO.write(sample_records(), source, "genbank")
    original = source.read_bytes()
    output = tmp_path / "out"
    output.mkdir()
    if force:
        for suffix in (".gbk", ".xlsx"):
            (output / f"input.kannotate{suffix}").write_bytes(b"old output")
    control = RunControl()

    def stop_in_second_writer(annotations, path, **kwargs):
        path.write_bytes(b"unpublished workbook")
        control.request_stop()

    monkeypatch.setattr("kannotate.batch.write_xlsx", stop_in_second_writer)
    with pytest.raises(RunCancelled):
        annotate_batch(
            [source],
            output,
            AnnotationEngine(NoHits(), NoHits()),
            genbank_output="both",
            force=force,
            control=control,
            table_format="xlsx",
        )
    assert not list(output.glob(".kannotate-*"))
    assert source.read_bytes() == original
    assert len(list(output.iterdir())) == (2 if force else 0)
    assert all(p.read_bytes() == b"old output" for p in output.iterdir())


def test_both_reports_mixed_input_collisions_even_with_force(tmp_path):
    gb, fa = tmp_path / "same.gbk", tmp_path / "same.faa"
    SeqIO.write(sample_records(), gb, "genbank")
    fa.write_text(">protein\nMACDEFG\n")
    for index, inputs in enumerate(([gb, fa], [fa, gb])):
        out = tmp_path / f"out{index}"
        results = annotate_batch(
            inputs,
            out,
            AnnotationEngine(NoHits(), NoHits()),
            genbank_output="both",
            force=True,
            table_format="xlsx",
        )
        assert not results[0].error and "collision" in results[1].error
        assert len(list(out.iterdir())) == (2 if index == 0 else 1)


def test_cli_defaults_folder_and_reports_both_outputs(tmp_path, monkeypatch, capsys):
    from kannotate.cli import main

    source = tmp_path / "cluster.gbk"
    SeqIO.write(sample_records(), source, "genbank")
    monkeypatch.setattr(
        "kannotate.cli.AnnotationEngine.from_config", lambda **kwargs: AnnotationEngine(NoHits(), NoHits())
    )
    assert main(["annotate", str(source), "--genbank-output", "both"]) == 0
    printed = capsys.readouterr().out
    assert "cluster.kannotate.gbk, cluster.kannotate.csv" in printed
    assert "1 succeeded, 0 failed" in printed
