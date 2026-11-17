import csv

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from kannotate import AnnotationEngine, Protein
from kannotate.batch import annotate_batch
from kannotate.cancellation import RunCancelled, RunControl
from kannotate.io.csv import GENBANK_CSV_HEADERS, PROVENANCE_HEADERS, write_csv
from kannotate.io.genbank import annotated_records, read_genbank
from kannotate.io.tabular import HEADERS
from kannotate.selection import select_annotation


class NoHits:
    def search(self, *args):
        return []


def csv_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.reader(handle))


def test_csv_quotes_full_text_and_self_contained_provenance(tmp_path):
    headers = ["000123", "=1+1", "2026-01-01", 'comma, "quotes"\r\nβ\ttab', "L" * 40000]
    annotations = [
        select_annotation(Protein(str(i), h, "MACDEFG"), [], []) for i, h in enumerate(headers)
    ]
    path = tmp_path / "table.csv"
    write_csv(annotations, path, database_versions={"pfam": "38.2", "cdd": "3.21"})
    assert path.read_bytes().startswith(b"\xef\xbb\xbf")
    rows = csv_rows(path)
    assert rows[0] == HEADERS + PROVENANCE_HEADERS
    assert [row[0] for row in rows[1:]] == headers
    assert all(len(row) == 28 and row[-2:] == ["38.2", "3.21"] for row in rows[1:])
    assert all(row[4:24] == [""] * 20 for row in rows[1:])


def test_default_csv_streams_fasta_and_keeps_only_counts(tmp_path):
    path = tmp_path / "many.faa"
    path.write_text("".join(f">duplicate\n{'MACDEFG' if i % 2 else 'X'}\n" for i in range(1601)))
    engine = AnnotationEngine(NoHits(), NoHits(), database_versions={"pfam": "38.2", "cdd": "3.21"})
    original = engine.annotate
    windows, progress = [], []

    def annotate(proteins, **kwargs):
        windows.append(len(proteins))
        return original(proteins, **kwargs)

    engine.annotate = annotate
    result = annotate_batch(
        [path],
        tmp_path / "out",
        engine,
        retain_annotations=False,
        on_progress=lambda p, n: progress.append(n),
    )[0]
    assert not result.error and result.output.suffix == ".csv"
    assert windows == [512, 512, 512, 65]
    assert result.run.protein_count == 1601 and result.run.annotations == []
    assert progress == [0, 512, 1024, 1536, 1601]
    assert len(csv_rows(result.output)) == 1602
    assert len(list(result.output.parent.iterdir())) == 1


def test_genbank_csv_lengths_use_proteins_and_keep_empty_cds_rows(tmp_path):
    record = SeqRecord(
        Seq("ATGAAATTTTAA"), id="record", annotations={"molecule_type": "DNA"}
    )
    record.features = [
        SeqFeature(
            SimpleLocation(0, 12, strand=1),
            type="CDS",
            qualifiers={} if translation is None else {"translation": [translation]},
        )
        for translation in ["MACDEFG*", "MA*CX", "", None, "ma cd!fg"]
    ]
    record.features.append(
        SeqFeature(SimpleLocation(0, 12, strand=1), type="CDS", qualifiers={"pseudo": [""]})
    )
    source = tmp_path / "input.gbk"
    SeqIO.write([record], source, "genbank")
    original = source.read_bytes()
    data = read_genbank(source)
    annotations = AnnotationEngine(NoHits(), NoHits()).annotate(data.proteins).annotations
    output = tmp_path / "table.csv"
    write_csv(
        annotations,
        output,
        genbank=data,
        input_name=source.name,
        database_versions={"pfam": "38.2", "cdd": "3.21"},
    )
    table = csv_rows(output)
    assert table[0][5:9] == ["Strand", "Length", "Original product", "Normalized/search sequence"]
    assert len(table[0]) == 31
    assert not set(["Location note"] + PROVENANCE_HEADERS).intersection(table[0])
    assert [row[6] for row in table[1:]] == ["7", "4", "0", "3", "4", "0"]
    assert [row[8] for row in table[1:]] == ["MACDEFG*", "MA*CX", "", "MKF", "MACD", ""]
    assert all(len(row) == 31 and row[11:] == [""] * 20 for row in table[1:])
    assert source.read_bytes() == original

    record.features.clear()
    SeqIO.write([record], source, "genbank")
    data = read_genbank(source)
    write_csv([], output, genbank=data, input_name=source.name)
    assert csv_rows(output) == [GENBANK_CSV_HEADERS]


def records(count):
    for i in range(count):
        record = SeqRecord(Seq("ATG" * 12), id=f"record_{i}", annotations={"molecule_type": "DNA"})
        location = CompoundLocation(
            [SimpleLocation(0, 6, strand=1), SimpleLocation(12, 18, strand=1)]
        )
        record.features = [
            SeqFeature(
                location,
                type="CDS",
                qualifiers={
                    "product": [f"original product {i}"],
                    "translation": ["MACDEFG"],
                    "locus_tag": [f"000{i}"],
                },
            )
        ]
        yield record


@pytest.mark.parametrize("mode", ["csv", "genbank", "both"])
@pytest.mark.parametrize("retain_annotations", [False, True])
def test_genbank_streaming_preserves_order_metadata_and_input(tmp_path, mode, retain_annotations):
    path = tmp_path / "many.gbff"
    SeqIO.write(records(601), path, "genbank")
    before = path.read_bytes()
    result = annotate_batch(
        [path],
        tmp_path / "out",
        AnnotationEngine(NoHits(), NoHits()),
        genbank_output=mode,
        retain_annotations=retain_annotations,
    )[0]
    assert not result.error and result.run.protein_count == 601
    assert len(result.run.annotations) == (601 if retain_annotations else 0)
    assert path.read_bytes() == before
    assert len(result.outputs) == (2 if mode == "both" else 1)
    if mode in ("csv", "both"):
        table = csv_rows(result.output.parent / "many.kannotate.csv")
        assert table[0] == GENBANK_CSV_HEADERS
        assert len(table) == 602
        assert [row[1] for row in table[1:]] == [f"record_{i}" for i in range(601)]
        assert [row[7] for row in table[1:]] == [f"original product {i}" for i in range(601)]
        assert all(len(row) == 31 and row[6] == "7" for row in table[1:])
    if mode in ("genbank", "both"):
        output = list(SeqIO.parse(result.output.parent / "many.kannotate.gbk", "genbank"))
        assert len(output) == 601
        assert (
            sum("Annotated with Kannotate" in r.annotations.get("comment", "") for r in output) == 1
        )
        assert all(r.features[0].qualifiers["product"] == ["hypothetical protein"] for r in output)


def test_same_count_wrong_genbank_annotation_order_fails(tmp_path):
    source = tmp_path / "input.gbk"
    SeqIO.write(records(2), source, "genbank")
    data = read_genbank(source)
    annotations = AnnotationEngine(NoHits(), NoHits()).annotate(data.proteins).annotations
    with pytest.raises(ValueError, match="order"):
        annotated_records(data, list(reversed(annotations)))


def test_stop_after_partial_stream_keeps_existing_output_and_source(tmp_path):
    path = tmp_path / "input.faa"
    path.write_text(">p\nMACDEFG\n" * 600)
    original = path.read_bytes()
    out = tmp_path / "out"
    out.mkdir()
    destination = out / "input.kannotate.csv"
    destination.write_bytes(b"existing output")
    control = RunControl()

    def stop(path, count):
        if count:
            control.request_stop()

    with pytest.raises(RunCancelled):
        annotate_batch(
            [path],
            out,
            AnnotationEngine(NoHits(), NoHits()),
            retain_annotations=False,
            force=True,
            control=control,
            on_progress=stop,
        )
    assert path.read_bytes() == original
    assert destination.read_bytes() == b"existing output"
    assert list(out.iterdir()) == [destination]
