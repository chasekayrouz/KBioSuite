import io
import tarfile

import pytest

from kannotate.databases import safe_extract, specification
from kannotate.dependencies import resolve_executable
from kannotate.downloads import fetch
from kannotate.search.hmmer import parse_domtbl
from kannotate.search.rpsblast import parse_rps


def test_hmmer_profile_collapse_and_query_coordinates(tmp_path):
    table = tmp_path / "hits"
    table.write_text(
        "# comment\n"
        + "AMP-binding PF00501.1 400 q - 2000 1e-60 210 0 1 2 1e-30 1e-30 100 0 1 380 101 480 90 490 0.99 AMP-binding enzyme\n"
        + "AMP-binding PF00501.1 400 q - 2000 1e-60 210 0 2 2 1e-30 1e-30 100 0 1 380 1001 1380 990 1390 0.99 AMP-binding enzyme\n"
    )
    hits = parse_domtbl(table)
    assert len(hits) == 1
    assert hits[0].spans == ((101, 480), (1001, 1380))
    assert hits[0].evalue == 1e-60
    assert hits[0].query_length == 2000
    assert hits[0].bitscore == 210  # full-sequence score, not each domain's score


def test_rps_descriptions_and_profile_collapse(tmp_path):
    table = tmp_path / "hits"
    table.write_text(
        "q\tgnl|CDD|123\t1e-20\t80\t4\t100\t200\t150\t97\tfull title\n"
        + "q\tgnl|CDD|123\t1e-10\t50\t110\t180\t200\t150\t71\tfull title\n"
    )
    hits = parse_rps(table, {"123": ("cd42", "enzyme", "Enzyme full description")})
    assert len(hits) == 1
    assert hits[0].profile_id == "cd42"
    assert (hits[0].start, hits[0].end) == (4, 180)
    assert hits[0].description == "Enzyme full description"
    assert hits[0].evalue == 1e-20 and hits[0].bitscore == 80
    assert hits[0].spans == ((4, 100), (110, 180))


def test_explicit_dependencies_do_not_silently_fall_back(tmp_path, monkeypatch):
    executable = tmp_path / "hmmscan"
    executable.write_text("#!/bin/sh\nexit 0\n")
    executable.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    assert resolve_executable("hmmscan", "/wrong/path") is None
    assert resolve_executable("hmmscan") == str(executable)
    explicit = tmp_path / "pinned" / "hmmscan"
    explicit.parent.mkdir()
    explicit.write_text("#!/bin/sh\nexit 0\n")
    explicit.chmod(0o755)
    assert resolve_executable("hmmscan", str(explicit)) == str(explicit)


def test_archive_traversal_rejected(tmp_path):
    archive = tmp_path / "bad.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        member = tarfile.TarInfo("../outside")
        member.size = 1
        tar.addfile(member, io.BytesIO(b"x"))
    with pytest.raises(RuntimeError, match="Unsafe"):
        safe_extract(archive, tmp_path / "out")
    assert not (tmp_path / "outside").exists()


def test_checksum_pin_refuses_changed_content(tmp_path):
    source = tmp_path / "source"
    source.write_text("changed")
    with pytest.raises(RuntimeError, match="Checksum mismatch"):
        fetch(source.as_uri(), tmp_path / "download", sha256="0" * 64)
    assert not (tmp_path / "download").exists()
    spec = specification()
    assert spec["cdd"]["version"] == "3.21"
    assert all(len(a["sha256"]) == 64 for e in (spec["pfam"], spec["cdd"]) for a in e["assets"])


def test_rps_gapped_alignment_coordinates_and_repeated_hsps(tmp_path):
    table = tmp_path / "aligned.tsv"
    table.write_text(
        "q\tgnl|CDD|123\t1e-20\t80\t4\t9\t200\t150\t7\tfull title\t11\t16\tAC-DEFG\tACD-EFG\n"
        "q\tgnl|CDD|123\t1e-10\t50\t110\t113\t200\t150\t4\tfull title\t20\t23\tAAAA\tAADA\n"
    )
    hit = parse_rps(table, include_alignments=True)[0]
    assert hit.spans == ((4, 9), (110, 113))
    assert hit.evalue == 1e-20 and hit.bitscore == 80
    assert len(hit.alignments) == 2
    a, b = hit.alignments
    assert (a.query_start, a.query_end, a.model_start, a.model_end) == (4, 9, 11, 16)
    assert a.query == "AC-DEFG" and a.model == "ACD-EFG"
    assert b.evalue == 1e-10 and b.bitscore == 50
    table.write_text(table.read_text().replace("AC-DEFG", "ACDEFG"))
    with pytest.raises(ValueError, match="alignment coordinates"):
        parse_rps(table, include_alignments=True)


def test_rps_alignment_output_is_opt_in_and_keeps_search_arguments(tmp_path, monkeypatch):
    from kannotate.search.rpsblast import RPS_ARGUMENTS, RpsBlastAdapter

    commands = []

    def run(command):
        commands.append(command)
        (tmp_path / "rps.tsv").write_text("")

    monkeypatch.setattr("kannotate.search.rpsblast.run_search", run)
    plain = RpsBlastAdapter("rpsblast", "/db/Cdd")
    aligned = RpsBlastAdapter("rpsblast", "/db/Cdd", include_alignments=True)
    for adapter in (plain, aligned):
        assert adapter.search(tmp_path / "q.faa", tmp_path, 4) == []
    assert commands[0][:-4] == commands[1][:-4]
    assert all(arg in commands[0] for arg in RPS_ARGUMENTS)
    assert "qseq sseq" not in plain.output_format
    assert aligned.output_format.endswith("sstart send qseq sseq")
