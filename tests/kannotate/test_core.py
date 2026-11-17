import math
from dataclasses import replace

import pytest

from kannotate import AnnotationEngine, Hit, Protein
from kannotate.architectures import classify_architecture
from kannotate.descriptions import normalize_description
from kannotate.hits import collapse_profiles
from kannotate.normalization import normalize_sequence
from kannotate.scoring import coverage, score_hit
from kannotate.search.common import SearchError
from kannotate.selection import select_annotation


def hit(
    engine="hmm", profile="PF00001.1", start=1, end=100, length=100, evalue=1e-20, desc="enzyme"
):
    return Hit(engine, "q", profile, profile, desc, evalue, start, end, length)


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("  ma c\nDef Bjzoux ", "MACDEFBJZOUX"),
        ("MKT!QQ", "MKT"),
        ("?MKT", ""),
        ("MKT123", "MKT"),
        ("MKT-QQ", "MKT"),
        ("MKT*QQ", "MKT*QQ"),
        ("MKTéQQ", "MKT"),
        ("", ""),
    ],
)
def test_normalization(raw, expected):
    assert normalize_sequence(raw).sequence == expected


def test_scoring_explicit_formula_and_engine_independence():
    a = score_hit(hit())
    b = score_hit(hit("rps"))
    assert a.score == b.score == 25
    assert score_hit(hit(end=25)).score == pytest.approx(16.25)
    assert score_hit(hit(evalue=0)).score == 305
    assert score_hit(hit(evalue=1)).score == 5
    assert score_hit(hit(evalue=1e-30)).score > a.score
    for evalue in (math.nan, math.inf, -1):
        with pytest.raises(ValueError):
            score_hit(hit(evalue=evalue))


def test_profile_collapse_does_not_fill_unaligned_gap():
    assert normalize_sequence("BJZOUX").searchable
    assert not normalize_sequence("XXXXXX").searchable
    rows = [hit(start=1, end=20), hit(start=80, end=100), hit(start=10, end=25)]
    collapsed = collapse_profiles(rows)
    assert len(collapsed) == 1
    assert collapsed[0].spans == ((1, 25), (80, 100))
    assert (collapsed[0].start, collapsed[0].end) == (1, 100)
    assert coverage(collapsed[0]) == 0.46
    assert collapse_profiles(list(reversed(rows))) == collapsed


def test_threshold_ties_and_hypothetical():
    protein = Protein("q", "header", "M" * 100)
    hmm = hit(profile="Z", desc="broad family")
    rps = hit("rps", profile="A", desc="specific family")
    assert select_annotation(protein, [hmm], [rps]).product == "specific family"
    assert select_annotation(protein, [hmm], [rps], 26).product == "hypothetical protein"
    assert select_annotation(protein, [], []).score == 0
    assert (
        select_annotation(protein, [replace(hmm, description="")], []).product
        == "hypothetical protein"
    )


def test_conservative_description_cleanup():
    assert (
        normalize_description("glycosyltransferase, activator-dependent family")
        == "glycosyltransferase, activator-dependent family"
    )
    assert (
        normalize_description("methyltransferase, class I family")
        == "methyltransferase, class I family"
    )
    assert normalize_description("ferredoxin [4Fe-4S]") == "ferredoxin [4Fe-4S]"
    assert (
        normalize_description("FAD synthase [Coenzyme transport and metabolism]. ")
        == "FAD synthase"
    )
    assert (
        normalize_description(
            "Threonine synthase catalyzes the final step of threonine biosynthesis."
        )
        == "Threonine synthase"
    )
    assert (
        normalize_description("aspartate ammonia-lyase; Provisional") == "aspartate ammonia-lyase"
    )
    assert (
        normalize_description("cd12345, Acylase, Acylase domain. This is a long explanation.")
        == "Acylase domain"
    )
    assert normalize_description("  ATP-dependent   helicase... ") == "ATP-dependent helicase"
    assert (
        normalize_description("Domain of unknown function (DUF42)")
        == "Domain of unknown function (DUF42)"
    )


@pytest.mark.parametrize(
    "profiles,expected",
    [
        (["PF00501", "PF00668"], "nonribosomal peptide synthetase"),
        (["PF00109", "PF00698"], "polyketide synthase"),
        (
            ["PF00501", "PF00668", "PF00109", "PF00698"],
            "hybrid nonribosomal peptide synthetase-polyketide synthase",
        ),
        (
            ["PF00668", "PF00109", "PF00698"],
            "hybrid nonribosomal peptide synthetase-polyketide synthase",
        ),
    ],
)
def test_architecture_rules(profiles, expected):
    hits = [
        score_hit(hit(profile=p, start=1 + i * 500, end=400 + i * 500, length=2500, evalue=1e-60))
        for i, p in enumerate(profiles)
    ]
    assert classify_architecture(hits, 2500, hits[0]) == expected
    assert classify_architecture([replace(h, score=5) for h in hits], 2500, hits[0]) is None
    assert (
        classify_architecture(
            [replace(h, start=1, end=400, spans=((1, 400),)) for h in hits], 2500, hits[0]
        )
        is None
    )
    assert classify_architecture(hits, 300, hits[0]) is None


def test_full_profile_blocks_architecture_and_fragments_do_not_classify():
    domains = [
        score_hit(hit(profile=p, start=1 + i * 400, end=300 + i * 400, length=2000, evalue=1e-60))
        for i, p in enumerate(("PF00501", "PF00668"))
    ]
    best = score_hit(
        hit(profile="Other", length=2000, end=2000, evalue=1e-60, desc="large helicase")
    )
    assert classify_architecture(domains + [best], 2000, best) is None
    assert (
        classify_architecture([replace(h, end=h.start + 15) for h in domains], 2000, domains[0])
        is None
    )


class FakeAdapter:
    def search(self, queries, workdir, threads):
        text = queries.read_text()
        if "QQQQQQ" in text:
            raise SearchError("unusable sequence")
        return []


def test_engine_retains_all_records_and_isolates_failures():
    engine = AnnotationEngine(FakeAdapter(), FakeAdapter())
    proteins = [
        Protein("duplicate", "first", "acdefg"),
        Protein("duplicate", "empty", "?"),
        Protein("bad", "fails search", "QQQQQQ"),
        Protein("x", "short", "M"),
    ]
    result = engine.annotate(proteins, threads=2)
    assert len(result.annotations) == 4
    assert [a.protein.id for a in result.annotations] == [p.id for p in proteins]
    assert result.annotations[0].protein.sequence == "ACDEFG"
    assert all(a.product == "hypothetical protein" for a in result.annotations)
    assert result.warnings == [
        "hmmscan: 1 proteins could not be searched",
        "rpsblast: 1 proteins could not be searched",
    ]
