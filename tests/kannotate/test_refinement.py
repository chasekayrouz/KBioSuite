from dataclasses import replace

import pytest

from kannotate import Hit, Protein
from kannotate.architectures import classify_architecture
from kannotate.products import profile_product
from kannotate.scoring import DEFAULT_MIN_SCORE, score_hit
from kannotate.search import rpsblast
from kannotate.selection import select_annotation


def domains(profiles):
    intervals = {
        "PF00668": (1, 350),
        "PF00501": (400, 750),
        "PF00109": (800, 1150),
        "PF08659": (1300, 1500),
        "PF00550": (1700, 1760),
    }
    return [
        score_hit(Hit("hmm", "q", p, p, "domain", 1e-60, *intervals[p], 2200)) for p in profiles
    ]


@pytest.mark.parametrize(
    "profiles,product",
    [
        (["PF00109", "PF08659", "PF00550"], "polyketide synthase"),
        (
            ["PF00668", "PF00109", "PF08659", "PF00550"],
            "hybrid nonribosomal peptide synthetase-polyketide synthase",
        ),
    ],
)
def test_trans_at_architectures(profiles, product):
    hits = domains(profiles)
    assert classify_architecture(hits, 2200, hits[0]) == product
    assert classify_architecture(list(reversed(hits)), 2200, hits[0]) == product
    assert classify_architecture(hits, 1499, hits[0]) is None
    weak = [replace(h, score=1) if h.profile_id == "PF00550" else h for h in hits]
    assert classify_architecture(weak, 2200, weak[0]) is None


def test_trans_at_fragments_and_order_do_not_classify():
    hits = domains(["PF00109", "PF08659", "PF00550"])
    fragment = [
        replace(h, end=h.start + 10, spans=((h.start, h.start + 10),))
        if h.profile_id == "PF00550"
        else h
        for h in hits
    ]
    assert classify_architecture(fragment, 2200, fragment[0]) is None
    reversed_order = [
        replace(h, start=1200, end=1260, spans=((1200, 1260),)) if h.profile_id == "PF00550" else h
        for h in hits
    ]
    assert classify_architecture(reversed_order, 2200, reversed_order[0]) is None


def test_explicit_rps_arguments(tmp_path, monkeypatch):
    captured = []

    def fake(command):
        captured.extend(command)
        (tmp_path / "rps.tsv").write_text("")

    monkeypatch.setattr(rpsblast, "run_search", fake)
    adapter = rpsblast.RpsBlastAdapter("rpsblast", "Cdd")
    assert adapter.search(tmp_path / "queries.faa", tmp_path, 3) == []
    for name, value in {
        "-comp_based_stats": "1",
        "-seg": "no",
        "-evalue": "0.01",
        "-max_target_seqs": "100",
        "-mt_mode": "1",
        "-num_threads": "3",
    }.items():
        assert captured.count(name) == 1
        assert captured[captured.index(name) + 1] == value


def test_functional_prefix_preserved_through_selection():
    protein = Protein("q", "q", "M" * 400)
    hit = Hit(
        "rps",
        "q",
        "COG0438",
        "RfbF",
        "glycosyltransferase, activator-dependent family",
        1e-70,
        1,
        400,
        400,
    )
    annotation = select_annotation(protein, [], [hit])
    assert annotation.product == "glycosyltransferase, activator-dependent family"
    assert annotation.rps_hits[0].description == hit.description


def test_carrier_lookup_is_broad_and_limited_to_short_supported_proteins():
    hit = score_hit(
        Hit(
            "hmm",
            "q",
            "PF00550.32",
            "PP-binding",
            "Phosphopantetheine attachment site",
            1e-7,
            5,
            70,
            90,
        )
    )
    assert profile_product(hit, 90) == "phosphopantetheine-binding carrier protein"
    assert profile_product(hit, 1000) is None
    assert profile_product(replace(hit, score=3), 90) is None
    assert profile_product(replace(hit, profile_id="OTHER", name="OTHER"), 90) is None
    assert profile_product(replace(hit, end=30, spans=((5, 30),)), 90) is None
    protein = Protein("q", "q", "M" * 90)
    result = select_annotation(protein, [hit], [])
    assert result.selection == "lookup"
    assert result.hmm_hits[0].description == "Phosphopantetheine attachment site"


def test_sensitive_default_and_threshold_override():
    assert DEFAULT_MIN_SCORE == 6
    protein = Protein("q", "q", "M" * 100)
    hit = Hit("hmm", "q", "profile", "profile", "carrier protein", 1e-2, 1, 100, 100)
    assert select_annotation(protein, [hit], []).product == "carrier protein"
    assert select_annotation(protein, [hit], [], min_score=8).product == "hypothetical protein"
