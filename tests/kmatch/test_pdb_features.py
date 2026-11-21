"""Public metadata, frame reconciliation, ligand filtering and residue correspondence."""

import copy
import http.client
import json
import re
import ssl
import urllib.error
from dataclasses import asdict
from email.message import Message
from pathlib import Path

import gemmi
import numpy as np
import pytest
from test_alignment import case as case
from test_alignment import make_structure

from kmatch import alignment, hit_structure, pdb_data
from kmatch.config import Config
from kmatch.gui import Job, create_app
from kmatch.process import Cancellation, Cancelled


def entry():
    return {
        "rcsb_id": "1ABC",
        "struct": {"title": "A public structure"},
        "polymer_entities": [
            {
                "rcsb_polymer_entity": {"pdbx_description": "Enzyme alpha"},
                "rcsb_polymer_entity_container_identifiers": {
                    "auth_asym_ids": ["B"],
                    "asym_ids": ["X"],
                    "entity_id": "1",
                },
            }
        ],
        "nonpolymer_entities": [
            {
                "pdbx_entity_nonpoly": {"comp_id": code, "name": name},
                "rcsb_nonpolymer_entity_container_identifiers": {
                    "auth_asym_ids": ["B"],
                    "asym_ids": ["Y"],
                },
                "nonpolymer_comp": {"chem_comp": {"id": code, "name": name, "formula": None}},
            }
            for code, name in [("ZN", "Zinc ion"), ("HOH", "Water"), ("ATP", "ATP")]
        ],
    }


def test_metadata_name_chain_ligands_cache_and_public_only_requests(tmp_path, monkeypatch):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    calls = []

    def fetch(url, *args):
        calls.append(url)
        return json.dumps({"data": {"entries": [entry()]}}).encode()

    monkeypatch.setattr(pdb_data, "fetch", fetch)
    target = "1abc-assembly1_B"
    result = pdb_data.load_details(cfg, [target], True, Cancellation())[target]
    assert result["protein_name"] == "Enzyme alpha"
    assert [ligand["id"] for ligand in result["ligands"]] == ["ZN", "ATP"]
    assert result["url"] == "https://www.rcsb.org/structure/1ABC"
    assert pdb_data.load_details(cfg, [target], False, Cancellation())[target] == result
    assert len(calls) == 1 and "1ABC" in calls[0] and "sequence" not in calls[0]
    assert (
        pdb_data.entry_for_hit(pdb_data.normalize_entry(entry()), "1abc-assembly1_B-2")[
            "protein_name"
        ]
        == "Enzyme alpha"
    )
    assert (
        pdb_data.entry_for_hit(pdb_data.normalize_entry(entry()), "1abc-assembly1_C")[
            "protein_name"
        ]
        is None
    )


@pytest.mark.parametrize(
    "target", ["../secret", "1abc_A?x=secret", "AF-Q1234", "https://evil/1abc_A"]
)
def test_identifiers_cannot_supply_arbitrary_urls(target):
    with pytest.raises(ValueError):
        pdb_data.pdb_url(target)


def refuse(reason):
    def urlopen(request, *, timeout):
        raise urllib.error.URLError(reason)

    return urlopen


def test_network_failure_is_unavailable_not_no_ligands(tmp_path, monkeypatch):
    cfg = Config(tmp_path, tmp_path / "cache")
    monkeypatch.setattr(pdb_data.tls, "urlopen", refuse(OSError("offline")))
    result = pdb_data.load_details(cfg, ["1abc_A"], True, Cancellation())["1abc_A"]
    assert result["status"] == "unavailable" and result["retryable"]
    assert result["message"].startswith("RCSB could not be reached.")
    assert "ligands" not in result
    cancel = Cancellation()
    cancel.cancel()
    with pytest.raises(Cancelled):
        pdb_data.load_details(cfg, ["1abc_A"], False, cancel)
    with pytest.raises(RuntimeError, match="not cached"):
        pdb_data.assembly_file(cfg, "1abc-assembly1_A", False, Cancellation(), lambda e: None)


def certificate_failure():
    return ssl.SSLCertVerificationError(
        1, "[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local "
        "issuer certificate (_ssl.c:1085)",
    )


def test_fetch_supplies_a_verified_truststore_context_to_urllib(monkeypatch):
    """Record the context passed to urllib's HTTPS handler before any connection."""
    import truststore

    seen = []

    class Recorder(pdb_data.urllib.request.HTTPSHandler):
        def __init__(self, context):
            seen.append(context)
            super().__init__(context=context)

        def https_open(self, request):
            raise urllib.error.URLError("no network in tests")

    monkeypatch.setattr(pdb_data.urllib.request, "HTTPSHandler", Recorder)
    with pytest.raises(pdb_data.RequestFailed):
        pdb_data.fetch("https://data.rcsb.org/graphql", Cancellation(), 1024)
    assert len(seen) == 1 and isinstance(seen[0], truststore.SSLContext)
    assert seen[0].verify_mode == ssl.CERT_REQUIRED and seen[0].check_hostname


def test_certificate_failure_is_explained_retryable_and_never_cached(tmp_path, monkeypatch):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    cached, missing = "1abc-assembly1_B", "2xyz_A"
    monkeypatch.setattr(
        pdb_data, "fetch", lambda *a: json.dumps({"data": {"entries": [entry()]}}).encode()
    )
    ready = pdb_data.load_details(cfg, [cached], True, Cancellation())[cached]
    monkeypatch.undo()

    monkeypatch.setattr(pdb_data.tls, "urlopen", refuse(certificate_failure()))
    result = pdb_data.load_details(cfg, [cached, missing], True, Cancellation())
    assert result[cached] == ready  # an earlier success survives a later failure
    failed = result[missing]
    assert failed["status"] == "unavailable" and failed["retryable"]
    assert failed["message"].startswith("macOS could not verify RCSB's security certificate.")
    assert "Reload details" in failed["message"]
    assert "urlopen" not in failed["message"] and "_ssl.c" not in failed["message"]
    # Technical detail returned separately from the user message.
    assert "CERTIFICATE_VERIFY_FAILED" in failed["detail"]
    assert sorted(p.name for p in (cfg.pdb_cache / "metadata").iterdir()) == ["1ABC.json"]
    monkeypatch.undo()

    # Reload details: only the entry that failed is requested again, and now it is cached.
    requested = []

    def fetch(url, *args):
        requested.append(url)
        return json.dumps({"data": {"entries": [{**entry(), "rcsb_id": "2XYZ"}]}}).encode()

    monkeypatch.setattr(pdb_data, "fetch", fetch)
    reloaded = pdb_data.load_details(cfg, [cached, missing], True, Cancellation())
    assert reloaded[missing]["status"] == "ready" and "retryable" not in reloaded[missing]
    assert len(requested) == 1 and "2XYZ" in requested[0] and "1ABC" not in requested[0]
    # Cached data only: no request at all.
    assert pdb_data.load_details(cfg, [cached, missing], False, Cancellation()) == reloaded
    assert len(requested) == 1


@pytest.mark.parametrize(
    "error, expected",
    [
        (urllib.error.HTTPError("https://data.rcsb.org/", 503, "busy", {}, None), "HTTP 503"),
        (TimeoutError("timed out"), "could not be reached"),
        (urllib.error.URLError(certificate_failure()), "security certificate"),
        (certificate_failure(), "security certificate"),
        (http.client.IncompleteRead(b"partial", 10), "interrupted"),
    ],
)
def test_failure_messages_include_the_expected_cause_and_omit_reload_advice(error, expected):
    message = pdb_data.failure_message(error)
    assert expected in message and "<" not in message and "Errno" not in message
    # Caller-specific retry instructions are checked in the download tests.
    assert "Reload details" not in message


class Body:
    """A response whose body ends as scripted: bytes are returned, exceptions raised."""

    def __init__(self, *reads, length=None, chunked=False):
        self.reads = list(reads)
        self.headers = Message()
        for value in [] if length is None else [length] if isinstance(length, str | int) else length:
            self.headers["Content-Length"] = str(value)
        if chunked:
            self.headers["Transfer-Encoding"] = "chunked"

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def read(self, amount):
        step = self.reads.pop(0) if self.reads else b""
        if isinstance(step, Exception):
            raise step
        return step


def serve(monkeypatch, make_body):
    monkeypatch.setattr(pdb_data.tls, "urlopen", lambda request, *, timeout: make_body())


@pytest.mark.parametrize(
    "make_body",
    [
        # Simulate IncompleteRead from an interrupted chunk read.
        lambda: Body(b'{"data": ', http.client.IncompleteRead(b"", 40)),
        lambda: Body(b'{"data": ', length=64),  # fixed length, early end of file
        lambda: Body(b'{"data": ', length="64, 64"),  # the same length joined by an intermediary
        lambda: Body(b'{"data": ', length=[64, 64]),  # ... or repeated
        lambda: Body(b'{"data": ', length="9, 64"),  # conflicting lengths, even if one fits
        lambda: Body(b'{"data": ', length="nine"),
    ],
)
def test_detected_interruption_is_a_retryable_failure_and_is_not_cached(
    tmp_path, monkeypatch, make_body
):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    serve(monkeypatch, make_body)
    with pytest.raises(pdb_data.RequestFailed, match="interrupted"):
        pdb_data.fetch("https://data.rcsb.org/graphql", Cancellation(), 1024)
    result = pdb_data.load_details(cfg, ["1abc_A"], True, Cancellation())["1abc_A"]
    assert result["retryable"] and result["message"] == (
        "RCSB's response was interrupted. Choose Reload details to try again."
    )
    with pytest.raises(RuntimeError, match="interrupted. Export the 3D alignment again"):
        pdb_data.assembly_file(cfg, "1abc-assembly1_A", True, Cancellation(), lambda e: None)
    assert not [p for p in cfg.pdb_cache.rglob("*") if p.is_file()]


@pytest.mark.parametrize(
    "framing",
    [
        {"length": 6},
        {"length": "6, 6"},
        {"length": [6, 6]},
        {},  # no declared length: nothing to compare
        {"length": 99, "chunked": True},  # chunked framing is http.client's to check
    ],
)
def test_complete_response_is_accepted(monkeypatch, framing):
    serve(monkeypatch, lambda: Body(b"abc", b"def", **framing))
    assert pdb_data.fetch("https://data.rcsb.org/graphql", Cancellation(), 1024) == b"abcdef"


def test_parseable_coordinate_prefix_shorter_than_declared_length_is_not_cached(
    tmp_path, monkeypatch
):
    cfg = Config(tmp_path, tmp_path / "cache")
    structure = make_structure("ACDKEFG", np.arange(21, dtype=float).reshape(7, 3))
    structure.name = "1ABC"
    whole = structure.make_mmcif_document().as_string().encode()
    cut = whole.rstrip().rsplit(b"\n", 2)[0] + b"\n"  # the last two atom rows never arrived
    # Control: the prefix alone would pass every content check, so only its length gives it away.
    block = gemmi.cif.read_string(cut.decode()).sole_block()
    assert block.find_value("_entry.id") == "1ABC" and len(cut) < len(whole)
    assert len(gemmi.make_structure_from_block(block)) == 1
    assert len(gemmi.make_structure_from_block(block)[0][0]) == 5
    serve(monkeypatch, lambda: Body(cut, length=f"{len(whole)}, {len(whole)}"))
    with pytest.raises(RuntimeError, match="interrupted"):
        pdb_data.assembly_file(cfg, "1abc_A", True, Cancellation(), lambda e: None)
    assert not [p for p in cfg.pdb_cache.rglob("*") if p.is_file()]
    # The complete file, with the same framing, is saved.
    serve(monkeypatch, lambda: Body(whole, length=f"{len(whole)}, {len(whole)}"))
    path, _ = pdb_data.assembly_file(cfg, "1abc_A", True, Cancellation(), lambda e: None)
    assert path.read_bytes() == whole


def test_cancellation_wins_over_an_interrupted_response(monkeypatch):
    cancel = Cancellation()

    class Cancelling(Body):
        def read(self, amount):
            cancel.cancel()
            raise http.client.IncompleteRead(b"", 10)

    serve(monkeypatch, Cancelling)
    with pytest.raises(Cancelled):
        pdb_data.fetch("https://data.rcsb.org/graphql", cancel, 1024)


def answer(*entries):
    return lambda *a: json.dumps({"data": {"entries": list(entries)}}).encode()


def test_null_author_chain_and_ligand_instance_lists_survive_cache_reload(tmp_path, monkeypatch):
    """Keep the label-chain protein name and an unknown ligand count across a cache read."""
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    raw = entry()
    identifiers = raw["polymer_entities"][0]["rcsb_polymer_entity_container_identifiers"]
    identifiers["auth_asym_ids"] = None
    raw["nonpolymer_entities"][0]["rcsb_nonpolymer_entity_container_identifiers"]["asym_ids"] = None
    monkeypatch.setattr(pdb_data, "fetch", answer(raw))
    result = pdb_data.load_details(cfg, ["1abc_X"], True, Cancellation())["1abc_X"]
    assert result["status"] == "ready" and result["protein_name"] == "Enzyme alpha"  # label chain
    # Zinc keeps its component ID and an unknown instance count.
    assert [(ligand["id"], ligand["copies"]) for ligand in result["ligands"]] == [
        ("ZN", None), ("ATP", 1),
    ]
    assert pdb_data.load_details(cfg, ["1abc_X"], False, Cancellation())["1abc_X"] == result
    saved = json.loads((cfg.pdb_cache / "metadata" / "1ABC.json").read_text())
    assert saved["ligands"][0]["copies"] is None and saved["ligands"][1]["copies"] == 1
    # A list RCSB does give is counted, including an empty one.
    raw["nonpolymer_entities"][0]["rcsb_nonpolymer_entity_container_identifiers"]["asym_ids"] = []
    assert pdb_data.normalize_entry(raw)["ligands"][0]["copies"] == 0


def selected_fields():
    """Build field paths and a matching answer for the current entries selection."""

    def brace(tokens):  # consume "{", parse to the matching "}"
        tokens.pop(0)
        return parse(tokens)

    def parse(tokens):
        fields = {}
        while tokens:
            token = tokens.pop(0)
            if token == "}":
                break
            fields[token] = brace(tokens) if tokens and tokens[0] == "{" else None
        return fields

    body = pdb_data.GRAPHQL.split("{", 2)[2]  # after "query { entries(entry_ids: IDS) {"
    selection = parse(re.findall(r"[{}]|\w+", body))
    lists = {"polymer_entities", "nonpolymer_entities"}
    values = {"rcsb_id": "1ABC", "auth_asym_ids": ["A"], "asym_ids": ["A"]}

    def build(fields):
        made = {}
        for name, inner in fields.items():
            value = values.get(name, name) if inner is None else build(inner)
            made[name] = [value] if name in lists else value
        return made

    def paths(fields, prefix=()):
        for name, inner in fields.items():
            yield (*prefix, name)
            if inner is not None:
                yield from paths(inner, (*prefix, name, 0) if name in lists else (*prefix, name))

    return list(paths(selection)), build(selection)


def test_every_selected_field_is_required_and_type_checked():
    """For each field in the current entries selection, check omission and a numeric value."""
    fields, complete = selected_fields()
    assert len(fields) == 22 and ("nonpolymer_entities", 0, "nonpolymer_comp", "chem_comp", "formula") in fields
    record = pdb_data.normalize_entry(complete)
    assert record["polymers"][0]["name"] == "pdbx_description"
    assert record["ligands"] == [
        {"id": "id", "name": "name", "formula": "formula", "author_chains": ["A"], "copies": 1}
    ]
    for path in fields:
        for change in ("omitted", "a number"):
            raw = copy.deepcopy(complete)
            parent = raw
            for step in path[:-1]:
                parent = parent[step]
            if change == "omitted":
                del parent[path[-1]]
            else:
                parent[path[-1]] = 7
            with pytest.raises(ValueError):
                pytest.fail(f"{path} {change}: accepted as {pdb_data.normalize_entry(raw)}")


def incomplete_answers():
    """Examples with missing fields, wrong types, or no usable ligand ID."""
    no_ligand_field = entry()
    del no_ligand_field["nonpolymer_entities"]
    text_instances = entry()
    text_instances["nonpolymer_entities"][0]["rcsb_nonpolymer_entity_container_identifiers"][
        "asym_ids"
    ] = "AA"  # Taking the length of this string would give a spurious count of two.
    no_polymer_field = entry()
    del no_polymer_field["polymer_entities"]
    no_component = entry()
    del no_component["nonpolymer_entities"][0]["nonpolymer_comp"]
    unnamed = entry()
    unnamed["nonpolymer_entities"][0]["nonpolymer_comp"] = None
    unnamed["nonpolymer_entities"][0]["pdbx_entity_nonpoly"] = None
    empty_component = entry()
    empty_component["nonpolymer_entities"][0]["nonpolymer_comp"] = {"chem_comp": {}}
    numeric_name = entry()
    numeric_name["polymer_entities"][0]["rcsb_polymer_entity"]["pdbx_description"] = 7
    return {
        "no ligand field": no_ligand_field,
        "text instance IDs": text_instances,
        "no polymer field": no_polymer_field,
        "no component field": no_component,
        "ligand without an ID": unnamed,
        "component without its fields": empty_component,
        "numeric name": numeric_name,
        "no title field": {k: v for k, v in entry().items() if k != "struct"},
    }


@pytest.mark.parametrize("case", incomplete_answers())
def test_incomplete_answer_is_never_saved_and_reload_gets_the_complete_one(
    tmp_path, monkeypatch, case
):
    """An answer without its ligand field must not become a saved "no ligands" result."""
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    with pytest.raises(ValueError):
        pdb_data.normalize_entry(incomplete_answers()[case])
    monkeypatch.setattr(pdb_data, "fetch", answer(incomplete_answers()[case]))
    first = pdb_data.load_details(cfg, ["1abc_B"], True, Cancellation())["1abc_B"]
    assert first["retryable"] and "ligands" not in first
    assert first["message"] == "RCSB returned an unexpected response. Choose Reload details to try again."
    assert not [p for p in cfg.pdb_cache.rglob("*") if p.is_file()]
    requested = []
    monkeypatch.setattr(pdb_data, "fetch", lambda url, *a: requested.append(url) or answer(entry())())
    reloaded = pdb_data.load_details(cfg, ["1abc_B"], True, Cancellation())["1abc_B"]
    assert len(requested) == 1 and [ligand["id"] for ligand in reloaded["ligands"]] == ["ZN", "ATP"]
    assert [ligand["copies"] for ligand in reloaded["ligands"]] == [1, 1]


def test_nullable_descriptions_and_null_or_empty_ligand_collections_are_accepted(
    tmp_path, monkeypatch
):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    raw = entry()
    raw.update(struct=None, nonpolymer_entities=None)
    raw["polymer_entities"][0]["rcsb_polymer_entity"] = None
    monkeypatch.setattr(pdb_data, "fetch", answer(raw, {**entry(), "rcsb_id": "2XYZ", "nonpolymer_entities": []}))
    result = pdb_data.load_details(cfg, ["1abc_B", "2xyz_B"], True, Cancellation())
    assert result["1abc_B"]["status"] == "ready" and result["1abc_B"]["ligands"] == []
    assert result["1abc_B"]["title"] is None and result["1abc_B"]["protein_name"] is None
    assert result["2xyz_B"]["ligands"] == [] and result["2xyz_B"]["protein_name"] == "Enzyme alpha"


def legacy_record(**polymer):
    """A record as Kmatch 0.6.1 saved it: RCSB's null chain lists were stored as null."""
    return {
        "pdb_id": "1ABC",
        "title": "A public structure",
        "polymers": [
            {"name": "Enzyme alpha", "author_chains": ["A"], "label_chains": None,
             "entity_id": "1", **polymer}
        ],
        "ligands": [
            {"id": "ZN", "name": "ZINC ION", "formula": "Zn", "author_chains": None, "copies": 1}
        ],
        "source": "https://data.rcsb.org/graphql",
        "retrieved_utc": "2026-11-01T00:00:00Z",
    }


def test_legacy_format_record_stays_visible_without_a_fetch(tmp_path, monkeypatch):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    path = cfg.pdb_cache / "metadata" / "1ABC.json"
    saved = json.dumps(legacy_record(), indent=2).encode()
    pdb_data.cache_write(path, saved)
    monkeypatch.setattr(pdb_data, "fetch", lambda *a: pytest.fail("a readable record was refetched"))
    for allow_network in (False, True):
        result = pdb_data.load_details(cfg, ["1abc_A"], allow_network, Cancellation())["1abc_A"]
        assert result["status"] == "ready" and result["protein_name"] == "Enzyme alpha"
        assert [ligand["id"] for ligand in result["ligands"]] == ["ZN"]
    assert path.read_bytes() == saved  # Cache contents still match the original bytes.


def test_unreadable_cached_record_is_requested_again(tmp_path, monkeypatch):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    path = cfg.pdb_cache / "metadata" / "1ABC.json"
    good = legacy_record()
    negative = copy.deepcopy(good)
    negative["ligands"][0]["copies"] = -1
    unreadable_records = (
        legacy_record(author_chains="A"),  # text where a list of chains belongs
        legacy_record(label_chains=[1]),
        negative,
        [good],
        {**good, "ligands": None},
        {**good, "pdb_id": "9XYZ"},
        {**good, "retrieved_utc": 20261101},
    )
    requested = []
    monkeypatch.setattr(
        pdb_data, "fetch", lambda url, *a: requested.append(url) or answer(entry())()
    )
    for unreadable in unreadable_records:
        pdb_data.cache_write(path, json.dumps(unreadable).encode())
        assert pdb_data.cached_entry(path, "1ABC") is None
        requested.clear()
        # Cached data only: reported as not cached rather than raising for every hit.
        offline = pdb_data.load_details(cfg, ["1abc_B"], False, Cancellation())["1abc_B"]
        assert offline["message"] == "PDB details are not cached." and not requested
        result = pdb_data.load_details(cfg, ["1abc_B"], True, Cancellation())["1abc_B"]
        assert len(requested) == 1 and result["protein_name"] == "Enzyme alpha"
        assert pdb_data.cached_entry(path, "1ABC") is not None
    # The cache reader accepts this record without a retrieval timestamp.
    untimed = {k: v for k, v in good.items() if k != "retrieved_utc"}
    pdb_data.cache_write(path, json.dumps(untimed).encode())
    assert pdb_data.cached_entry(path, "1ABC") is not None


def named(identifier):
    return {**entry(), "rcsb_id": identifier}


@pytest.mark.parametrize(
    "bad",
    [
        {**entry(), "rcsb_id": ["2XYZ"]},  # an ID that cannot even be looked up
        {**entry(), "rcsb_id": None},
        {**entry(), "rcsb_id": 2},
        "text",
        {"unexpected": True},
    ],
)
def test_record_without_a_usable_id_never_costs_its_neighbours(tmp_path, monkeypatch, bad):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    monkeypatch.setattr(pdb_data, "fetch", answer(named("1ABC"), bad, named("3DEF")))
    result = pdb_data.load_details(cfg, ["1abc_B", "2xyz_B", "3def_B"], True, Cancellation())
    assert result["1abc_B"]["status"] == result["3def_B"]["status"] == "ready"
    # The record that could not be attributed may have been the one still outstanding.
    assert result["2xyz_B"]["retryable"] and "unexpected response" in result["2xyz_B"]["message"]
    assert sorted(p.name for p in (cfg.pdb_cache / "metadata").iterdir()) == [
        "1ABC.json", "3DEF.json",
    ]


def test_known_bad_record_is_retryable_and_a_genuine_omission_is_not(tmp_path, monkeypatch):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    broken = named("2XYZ")
    del broken["nonpolymer_entities"]
    monkeypatch.setattr(pdb_data, "fetch", answer(broken, None, named("1ABC")))
    result = pdb_data.load_details(cfg, ["1abc_B", "2xyz_B", "4ghi_A"], True, Cancellation())
    assert result["1abc_B"]["status"] == "ready"
    assert result["2xyz_B"]["retryable"] and "nonpolymer_entities" in result["2xyz_B"]["detail"]
    assert result["4ghi_A"]["message"] == "This entry was not returned by RCSB."
    assert "retryable" not in result["4ghi_A"]
    assert sorted(p.name for p in (cfg.pdb_cache / "metadata").iterdir()) == ["1ABC.json"]


def test_one_unwritable_cache_file_does_not_stop_the_others(tmp_path, monkeypatch):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    (cfg.pdb_cache / "metadata" / "1ABC.json").mkdir(parents=True)  # cannot be replaced by a file
    monkeypatch.setattr(pdb_data, "fetch", answer(named("1ABC"), named("3DEF")))
    result = pdb_data.load_details(cfg, ["1abc_B", "3def_B"], True, Cancellation())
    assert result["3def_B"]["status"] == "ready"
    assert result["1abc_B"]["retryable"]
    assert result["1abc_B"]["message"] == "PDB details for 1ABC could not be saved to Kmatch's cache."
    assert result["1abc_B"]["detail"]  # A technical detail accompanies the cache-save message.
    assert (cfg.pdb_cache / "metadata" / "3DEF.json").is_file()
    assert not list((cfg.pdb_cache / "metadata").glob(".pdb-*"))  # no temporary file left behind


@pytest.mark.parametrize("body", [b"[]", b'"text"', b'{"data": null}', b'{"data": {"entries": {}}}'])
def test_returned_hits_report_an_unexpected_answer_as_retryable(tmp_path, monkeypatch, body):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    monkeypatch.setattr(pdb_data, "fetch", lambda *a: body)
    result = pdb_data.load_details(cfg, ["1abc_A", "2xyz_A"], True, Cancellation())
    assert all(r["retryable"] and "unexpected response" in r["message"] for r in result.values())


def test_entry_missing_from_an_answered_request_is_not_retryable(tmp_path, monkeypatch):
    cfg = Config(tmp_path / "home", tmp_path / "cache")
    monkeypatch.setattr(pdb_data, "fetch", lambda *a: b'{"data": {"entries": []}}')
    result = pdb_data.load_details(cfg, ["1abc_A"], True, Cancellation())["1abc_A"]
    assert result["message"] == "This entry was not returned by RCSB."
    assert "retryable" not in result
    monkeypatch.setattr(pdb_data, "fetch", lambda *a: b"<html>maintenance</html>")
    result = pdb_data.load_details(cfg, ["1abc_A"], True, Cancellation())["1abc_A"]
    assert result["retryable"] and "unexpected response" in result["message"]


def test_coordinate_download_failure_is_explained(tmp_path, monkeypatch):
    cfg = Config(tmp_path, tmp_path / "cache")
    monkeypatch.setattr(pdb_data.tls, "urlopen", refuse(certificate_failure()))
    progress = []
    with pytest.raises(RuntimeError, match="could not verify RCSB's security certificate") as e:
        pdb_data.assembly_file(cfg, "1abc-assembly1_A", True, Cancellation(), progress.append)
    # Coordinates are requested by the export, so that is the action to repeat.
    assert "Export the 3D alignment again to retry" in str(e.value)
    assert "Reload details" not in str(e.value)
    assert "cached" not in progress[0].message  # nothing has been downloaded yet
    assert not list(cfg.pdb_cache.rglob("*.cif"))


def test_invalid_coordinate_response_is_not_cached(tmp_path, monkeypatch):
    cfg = Config(tmp_path, tmp_path / "cache")
    monkeypatch.setattr(pdb_data, "fetch", lambda *a: b"data_wrong\n_entry.id 9XYZ\n")
    with pytest.raises(RuntimeError, match="does not identify"):
        pdb_data.assembly_file(cfg, "1abc-assembly1_A", True, Cancellation(), lambda e: None)
    assert not list(cfg.pdb_cache.rglob("*.cif"))


def add_ligand(model, name, element, position, entity_type, number):
    chain = gemmi.Chain("B")
    residue = gemmi.Residue()
    residue.name = name
    residue.seqid = gemmi.SeqId(number, " ")
    residue.het_flag = "H"
    residue.entity_type = entity_type
    residue.subchain = name
    atom = gemmi.Atom()
    atom.name = element
    atom.element = gemmi.Element(element)
    atom.pos = list(position)
    atom.occ = 1
    residue.add_atom(atom)
    chain.add_residue(residue)
    model.add_chain(chain)


@pytest.fixture
def full_case(case, monkeypatch, tmp_path):
    cfg, run, query, trace, hit, _ = case
    residues = hit_structure.ca_residues(trace[0][0])
    # Independent source→database transform in addition to nonidentity Foldseek u/t.
    source_u = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    source_t = np.array([40, 25, -15])
    source_xyz = (hit_structure.xyz(residues) - source_t) @ source_u
    full = make_structure("ACDKEFG", source_xyz, full=True)
    full[0][0].name = "B"
    for residue in full[0][0]:
        residue.entity_type = gemmi.EntityType.Polymer
        residue.subchain = "B"
    # Close to the first target residue only; the alignment contains an insertion later.
    contact = source_xyz[0] + [0, 0, 2]
    add_ligand(full[0], "ZN", "Zn", contact, gemmi.EntityType.NonPolymer, 200)
    add_ligand(full[0], "HOH", "O", contact, gemmi.EntityType.Water, 201)
    add_ligand(full[0], "NAG", "C", contact, gemmi.EntityType.Branched, 202)
    add_ligand(full[0], "SO4", "S", [999, 999, 999], gemmi.EntityType.NonPolymer, 203)
    # Supply parsed source to isolate scientific logic; fixture never changes normal cached data.
    path = tmp_path / "assembly.cif"
    path.touch()
    monkeypatch.setattr(hit_structure, "assembly_file", lambda *a: (path, {"sha256": "fixture"}))
    original_reader = gemmi.read_structure
    monkeypatch.setattr(
        gemmi, "read_structure", lambda p: full.clone() if Path(p) == path else original_reader(p)
    )
    return cfg, run, query, trace, hit, full


def test_full_atom_frame_and_ligand_contacts_follow_foldseek_alignment(full_case):
    cfg, _, query, trace, hit, full = full_case
    original = full.make_pdb_string()
    combined, _, _ = alignment.superpose(query, trace, hit)
    info = hit_structure.full_atom_superposition(
        cfg, combined, trace, hit, False, Cancellation(), lambda e: None
    )
    assert info["registration"]["ca_rmsd_angstrom"] < 1e-9
    assert [x["component_id"] for x in info["ligands"]] == ["ZN"]
    assert info["aligned_query_contacts"][0]["query_position"] == 1
    assert info["aligned_query_contacts"][0]["database_target_position"] == 1
    assert info["aligned_query_contacts"][0]["hit_residue"]["number"] == 10
    assert [c.name for c in combined[0]] == ["Q", "T", "L"]
    q = hit_structure.xyz(hit_structure.ca_residues(combined[0]["Q"]))
    t = hit_structure.xyz(hit_structure.ca_residues(combined[0]["T"]))
    pairs = alignment.aligned_pairs(hit, "ACDEFG", "ACDKEFG")
    assert np.allclose(q[pairs[:, 0]], t[pairs[:, 1]])
    assert full.make_pdb_string() == original


@pytest.mark.parametrize("fmt", ["pdb", "cif"])
def test_full_atom_coordinate_export_preserves_ligands_and_provenance(full_case, fmt):
    cfg, run, _, _, _, _ = full_case
    result = alignment.export_alignment(cfg, run, 0, fmt, full_atom=True)
    structure = gemmi.read_structure(str(result.artifact))
    assert [c.name for c in structure[0]] == ["Q", "T", "L"]
    assert structure[0]["L"][0].name == "ZN"
    info = json.loads((result.directory / "provenance.json").read_text())
    assert info["full_atom"]["contact_cutoff_angstrom"] == 4
    assert info["source_chain_names"]["T"] == "B"
    assert info["refitted"] is False


def test_wrong_conformation_rejected_instead_of_moving_ligands(full_case):
    cfg, _, query, trace, hit, full = full_case
    full[0][0][2]["CA"][0].pos.x += 9
    combined, _, _ = alignment.superpose(query, trace, hit)
    with pytest.raises(RuntimeError, match="does not reproduce"):
        hit_structure.full_atom_superposition(
            cfg, combined, trace, hit, False, Cancellation(), lambda e: None
        )


def test_mapping_extra_residues_missing_residues_and_alternate_conformer():
    assert hit_structure.trace_mapping("ACDE", "XACDE") == [1, 2, 3, 4]
    with pytest.raises(ValueError):
        hit_structure.trace_mapping("ACDE", "ACE")
    with pytest.raises(ValueError):
        hit_structure.trace_mapping("ACDE", "ACKE")
    r = gemmi.Residue()
    for label, occ in [("A", 0.2), ("B", 0.8)]:
        a = gemmi.Atom()
        a.name = "CA"
        a.element = gemmi.Element("C")
        a.altloc = label
        a.occ = occ
        r.add_atom(a)
    best = hit_structure.major_conformer(r)
    assert len(best) == 1 and best[0].occ == pytest.approx(0.8) and not best[0].has_altloc()


def test_metadata_job_is_independent_and_error_does_not_destroy_hits(case, monkeypatch):
    cfg, run, _, _, hit, _ = case
    job = Job(cfg)
    job.update(run_id=run.name, directory=str(run), prediction={"ready": True}, hits=[asdict(hit)])

    def fail(*args):
        raise OSError("fixture offline")

    monkeypatch.setattr("kmatch.gui.load_details", fail)
    job.load_pdb_details(run.name, False)
    job.details_thread.join(2)
    assert not job.snapshot()["busy"] and len(job.snapshot()["hits"]) == 1
    assert "fixture offline" in job.snapshot()["pdb_details_status"]
    # The certificate error appears in the summary and the hit is retryable.
    monkeypatch.undo()
    monkeypatch.setattr(pdb_data.tls, "urlopen", refuse(certificate_failure()))
    job.load_pdb_details(run.name, True)
    job.details_thread.join(2)
    status = job.snapshot()["pdb_details_status"]
    assert status.startswith("PDB details: 0/1 hits available. 1 not loaded. macOS could not")
    assert status.endswith("Choose Reload details to try again.")
    assert job.snapshot()["pdb_details"][hit.target]["retryable"]
    monkeypatch.undo()
    monkeypatch.setattr(pdb_data, "fetch", lambda *a: b'{"data": {"entries": []}}')
    job.load_pdb_details(run.name, True)
    job.details_thread.join(2)
    # An empty entries response yields only the available-hit count in the summary.
    assert job.snapshot()["pdb_details_status"] == "PDB details: 0/1 hits available."
    with pytest.raises(ValueError):
        job.load_pdb_details("stale", True)
    assert (
        create_app(job)
        .test_client()
        .post("/api/pdb-details", json={"run_id": run.name})
        .status_code
        == 403
    )


def test_cancellation_during_network_error_remains_cancellation(monkeypatch):
    cancel = Cancellation()

    def timeout(*args, **kwargs):
        cancel.cancel()
        raise TimeoutError("fixture request timed out")

    monkeypatch.setattr(pdb_data.tls, "urlopen", timeout)
    with pytest.raises(Cancelled):
        pdb_data.fetch("https://data.rcsb.org/graphql", cancel, 1024)


def test_query_contact_mapping_accounts_for_target_insertion(full_case):
    cfg, _, query, trace, hit, full = full_case
    position = list(full[0][0][4]["CA"][0].pos)
    position[2] += 1
    add_ligand(full[0], "LIG", "C", position, gemmi.EntityType.NonPolymer, 204)
    combined, _, _ = alignment.superpose(query, trace, hit)
    info = hit_structure.full_atom_superposition(
        cfg, combined, trace, hit, False, Cancellation(), lambda e: None
    )
    pair = next(p for p in info["aligned_query_contacts"] if p["database_target_position"] == 5)
    assert pair["query_position"] == 4
    assert pair["query_residue"]["number"] == 13
    assert pair["hit_residue"]["number"] == 14


def test_native_export_forwards_full_atom_and_network_options(case, tmp_path, monkeypatch):
    from kmatch.gui import DesktopApi

    cfg, run, _, _, hit, _ = case
    job = Job(cfg)
    job.update(run_id=run.name, directory=str(run), prediction={"ready": True}, hits=[asdict(hit)])
    api = DesktopApi(job)

    class Window:
        def create_file_dialog(self, *args, **kwargs):
            return str(tmp_path / "alignment.pse")

    api._window = Window()
    calls = []
    monkeypatch.setattr(job, "export_alignment", lambda *a, **kw: calls.append((a, kw)))
    assert api.export_alignment(run.name, 0, "pse", True, False)
    assert calls[0][1] == {"full_atom": True, "allow_network": False}


def test_backbone_alternate_follows_database_not_occupancy(full_case):
    cfg, _, query, trace, hit, full = full_case
    residue = full[0][0][2]
    original = residue["CA"][0]
    original.altloc = "A"
    original.occ = 0.2
    alternate = original.clone()
    alternate.altloc = "B"
    alternate.occ = 0.8
    alternate.pos.x += 2
    residue.add_atom(alternate)
    combined, _, _ = alignment.superpose(query, trace, hit)
    info = hit_structure.full_atom_superposition(
        cfg, combined, trace, hit, False, Cancellation(), lambda e: None
    )
    assert info["registration"]["ca_rmsd_angstrom"] < 1e-8
    assert len([a for a in combined[0]["T"][2] if a.name == "CA"]) == 1
