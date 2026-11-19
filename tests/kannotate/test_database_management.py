import gzip
import io
import json
import tarfile
from pathlib import Path
from types import SimpleNamespace

import pytest

import kannotate.databases as database_manager
from kannotate.config import load_config
from kannotate.downloads import digest


def release_note(path, version, profiles):
    path.write_text(
        f"cdd version {version}\n{profiles} models in the default search set\n"
        "1009 models from SMART v6.0\n"
    )


@pytest.fixture
def runtime(tmp_path, monkeypatch):
    """Local sources standing in for EBI/NCBI, with search tools stubbed out."""
    sources = tmp_path / "sources"
    sources.mkdir()
    hmm = sources / "Pfam-A.hmm.gz"
    hmm.write_bytes(gzip.compress(b"HMMER3/f\nNAME Test\n//\n"))
    archive = sources / "Cdd_LE.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for name in ("Cdd.pal", "Cdd.00.rps"):
            member = tarfile.TarInfo(name)
            member.size = 3
            tar.addfile(member, io.BytesIO(b"CDD"))
    descriptions = sources / "cddid.tbl.gz"
    descriptions.write_bytes(gzip.compress(b"1\tcd1\tTest\tTest enzyme\t3\n"))
    note = sources / "cdd.info"
    release_note(note, "3.22", 1)
    spec = {
        "pfam": {
            "version": "test",
            "expected_profiles": 1,
            "database_prefix": "Pfam-A.hmm",
            "assets": [{"filename": hmm.name, "url": hmm.as_uri(), "sha256": digest(hmm)}],
        },
        "cdd": {
            "release_info": note.as_uri(),
            "database_prefix": "Cdd",
            "assets": [{"filename": p.name, "url": p.as_uri()} for p in (archive, descriptions)],
        },
    }
    monkeypatch.setattr(database_manager, "specification", lambda: spec)
    monkeypatch.setattr(
        database_manager,
        "resolve_dependencies",
        lambda _: {"hmmscan": "test", "hmmpress": "test", "rpsblast": "test", "blastdbcmd": "test"},
    )
    monkeypatch.setattr(database_manager, "executable_version", lambda _: "test")

    def press(command, **kwargs):
        for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
            Path(command[1] + suffix).write_text("test index")
        return SimpleNamespace(returncode=0, stderr="")

    monkeypatch.setattr(database_manager.subprocess, "run", press)
    verified = []
    monkeypatch.setattr(
        database_manager, "verify_databases", lambda c: verified.append(c) or {"passed": True}
    )
    return SimpleNamespace(
        root=tmp_path / "runtime", note=note, spec=spec, archive=archive, verified=verified
    )


def setup(runtime):
    return database_manager.setup(runtime.root, report=lambda _: None)


def test_setup_installs_current_cdd_release_and_reuses_it_offline(runtime):
    config = setup(runtime)
    storage = runtime.root / "data/databases"
    assert config["databases"]["cdd"]["version"] == "3.22"
    assert config["databases"]["cdd"]["path"] == str(storage / "cdd-3.22/Cdd")
    metadata = json.loads((storage / "cdd-3.22/installation.json").read_text())
    assert metadata["release"]["version"] == "3.22" and metadata["release"]["profiles"] == 1
    assert metadata["release"]["assets"]["Cdd_LE.tar.gz"] == digest(runtime.archive)
    assert set(metadata["files"]) == {"Cdd.pal", "Cdd.00.rps", "cddid.tbl"}
    assert sorted(p.name for p in storage.iterdir()) == [".downloads", "cdd-3.22", "pfam-test"]
    assert not list((storage / ".downloads").iterdir())
    # The fresh CDD was checked in staging before it was published.
    assert ".cdd-" in runtime.verified[0]["databases"]["cdd"]["path"]
    pfam = json.loads((storage / "pfam-test/installation.json").read_text())
    assert pfam["specification"] == runtime.spec["pfam"]

    # A configured installation is kept without contacting NCBI, even after a new release.
    runtime.note.unlink()
    assert setup(runtime)["databases"] == config["databases"]
    assert database_manager.status(runtime.root)["databases"]["cdd"] == {
        "version": "3.22",
        "ready": True,
    }
    Path(config["databases"]["pfam"]["path"] + ".h3f").write_text("changed")
    with pytest.raises(RuntimeError, match="missing or changed"):
        setup(runtime)


def test_release_change_during_download_publishes_nothing(runtime, monkeypatch):
    fetch = database_manager.fetch

    def rollover(url, path, **kwargs):
        release_note(runtime.note, "3.23", 2)
        return fetch(url, path, **kwargs)

    monkeypatch.setattr(database_manager, "fetch", rollover)
    with pytest.raises(RuntimeError, match="new CDD release during the download"):
        setup(runtime)
    assert not list((runtime.root / "data/databases").glob("*cdd*"))


def test_failed_cdd_check_leaves_nothing_to_reuse(runtime, monkeypatch):
    def reject(config):
        raise RuntimeError("CDD profile count differs")

    monkeypatch.setattr(database_manager, "verify_databases", reject)
    with pytest.raises(RuntimeError, match="profile count"):
        setup(runtime)
    assert not list((runtime.root / "data/databases").glob("*cdd*"))


def test_installations_pinned_by_earlier_versions_stay_ready(runtime):
    setup(runtime)
    storage = runtime.root / "data/databases"
    # Rewrite the receipt the way Kannotate 0.6.0 recorded its single pinned release.
    legacy = storage / "cdd-3.21"
    (storage / "cdd-3.22").rename(legacy)
    receipt = json.loads((legacy / "installation.json").read_text())
    del receipt["release"]
    receipt["specification"] = {
        "version": "3.21",
        "expected_profiles": 62456,
        "database_prefix": "Cdd",
        "assets": [{"filename": "Cdd_LE.tar.gz", "url": "u", "sha256": "0" * 64}],
    }
    (legacy / "installation.json").write_text(json.dumps(receipt))
    path = runtime.root / ".kannotate/config.json"
    config = json.loads(path.read_text())
    config["databases"]["cdd"].update(
        version="stale", path=str(legacy / "Cdd"), descriptions=str(legacy / "cddid.tbl")
    )
    path.write_text(json.dumps(config))

    assert load_config(runtime.root)["databases"]["cdd"]["version"] == "3.21"
    assert database_manager.status(runtime.root)["databases"]["cdd"] == {
        "version": "3.21",
        "ready": True,
    }
    runtime.note.unlink()  # proves setup keeps the old release without contacting NCBI
    assert setup(runtime)["databases"]["cdd"]["path"] == str(legacy / "Cdd")


def test_release_information_is_parsed_or_reported(tmp_path):
    note = tmp_path / "cdd.info"
    release_note(note, "3.22", 71136)
    assert database_manager.read_release_info(note.as_uri()) == {
        "version": "3.22",
        "profiles": 71136,
    }
    note.write_text("<html>maintenance</html>")
    with pytest.raises(RuntimeError, match="Could not read"):
        database_manager.read_release_info(note.as_uri())
    with pytest.raises(RuntimeError, match="Could not reach"):
        database_manager.read_release_info((tmp_path / "missing").as_uri())


@pytest.mark.parametrize(
    ("count", "title", "described", "error"),
    [
        ("200", "CDD.v3.22", 200, None),
        ("200", "CDD (renamed)", 199, None),  # NCBI's 3.21 table also lacked one model
        ("1,200", "CDD.v3.22", 200, "profile count"),
        ("200", "CDD.v3.21", 200, "is not release 3.22"),
        ("200", "CDD.v3.22", 197, "description table"),
    ],
)
def test_cdd_check_matches_the_recorded_release(tmp_path, monkeypatch, count, title, described, error):
    table = tmp_path / f"cddid-{described}.tbl"
    table.write_text("".join(f"{i}\tcd{i}\tname\tdescription\t1\n" for i in range(described)))
    outputs = {
        "-info": f"Database: {title}\n\t{count} sequences; 1 total residues\n",
        "-entry": "".join(f"gnl|CDD|{i}\n" for i in range(200)),
    }
    monkeypatch.setattr(
        database_manager.subprocess,
        "run",
        lambda command, **_: SimpleNamespace(returncode=0, stdout=outputs[command[3]], stderr=""),
    )
    release = {"version": "3.22", "profiles": 200}
    if error:
        with pytest.raises(RuntimeError, match=error):
            database_manager.verify_cdd("blastdbcmd", "/db/Cdd", str(table), release)
    else:
        assert database_manager.verify_cdd("blastdbcmd", "/db/Cdd", str(table), release)


def test_configured_cdd_without_receipt_is_reported_not_replaced(runtime):
    config = setup(runtime)
    (Path(config["databases"]["cdd"]["path"]).parent / "installation.json").unlink()
    runtime.note.unlink()  # a newer NCBI release must not be consulted or substituted
    with pytest.raises(RuntimeError, match="Incomplete CDD installation"):
        setup(runtime)


def test_malformed_receipt_is_a_handled_error(runtime):
    config = setup(runtime)
    receipt = Path(config["databases"]["cdd"]["path"]).parent / "installation.json"
    metadata = json.loads(receipt.read_text())
    del metadata["release"]
    receipt.write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="does not name its release"):
        load_config(runtime.root)
    assert database_manager.status(runtime.root)["databases"] == {}


def test_release_name_must_be_a_plain_directory_name(tmp_path):
    note = tmp_path / "cdd.info"
    release_note(note, "../../escape", 1)
    with pytest.raises(RuntimeError, match="Could not read"):
        database_manager.read_release_info(note.as_uri())

