import gzip
import io
import json
import tarfile
from pathlib import Path
from types import SimpleNamespace

import pytest

import kannotate.databases as database_manager
from kannotate.downloads import digest


def test_setup_reuses_installs_removes_archives_and_detects_changed_index(tmp_path, monkeypatch):
    sources = tmp_path / "sources"
    sources.mkdir()
    hmm = sources / "Pfam-A.hmm.gz"
    hmm.write_bytes(gzip.compress(b"HMMER3/f\nNAME Test\n//\n"))
    archive = sources / "Cdd.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for name in ("Cdd.pal", "Cdd.00.rps"):
            member = tarfile.TarInfo(name)
            member.size = 3
            tar.addfile(member, io.BytesIO(b"CDD"))
    descriptions = sources / "cddid.tbl.gz"
    descriptions.write_bytes(gzip.compress(b"1\tcd1\tTest\tTest enzyme\t3\n"))

    def asset(path):
        return {"filename": path.name, "url": path.as_uri(), "sha256": digest(path)}

    spec = {
        "pfam": {
            "version": "test",
            "expected_profiles": 1,
            "database_prefix": "Pfam-A.hmm",
            "assets": [asset(hmm)],
        },
        "cdd": {
            "version": "test",
            "expected_profiles": 1,
            "database_prefix": "Cdd",
            "assets": [asset(archive), asset(descriptions)],
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
    monkeypatch.setattr(database_manager, "verify_databases", lambda config: {"passed": True})
    runtime = tmp_path / "runtime"
    config = database_manager.setup(runtime, report=lambda _: None)
    assert config["validation"]["passed"]
    assert not list((runtime / "data/databases/.downloads").iterdir())
    metadata = json.loads((runtime / "data/databases/pfam-test/installation.json").read_text())
    assert metadata["specification"] == spec["pfam"]
    assert (
        database_manager.setup(runtime, report=lambda _: None)["databases"] == config["databases"]
    )
    Path(config["databases"]["pfam"]["path"] + ".h3f").write_text("changed")
    with pytest.raises(RuntimeError, match="missing or changed"):
        database_manager.setup(runtime, report=lambda _: None)
