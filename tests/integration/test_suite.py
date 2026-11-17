import importlib
import importlib.metadata
import subprocess
from importlib.resources import files

from kbiosuite.version import __version__


def test_independent_namespaces_share_one_version():
    for name in ("kbiosuite", "kblast", "kannotate", "kstack", "kmatch"):
        assert importlib.import_module(name).__version__ == __version__
    assert importlib.metadata.version("kbiosuite") == __version__
    assert importlib.util.find_spec("cluster" + "stack") is None
    assert importlib.util.find_spec("kbiosuite.pipeline") is None


def test_developer_entry_points_and_packaged_resources():
    scripts = importlib.metadata.distribution("kbiosuite").entry_points
    for name in ("kblast", "kannotate", "kstack", "kmatch", "kmatch-gui"):
        entry = next(p for p in scripts if p.group == "console_scripts" and p.name == name)
        assert callable(entry.load())
        result = subprocess.run([name, "--version"], capture_output=True, text=True, check=True)
        assert __version__ in result.stdout
    for name, resource in (
        ("kblast", "static/index.html"),
        ("kannotate", "resources/scoring.json"),
        ("kannotate", "static/index.html"),
        ("kstack", "resources/renderer/template.html"),
        ("kstack", "resources/shared/project.schema.json"),
        ("kmatch", "static/index.html"),
        ("kmatch", "static/bridge.js"),
        ("kmatch", "static/vendor/3Dmol-min.js"),
        ("kmatch", "static/vendor/3Dmol-LICENSE"),
        ("kmatch", "static/vendor/3Dmol-min.js.LICENSE.txt"),
        ("kmatch", "resources/boltz-requirements.txt"),
        # External engine and PyMOL interpreters execute these package files directly.
        ("kmatch", "engine_runner.py"),
        ("kmatch", "pymol_runner.py"),
    ):
        assert files(name).joinpath(resource).is_file()


def test_imports_clis_and_resources_from_unrelated_working_directory(tmp_path):
    import os
    import sys

    environment = {
        "HOME": str(tmp_path),
        "PATH": os.environ["PATH"],
        "PYTHONDONTWRITEBYTECODE": "1",
    }
    probe = """
from importlib.resources import files
import kbiosuite, kblast, kannotate, kstack, kmatch
from kblast.config import Config
assert not Config.discover().api_key
for package, resource in [('kblast','static/index.html'), ('kannotate','static/bridge.js'),
                          ('kstack','resources/renderer/template.html'),
                          ('kmatch','static/index.html')]:
    assert files(package).joinpath(resource).is_file()
print('outside-checkout imports/resources passed')
"""
    result = subprocess.run(
        [sys.executable, "-B", "-c", probe],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    for module in ["kblast.app", "kannotate.cli", "kstack.cli", "kmatch.cli"]:
        result = subprocess.run(
            [sys.executable, "-B", "-m", module, "--version"],
            cwd=tmp_path,
            env=environment,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0 and __version__ in result.stdout
