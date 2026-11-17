import json
import os
from pathlib import Path

import psutil


def default_threads() -> int:
    physical = psutil.cpu_count(logical=False) or os.cpu_count() or 1
    return max(1, physical - 2)


def default_root() -> Path:
    if os.environ.get("KANNOTATE_HOME"):
        return Path(os.environ["KANNOTATE_HOME"]).expanduser().resolve()
    if __import__("sys").platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "KBioSuite" / "Kannotate"
    return Path.home() / ".local" / "share" / "kbiosuite" / "kannotate"


def load_config(root: Path | None = None) -> dict:
    root = root or default_root()
    path = root / ".kannotate" / "config.json"
    if not path.is_file():
        raise RuntimeError("Databases are not configured. Run kannotate setup first.")
    config = json.loads(path.read_text())
    if config.get("shared_resources"):
        from kbiosuite.resource_client import bindings

        try:
            paths = bindings(config["shared_resources"])
            for name in ("pfam", "cdd"):
                config["databases"][name]["path"] = paths[name]
            config["databases"]["cdd"]["descriptions"] = paths["cdd_descriptions"]
            config["executables"] = {name: paths[name] for name in config["executables"]}
        except (OSError, KeyError, ValueError) as exc:
            raise RuntimeError(f"Shared resources need attention: {exc}") from exc
    return config


def save_config(config: dict, root: Path):
    directory = root / ".kannotate"
    directory.mkdir(parents=True, exist_ok=True)
    temporary = directory / "config.json.tmp"
    temporary.write_text(json.dumps(config, indent=2) + "\n")
    temporary.replace(directory / "config.json")
