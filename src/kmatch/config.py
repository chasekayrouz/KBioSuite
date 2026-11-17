"""One configuration boundary for the Kmatch desktop app, CLI and library callers."""

import json
import os
import stat
from dataclasses import dataclass
from pathlib import Path

SUITE_HOME = "Library/Application Support/KBioSuite/Kmatch"
# Consumer whose shared Foldseek/PDB100 pins a home uses when config.json chooses nothing.
DEFAULT_SHARED_RESOURCES = "kmatch"
# What a missing root/record or a malformed catalog/manifest raises from the read-only client.
SHARED_FAILURES = (OSError, KeyError, TypeError, ValueError, AttributeError)


# kbiosuite.resource_client.resource_root()'s default, before symlinks are resolved.
RESOURCES_HOME = "Library/Application Support/KBioSuite/Resources"


def shared_record_registered(name):
    """Whether consumer NAME may be registered in the shared library (existence only).

    Only a library root that does not exist at all, or a root without that consumer record,
    counts as unregistered. A root that cannot be inspected (a dangling symlink such as an
    unmounted volume, a file, or an unreadable folder) counts as registered, so binding
    reports it instead of setup filling the gap with a private download.
    """
    from kbiosuite.resource_client import resource_root

    configured = Path(
        os.environ.get("KBIOSUITE_RESOURCES", Path.home() / RESOURCES_HOME)
    ).expanduser()
    try:
        os.lstat(configured)
    except FileNotFoundError:
        return False
    except OSError:
        return True
    try:
        if not stat.S_ISDIR(os.stat(configured).st_mode):
            return True
    except OSError:
        return True
    try:
        os.lstat(resource_root() / "consumers" / f"{name}.json")
    except FileNotFoundError:
        return False
    except OSError:
        return True
    return True


def shared_paths(name, *, default=False):
    """Foldseek executable and PDB100 prefix pinned for consumer NAME (read-only).

    Any unresolvable pin is a startup error, never a fallback that setup could fill with a
    private download.
    """
    from kbiosuite.resource_client import bindings

    try:
        paths = bindings(name)
        return {"foldseek_override": paths["foldseek"], "pdb_prefix_override": paths["pdb100"]}
    except SHARED_FAILURES as exc:
        if not default:
            raise RuntimeError(f"Kmatch shared resources ({name!r}) need attention: {exc}") from exc
        raise RuntimeError(
            f"Kmatch default shared resources ({name!r}) need attention: {exc}. config.json "
            f"sets no shared_resources, so Kmatch uses the registered {name!r} pins by default. "
            'Repair them, or set "shared_resources": null in config.json to use '
            "Kmatch-managed Foldseek and PDB100 instead."
        ) from exc


def default_home():
    """Persistent Kmatch data: configuration, engines, models, runs and evidence."""
    fallback = Path.home() / SUITE_HOME
    return Path(os.environ.get("KMATCH_HOME", fallback)).expanduser().resolve()


def default_cache():
    """Purgeable Kmatch scratch, MSA and public-PDB caches."""
    fallback = Path.home() / "Library/Caches/KBioSuite/Kmatch"
    return Path(os.environ.get("KMATCH_CACHE", fallback)).expanduser().resolve()


@dataclass(frozen=True)
class Config:
    home: Path
    cache: Path
    predictor: str = "boltz2"
    threads: int = 8
    diffusion_samples: int = 3
    seed: int = 42
    boltz_python_override: str | None = None
    boltz_models_override: str | None = None
    af3_python_override: str | None = None
    af3_source_override: str | None = None
    af3_models_override: str | None = None
    shared_resources: str | None = None
    foldseek_override: str | None = None
    pdb_prefix_override: str | None = None
    pymol_override: str | None = None
    # Caller-selected run root (e.g. an isolated acceptance directory). Never a config.json key.
    runs_root: Path | None = None

    @classmethod
    def load(cls):
        home, cache = default_home(), default_cache()
        path = home / "config.json"
        try:
            settings = json.loads(path.read_text()) if path.exists() else {}
        except json.JSONDecodeError as exc:
            raise ValueError(f"{path} is not valid JSON: {exc}") from exc
        if not isinstance(settings, dict):
            raise ValueError("Kmatch config.json must contain a JSON object.")
        allowed = set(cls.__dataclass_fields__) - {"home", "cache", "runs_root"}
        unknown = set(settings) - allowed
        if unknown:
            raise ValueError("Unknown Kmatch configuration: " + ", ".join(sorted(unknown)))
        shared = settings.get("shared_resources")
        if shared is not None and shared is not False:
            if not isinstance(shared, str) or not shared.strip():
                raise ValueError(
                    "shared_resources must name a Resources consumer, or be null or false to "
                    "use Kmatch-managed Foldseek and PDB100."
                )
        # Foldseek/PDB100 precedence: (1) an explicit shared_resources consumer binds its pins,
        # while null/false opts out; (2) with the key absent and neither path override set
        # (unset or null), a registered default consumer binds its pins; (3) otherwise the
        # overrides, then Kmatch-managed home paths. Only a missing library or default record
        # reaches (3) silently; an unresolvable pin in (1) or (2) is a startup error.
        if settings.get("shared_resources"):
            settings.update(shared_paths(settings["shared_resources"]))
        elif (
            "shared_resources" not in settings
            and all(
                settings.get(key) is None for key in ("foldseek_override", "pdb_prefix_override")
            )
            and shared_record_registered(DEFAULT_SHARED_RESOURCES)
        ):
            settings.update(
                shared_paths(DEFAULT_SHARED_RESOURCES, default=True),
                shared_resources=DEFAULT_SHARED_RESOURCES,
            )
        config = cls(home, cache, **settings)
        from .predictors import BACKENDS

        if not isinstance(config.predictor, str) or config.predictor not in BACKENDS:
            raise ValueError("Unknown prediction backend in config.json.")
        if (
            type(config.threads) is not int
            or type(config.diffusion_samples) is not int
            or not 1 <= config.threads <= 32
            or not 1 <= config.diffusion_samples <= 5
        ):
            raise ValueError("Use 1–32 threads and 1–5 diffusion samples.")
        if type(config.seed) is not int or not 0 <= config.seed < 2**32:
            raise ValueError("Use an integer seed from 0 through 4294967295.")
        for name in allowed:
            if name.endswith("_override"):
                value = getattr(config, name)
                if value is not None and (not isinstance(value, str) or not value.strip()):
                    raise ValueError(f"{name} must be a nonempty path string or null.")
        return config

    def path(self, override, fallback):
        return Path(override).expanduser().resolve() if override else fallback

    @property
    def boltz_python(self):
        return self.path(
            self.boltz_python_override, self.home / "backends/boltz-community/.venv/bin/python"
        )

    @property
    def boltz_models(self):
        return self.path(self.boltz_models_override, self.home / "models/boltz2")

    @property
    def af3_source(self):
        return self.path(self.af3_source_override, self.home / "backends/alphafold3")

    @property
    def af3_python(self):
        return self.path(self.af3_python_override, self.af3_source / ".venv/bin/python")

    @property
    def af3_models(self):
        return self.path(self.af3_models_override, self.home / "models/af3")

    @property
    def foldseek(self):
        return self.path(self.foldseek_override, self.home / "tools/foldseek/bin/foldseek")

    @property
    def pdb_prefix(self):
        return self.path(self.pdb_prefix_override, self.home / "databases/pdb100/current/pdb")

    @property
    def runs(self):
        return Path(self.runs_root) if self.runs_root else self.home / "runs"

    @property
    def pdb_cache(self):
        return self.cache / "pdb"
