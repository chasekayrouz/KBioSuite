"""Select only the reviewed docstring-trimmed source for this frozen module."""

import importlib.util
from pathlib import Path


def pre_find_module_path(api):
    path = Path(__file__).resolve().parents[1] / "frozen_docstrings.py"
    spec = importlib.util.spec_from_file_location("kbiosuite_frozen_docstrings", path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    helper.pre_find_module_path(api)
