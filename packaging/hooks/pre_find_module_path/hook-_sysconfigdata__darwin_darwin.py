"""Freeze relocatable installation metadata through PyInstaller's pre-find hook API.

uv relocates python-build-standalone's installation fields to the current user's managed
interpreter directory. Those fields must describe the frozen installation instead. Build/ABI
facts remain byte-for-byte values from the selected interpreter; no installed files are edited.
See PreFindModulePathAPI in PyInstaller/depend/imphookapi.py (search_dirs is the public API).
"""

import importlib
import os
import pprint
import re
import sys
from pathlib import Path

from PyInstaller.config import CONF

TOKEN = "__KBIOSUITE_FROZEN_PREFIX__"
# PyInstaller puts Python's shared library, stdlib archive and extension modules under
# sys._MEIPASS, rather than CPython's usual <prefix>/lib/pythonX.Y installation tree.
RUNTIME_ROOT_FIELDS = ("BINLIBDEST", "DESTLIB", "DESTSHARED", "LIBDEST", "LIBDIR", "MACHDESTLIB")


def module_source(values, prefix):
    """Rebase installation strings only; keep compiler flags, ABI tags and feature facts.

    prefix is the selected interpreter's actual base installation, never a synthetic default.
    Refuse another personal path rather than hiding unrecognized metadata with broad rewriting.
    """
    prefix = os.fspath(prefix)
    if not os.path.isabs(prefix) or prefix == os.path.sep:
        raise ValueError("Expected an absolute interpreter installation prefix")
    installation = re.compile(re.escape(prefix) + r"(?=/|$|[\s\"'])")
    normalized = {}
    for key, value in values.items():
        if isinstance(value, str):
            if TOKEN in value:
                raise ValueError("Interpreter metadata collides with frozen prefix token")
            value = installation.sub(TOKEN, value)
            # A new non-installation private path needs explicit semantic review.
            pieces = value.replace("\\", "/").split("/")
            if any(piece in {"Users", "home"} for piece in pieces[:-1]):
                raise ValueError("Unrecognized personal path in interpreter metadata")
        normalized[key] = value
    return (
        '"""Selected CPython build facts with frozen-runtime installation paths."""\n'
        "import os as _os\nimport sys as _sys\n"
        "_root = _os.fspath(getattr(_sys, '_MEIPASS', _sys.prefix))\n"
        f"build_time_vars = {pprint.pformat(normalized, sort_dicts=True)}\n"
        f"_token = {TOKEN!r}\n"
        "for _key, _value in tuple(build_time_vars.items()):\n"
        "    if isinstance(_value, str):\n"
        "        build_time_vars[_key] = _value.replace(_token, _root)\n"
        f"for _key in {RUNTIME_ROOT_FIELDS!r}:\n"
        "    if _key in build_time_vars:\n"
        "        build_time_vars[_key] = _root\n"
        "if 'BINDIR' in build_time_vars:\n"
        "    build_time_vars['BINDIR'] = _os.path.dirname(_sys.executable)\n"
        "if 'EXENAME' in build_time_vars:\n"
        "    build_time_vars['EXENAME'] = _sys.executable\n"
        "del _key, _value, _token, _root, _os, _sys\n"
    )


def pre_find_module_path(api):
    # PyInstaller's normal sysconfig hook supplies this platform-specific hidden import.
    # Import the real module in the build interpreter, before redirecting analysis to our copy.
    values = importlib.import_module(api.module_name).build_time_vars
    source = module_source(values, sys.base_prefix)
    directory = Path(CONF["workpath"]) / "kbiosuite-relocatable-sysconfig"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / (api.module_name + ".py")).write_text(source, encoding="utf-8")
    api.search_dirs = [str(directory)]
