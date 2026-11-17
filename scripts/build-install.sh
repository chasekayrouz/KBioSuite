#!/bin/zsh
set -euo pipefail
cd "${0:A:h:h}"
export UV_PYTHON_PREFERENCE=only-managed
uv sync --locked
exec uv run --locked python scripts/build_apps.py --install "$@"
