#!/bin/zsh
set -euo pipefail
cd "${0:A:h:h}"
# Checks name directories, never individual files, so the same script passes in any
# checkout of this repository, whichever optional files it contains.
uv run --locked --offline python scripts/check_versions.py
uv run --locked --offline python scripts/check_naming.py src apps tests docs scripts packaging resources licenses examples README.md
uv run --locked --offline pytest -q
# Kstack's Python and its helper scripts keep their compact upstream style and are
# not linted, as before; validate_live.py and validate_mapping.py are not linted yet.
lint_excluded=(src/kstack tests/kstack scripts/kstack scripts/kblast/validate_live.py scripts/kblast/validate_mapping.py)
uv run --locked --offline ruff check src tests scripts packaging --extend-exclude "${(j:,:)lint_excluded}"
for file in src/kannotate/static/*.js; do node --check "$file"; done
# Kmatch vendor/ holds the upstream minified 3Dmol.js build.
for file in src/kmatch/static/*.js; do node --check "$file"; done
for file in src/kstack/resources/{renderer,export}/*.js; do node --check "$file"; done
node --test tests/integration/gui/*.test.cjs
node --test tests/kstack/renderer/*.test.cjs
KSTACK_PYTHON="$PWD/.venv/bin/python" node tests/kstack/renderer/e2e.cjs
swift test --package-path apps/kstack-macos
