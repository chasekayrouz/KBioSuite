"""Reject the retired identity in active source and built/installed suite resources."""

import argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RETIRED = ("cluster" + "stack", "cluster" + " stack")


def scan(directory):
    matches = []
    paths = directory.rglob("*") if directory.is_dir() else [directory]
    for path in paths:
        relative = path.relative_to(directory)
        # Git history, dependency environments and preserved ignored review archives
        # are not application source. Build/install roots are scanned explicitly too.
        if any(
            p
            in {
                ".git",
                ".venv",
                "__pycache__",
                ".pytest_cache",
                ".ruff_cache",
                "review-bundle",
                "cache",
                "results",
            }
            for p in relative.parts
        ):
            continue
        if any(s in str(relative).lower() for s in RETIRED):
            matches.append(str(relative))
        if path.is_file() and not path.is_symlink():
            data = path.read_bytes().lower()
            if any(s.encode() in data for s in RETIRED):
                matches.append(str(relative))
    return sorted(set(matches))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", type=Path, nargs="*", default=[ROOT])
    args = parser.parse_args()
    failures = {}
    for directory in args.directories:
        found = scan(directory.resolve())
        if found:
            failures[str(directory)] = found
    if failures:
        for directory, files in failures.items():
            print(directory)
            for name in files:
                print(f"  {name}")
        raise SystemExit(1)
    print("Retired-name audit: zero active path/content references")


if __name__ == "__main__":
    main()
