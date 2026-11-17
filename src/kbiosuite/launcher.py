"""Developer entry point for the installed lightweight native launcher."""

import argparse
import subprocess
import sys
from pathlib import Path

from .version import __version__


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "resources":
        from .resources import main as resources_main
        return resources_main(sys.argv[2:])
    parser = argparse.ArgumentParser(description="Open the KBioSuite launcher")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()
    del args
    app = Path.home() / "Applications" / "KBioSuite" / "KBioSuite.app"
    if not app.is_dir():
        parser.exit(1, "Install the apps with ./scripts/build-install.sh first.\n")
    return subprocess.call(["/usr/bin/open", str(app)])
