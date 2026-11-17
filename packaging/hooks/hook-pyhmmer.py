"""Keep PyHMMER's runtime and metadata without its upstream test corpus."""

from PyInstaller.utils.hooks import collect_all

datas, binaries, hiddenimports = collect_all(
    "pyhmmer",
    exclude_datas=["tests/**"],
    filter_submodules=lambda name: name != "pyhmmer.tests"
    and not name.startswith("pyhmmer.tests."),
)
