import json
from importlib.resources import files
from pathlib import Path

from ..cancellation import check_cancelled
from ..hits import collapse_profiles
from ..models import Hit
from .common import run_search

HMM_ARGUMENTS = tuple(
    json.loads(files("kannotate.resources").joinpath("search.json").read_text())["hmmscan"][
        "arguments"
    ]
)


def parse_domtbl(path: Path) -> list[Hit]:
    hits = []
    with path.open() as handle:
        for line in handle:
            check_cancelled()
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split(maxsplit=22)
            if len(fields) < 22:
                raise ValueError("Malformed HMMER domain table")
            hits.append(
                Hit(
                    "hmm",
                    fields[3],
                    fields[1] if fields[1] != "-" else fields[0],
                    fields[0],
                    fields[22].strip() if len(fields) > 22 else fields[0],
                    float(fields[6]),
                    int(fields[17]),
                    int(fields[18]),
                    int(fields[5]),
                    int(fields[2]),
                    float(fields[7]),
                )
            )
    return collapse_profiles(hits)


class HmmerAdapter:
    def __init__(self, executable: str, database: str):
        self.executable = executable
        self.database = database

    def search(self, queries: Path, workdir: Path, threads: int):
        from .scheduling import hmmer_search

        return hmmer_search(self._search, queries, workdir, threads)

    def _search(self, queries: Path, workdir: Path, threads: int):
        table = workdir / "hmmer.domtbl"
        run_search(
            [
                self.executable,
                "--cpu",
                str(max(0, threads - 1)),
                "--domtblout",
                str(table),
                "-o",
                str(workdir / "hmmer.log"),
                *HMM_ARGUMENTS,
                self.database,
                str(queries),
            ]
        )
        return parse_domtbl(table)
