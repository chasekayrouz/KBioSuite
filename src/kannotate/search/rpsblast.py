import json
import tempfile
from contextlib import contextmanager
from functools import lru_cache
from importlib.resources import files
from pathlib import Path
from types import MappingProxyType

from ..cancellation import check_cancelled
from ..hits import collapse_profiles
from ..models import Alignment, Hit
from .common import run_search
from .identity import file_identity

SEARCH_PARAMETERS = json.loads(files("kannotate.resources").joinpath("search.json").read_text())
RPS_ARGUMENTS = tuple(SEARCH_PARAMETERS["rpsblast"]["arguments"])
RPS_COLUMNS = "6 qseqid sseqid evalue bitscore qstart qend qlen slen length stitle"
ALIGNMENT_COLUMNS = " sstart send qseq sseq"


def load_descriptions(path: str):
    check_cancelled()
    path = Path(path).resolve()
    stat = path.stat()
    return _descriptions((str(path), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns))


@lru_cache(maxsize=4)
def _descriptions(identity):
    mapping = {}
    with open(identity[0]) as handle:
        for line in handle:
            check_cancelled()
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4:
                mapping[fields[0]] = (fields[1], fields[2], fields[3])
    return MappingProxyType(mapping)


def parse_rps(path: Path, descriptions=None, *, include_alignments=False) -> list[Hit]:
    hits = []
    descriptions = descriptions or {}
    with path.open() as handle:
        for line in handle:
            check_cancelled()
            if include_alignments:
                prefix, *alignment_fields = line.rstrip("\n").rsplit("\t", 4)
                fields = prefix.split("\t", 9) + alignment_fields
            else:
                fields = line.rstrip("\n").split("\t", 9)
            if len(fields) != (14 if include_alignments else 10):
                raise ValueError("Malformed RPS-BLAST table")
            query, subject, evalue, bits, start, end, qlen, slen, length, title = fields[:10]
            alignments = ()
            if include_alignments:
                model_start, model_end = map(int, fields[10:12])
                aligned_query, aligned_model = fields[12:14]
                if (
                    len(aligned_query) != len(aligned_model)
                    or len(aligned_query) != int(length)
                    or len(aligned_query.replace("-", "")) != int(end) - int(start) + 1
                    or len(aligned_model.replace("-", "")) != model_end - model_start + 1
                ):
                    raise ValueError("Malformed RPS-BLAST alignment coordinates")
                alignments = (
                    Alignment(
                        int(start),
                        int(end),
                        model_start,
                        model_end,
                        aligned_query,
                        aligned_model,
                        float(evalue),
                        float(bits),
                    ),
                )
            uid = subject.rsplit("|", 1)[-1]
            if uid in descriptions:
                accession, name, description = descriptions[uid]
            else:
                parts = title.split(", ", 2)
                accession = parts[0].split()[-1] if parts else subject
                name = parts[1] if len(parts) > 1 else accession
                description = parts[2] if len(parts) > 2 else title
            a, b = sorted((int(start), int(end)))
            hits.append(
                Hit(
                    "rps",
                    query,
                    accession,
                    name,
                    description,
                    float(evalue),
                    a,
                    b,
                    int(qlen),
                    int(slen),
                    float(bits),
                    alignments=alignments,
                )
            )
    return collapse_profiles(hits)


@contextmanager
def database_path(path):
    """Expose the same CDD files without spaces for RPS-BLAST's worker commands.

    Directory symlinks are canonicalized by BLAST, so link individual files in a
    temporary directory instead. No database bytes or search parameters change.
    """
    if not any(c.isspace() for c in path) and '"' not in path:
        yield path
        return
    prefix = Path(path).resolve()
    with tempfile.TemporaryDirectory(prefix="kannotate-cdd-", dir="/tmp") as temporary:
        directory = Path(temporary)
        for asset in prefix.parent.iterdir():
            if asset.is_file():
                (directory / asset.name).symlink_to(asset)
        yield str(directory / prefix.name)


class RpsBlastAdapter:
    def __init__(
        self,
        executable: str,
        database: str,
        descriptions: str | None = None,
        *,
        include_alignments=False,
    ):
        self.executable = executable
        self.database = database
        self.description_identity = file_identity(descriptions) if descriptions else None
        self.descriptions = load_descriptions(descriptions) if descriptions else {}
        self.include_alignments = include_alignments
        self.output_format = RPS_COLUMNS + (ALIGNMENT_COLUMNS if include_alignments else "")

    def search(self, queries: Path, workdir: Path, threads: int):
        table = workdir / "rps.tsv"
        with database_path(self.database) as database:
            run_search(
                [
                    self.executable,
                    "-query",
                    str(queries),
                    "-db",
                    database,
                    "-num_threads",
                    str(threads),
                    *RPS_ARGUMENTS,
                    "-outfmt",
                    self.output_format,
                    "-out",
                    str(table),
                ]
            )
        return parse_rps(table, self.descriptions, include_alignments=self.include_alignments)
