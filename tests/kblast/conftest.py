import io
import threading
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from kblast.models import Run, Settings
from kblast.network import RemoteError
from kblast.pipeline import Pipeline
from kblast.sequence import parse_query
from kblast.storage import RunStore

PROTEIN = "MAVKDLGATNRELVKATGQV"
HEADER = (
    "Id\tSource\tNucleotide Accession\tStart\tStop\tStrand\tProtein\t"
    "Protein Name\tOrganism\tStrain\tAssembly\n"
)


def hsp_xml(bit_score=120, start=1, end=20, identities=18, evalue="1e-30"):
    return f"""<Hsp><num>1</num><bit-score>{bit_score}</bit-score><score>250</score>
    <evalue>{evalue}</evalue><query-from>{start}</query-from><query-to>{end}</query-to>
    <hit-from>1</hit-from><hit-to>{end - start + 1}</hit-to><identity>{identities}</identity>
    <positive>20</positive><gaps>0</gaps><align-len>{end - start + 1}</align-len></Hsp>"""


def blast_xml(database="nr"):
    hits = []
    for i, accession in enumerate(["PROTA.1", "PROTB.1", "PROTC.1"], 1):
        hsps = hsp_xml(100 - i)
        if i == 1:
            hsps += hsp_xml(200, 5, 14, 8)
        hits.append(f"""<Hit><num>{i}</num><description><HitDescr>
        <id>ref|{accession}|</id><accession>{accession}</accession>
        <title>Test synthase [Example species]</title><taxid>1</taxid>
        <sciname>Example species</sciname>
        </HitDescr><HitDescr><id>gb|ALIAS{i}.1|</id><accession>ALIAS{i}.1</accession>
        <title>Alternative defline</title></HitDescr></description><len>20</len>
        <hsps>{hsps}</hsps></Hit>""")
    return f"""<?xml version="1.0"?><BlastXML2 xmlns="http://www.ncbi.nlm.nih.gov">
    <BlastOutput2><report><Report><program>blastp</program><version>BLASTP test</version>
    <search-target><Target><db>{database}</db></Target></search-target>
    <params><Parameters><matrix>BLOSUM62</matrix><expect>0.05</expect></Parameters></params>
    <results><Results><search><Search><query-id>Query_1</query-id><query-len>20</query-len>
    <hits>{"".join(hits)}</hits><stat><Statistics><db-num>123</db-num><db-len>456</db-len>
    </Statistics></stat></Search></search></Results></results></Report></report></BlastOutput2>
    </BlastXML2>"""


def ipg_text():
    return (
        HEADER
        + "\n".join(
            [
                "1\tRefSeq\tNC_000001.1\t10\t69\t+\tPROTA.1\tsynthase\t"
                "Example species\tstrain A\tGCF_1.1",
                "1\tRefSeq\tNC_000001.1\t10\t69\t+\tWP_000001.1\tsynthase\t"
                "Example species\tstrain A\tGCF_1.1",
                "1\tINSDC\tNC_000002.1\t959\t900\t-\tPROTA.1\tsynthase\t"
                "Example species\tstrain B\tGCA_2.1",
                "1\tSwissProt\t\t\t\t\tPROTA.1\tsynthase\tExample species\t\t",
                "1\tRefSeq\tNC_000003.1\t250\t309\t+\tPROTA.1\tsynthase\t"
                "Example species\tstrain C\t",
                "1\tRefSeq\tNC_000004.1\t980\t1040\t+\tPROTA.1\tsynthase\t"
                "Example species\tstrain D\t",
            ]
        )
        + "\n"
    )


def genomic_region(accession, start, end):
    sequence = ("ACGTTAGC" * 125)[start - 1 : end]
    record = SeqRecord(
        Seq(sequence),
        id=accession,
        name=accession.split(".")[0],
        description="Synthetic annotated test fixture",
    )
    record.annotations = {
        "molecule_type": "DNA",
        "organism": "Example species",
        "topology": "linear",
    }
    record.features.append(SeqFeature(SimpleLocation(0, len(sequence)), type="source"))
    target_low, target_high, strand = {
        "NC_000001.1": (10, 69, 1),
        "NC_000003.1": (250, 309, 1),
    }.get(accession, (900, 959, -1))
    record.features.append(
        SeqFeature(
            SimpleLocation(target_low - start, target_high - start + 1, strand=strand),
            type="CDS",
            qualifiers={
                "gene": ["synthase"],
                "locus_tag": ["TEST_001"],
                "protein_id": ["PROTA.1"],
                "product": ["test synthase"],
                "translation": [PROTEIN],
            },
        )
    )
    if start == 1:
        record.features.append(
            SeqFeature(
                SimpleLocation(99, 140, strand=-1),
                type="CDS",
                qualifiers={"locus_tag": ["TEST_002"], "product": ["neighbor"]},
            )
        )
    output = io.StringIO()
    SeqIO.write(record, output, "genbank")
    return output.getvalue(), record


class FakeHttp:
    def __init__(self):
        self.cancel = threading.Event()

    def safe_message(self, text):
        return text

    def close(self):
        pass


class FixtureClient:
    def __init__(self):
        self.http = FakeHttp()
        self.calls = []
        self.fail_context = True
        self.cancel_after_context = False

    def submit_blast(self, run):
        self.calls.append(("submit",))
        return "FIXTURE-RID", 1

    def prefetch_hits(self, hits):
        pass

    def linked_occurrences(self, hit, flank, members, **kwargs):
        return []

    def prefetch_nucleotide_info(self, accessions):
        pass

    def context_results(self, requests, whole_records, unindexed=()):
        for key in requests:
            try:
                yield key, self.fetch_context(*key)
            except RemoteError as exc:
                yield key, exc

    def wait_for_blast(self, run, checkpoint):
        self.calls.append(("blast",))
        return blast_xml()

    def fetch_protein(self, accession, rank):
        self.calls.append(("protein", accession))
        if rank == 2:
            raise RemoteError("Synthetic protein fetch failure")
        return PROTEIN, accession

    def fetch_ipg(self, accession, rank):
        self.calls.append(("ipg", accession))
        return ipg_text() if rank == 1 else HEADER

    def nucleotide_info(self, accession, **kwargs):
        self.calls.append(("info", accession))
        return {"accession": accession, "length": 1000}

    def fetch_context(self, accession, start, end):
        self.calls.append(("context", accession, start, end))
        if accession == "NC_000003.1" and self.fail_context:
            raise RemoteError("Synthetic context fetch failure")
        if self.cancel_after_context:
            self.http.cancel.set()
        text, record = genomic_region(accession, start, end)
        return text, record, "2026-01-01T00:00:00+00:00"


@pytest.fixture
def run_fixture(tmp_path: Path):
    store = RunStore(tmp_path / "runs")
    run = Run("a" * 32, parse_query(PROTEIN), Settings(2, 100))
    store.save(run)
    client = FixtureClient()
    pipeline = Pipeline(store, run, client)
    return store, run, client, pipeline
