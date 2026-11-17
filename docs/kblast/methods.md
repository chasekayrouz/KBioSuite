# Kblast methods

How Kblast searches, maps hits to genomes, retrieves contexts and reports results.
For everyday use, see the [user guide](../user-guide.md#find-genomic-neighborhoods-kblast).

Kblast works in three stages: **BLAST hits** (protein search), **occurrences** (every
genome location that encodes each hit) and **contexts** (the GenBank region around
each occurrence). Every stage saves its raw NCBI responses before parsing them.

## Query

The query is one protein: pasted amino acids or a one-record FASTA file. Whitespace
is removed, letters are uppercased and one terminal `*` is removed. Invalid characters
and multiple records are rejected. The normalized query is exported as `query.fasta`.

## BLAST search

New searches are submitted once through the
[NCBI BLAST URL API](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html) with:

| Parameter | Value |
| --- | --- |
| `PROGRAM` | `blastp` |
| `DATABASE` | `nr` |
| `HITLIST_SIZE` | the requested number of hits, N |
| `EXPECT` | `0.05` |
| `WORD_SIZE` | `5` |
| `GAPCOSTS` | `11 1` |
| `SHORT_QUERY_ADJUST` | `true` |

These values match the standard Protein BLAST web form for nr; they are deliberate
overrides of the URL API defaults (expect 10, word size 3, gap costs 11 2, no short-query
adjustment). Matrix, composition-based statistics and filtering keep the API defaults
(BLOSUM62, mode 2, no filter). The parameters NCBI actually used, including any
short-query adjustments, are read from the result and saved. Quick BLASTP and
`blastp-fast` are different, accelerated methods; Kblast does not use them.

The first N protein subjects are kept in NCBI's order. A subject can carry several
identical-sequence deflines: the first identifies it and the others are kept in
provenance. The hit table summarizes each subject across all its HSPs:

- **query coverage**: the union of aligned query positions divided by query length;
  overlapping HSPs count once;
- **percent identity**: the highest HSP identity divided by its alignment length;
- **bit score**: the highest HSP bit score;
- **E-value**: the lowest HSP E-value.

The original HSPs are kept unchanged in provenance. NCBI's databases change, so a
later search of the same protein can return different hits. Processing is
deterministic for a given result.

## NCBI pacing and retries

Kblast follows NCBI's [BLAST usage guidance](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html)
and [E-utilities usage rules](https://www.ncbi.nlm.nih.gov/books/NBK25497/):

- **BLAST:** at least 10 seconds between BLAST requests, and each search's status is
  checked at most once per minute. NCBI's estimated time to completion is advisory.
  Polling stops after six hours with a resumable status.
- **E-utilities:** fewer than 3 requests per second without an API key and fewer than
  10 with one, shared by all of Kblast's retrieval workers (at most three requests in
  flight).
- **Retries:** connection failures, interrupted responses and HTTP 408, 429 and common
  5xx responses are retried up to three times with backoff, honoring `Retry-After`.
  A very long server-requested delay stops that item for a later resume.
- **Submissions are never repeated automatically**, because a timeout can hide a
  successful submission.

Every NCBI request carries `tool=Kblast` and, if configured, your contact email. An
API key is sent only to E-utilities, never to BLAST; it raises the E-utilities rate
limit and gives no BLAST priority. Keys never appear in runs, exports or logs.

## BLAST statuses

| NCBI status | What Kblast does |
| --- | --- |
| `WAITING` | Keeps polling once per minute. A waiting search is pending at NCBI, not a client error. |
| `READY` | Downloads the XML2 result and parses it. |
| `FAILED`, `UNKNOWN` | Stops with an explicit message; `UNKNOWN` usually means the RID has expired. |
| Unexpected page | Saved and retried at the polling interval; three in a row stop the run. |

The returned program and database must be blastp and the run's database; otherwise
Kblast stops before mapping and keeps the report. A completed search with no hits
ends with status `no_hits` and makes no further requests.

## Existing and browser RIDs

A run can use the RID of a search started elsewhere (for example on NCBI's website)
instead of submitting its own. An unfinished run can also switch from its own waiting
RID to another one; Kblast stops its worker first, and the original RID can be
reactivated until a result has been parsed. The original NCBI search keeps running at
NCBI either way.

An imported result must contain exactly one report and query, program blastp,
database nr, the same query length and valid hits. The aligned query residues NCBI
reports are compared with your sequence at their positions. Masked, ambiguous or
unaligned positions can't be checked, so a different sequence of the same length can
go undetected when NCBI reports little alignment (for example, with zero hits); the
checked positions are recorded. Query titles alone are never taken as proof of
identity.

Imported results may use other valid search parameters. Kblast labels them as
externally supplied and shows the program, database, version, parameters and hit
counts that NCBI reports; it does not invent values that the XML2 report doesn't
contain. Fewer returned hits than requested are used as available; extra hits are
truncated in NCBI's order.

## From hits to occurrences

1. **Identical Protein Groups (IPG).** For each hit, Kblast fetches the IPG report
   (`efetch db=protein rettype=ipg`), batched up to 100 hits and matched by versioned
   accession. Every row with a nucleotide accession and CDS coordinates becomes an
   occurrence, including repeated rows. Protein-only rows stay in provenance. If a
   batch fails, reports are fetched one at a time.
2. **Linked nucleotide records.** When a valid IPG report has no genomic rows, or a
   group lacks RefSeq/INSDC nucleotide entries, Kblast follows protein → nuccore
   links (ELink with history), pages through every linked record in groups of 10 and
   keeps every CDS whose versioned `protein_id` matches the hit. These occurrences are
   marked `ELink + GenBank CDS` with `ipg_row=0`. A linked CDS that duplicates an IPG
   occurrence (same protein, accession, coordinates and strand) adds evidence rather
   than a second copy.

There is no deduplication, clustering, representative selection or taxonomic
filtering, and no cap on occurrences. Identical locations under different proteins,
records or hits stay separate.

### Hit statuses

| Status | Meaning |
| --- | --- |
| `no_mapping` | Both a valid IPG report and a valid, empty link response found no genome location. A normal outcome. |
| `retrieval_failure` | A mapping request failed. It can be retried on resume. |
| `parse_failure` | NCBI returned an invalid or malformed response (for example, an IPG row without usable coordinates). |
| `unresolved_mapping` | Linked records exist but none establishes the target CDS, or the CDS location is unsupported. |
| `integrity_failure` | A downloaded context doesn't contain the target CDS at the expected coordinates and strand. |

Empty, malformed, truncated or error responses never count as biological absence.
A run whose items did not all succeed reports `partial_success`; read the statuses
and the list of missing protein ranks before treating results as complete.

## Contexts and coordinates

All coordinates are **1-based and inclusive**, in the native nucleotide record.

- Requested start = CDS low − flank; requested end = CDS high + flank.
- Regions are clipped at the record's ends and never wrap around the origin of a
  circular record. Both clips are recorded.
- Record lengths come from ESummary (batched up to 200 accessions) when available;
  otherwise Kblast requests the region directly and validates its actual range from
  the returned record, rather than downloading a whole chromosome to learn its
  length.
- Contexts are fetched with `rettype=gbwithparts`, `seq_start`, `seq_stop` and
  `strand=1`, so the record keeps its native orientation; a target on the minus strand
  stays a minus-strand gene.

Inside an exported GenBank file, feature coordinates are relative to the extracted
region. The table and the arrow view use native record coordinates. A context counts
as successful only if it contains a CDS matching the target's coordinates and strand;
otherwise it gets `integrity_failure` and its downloaded text is kept locally for
inspection but not exported.

Full protein sequences are fetched in batches of up to 100; complete nucleotide
records in batches of up to 20 records and 2 Mbp. Failed batch members are retried
individually. A failed protein sequence fetch doesn't affect mapping; missing ranks
are listed in the export's provenance.

## Export folder

An export is a new folder (a ZIP in browser mode) containing:

| File | Contents |
| --- | --- |
| `contexts/*.gb` | One GenBank file per selected successful occurrence, exactly as returned by NCBI. Names include the occurrence ID, nucleotide version and native range. |
| `contexts.gb` | The same records in BLAST rank and occurrence order; repeated record IDs are intentional. |
| `blast_hits.fasta` | Full sequences of the hit proteins in rank order (never HSP fragments). |
| `metadata.tsv` | One row per hit and per occurrence: RID and its origin, mapping method, statuses and errors, requested and actual coordinates, strand, clipping, assembly and export selection. |
| `query.fasta` | The normalized query. |
| `provenance/run.json` | Settings, all HSPs and deflines, mappings, statuses, timestamps, checksums and the selected contexts. |
| `provenance/` | Raw BLAST and IPG XML, per-hit IPG tables, link and linked-record evidence, length lookups and a credential-free request log. |

The selection changes only the context files; `blast_hits.fasta` and `metadata.tsv`
keep every hit and occurrence. An export made during a run reflects the saved state at
that moment.

## Run records

Each run is saved under `~/Library/Application Support/KBioSuite/Kblast/runs/`.
`run.json` holds the settings, results, statuses and RID history (the original RID,
the active RID and how it was obtained). `blast-lifecycle.jsonl` records submission,
status checks, transitions, result retrieval, parsing, cancellation and resumption.
`requests.jsonl` records each request without credentials, and raw status pages are
kept under `raw/`. Only one Kblast process can use a run store at a time.

Cancelling waits for requests already in flight and keeps everything saved;
resuming reuses saved responses and continues polling the same RID. NCBI deletes
results after a while, so an old RID may no longer be retrievable.

## Completeness

Kblast keeps every occurrence returned by the mapping routes it uses. NCBI's indexes
can lag or omit links, and
[Entrez link limits](https://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html)
apply, so Kblast cannot prove that no other locus exists. The link route reads whole
annotated records and can be slow for large or numerous linked chromosomes.
