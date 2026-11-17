# Public NCBI protocol fixtures

`biopython_blastp.xml` is the upstream `xml2_21500_blastp_001.xml` test fixture from the Biopython project, retrieved on 2026-10-03. It contains an NCBI BLASTP 2.15.0+ search against nr (239-aa query, 10 subject hits). It is used only for format compatibility; it is not a Kblast live search result.

Source: https://github.com/biopython/biopython/blob/master/Tests/Blast/xml2_21500_blastp_001.xml

Biopython’s license is included as `BIOPYTHON_LICENSE.txt`. Other tests generate clearly labeled synthetic records in temporary directories.

## WP_320672066.1 regression

The four `wp_320672066.*` fixtures were captured through NCBI E-utilities on 2026-10-04. The IPG TSV and XML both contain only a RefSeq `N/A` row (IPG 722037354, zero nucleotide entries). `linked.gb.gz` is the complete public nucleotide record `NZ_JAWZSQ010000020.1` (107,044 bp), with the exact protein CDS at `complement(24826..25701)`. `context.gb.gz` is its unmodified native slice 4826..45701. Gzip is only storage compression; decompression reproduces the captured text.

Requests: EFetch `db=protein`, `id=WP_320672066.1`, `rettype=ipg`, `retmode=text/xml`; ELink `dbfrom=protein`, `db=nuccore` resolves nucleotide UID 2626944563; EFetch `db=nuccore`, `rettype=gbwithparts`, `strand=1`, with the full record or the context range. These deliberate public response fixtures contain no local/private query or credential. History envelope/count fixtures in tests are synthetic. Normal tests make no network calls.
