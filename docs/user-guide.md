# KBioSuite user guide

This guide walks through each app step by step. To install KBioSuite and set up
Kannotate and Kmatch, follow the [README](../README.md#installation).

- [Before you start](#before-you-start)
- [Find genomic neighborhoods: Kblast](#find-genomic-neighborhoods-kblast)
- [Annotate proteins: Kannotate](#annotate-proteins-kannotate)
- [Compare gene clusters: Kstack](#compare-gene-clusters-kstack)
- [Predict a structure and find structural homologs: Kmatch](#predict-a-structure-and-find-structural-homologs-kmatch)
- [Command-line use](#command-line-use)
- [Troubleshooting](#troubleshooting)
- [Limitations](#limitations)

## Before you start

**Opening an app.** Press **Cmd+Space**, type the app name and press **Return**. Or
open **KBioSuite** and click a tool. Closing the launcher leaves the tools open.

**Your files.** Each app reads ordinary files and writes new ones. None of them
changes your original files, and nothing is passed from one app to another
automatically.

| App | Reads | Writes |
| --- | --- | --- |
| Kblast | A protein sequence or one-record FASTA | GenBank contexts, hit FASTA, tables and provenance |
| Kannotate | Protein FASTA and GenBank files | Annotated GenBank copies and CSV tables |
| Kstack | GenBank files, one record each | A portable HTML project, SVG figures and protein FASTA |
| Kmatch | A protein sequence or one-record FASTA, optionally an A3M alignment | Predicted structures, homolog tables and superposed 3D alignments |

**Where the apps keep their data.** Everything lives in your Library folder, not in
the apps, so updating the apps keeps it:

| Location | Contents |
| --- | --- |
| `~/Library/Application Support/KBioSuite/Kblast/runs/` | Saved Kblast runs |
| `~/Library/Application Support/KBioSuite/Kannotate/` | Kannotate databases, Quick Search cache and checkpoints |
| `~/Library/Application Support/KBioSuite/Kmatch/` | Kmatch engines, models, PDB100 database, settings and runs |
| `~/Library/Caches/KBioSuite/Kmatch/` | Kmatch sequence alignments and downloaded public PDB data |

**Shared resources.** The launcher's **Shared resources…** button shows an inventory
of a shared reference library, if you have set one up. Most people don't need one:
Kannotate and Kmatch install their own databases. The library's protocol is
described in `src/kbiosuite/resource_guide.md` in the source code.

## Find genomic neighborhoods: Kblast

Kblast searches one protein against NCBI's nr protein database with BLASTP, finds
every genome location that encodes each hit and downloads the surrounding region.

### Run a search

1. Paste the amino acids, or click **Import FASTA** and choose a file with one record.
   Spaces and line breaks are removed, letters are uppercased and one final `*` is
   dropped. Other characters and files with several records are rejected.
   **Use public example** fills in a public protein (CAC44203.1); nothing is sent
   until you run.
2. Choose **NR protein hits**, the number of BLAST hits to keep (default 10), and
   **Flank on each side**, the length of DNA to include on each side of each gene
   (default 20,000 bp).
3. Enter your **NCBI contact email** (recommended; see
   [below](#ncbi-contact-email-and-api-key)).
4. Click **Run search**. Kblast sends the protein to NCBI and waits for the result.

NCBI searches can take minutes or, when NCBI is busy, hours. The progress area shows
the current stage, elapsed time and NCBI's latest status. NCBI gives each search a
request ID (RID); click it to view the search on NCBI's website.

### While you wait

- **Cancel run** stops Kblast's work and saves the search and everything retrieved so
  far. It does not cancel the search at NCBI.
- **Saved runs** lists earlier runs. Open one to inspect it, or **Resume saved run**
  to continue an unfinished one with the same NCBI search, without searching again.
- NCBI keeps results for a limited time, so a very old unfinished run may no longer
  be resumable.

### Use an existing NCBI search

Scripted searches sometimes wait much longer in NCBI's queue than the same search
started on NCBI's website. You can run the search in your browser and give Kblast
its RID:

1. Click **Open NCBI Protein BLAST**. **Copy query** copies your sequence.
2. On the NCBI page, paste the sequence, choose the **nr** database and the standard
   **blastp** algorithm, and start the search.
3. Copy the RID shown by NCBI. In Kblast, click **Use existing NCBI RID…**, paste it
   and click **Retrieve existing result**.

Kblast checks that the result is a single protein search against nr with a query of
the same length, and compares the aligned query residues with your sequence. It
labels the result as externally supplied and shows the settings NCBI reports.

If a search has been waiting in NCBI's queue for five minutes, Kblast offers this
route. **Use / switch to existing RID…** switches the run to your browser search;
**Resume original RID** switches back, as long as no result has been read yet.

### Read the results

Kblast works through three stages, each with its own tab:

1. **BLAST hits**: the protein hits in NCBI's order, with coverage, identity, bit
   score and E-value.
2. **Occurrences**: every genome location that encodes each hit, found through NCBI's
   Identical Protein Groups and, when those list none, through linked nucleotide
   records. Repeated locations are all kept.
3. **Contexts**: the downloaded GenBank region around each occurrence, drawn as gene
   arrows. Each context must contain the expected gene at the expected coordinates.

Coordinates are 1-based and inclusive. A region stops at the end of its record.
Kblast does not remove duplicates, filter by taxonomy or predict gene clusters. Each
hit and context shows its status; a run can succeed partly, so check the statuses
before assuming the results are complete. [Kblast methods](kblast/methods.md)
explains the settings, statuses and coordinates in full.

### Export

**Export all** saves every successful context; **Export selected** saves only the
contexts you tick (all successful contexts are ticked at first). The export is a new
folder with:

- `contexts/`: one GenBank file per context, exactly as downloaded from NCBI;
- `contexts.gb`: the same records in one file;
- `blast_hits.fasta`: the full sequences of the hit proteins;
- `metadata.tsv`: one row per hit and per occurrence, with statuses and coordinates;
- `query.fasta` and `provenance/`: your query, settings and NCBI's raw responses.

The selection changes only the context files; the table and FASTA keep every hit. See
[Kblast methods](kblast/methods.md#export-folder) for every column and file.

### NCBI contact email and API key

The contact email tells NCBI who to contact if your searches cause a problem. It is
not a My NCBI sign-in, and searches don't appear in your NCBI search history. Type it
into Kblast, or save it once:

```sh
mkdir -p ~/.config/kblast
echo 'NCBI_EMAIL=you@example.org' > ~/.config/kblast/ncbi.env
chmod 600 ~/.config/kblast/ncbi.env
```

An NCBI API key lets Kblast fetch genomic records faster (up to 10 instead of 3
requests per second). It does not make the BLAST search itself faster. Create a key in your NCBI
account settings, save it as the only content of `~/.ncbi/api_key` (or as
`NCBI_API_KEY=...` in `~/.config/kblast/ncbi.env`), and protect it:

```sh
chmod 600 ~/.ncbi/api_key
```

Kblast never writes the key into runs, exports or logs; the app only shows whether a
key is configured. Restart Kblast after changing these files.

## Annotate proteins: Kannotate

Kannotate searches proteins against Pfam (with HMMER) and CDD (with RPS-BLAST) on
your Mac and suggests a product name for each protein. Set up the databases first
(see the [README](../README.md#kannotate)).

### Quick Search

1. Open **Quick Search** and paste 1 to 20 protein FASTA records, or one bare amino-acid
   sequence (it is named `query_1`). **Load protein FASTA or GenBank** reads a file
   instead; GenBank genes are converted to protein FASTA you can edit.
2. Click **Search** (or press **Cmd+Return**).
3. Choose a protein from the list (or press **Up/Down** in the list). The **Domains**
   diagram shows where each match aligns; overlapping matches stack on extra rows.
   Click a block or a table row to highlight both.
4. Use **All**, **CDD** or **Pfam** to filter, sort the table columns, and turn on
   **Show weaker hits** to see matches below Kannotate's score cutoff.
5. Select a CDD hit and open **Details / alignment** to see your protein aligned with
   the CDD model's consensus (`|` marks identical residues). Pfam alignments are not
   kept; their coordinates and statistics are.

**Edit FASTA** returns to the input; **Clear / reset** starts over. Quick Search
needs no output folder and saves nothing except a local cache that makes repeated
searches fast. **Search provenance** shows the tool and database versions used.

### Annotate files

1. Check that **Local databases** reads *ready*.
2. Choose one or more protein FASTA and/or GenBank files.
3. Choose the output folder. It follows the first input's folder until you pick
   another.
4. For GenBank inputs, choose **Reannotated GenBank only** (default), **CSV table
   only** or **Both GenBank and CSV**. Protein FASTA inputs always give a CSV table.
5. For large batches, tick **Keep / resume batch checkpoints**. After a stop or crash,
   run the same inputs, folder, threads and output choice again with it ticked to
   continue where you left off. Checkpoints stay until you delete them.
6. Click **Run annotation**. Each finished file appears with a **Show in Finder**
   button; a problem with one file doesn't stop the others.

Kannotate never overwrites existing outputs: choose another folder to run again.

**Stopping.** **Stop…** asks for confirmation. Stopping discards the unfinished file
and skips the remaining ones; finished files, your inputs and the databases are kept.

### Output files

- Protein FASTA gives `<name>.kannotate.csv`.
- GenBank gives `<name>.kannotate.gbk`, `<name>.kannotate.csv` or both.

The annotated GenBank keeps every gene location, locus tag and translation but
**replaces the functional annotations** (products, old protein IDs and similar) in the
copy. Keep the original file if you need those.

CSV tables are UTF-8 and open in Excel. To keep identifiers that look like numbers or
dates intact, import with **Data → From Text/CSV** and set those columns to **Text**.

- FASTA tables have one row per protein: header, sequence, Kannotate product and
  score, the top two Pfam and top two CDD matches (ID, name, description, E-value,
  residues and score), and the Kannotate, scoring, Pfam and CDD versions.
- GenBank tables have one row per gene: file, record, locus tag, **From**, **To**,
  **Strand**, **Length** (amino acids), original product, sequence, Kannotate product
  and score, and the same Pfam and CDD columns. From and To are 1-based and
  inclusive; joined genes show their outer bounds.

Kannotate annotates existing genes only: it doesn't find genes or translate
nucleotide FASTA.

### How hits are scored

Each Pfam or CDD match gets a score from its E-value and how much of the protein it
covers:

```text
L = min(300, -log10(E-value)),  c = fraction of the protein covered by the alignment
score = L × (0.5 + 0.5 × √c) + 5 × c
```

The best-scoring match names the protein if its score is at least **6**; otherwise
the protein is called a hypothetical protein. A few rules recognize multi-domain
enzymes (NRPS, PKS and hybrids) from several strong domains. **Show weaker hits**
also lists other matches the search programs report at E-value 0.01 or lower.

How to read the results:

- Pfam and CDD are independent evidence, even when a CDD model comes from Pfam.
- Coordinates are 1-based and inclusive in your protein. When a model matches
  several regions, the table shows their outer span; the diagram shows the gaps.
- HMMER (Pfam) E-values are for the whole protein; CDD values are for the best
  aligned region. Compare bit scores only within one database.
- A smaller E-value is stronger statistical evidence, not proof of function. The
  score ranks evidence; it is not a probability. Treat products as family-level
  hypotheses, and review generic, partial or very specific names.

## Compare gene clusters: Kstack

Kstack lines up gene clusters from GenBank files, groups similar proteins with
MMseqs2 (built in) and draws them as an editable figure.

### Build a stack

1. Drag GenBank files onto the window, or click **Browse / Add Files**. Use 1 to 50
   files, each with one record of up to 100,000 bp. Check the warnings, names and
   order.
2. For circular records, choose **Linearize at Origin** or **Omit**.
3. The project is saved as `Kstack.html` beside the first file (with a number added if
   that name exists). Click **Choose…** to save it elsewhere.
4. Adjust the grouping thresholds if needed (default 30% identity and 70% coverage of
   the shorter protein), the order and orientation, and the connectors.
5. Click **Generate Stack**.

### Edit the figure

- Drag rows to reorder or shift them; reverse rows; set gene and group colors; show
  or hide labels and connectors; filter connectors by identity; zoom and pan.
- Width, height and row spacing range from 25% to 400%.
- Right-click a gene and choose **Auto-align** to line up its group vertically under
  that gene, whose row stays fixed.
- Click a gene to see its details and protein. Rows with a **!** list genes that
  couldn't be drawn.

Figure edits never change the biology or the genome coordinates.

### Save, reopen and recompute

- **Cmd+S** saves and **Shift+Cmd+S** saves a copy. There is no autosave; Kstack
  warns before closing unsaved work.
- **Cmd+O** or **Recent Stacks** reopens a project.
- **Homology** recomputes the groups at new thresholds, keeping your layout.

The saved HTML is the whole project. It opens in any web browser for viewing,
selecting and exporting, with no other files needed. Its format is described in
[Kstack project format](kstack/project-format.md).

### Export SVG and FASTA

1. **Export SVG** opens a preview. Choose which clusters to include, drag rows to
   reorder or shift them, and adjust size, labels and connectors. These changes
   affect only the exported figure.
2. Click **Export SVG…** to save. The preview stays open, so cancelling loses nothing;
   **Reset to stack** starts the preview over.
3. To export proteins, add genes or whole clusters to the FASTA selection, then choose
   **Export Protein FASTA…**.

## Predict a structure and find structural homologs: Kmatch

Kmatch predicts one protein's structure with Boltz-2 on your Mac's GPU, then searches
it against the PDB100 database with Foldseek. Set it up first (see the
[README](../README.md#kmatch)).

### Predict and search

1. Paste a protein sequence or one FASTA record. Use only the 20 standard amino acids,
   at least 10 residues. Ambiguous letters (X, B, Z, J, U, O), digits, gaps and `*`
   are rejected, not replaced.
2. Choose a **Sequence alignment** mode:
   - **Local only · no MSA** (the default): nothing leaves your Mac.
   - **Local A3M file**: click **Choose A3M…** and select an alignment whose first
     sequence is your protein.
   - **Network-assisted MSA · recommended for accuracy**: **sends your sequence to
     the ColabFold MSA service** (`api.colabfold.com`), keeps the alignment in the
     cache and folds on your Mac. Use it only for sequences you may share. Running
     the same sequence again reuses the cached alignment.
3. Click **Run Kmatch**. The prediction runs first, then the search.

An alignment can improve difficult predictions considerably. Proteins longer than
about 400 residues can need much more time and memory.

**Cancel** stops the run and keeps whatever has finished. You can run again without
restarting the app. The window shows only the current run; earlier runs stay in
their folders.

### Read the structure

Drag to rotate and scroll to zoom; **Reset view** recenters. **Color by confidence**
colors residues by pLDDT (blue ≥ 90, light blue 70 to 90, yellow 50 to 70, orange
< 50). The summary shows mean pLDDT, pTM and runtime. A prediction made without an
alignment that has mean pLDDT below 70 or pTM below 0.5 shows a warning suggesting an
alignment. Confidence is the model's own estimate, not experimental validation.

### Structural homologs

The table shows up to 50 hits:

- **Structural similarity / qTM** (default) ranks by TM-score normalized to your
  protein's length.
- **Homology significance** ranks by Foldseek's probability and E-value.

Switching doesn't rerun the search; all hits are saved in the run folder. Names such
as `1hza-assembly1_B` give the PDB entry, biological assembly and chain. Select a row
for more scores, the alignment and Foldseek's rotation and translation.

With **Fetch public PDB details and coordinates (cached)** on, Kmatch shows protein
names and the entry's small molecules and ions, fetched from RCSB using only the PDB
IDs. These describe the whole PDB entry, not your protein's function or binding.
Switch it off to use cached details only.

### Save structures and 3D alignments

- **Save structure…** saves the predicted structure (PDB format).
- **Open run folder** shows all files of the run.
- **Export 3D alignment…** saves the selected hit superposed on your prediction, as
  a PyMOL session (`.pse`) or as coordinates (`.cif` or `.pdb`).

With **Include full-atom hit and contacting small molecules / ions** on (the
default), Kmatch downloads the hit's public structure from RCSB (when public PDB
access is on), checks it against the search result and keeps small molecules and
ions within 4 Å of the hit chain. These are neighbors of the hit, not predicted
ligands of your protein. Turn it off for a simpler export with only the hit's Cα
trace, made without downloading anything.

Exported files contain your protein as chain `Q`, the hit as `T` and retained ligands
as `L`. Nothing is refitted: the superposition is Foldseek's. PyMOL sessions need
PyMOL, which Kmatch doesn't install; it looks in `/Applications`, `~/Applications`
and on `PATH`. Coordinate files need nothing extra.

### Kmatch files and settings

Kmatch keeps its data in `~/Library/Application Support/KBioSuite/Kmatch/`:

| Folder or file | Contents |
| --- | --- |
| `config.json` | Optional settings (see below) |
| `backends/` | The Boltz-2 Python environment (and AlphaFold 3 code, if installed) |
| `models/` | Model files |
| `tools/` | Foldseek and the pinned uv used by setup |
| `databases/pdb100/` | The PDB100 database |
| `runs/` | One folder per run |

Alignments and downloaded public PDB data are kept in
`~/Library/Caches/KBioSuite/Kmatch/` (`msa/` and `pdb/`). Kmatch never deletes runs.

`config.json` is optional. Without it, Kmatch uses:

```json
{"predictor": "boltz2", "threads": 8, "diffusion_samples": 3, "seed": 42}
```

`threads` can be 1 to 32, `diffusion_samples` 1 to 5 and `seed` 0 to 4294967295.
Optional keys point Kmatch at existing installations by absolute path:
`boltz_python_override`, `boltz_models_override`, `foldseek_override`,
`pdb_prefix_override` and `pymol_override` (plus `af3_python_override`,
`af3_source_override` and `af3_models_override` for AlphaFold 3). Setup never
overwrites a path you configured. `shared_resources` names a record in a shared
Resources library; set it to `null` to use Kmatch's own Foldseek and PDB100. Unknown
keys are rejected. Restart Kmatch after editing the file.

### Run folders

Each run is a folder `runs/<date-time>-<id>/` containing:

- `sequence.fasta` and `run.json` (status and settings);
- `prediction/`: the predicted structure, `prediction.json` (confidence and
  provenance), the engine's inputs and raw outputs, and `predictor.log`;
- `search/`: `foldseek.tsv`, `search.json` with every hit, `foldseek.log` and
  `target-traces/` (hit coordinates kept for later exports);
- `exports/`: one folder per 3D alignment export, with its status and provenance.

Failed and cancelled runs keep their partial files and logs. Because hit coordinates
are kept with the run, exports from both top-50 lists keep working after PDB100 is
updated.

## Command-line use

A [source build](../README.md#build-from-source) also provides command-line tools.
Run them from the `KBioSuite` folder:

```sh
uv run kblast                                   # Kblast desktop app
uv run kblast --data-dir /path/to/runs          # open another run store
uv run kannotate status
uv run kannotate setup
uv run kannotate annotate proteins.faa cluster.gbk -o results/
uv run kannotate annotate cluster.gbk --genbank-output both --threads 4
uv run kannotate annotate cluster.gbk -o results/ --checkpoint-dir batch-state/ --resume
uv run kmatch status
uv run kmatch setup boltz                       # also: foldseek, pdb100
uv run kmatch run query.fasta
uv run kmatch run query.fasta --msa query.a3m
uv run kmatch run query.fasta --network-msa     # sends the sequence to api.colabfold.com
printf '%s\n' '{"action":"inspect","paths":["cluster.gbk"]}' | uv run kstack
```

`kannotate annotate` options: `--genbank-output genbank|csv|both|xlsx`,
`--table-format csv|xlsx`, `--min-score` (default 6), `--threads` and `--force` (to
replace existing outputs). `--checkpoint-dir` keeps resumable checkpoints; repeat the
same command with `--resume` to continue. Add `--help` to any command for details.

## Troubleshooting

**KBioSuite says it can't find the other apps, or asks you to move it to
Applications.** Quit, drag the whole KBioSuite folder into Applications and open it
from there. Keep all five apps together in that folder.

**macOS asks whether to open an app downloaded from the internet.** This happens once
for each app after installing or updating. Click **Open**. If macOS says an app is
damaged or can't be checked, delete the KBioSuite folder and install a fresh
download; don't turn off macOS security settings.

**Kannotate says tools are missing.** Install HMMER and BLAST+ (`brew install hmmer
blast`), then click **Refresh status**.

**Kmatch setup asks for the Command Line Tools.** Run `xcode-select --install`, wait
for the installation to finish, then click **Install** again.

**A setup download stops.** Check your internet connection and click the setup button
again. Kmatch downloads resume where they stopped. Kannotate keeps files it has
already finished and downloads an interrupted file again from the start.

**Kmatch shows "Not loaded" for protein names.** Kmatch couldn't get the details
from RCSB or couldn't save them. Select the hit to read the reason for that hit under
**Small-molecule ligands in this PDB entry**; the line above the table counts the hits
that didn't load and gives the first reason. For a connection problem, check your
internet connection and click **Reload details**. Your
prediction, the hits and the scores are not affected, and the **RCSB** link under
each hit opens the entry in your browser. If the message mentions a security
certificate, check that your Mac's date and time are correct. On a managed network,
ask your administrator whether the network inspects HTTPS traffic.

**Kmatch 0.6.1 or earlier shows "Unavailable" for protein names, with a certificate
error in the ligand section.** RCSB changed its certificate authority in November
2026 and those versions can't verify it. Update to 0.6.2 or later. Details that were
already saved on your Mac still appear.

**A Kblast search waits for a long time.** NCBI queues scripted searches behind
interactive ones. Keep waiting, cancel and resume later, or
[use an existing NCBI search](#use-an-existing-ncbi-search).

**An app doesn't appear in Spotlight.** Open the KBioSuite folder in Applications in
Finder and double-click the app. Spotlight may need a few minutes to find newly
installed apps.

## Limitations

- **Kblast** depends on NCBI's live databases and queues: later searches can return
  different hits, and NCBI's links between proteins and genomes can be incomplete.
  It does not remove duplicates, filter by taxonomy or predict clusters.
- **Kannotate** annotates existing protein-coding genes only. Products are
  family-level hypotheses from Pfam and CDD, not experimental validation. The
  multi-domain rules are provisional.
- **Kstack** accepts up to 50 single-record GenBank files of up to 100,000 bp each.
- **Kmatch** predicts one protein chain per run: no complexes, ligands, templates or
  batches. Boltz-2 is the only predictor; AlphaFold 3 is not offered (see
  [AlphaFold 3 notes](kmatch/af3.md)). It has been checked on public proteins of
  67 to 394 residues; larger proteins are untested. Similar structures can have
  different functions.
- All apps run only on Macs with Apple silicon and macOS 26 or later.
