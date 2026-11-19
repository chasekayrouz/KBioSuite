# KBioSuite

Protein bioinformatics tools for Macs with Apple silicon. KBioSuite finds the genomic
neighborhoods of a protein's relatives, annotates protein domains, compares gene
clusters and predicts a protein's structure to search for structural relatives.

| App | What it does |
| --- | --- |
| **Kblast** | Searches a protein against NCBI's nr database with BLASTP and retrieves the genomic neighborhood of every hit |
| **Kannotate** | Annotates proteins and GenBank genes with Pfam and CDD domains, entirely on your Mac |
| **Kstack** | Compares and draws gene clusters from GenBank files; saves portable HTML projects and SVG figures |
| **Kmatch** | Predicts one protein's structure on your Mac and searches it against the PDB for structural homologs |
| **KBioSuite** | Launcher for the four apps above |

Each app runs on its own and works with ordinary files. Nothing is passed between
the apps automatically.

> **Status:** research software under active development. Domain annotations,
> predicted structures and homolog matches are hypotheses to check, not experimental
> evidence. See [Limitations](docs/user-guide.md#limitations).

## Requirements

- A Mac with Apple silicon (M1 or newer). Intel Macs are not supported.
- macOS 26 or later.
- About 300 MB for the apps. Kannotate and Kmatch also need large reference data,
  downloaded when you set them up (see [First-time setup](#first-time-setup)).
- An internet connection for Kblast searches and for setup downloads.

## Installation

**1. Download.** Open the [latest release](https://github.com/chasekayrouz/KBioSuite/releases/latest)
and click the file whose name ends in `-macOS-arm64.dmg` to download it.

**2. Install.** Double-click the downloaded file. In the window that opens, drag the
**KBioSuite** folder onto the **Applications** folder. Keep the five apps together in
that folder: the launcher opens the tools that sit beside it. Then eject the disk image
(the eject button next to KBioSuite in the Finder sidebar) and delete the download.

If you can't install into Applications (some managed Macs don't allow it), drag the
KBioSuite folder into the **Applications** folder inside your home folder instead.
Create that folder if it doesn't exist.

**3. Open an app.** Press **Cmd+Space**, type **KBioSuite** (or any app name) and press
**Return**. The first time you open each app, macOS asks whether to open an app
downloaded from the internet. Click **Open**.

If you open KBioSuite directly from the disk image, it can't find the other apps and
asks you to move the KBioSuite folder to Applications first.

## First-time setup

Kblast and Kstack are ready to use. Kannotate and Kmatch need tools and reference
data that you install once. Opening an app never installs tools or databases;
their downloads start only through explicit setup actions.

| App | Before first use | Download | Disk space |
| --- | --- | --- | --- |
| Kblast | Nothing (an NCBI contact email is recommended) | None | Saved runs only |
| Kstack | Nothing | None | None |
| Kannotate | Install HMMER and BLAST+, then set up the databases | About 2.2 GB | About 9 GB (12 GiB free during setup) |
| Kmatch | Install Apple's Command Line Tools, then set up Boltz-2, Foldseek and PDB100 | About 9 GB | About 18 GB (allow 25 GB during setup) |

### Kannotate

Kannotate runs two standard search programs that you install yourself: **HMMER**
(`hmmscan`, `hmmpress`) and **NCBI BLAST+** (`rpsblast`, `blastdbcmd`).

1. Install [Homebrew](https://brew.sh) if you don't have it.
2. In **Terminal**, run:

   ```sh
   brew install hmmer blast
   ```

3. Open **Kannotate**. Under **Local databases**, open **Database setup & location**
   and click **Set up / verify databases**.

Setup downloads Pfam-A 38.2 from EBI and NCBI's current CDD release, checks Pfam
against its pinned checksum and CDD against NCBI's release notes, prepares them and
runs a public test protein through both. When it finishes, the status line reads
*ready* and names both releases; results record them too. An installed CDD keeps its
release when NCBI publishes a newer one. If tools are missing, the status line names
them.

Kannotate looks for the tools on its `PATH` and in `/opt/homebrew/bin`,
`/usr/local/bin` and `/usr/bin`. Apps opened from Finder don't see a conda
environment, so Homebrew is the simplest route. If you build from source and prefer
[Bioconda](https://bioconda.github.io) (`conda install -c conda-forge -c bioconda hmmer blast`),
run `uv run kannotate setup` once in a Terminal where that environment is active:
setup records where it found the tools.

### Kmatch

1. Install Apple's Command Line Tools (skip this if you have Xcode). In **Terminal**, run:

   ```sh
   xcode-select --install
   ```

   and confirm in the window that appears. Boltz-2 setup checks for these tools first
   and downloads nothing without them.
2. Open **Kmatch**. The **Dependencies & setup** panel opens by itself while anything
   is missing.
3. Click **Install** for **Boltz-2**, **Foldseek** and **PDB100**, one at a time:
   - **Boltz-2** (the structure predictor): downloads a pinned copy of the
     [uv](https://docs.astral.sh/uv/) installer tool (about 20 MB), Python 3.12 and
     the Boltz-2 packages (about 1 GB), and the model files (about 6 GB). About 10 GB
     on disk.
   - **Foldseek** (the structure search program): about 20 MB.
   - **PDB100** (the PDB structure database): about 2.3 GB to download and 8.5 GB on
     disk. Needs 12 GB free while it is set up.

An interrupted download resumes when you click the button again.

The panel also shows an **AlphaFold 3** card. It is optional and untested in Kmatch:
**Install AF3 code** installs AlphaFold 3's code only, without its model parameters,
and Kmatch doesn't offer AlphaFold 3 as a predictor. See
[AlphaFold 3 notes](docs/kmatch/af3.md).

### Kblast

No setup is needed. Enter an email address in **NCBI contact email** so NCBI can
contact you if your searches cause a problem. An optional NCBI API key speeds up
the retrieval of genomic records; see the
[user guide](docs/user-guide.md#ncbi-contact-email-and-api-key).

## Update

Quit all KBioSuite apps, download the new version and drag the KBioSuite folder onto
Applications again. Choose **Replace** when asked. Your saved runs, settings and
downloaded databases are kept, because they are stored in your Library folder, not in
the apps.

## Uninstall

Delete these items in Finder (press **Cmd+Shift+G** to type a path):

- `/Applications/KBioSuite` (or `~/Applications/KBioSuite`): the apps.
- `~/Library/Application Support/KBioSuite`: saved Kblast runs; Kannotate databases,
  search cache and checkpoints; Kmatch engines, models, PDB100 database and runs; and
  a shared Resources library, if you set one up. This folder can be very large.
- `~/Library/Caches/KBioSuite`: Kmatch's sequence alignments and downloaded public PDB
  data.
- Items whose names start with `org.kbiosuite.` in `~/Library/Preferences`,
  `~/Library/Caches` and `~/Library/WebKit`: app settings and web view data.
- `~/.config/kblast`, if you created it. `~/.ncbi/api_key` only if no other software
  uses it.

Kmatch's Boltz-2 setup uses uv, which keeps Python 3.12 in `~/.local/share/uv/python`
and downloaded packages in `~/.cache/uv`. Delete those two folders only if you don't
use uv for anything else. If you installed HMMER and BLAST+ only for Kannotate, run
`brew uninstall hmmer blast`.

## Build from source

Building the apps yourself gives the same apps, installed in `~/Applications/KBioSuite/`.
It needs about 3 GB of free disk space and an internet connection for the first build,
which downloads Python 3.14 and the packages KBioSuite needs. You will install Apple's
Command Line Tools (or Xcode) and uv. You don't need to install Python yourself.

You type these commands in the **Terminal** app (press **Cmd+Space**, type *Terminal*,
press **Return**). Copy one command at a time, paste it and press **Return**.

**1. Install Apple's Command Line Tools.** These provide Swift, Git and a C compiler.
Skip this step if you have Xcode.

```sh
xcode-select --install
```

**2. Install uv.** [uv](https://docs.astral.sh/uv/) downloads Python and the packages
KBioSuite needs. Run its official installer, then open a new Terminal window:

```sh
curl -LsSf https://astral.sh/uv/install.sh | sh
```

If you use Homebrew, `brew install uv` works too.

**3. Download the source.** This creates a `KBioSuite` folder in your home folder.

```sh
cd ~
git clone https://github.com/chasekayrouz/KBioSuite.git
cd KBioSuite
```

**4. Download MMseqs2 for Kstack.** This fetches the pinned MMseqs2 release from
GitHub, checks its checksum and keeps it inside the `KBioSuite` folder.

```sh
./scripts/kstack/setup-mmseqs.sh
```

**5. Build and install.** Quit any KBioSuite apps that are open, then run:

```sh
./scripts/build-install.sh
```

The first build downloads Python 3.14 (managed by uv, without replacing any other
Python on your Mac), creates a private Python environment inside the `KBioSuite`
folder, compiles and signs the five apps for this Mac, checks them and installs them
in `~/Applications/KBioSuite/`. A local build is for your own use: it accepts any
Python 3.14 release that uv provides, while the published download is built with one
exact Python release. If a check fails, your previous installation is left
untouched. If you also installed the download in `/Applications/KBioSuite`, delete
one of the two copies so that Spotlight finds only one.

**6. Run the tests (optional).** The full test suite also needs
[Node.js](https://nodejs.org) (`brew install node`) and Xcode, because the Kstack
native tests use XCTest, which the Command Line Tools don't include.

```sh
./scripts/test.sh
```

To update, run `git pull` and then `./scripts/build-install.sh` again. The builder
installs only from a clean checkout of `main` and refuses to replace a newer
installed version with an older one. `src/kbiosuite/version.py` is the sole suite
version source: **0.6.1**.

A source build also gives command-line versions of the tools (`uv run kblast`,
`uv run kannotate`, `uv run kstack`, `uv run kmatch`); see the
[user guide](docs/user-guide.md#command-line-use).

## Using the apps

The [user guide](docs/user-guide.md) walks through each app step by step. The
[Kblast methods](docs/kblast/methods.md) describe its search settings and output
files, and the [Kstack project format](docs/kstack/project-format.md) documents its
saved HTML projects.

## Privacy

Kannotate and Kstack keep your sequences on your Mac. Kblast's purpose is an NCBI
search, so it sends your protein to NCBI. Kmatch folds and searches on your Mac and
sends your sequence elsewhere only if you choose a network alignment. No app sends
analytics or usage statistics, and the launcher never connects to the internet.

**Kblast**

- **Run search** sends your protein sequence to NCBI BLAST (`blast.ncbi.nlm.nih.gov`).
- It then requests records for the hits from NCBI E-utilities
  (`eutils.ncbi.nlm.nih.gov`). Your contact email goes with every NCBI request; an API
  key, if configured, goes only to E-utilities.
- **Open NCBI Protein BLAST** and the RID links open NCBI's website in your browser.

**Kannotate**

- Searches run on your Mac; sequences are never sent anywhere.
- Database setup downloads Pfam from EBI (`ftp.ebi.ac.uk`) and CDD from NCBI
  (`ftp.ncbi.nlm.nih.gov`). Kannotate turns off BLAST+'s usage reporting.

**Kstack** never connects to the internet.

**Kmatch**

- Prediction and structure search run on your Mac; the predictor runs with network
  access blocked.
- **Network-assisted MSA** is off unless you choose it for a run. It sends your
  sequence to the ColabFold MSA service (`api.colabfold.com`).
- **Fetch public PDB details and coordinates** (on by default; you can switch it off)
  sends only the PDB IDs of your hits to RCSB (`data.rcsb.org`) to show protein and
  ligand names. Full-atom 3D alignment export downloads public structure files from
  `files.rcsb.org`. Your sequence and predicted structure are never sent. **RCSB ↗**
  links open `www.rcsb.org` in your browser.
- Setup downloads, only when you click a setup button:
  - uv from GitHub (`github.com`);
  - Python 3.12 and the Boltz-2 packages through uv, from GitHub and the Python
    Package Index (`pypi.org`, `files.pythonhosted.org`);
  - Boltz-2 model files from Hugging Face (`huggingface.co`) and Boltz
    (`model-gateway.boltz.bio`);
  - Foldseek from `mmseqs.com` and PDB100 from `foldseek.steineggerlab.workers.dev`;
  - the optional AlphaFold 3 code and its packages from GitHub and the Python Package
    Index.

Building from source downloads Python and packages through uv (`github.com`,
`pypi.org`, `files.pythonhosted.org`) and MMseqs2 from `github.com`.

## License

KBioSuite is released under the [MIT License](LICENSE). Bundled third-party components
keep their own licenses; see [third-party notices](licenses/README.md).
