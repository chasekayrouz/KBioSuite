# Shared bioinformatics resources

A Resources library lets several applications use one verified copy of large
reference databases and search tools. It is optional: Kannotate and Kmatch install
their own databases when no library is configured.

The library lives in `~/Library/Application Support/KBioSuite/Resources/`.
`KBIOSUITE_RESOURCES` selects another root for one process; no shell configuration is
needed. KBioSuite manages the library; other applications read it through the
documented JSON protocol below. Opening an app never downloads, switches or deletes a
resource.

**Shared resources…** in the KBioSuite launcher opens `index.html`, an offline
inventory of installed versions and the applications that use them. It is regenerated
by adoption, registration and `resources view`. Its status column is a quick
file-presence check made at that time, not a checksum verification, and it does not
check for newer upstream releases.

## Commands

From a KBioSuite source checkout with its locked environment:

```sh
uv run --offline --locked kbiosuite resources list
uv run --offline --locked kbiosuite resources resolve 'RESOURCE_ID' --entrypoint database
uv run --offline --locked kbiosuite resources verify
uv run --offline --locked kbiosuite resources adopt /path/to/recipe.json
uv run --offline --locked kbiosuite resources register consumer-name /path/to/pins.json
uv run --offline --locked kbiosuite resources view
```

These commands need no network or background service. `adopt` copies files that
already exist on disk; it never downloads. To add a new release, install it with the
official, version-pinned installer in a separate staging folder, validate it, then
adopt it. `resolve` prints a JSON string holding an absolute path. A BLAST database
entrypoint is a prefix, which need not itself be a file.

## Protocol v1

- `catalog.json`: `schema_version: 1` and `resources`, keyed by immutable resource ID.
  Each entry has a relative `path`, content `fingerprint`, `kind`, `name`, `version`
  and apparent file size in `bytes`.
- `<path>/resource.json`: schema, ID, fingerprint, entrypoints, the complete file map,
  provenance and adoption time. Each file has `bytes`, `sha256` and `executable`. The
  fingerprint is the SHA-256 of compact, sorted-key JSON containing `kind`, `name`,
  `version`, `entrypoints` and `files`. The resource ID is `name@version:first16hash`.
- `consumers/<name>.json`: `schema_version: 1` and `bindings`, mapping an
  application's roles to `{"resource": "RESOURCE_ID", "entrypoint": "database"}`.
  Each binding pins one generation, never a mutable `latest` link. Optional
  `path_mappings` let a consumer translate its own historical relative paths: each
  maps a `prefix` to a `resource`, and may add an `entrypoint` to map one exact file.
  Other optional fields, such as `owner` and `retention`, give context.
- `history/`: earlier consumer registrations, kept on every change.
- `staging/`: unpublished copies. Interrupted work is kept for inspection.
- `migrations/`: optional local receipts for explicit configuration changes.

The read-only reference client is `kbiosuite/resource_client.py`, a small
standard-library module. Applications in KBioSuite import it directly; an
independent application can carry a copy that implements protocol v1 without
importing KBioSuite or sharing its Python environment. Keep such copies in step with
the protocol, and introduce a new schema version for any incompatible change.

Consumers validate the schema, path containment, identity, file presence and sizes.
They resolve pins once per configuration load or job and record the actual
generation and tool version in their normal provenance. A missing pinned resource is
an error, never permission to search elsewhere or silently pick another release.
Explicitly configured tool paths take precedence over `PATH`; `PATH` discovery
applies only when no executable was configured.

## Adoption and validation

A recipe names a source folder and the exact files to adopt:

```json
{
  "kind": "databases",
  "name": "example",
  "version": "1.0",
  "source": "/absolute/source/directory",
  "files": {"database.bin": "EXPECTED_SHA256"},
  "entrypoints": {"database": "database.bin"},
  "provenance": {"source": "official release and retrieval record"}
}
```

`kind` is `databases`, `tools` or `models`. An expected hash may be `null` for an
initial local inventory: adoption then computes and records it, but that alone does
not prove where the file came from. Supply expected hashes from the release's own
manifest whenever one exists. Hidden files, path traversal and symbolic links that
leave the source folder are rejected.

Adoption checks the source hashes, makes independent copies (APFS clones when
possible), checks the copies, then publishes the complete generation and updates the
catalog atomically under a writer lock. A published generation is never modified. If
a crash leaves an unregistered generation or staging folder, inspect it; a retry
never overwrites it. Registration validates every requested pin, archives the
previous registration and replaces the consumer record atomically.

`verify` hashes every file and checks the manifest fingerprint. File integrity is
separate from scientific qualification: before switching a consumer to a different
build, compare its results on known public controls. Keep several versions of a tool
when different consumers need them. Scoring, thresholds and scientific
interpretation stay with each application.

## Storage and retirement

There is deliberately no automatic deletion and no unqualified update command.
Retiring a resource is a separate, explicit decision that considers consumer pins,
earlier registrations, research provenance and running jobs. A registered consumer is
a declaration, not proof that every script using a resource is known. Unresolved or
legacy material stays in place. `inventory.json`, when present, lists other
discovered collections and open ownership questions; it is not a deletion list.

Databases, executable tools and optional public model files may be shared. Research
originals, private sequence searches and results, credentials and application Python
environments stay with their owners; receipts must not contain secrets or private
research inputs. On APFS, file sizes are not exact reclaimable space: clones and
retained original copies share blocks and must be accounted for separately.

## KBioSuite consumers

- **Kannotate** uses the library when its configuration contains
  `"shared_resources": "kannotate"`. The bindings supply its Pfam and CDD databases and
  its HMMER and BLAST+ executables. Setup then verifies the pinned files and runs its
  public search controls instead of downloading.
- **Kmatch** uses the `kmatch` record for Foldseek and PDB100 when its `config.json`
  sets `"shared_resources": "kmatch"`, or by default when that key is absent and
  neither `foldseek_override` nor `pdb_prefix_override` is set. If the library or the
  `kmatch` record doesn't exist, Kmatch uses its own managed copies. A record or root
  that exists but doesn't resolve stops Kmatch with an explanation rather than falling
  back. `"shared_resources": null` (or `false`) opts out. While bound, Foldseek and
  PDB100 setup only verify the pinned files; they never download, switch or retire a
  generation.

Before changing any consumer's configuration, keep a copy of its previous
configuration and test its existing outputs.
