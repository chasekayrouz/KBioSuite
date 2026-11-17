# HTML project format: schema 1

The saved HTML is the source of truth. There is no `.kstack` bundle, hidden copy, database, automatically exported FASTA, or automatically exported SVG.

A document contains exactly one marker:

```html
<meta name="kstack-format" content="Kstack; schema=1; app=0.1.0">
<script id="kstack-project" type="application/json" data-schema="1">…</script>
```

The JSON has `format:"Kstack"`, integer `schema_version:1`, and `app_version` from the suite version source. [The JSON schema](../../src/kstack/resources/shared/project.schema.json) describes structural fields; `src/kstack/model.py` owns runtime validation including cross-references, unique IDs, lengths, order permutation, coordinate bounds, colors, and protein/alignment consistency. Kstack rejects unknown versions with an actionable error rather than attempting a lossy interpretation.

| Field | Meaning |
| --- | --- |
| `metadata` | UTC creation/modification timestamps, initial title, exact analysis provenance |
| `analysis` | identity and coverage fractions, preserve-input option |
| `visualization` | connector visibility, display identity fraction, label state, optional `width_scale`, `height_scale`, and `row_spacing_scale` presentation factors |
| `clusters` | stable ID, imported display name, source filename (no original absolute path), accession/description, record length/topology, useful record annotations, warnings, CDSs |
| `clusters[].cds` | stable ID, cluster reference, gene/product/locus tag, compound segments, original location text, strand, CDS length, complete protein/protein length, translation source, qualifiers, partial/problem flags, origin-spanning flag, homolog group |
| `homolog_groups` | connected components, stable component ID, member IDs, generated RGB color |
| `pairwise_hits` | canonical undirected query/target IDs, identity, aligned length, both protein lengths, `query_coverage`/`target_coverage` from MMseqs2 `qcov`/`tcov`, derived `shorter_coverage`, bit score/E-value; includes below-threshold candidates |
| `layout` | cluster order, row orientation ±1, horizontal offset in bp |
| `color_overrides` | gene/group ID → `#rrggbb` |
| `fasta_selection` | unique selected CDS IDs |

Coordinates are zero-based, half-open genomic segments. Inspector presentation converts them to one-based inclusive text. Orientation and offset are figure properties, never rewritten genomic coordinates. Cluster annotations preserve the original circular topology; `linearized_at_origin:true` records the explicit import choice without rotation. A CDS without usable coordinates or protein still has one stable biological entry and warnings. Nucleotides are not retained.

The HTML embeds CSS plus the pure model, standalone renderer, SVG exporter, and FASTA exporter as inline resources. JSON serialization escapes `<`, `>`, `&`, U+2028 and U+2029; annotations cannot terminate the data script. User text enters the live DOM through textContent and enters SVG through XML escaping. A restrictive CSP forbids network/external resources. Opening a file in the desktop app extracts and validates only the JSON; the app regenerates its own trusted renderer instead of executing arbitrary scripts from the selected HTML.

Save is explicit and atomic. The initial Generate action writes the chosen HTML. Subsequent persistent edits update the working state/dirty indicator until ⌘S or Save As regenerates the HTML. Width, height, and row-spacing changes are persistent visual edits. All three factors are finite numbers in [0.25, 4]; omitted fields mean 1, so existing schema-1 documents require no migration or dirtying on open. Width projects genomic X positions and lengths; Height sets intrinsic glyph size (native font size, arrow height/head length, strokes) and figure margins. Row-center distance is `90 × height_scale × row_spacing_scale`; spacing changes no glyph size or X geometry. Text has no anisotropic transform. They never rewrite coordinates, alignment offsets, or analysis. Viewport zoom/pan and inspector focus are temporary and do not dirty or alter exported geometry. In a standalone browser, visualization/selection changes last for that session; the original file remains intact. The app's Recent Stacks stores only references/bookmarks, with no duplicate HTML contents.

Main-window horizontal drag and gene Auto-align write the existing layout row
`offset` fields, saved through the normal explicit Save flow. Auto-align changes
neither orientation nor genomic segments; no new schema fields are required.

SVG export-preview visibility, ordering, alignment offsets and presentation settings are temporary,
separate from these persistent fields. Opening the preview copies the current
stack; canceling or closing discards that copy. Export writes only the selected
figure and does not alter the main stack or the HTML's JSON marker. Reopening
starts from the main stack again. This requires no schema migration. Standalone
HTML saved with the new renderer embeds the same preview controls; older HTML
receives them when reopened and explicitly saved in the updated desktop app.

Threshold recomputation runs a fresh MMseqs2 search from persisted protein sequences. It does not require source GenBanks, and does not rely on the original candidate list being complete for newly lowered thresholds. It preserves row layout, gene overrides, selections and visualization settings. Group overrides survive only when the exact member set remains the same.

Sequence coverage is bounded to [0,1] and validated against the stored query/target coverage: use query coverage for a shorter query, target coverage for a shorter target, and their minimum for equal lengths. Gapped alignment length is not a coverage denominator. Earlier pre-baseline development HTML without these sequence-coverage fields must be regenerated from its stored proteins or source GenBanks; it is not silently interpreted using the old definition.
