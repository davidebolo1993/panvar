# Call Module (Module 3)

Date: 2026-04-28

CLI entrypoint:

- `panvar call`

## What it does

`call` consumes:

- module-1 bubbles (`panvar bubble`)
- module-2 allele clusters + assignments (`panvar allele`)

For each bubble, cluster representatives are compared to the reference-cluster representative, and only structural events are emitted:

- `INS`
- `DEL`
- `INV`

Optional INS subtyping/copy-number annotation can be enabled with `--classify-ins`.
Optional gene copy-delta annotation can be enabled with `--pangene-bed`.

## Required inputs

- `--gfa <graph.gfa>`
- `--bubbles-csv-in <module1.bubbles.csv>`
- `--clusters-csv-in <module2.allele_clusters.csv>`
- `--assignments-csv-in <module2.allele_assignments.csv>`
- `--reference-path <path>`

`call` does not re-cluster alleles.

## Calling algorithm

For each bubble:

1. Identify the reference cluster using `--reference-path`.
2. For each cluster representative, map to the reference-cluster sequence with minimap2 C API (`--minimap-preset`, default `asm20`; `-c`, `--eqx` behavior).
3. Choose the best orientation (`+`/`-`) from the best minimap2 hit and keep all per-cluster hits/chains for evidence.
4. Extract atomic `INS`/`DEL` events from CIGAR (`eqx`) on the best hit.
   - CIGAR `I` -> `INS` candidate
   - CIGAR `D` -> `DEL` candidate
5. Extract `INS`/`DEL` events from split-read/split-chain junction deviations (same orientation), by default.
   - split `INS` length (`SVLEN`) uses geometric deviation by default (SVIM-like),
     configurable with `--split-ins-svlen-mode`.
6. Union the candidate event list (CIGAR-derived + split-derived).
7. If no SV-sized indel is found but there is a large net length shift, recover one net `INS` or `DEL` from prefix/suffix decomposition.
8. Add `INV` when reverse orientation is selected with sufficiently low normalized edit distance.
9. Merge nearby concordant events within the same cluster.
10. Optional (`--classify-ins`): classify `INS` as `NOVEL` vs DUP-like subtype and estimate copy number fields.
    - duplicate-source search tries exact forward and reverse-complement matches first
    - if exact matching fails, bounded approximate matching is used in both orientations
    - approximate mode uses DP refinement for smaller inserts and sketch-based scoring for large inserts
      (large trigger: inserted length > 5000 bp or > 50% of bubble reference span)
11. Optional (`--pangene-bed`): annotate event clusters with gene copy deltas
    (representative haplotype vs reference haplotype) and optionally tune INS subtype.
    - `--pangene-gene-match <expr>` can restrict annotation to genes of interest
    - `--pangene-tune-ins` sets `INS_SUBTYPE=DUP_PANGENE` when INS has copy-gain support
      and no stronger local DUP evidence
12. Keep events passing `--min-sv-bp` (default `50`).
13. Merge equivalent events across clusters before writing VCF (same type, nearby coordinates, similar sequence).

## Within-cluster merge rule (important)

The per-cluster merge is intended to compact fragmented evidence of the same event, not to erase distinct events.

Current behavior:

- merge candidates must be same event type (`INS` with `INS`, `DEL` with `DEL`)
- both reference and alternate coordinates must be near (`merge_gap = max(1, min_sv_bp / 5)`)
- when two candidates come from different evidence sources (CIGAR vs split):
  - they are merged only if their sizes are similar enough (size similarity >= `0.70`)
  - this avoids collapsing a short CIGAR `DEL` into a much larger split-derived `DEL`

## Cross-cluster merge in VCF

When cluster A and cluster B carry the same event in the same context, one VCF record is emitted.

Default merge criteria:

- same `EVENT` type
- coordinate proximity within `20 bp` (`--vcf-merge-window-bp`)
- sequence similarity >= `0.80` (`--vcf-merge-min-seq-sim`)
- bounded comparison edit fraction <= `0.35` (`--vcf-merge-max-edit-frac`)

- `INFO/CLUSTERS`: cluster IDs carrying the merged event
- `INFO/MERGED_EVENTS`: source events merged into this record (`C<cluster>E<event>`)

Per-sample genotype for merged events:

- `GT=1` carrier cluster
- `GT=0` assigned non-carrier
- `GT=.` sample has no assignment at that bubble

## Outputs

### Main output

- region-level VCF (`--vcf-out`, default `<out>.region.vcf`)
- optional pangene copy table (`--pangene-copy-tsv`, default `<out>.pangene_copy.tsv` when `--pangene-bed` is enabled)

### Debug output (`--debug` or `--debug-out-dir <dir>`)

If debug is enabled, the tool writes per-bubble/per-cluster alignment artifacts only.

Layout:

- `<debug>/bubble_<id>/cluster_<cluster_id>/reference.fa`
- `<debug>/bubble_<id>/cluster_<cluster_id>/representative.fa`
- `<debug>/bubble_<id>/cluster_<cluster_id>/alignment.paf`
- `<debug>/bubble_<id>/cluster_<cluster_id>/dotplot.svg`
- `<debug>/bubble_<id>/cluster_<cluster_id>/cluster_vs_reference.vcf`

Notes:

- no decision traces and no CIGAR interpretation TSVs are emitted in debug mode
- `alignment.paf` is emitted from minimap2 C API hits (preset from `--minimap-preset`, default `asm20`; `-c`, `--eqx`)
- dotplot is drawn from alignment segments (CIGAR-aware), so indels appear as breaks rather than a single forced diagonal
- called variant breakpoints are overlaid as grey dashed guides:
  - `DEL`: vertical lines at reference breakpoints
  - `INS`: single vertical insertion breakpoint on the reference axis
  - `INV`: both vertical (reference) and horizontal (assembly/query) breakpoints
- reference coordinate ticks are drawn on the X axis with sparse, equally spaced labels
- title is kept compact (`Bubble <id> Cluster <id> Dotplot (PAF)`)
- X axis label is the reference name; Y axis label is the representative allele name
- primary alignment segments are solid; supplementary/secondary segments are dashed (no in-plot legend)
- optional gene overlay:
  - `--dotplot-gtf <path>` enables gene lookup
  - `--dotplot-gene-match <expr>` filters genes by case-insensitive regex/term (repeatable flag)
  - matched genes are drawn directly on the alignment line segments (not as axis bars)
  - a gene-color legend is rendered below the plot area (outside the track)

## VCF INFO fields (key)

- `SVTYPE`, `SVLEN`, `END`
- `EVENT`, `INS_SUBTYPE`, `ORIENT`, `BEST_NORM_ED`
- `INSSEQ` / `DELSEQ` / `INVSEQ`
- `DUP_SIM`, `DUP_REF_START`, `DUP_REF_END`, `DUP_ORIENT`, `DUP_UNIT_BP`
- `DUP_REF_CN`, `DUP_ALT_CN`, `DUP_ADDED`, `DUP_COPY_RATIO`
- `PANGENE_CN_DELTA`, `PANGENE_GAIN_GENES`, `PANGENE_LOSS_GENES`
- `PANGENE_GAIN_COPIES`, `PANGENE_LOSS_COPIES`

Notes:

- `ORIENT` describes global cluster-vs-reference orientation for the called allele
- `DUP_ORIENT` describes local duplicated-source orientation for INS DUP-like matching
  (they can differ, e.g. `ORIENT=+` with `DUP_ORIENT=-`)
- `INS_SUBTYPE=DUP_PANGENE` means INS tuning used pangene copy-gain support
  (without overriding `SVTYPE=INS`)
- `CLUSTERS`, `MERGED_EVENTS`
- `BUBBLE_ID`, `REF_CLUSTER_ID`, `BUBBLE_SOURCE`, `BUBBLE_SINK`

## Key options

- `--min-sv-bp <N>`
- `--vcf-out <path>`
- `--debug`
- `--debug-out-dir <dir>`
- `--dotplot-gtf <path>`
- `--dotplot-gene-match <expr>` (repeatable)
- `--minimap-preset <asm5|asm10|asm20>` (default `asm20`)
- `--minimap-best-n <N>` (default `8`)
- `--minimap-no-secondary`
- `--split-ins-svlen-mode <query-span|geometric>` (default `geometric`)
- `--classify-ins` (default off; uses bounded approximate DUP search for large INS via seeded candidates and sketch similarity to preserve sensitivity without pathological runtime)
- `--pangene-bed <path>` (optional pangene BED/BED.GZ copy-count annotation layer)
- `--pangene-gene-match <expr>` (repeatable gene filter for pangene annotation)
- `--pangene-tune-ins` (optional INS subtype tuning from pangene copy gains)
- `--pangene-copy-tsv <path>` (optional per-bubble/cluster copy-count table)
- `--vcf-merge-window-bp <N>` (default `20`)
- `--vcf-merge-mode <strict|lenient>` (default `strict`)
- `--vcf-merge-lenient-window-bp <N>` (default `100`)
- `--vcf-merge-lenient-min-ref-jaccard <X>` (default `0.30`)
- `--vcf-merge-min-seq-sim <X>` (default `0.80`)
- `--vcf-merge-max-edit-frac <X>` (default `0.35`)
  - effective cap is `min(--vcf-merge-max-edit-frac, 1 - --vcf-merge-min-seq-sim)`

Merge behavior:

- `strict` (default): requires start proximity within `--vcf-merge-window-bp`, and for non-INS also end proximity within the same window
- `lenient`: if strict checks fail, allows:
  - `INS` merge within `--vcf-merge-lenient-window-bp`
  - non-INS merge when starts are within `--vcf-merge-lenient-window-bp` and reference-interval Jaccard overlap is at least `--vcf-merge-lenient-min-ref-jaccard`

## Example (C4)

```bash
./build/panvar call \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/call \
  --bubbles-csv-in tests/results/c4/bubble.bubbles.csv \
  --clusters-csv-in tests/results/c4/allele.allele_clusters.csv \
  --assignments-csv-in tests/results/c4/allele.allele_assignments.csv \
  --reference-path grch38#1#chr6:31891045-32123783 \
  --dotplot-gtf /path/to/gencode.annotation.gtf.gz \
  --dotplot-gene-match C4 \
  --pangene-bed tests/real_data/c4.pangene.bed.gz \
  --pangene-gene-match C4 \
  --pangene-tune-ins \
  --debug
```
