# Module `call`

CLI: `panvar call`

## What it does

Types structural variants on the pangenome graph into a multi-sample VCF (Variant Call Format). For each bubble it:
- diffs every haplotype's `source → sink` walk against a reference walk and types each difference (`DEL`/`INS`/`INV`/`DUP`, see the [table](#event-types) below)
- reduces fragmentation by coalescing events within a haplotype and clustering equivalent events across haplotypes
- recovers sub-threshold carriers by force-calling them against the merged record
- reads copy number (`CN`) directly off the walk (`--cn`)

Algorithm and worked trace: [algorithms/call.md](../algorithms/call.md).

The recommended path is `bubble → panphorte → [refine] → call`, so `call` runs on the `panphorte` graph, or on the `refine` graph when the opt-in `refine` step is used. `panphorte` folds genuine tandem repeats into `REP` self-loops — which a tandem needs for an exact count — and leaves paralog clusters and rare duplications untouched; `refine` then POA-realigns the bubble interiors to remove graph-builder artifacts while holding those `REP` self-loops fixed, so it is copy-number-safe and the substrate `call` sees is the same folding with cleaner interiors. Calling directly on the `bubble` graph also works for a non-tandem paralog cluster and gives identical results there, since `panphorte` leaves such clusters untouched, but it misses the tandem folding.

#### Event types

| type | meaning |
|------|---------|
| `DEL` | reference-only nodes (deleted from the haplotype) |
| `INS` | haplotype-only nodes (inserted vs reference); `--classify-ins` adds `INS_SUBTYPE=NOVEL\|DUP` |
| `INV` | a haplotype run that is the reverse-complement node-walk of a reference run |
| `DUP` | a copy-number gain/loss (`--cn`) |

## Required inputs

- `-i, --gfa <graph.gfa>` — the call substrate from `panphorte`/`refine`/`bubble` (in any case, node ids should match the CSV to `-b`)
- one of `-b, --bubble-prefix-in <prefix>` or `-c, --bubbles-csv <path>`
- `-r, --reference-path <name>` — the diff baseline (full name or unique case-insensitive substring).
- `-o, --out-prefix <prefix>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--min-sv-bp <N>` | minimum size of a reported (merged) event | `50` |
| `--rescue-min-bp <N>` | floor for sub-threshold events kept for merge/rescue | `min-sv-bp/2` |
| `--merge-distance-bp <N>` | coalesce nearby same-type events (reference or haplotype space) | `100` |
| `--merge-jaccard <X>` | node-set [Jaccard](../algorithms/call.md#merge-keys--jaccard-vs-sequence-identity) to merge events | `0.80` |
| `--merge-seq-identity <X>` | event-sequence identity to merge | `0.80` |
| `--merge-size-ratio <X>` | length-ratio floor for the sequence merge (lower to merge more differing lengths) | `0` (off) |
| `--min-haplotypes <N>` / `--min-maf <X>` | drop records below N carriers / carrier frequency `AF=AC/AN` | `1` / `0` |
| `--cn` | emit the total copy number of the folded module; with `--gtf` it is reassigned to per-gene `CN` where the paralogs are separable. See [algorithms/call.md](../algorithms/call.md#copy-number) for the routes and [calling without `--cn`](../algorithms/call.md#copy-number) for what surfaces without it | off |
| `--classify-ins` | refine `INS` subtype `NOVEL`/`DUP` via minimap2 | off |
| `--multiallelic-loci` | collapse a bounded locus into one multiallelic record; `--multiallelic-max-bp` (5000) bounds it | off |
| `--gtf <path>` | gene annotation (needs a PanSN (Pangenome Sequence Naming) `--reference-path`); see [below](#gene-annotation---gtf) | — |
| `--bubble-id <N>` / `--no-per-bubble-vcf` / `--no-variant-nodes` / `-q, --quiet` | scope and output toggles | — |


## Outputs

| file | contents |
|------|----------|
| `<prefix>.region.vcf` | all records, coordinate-sorted, unique IDs (`bgzip`+`tabix`-able) |
| `<prefix>.bubble_<id>.vcf` | per-bubble VCF (unless `--no-per-bubble-vcf`) |
| `<prefix>.variant_nodes.tsv` | per-variant node set — the `describe --variant-nodes` handoff and the `benchmark` input (unless `--no-variant-nodes`) |
| `<prefix>.node_genes.tsv`, `<prefix>.dup_gene_cn.tsv` | with `--gtf` (see below) |

VCF 4.2; samples = haplotypes. `FORMAT`: `GT` (`1` carrier / `0` ref-like / `.` doesn't traverse), `CN` (per-sample copy number on `DUP`), `CNBP` (`DUP` only: the actual linear bp this haplotype gains `+` / loses `−` across the module vs the reference: `Σ node_length × traversal_multiplicity` from the bubble `source→sink`, minus the reference's). Key `INFO`:

| field | meaning |
|-------|---------|
| `SVTYPE`, `SVLEN` | event type and length difference (ALT − REF) |
| `END`, `BUBBLE_ID`, `START_NODE`, `END_NODE`, `EVENT_NODES` | locus span + graph anchors + the variant's node set |
| `AN`/`AC`/`AF` | haplotypes traversing the bubble/carriers/carrier freq `AC/AN` (ALT freq, not folded) |
| `NMERGED`, `SVLEN_RANGE`, `MERGE_JACCARD`/`MERGE_SEQID`/`MERGE_SIZE_RATIO` | merge count, size span, merge evidence (≥2-event records only) |
| `EVENTID` | links a co-located `DEL`+`INS` substitution |
| `INS_SUBTYPE` | `NOVEL`/`DUP` (`--classify-ins`) |
| `REF_CN` / `RU_LEN` | `DUP` only: reference copy number / one-copy repeat-unit length (the folded unit; for the real per-haplotype linear size see `FORMAT:CNBP`) |
| `GENES` | `--gtf`: gene(s) overlapped (whole folded module for a `DUP`) |
| `NALLELES` | `--multiallelic-loci`: number of alleles |
| `INSSEQ`/`DELSEQ`/`INVSEQ` | event sequence (omitted when very long) |


## Plotting

`scripts/plot_vcf_map.R` draws the headline oncoprint-style map (rows = haplotypes, columns = variants grouped by bubble; `DEL` red, `INS`-`NOVEL` green, `INS`-`DUP` purple, `INV` orange, multiallelic yellow-amber, `DUP` shaded blue by `FORMAT:CN`). It needs only the region VCF (`Rscript` + `ggplot2`):

```bash
Rscript scripts/plot_vcf_map.R \
  --vcf <prefix>.region.vcf \
  --out <prefix>.vcf_map
```

Flags:

- `--vcf <file>` — the `call` region VCF (required).
- `--out <prefix>` — output prefix (required).
- `--clusters` / `--cluster-by <tsv>` — plot only cluster representatives / group and order rows by cluster (mirrors the inspect heatmaps).
- `--max-paths <N>` — cap the number of haplotype rows.
- `--flip` — transpose (variants as rows).
- `--scale` / `--scale-transform <raw|sqrt|log1p>` — size each cell by `|SVLEN|` (`RU_LEN` for a `DUP`), with an optional transform.
- `--reference-path <name>` — optionally pin a haplotype (substring match) as the top (or leftmost, with `--flip`) row.
- `--width` / `--height` / `--dpi` — figure size (inches) and PNG resolution (default 300).

### Interactive node-coverage viewer

For a per-node view — how each haplotype actually traverses the graph under every call — `scripts/variant_node_heatmap_app.R` is a Shiny + plotly app. The top panel is a coverage heatmap (rows = haplotypes, columns = variant-affected nodes ordered by genomic position, width ∝ node length): white = not traversed, grey = ×1, red-gradient = ×2+ (so a `DUP` reads as a copy-number gradient). The bottom panel marks each `variant_id`'s node set. Hovering a node shows its gene, coverage, and — on a variant node — that haplotype's `GT` (+ `CN`/`CNBP` for `DUP`). It opens on a representative subset (the reference plus one carrier per variant / per distinct `DUP` CN); pick a bubble to zoom, clear the box to show all haplotypes, or click a `variant_id` to select its carriers.

First assemble the bundle with `scripts/build_variant_node_data.R` (needs `data.table`), then launch the app (`shiny`, `plotly`, `data.table`, `DT`):

```bash
Rscript scripts/build_variant_node_data.R \
  --gfa <call-graph.normalized.sorted.gfa> \   # the panphorte or refine graph call ran on
  --variant-nodes <prefix>.variant_nodes.tsv \
  --vcf <prefix>.region.vcf \
  --bubbles <panphorte-prefix>.bubbles.csv \
  --node-genes <prefix>.node_genes.tsv \
  --out variant_nodes.rds

VN_RDS=variant_nodes.rds Rscript scripts/variant_node_heatmap_app.R
```

`--node-genes` is optional (drops the gene hover if omitted). The `--gfa` is the graph `call` ran on and `--bubbles` its panphorte `bubbles.csv`, so node order and bubble spans match the VCF.


## Gene annotation (`--gtf`)

A reference-coordinate GTF (Gene Transfer Format; Ensembl/GENCODE) projected onto the graph via a PanSN reference (`sample#hap#chrom:start-end`); else skipped with a warning. It adds: `INFO=GENES` per record; `<prefix>.node_genes.tsv` (`node_id → gene(s)`); and `<prefix>.dup_gene_cn.tsv` — the per-haplotype per-gene copy number split of a folded paralog cluster, resolved from the GTF alone by k-mer dosage. Each gene's coding sequence (its gene span if it has no CDS) supplies a set of k-mers private to it vs its paralogs; a haplotype's per-copy count of those k-mers is the gene's copy number. Near-identical paralog pairs — where a plain count blurs under gene conversion — are split instead by per-site consensus over their divergent coding columns. See [algorithms/call.md](../algorithms/call.md#gene-annotation---gtf).

## Example

See the [LPA walkthrough](../walkthrough.md) for this module in a full end-to-end run.
