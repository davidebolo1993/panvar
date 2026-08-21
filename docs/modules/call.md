# Module `call`

CLI: `panvar call`

## What it does

Types structural variants on the pangenome graph into a multi-sample VCF (Variant Call Format). For each bubble it:
- diffs every haplotype's `source → sink` walk against a reference walk and types each difference (`DEL`/`INS`/`INV`/`DUP`, see the [table](#event-types) below)
- reduces fragmentation by coalescing events within a haplotype and clustering equivalent events across haplotypes
- recovers sub-threshold carriers by force-calling them against the merged record
- reads copy number (`CN`) directly off the walk (`--cn`)

Algorithm and worked trace: [algorithms/call.md](../algorithms/call.md).

The recommended path is `bubble → panphorte → [refine] → call`, so `call` runs on the `panphorte` graph, or on the `refine` graph when the opt-in `refine` step is used. `panphorte` folds genuine tandem repeats into repeat-unit self-loops, which a tandem needs for an exact count, and leaves paralog clusters and rare duplications untouched; `refine` then POA-realigns the bubble interiors to remove graph-builder artifacts while holding those self-loops fixed, so it is copy-number-safe and the substrate `call` sees is the same folding with cleaner interiors. Calling directly on the `bubble` graph also works for a non-tandem paralog cluster and gives identical results there, since `panphorte` leaves such clusters untouched, but it misses the tandem folding.

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

### Input contract

Checked before anything is read, and refused with the offender named rather than worked around. Each of these used to exit 0 and write a header-only VCF, which is indistinguishable from a locus with no variation:

| refused | why it cannot be tolerated |
|---|---|
| a path step naming a node with no `S` line | spelling skips it, so every sequence, coordinate and `SVLEN` from that walk is silently short |
| a duplicate path name | which haplotype a genotype column refers to would be undefined |
| a non-zero `L` overlap, or `*` | `call` spells by concatenating whole nodes, so an overlap is double-counted |
| a bubbles CSV naming nodes absent from the graph | the CSV belongs to a different graph; every bubble would be skipped |
| a duplicate `bubble_id` | the same site would be called twice |
| `--bubble-id` naming an id not in the CSV | a typo, not an empty result |
| an output path that is also an input | the run would overwrite the data it is reading |

Outputs are staged and renamed in only once the whole run succeeds, so a failure part-way through no longer leaves a complete-looking region VCF beside a non-zero exit. A narrowed rerun (`--bubble-id`) also clears per-bubble VCFs left by an earlier wider run at the same prefix, which would otherwise be indistinguishable from current output.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--min-sv-bp <N>` | minimum size of a reported (merged) event | `50` |
| `--rescue-min-bp <N>` | floor for sub-threshold events kept for merge/rescue | `min-sv-bp/2` |
| `--merge-distance-bp <N>` | coalesce nearby same-type events (reference or haplotype space) | `100` |
| `--merge-jaccard <X>` | node-set [Jaccard](../algorithms/call.md#merge-keys--jaccard-vs-sequence-identity) to merge events | `0.80` |
| `--merge-seq-identity <X>` | event-sequence identity to merge | `0.80` |
| `--merge-size-ratio <X>` | length-ratio floor for the sequence merge (lower to merge more differing lengths) | `0` (off) |
| `--min-haplotypes <N>` / `--min-alt-af <X>` | drop records below N carriers, or below that carrier frequency | `1` / `0` |
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
| `REF_CN` | `DUP` only: reference copy number |
| `CN_METHOD` / `CN_SCOPE` | `DUP` only: which of the three CN routes measured this — `REP` (exact traversal count of a panphorte self-loop), `MODULE_BP` (folded-node bp over a reference-calibrated unit), `PEAK` (highest interior-node multiplicity) — and whether a copy is a copy of a `REPEAT_UNIT` or of a `COLLAPSED_MODULE` |
| `RU_LEN` | `CN_METHOD=REP` only: the literal repeat unit, where `(CN − REF_CN) × RU_LEN` is the haplotype's size |
| `CN_UNIT_BP` | `CN_METHOD=MODULE_BP` only: the calibrated unit CN was divided by. This is the shared per-copy content, not a whole copy, so it understates a carrier's true gain or loss |
| `CN_SHARED_BP` / `CN_REF_FOLD` / `CN_MODULE_REF_BP` | the calibration inputs (`CN_UNIT_BP = CN_SHARED_BP / CN_REF_FOLD`) and the reference's total bp across the module |
| `GENES` | `--gtf`: gene(s) overlapped (whole folded module for a `DUP`) |
| `NALLELES` | `--multiallelic-loci`: number of alleles |
| `INSSEQ`/`DELSEQ`/`INVSEQ` | event sequence (omitted when very long) |


## Limitations

- A merged record carries one representative sequence for every carrier. `SVLEN_RANGE` and `MERGE_DIAMETER` show how far that reaches, but the per-carrier allele is only in the allele VCF.
- `CN_METHOD=PEAK` takes copy number from the highest interior-node multiplicity, which is a heuristic; those records carry `CN_CONFIDENCE=HEURISTIC`.
- Under `CN_METHOD=MODULE_BP`, `CN_UNIT_BP` is a calibration constant rather than a repeat-unit length, so `(CN − REF_CN) × CN_UNIT_BP` understates a carrier. `FORMAT:CNBP` is the per-haplotype size on every route.
- Where a module boundary is visited more than once, the span used is the widest one, from the first source occurrence to the last sink. It can enclose sequence lying between separate visits, and those records carry `CN_SPAN_AMBIGUOUS`.
- The event sequence is omitted from very large records, so a consumer of `INSSEQ`/`DELSEQ`/`INVSEQ` alone cannot reconstruct those; the allele VCF carries them.

## Lossless companion: the allele VCF

The region VCF is the interpreted output: records are merged so a reader can see what varies. `--allele-vcf` additionally writes `<prefix>.alleles.vcf`, one record per bubble whose ALT column carries every distinct allele at that site spelled out, with each haplotype's `GT` indexing its own. Nothing is merged and no sequence is summarized, so reconstructing a haplotype from it reproduces the input exactly. Use the region VCF to interpret a locus and the allele VCF when the exact sequence a given haplotype carries matters.

## Plotting

`scripts/plot_vcf_map.R` draws an oncoprint-style map of the region VCF, rows being haplotypes and columns variants grouped by bubble, with duplications shaded by copy number. `scripts/build_variant_node_data.R` and `scripts/variant_node_heatmap_app.R` build and serve an interactive per-node coverage view of the same calls. Each script documents its own options under `--help`.

## Gene annotation

With `--gtf` and a PanSN (Pangenome Sequence Naming) reference path, genes overlapping each record are named in `INFO/GENES`, `<prefix>.node_genes.tsv` maps every node to the genes covering it, and for a duplication whose paralogs are separable `<prefix>.dup_gene_cn.tsv` splits the module's copy number per gene.

## Example

See the [walkthrough](../walkthrough.md) for this module in a full end-to-end run.
