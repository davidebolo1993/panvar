# Module `call`

CLI: `panvar call`

## What it does

Types structural variants on the pangenome graph into a multi-sample VCF (Variant Call Format). For each bubble it:
- diffs every haplotype's `source → sink` walk against a reference walk and types each difference (`DEL`/`INS`/`INV`/`DUP`, see the [table](#event-types) below)
- reduces fragmentation by coalescing events within a haplotype and clustering equivalent events across haplotypes
- recovers sub-threshold carriers by force-calling them against the merged record
- reads exact tandem-repeat copy number directly from self-loop traversal counts, and optionally infers copy number for other repeated modules (`--cn`)

Algorithm and worked trace: [algorithms/call.md](../algorithms/call.md).

The recommended path is `bubble → panphorte → [refine] → call`, so `call` runs on the `panphorte` graph, or on the `refine` graph when the opt-in step is used. `panphorte` folds supported tandem repeats into repeat-unit self-loops, which provide exact traversal counts; `refine` can then re-align residual interiors while holding those loops fixed. `call` can also run directly on a `bubble` graph, but unfolded tandems lack the explicit repeat-unit representation.

#### Event types

| type | meaning |
|------|---------|
| `DEL` | reference-only nodes (deleted from the haplotype) |
| `INS` | haplotype-only nodes (inserted vs reference); `--classify-ins` adds `INS_SUBTYPE=NOVEL\|DUP` |
| `INV` | a haplotype run that is the reverse-complement node-walk of a reference run |
| `DUP` | a copy-number gain/loss; exact self-loop counts are always active, while `--cn` enables the inferred module routes |

## Required inputs

- `-i, --gfa <graph.gfa>` — the call substrate from `panphorte`/`refine`/`bubble` (in any case, node ids should match the CSV to `-b`)
- one of `-b, --bubble-prefix-in <prefix>` or `-c, --bubbles-csv-in <path>`
- `-r, --reference-path <name>` — the diff baseline (full name or unique case-insensitive substring).
- `-o, --out-prefix <prefix>`.

The graph and the bubbles CSV must describe each other, and both are checked before anything is read. A path step naming a node with no segment, a duplicate path name, a non-zero or unspecified link overlap, a bubbles CSV naming nodes the graph lacks, a duplicate bubble id, a `--bubble-id` that is not in the CSV, or an output path that is also an input: each is refused with the offender named. Every one of them would otherwise produce a plausible-looking VCF rather than an error.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--min-sv-bp <N>` | minimum size of a reported (merged) event | `50` |
| `--rescue-min-bp <N>` | floor for sub-threshold events kept for merge/rescue | `min-sv-bp/2` |
| `--merge-distance-bp <N>` | coalesce nearby same-type events (reference or haplotype space) | `100` |
| `--merge-jaccard <X>` | node-set [Jaccard](../algorithms/call.md#5-merge-across-haplotypes) to merge events | `0.80` |
| `--merge-seq-identity <X>` | event-sequence identity to merge | `0.80` |
| `--merge-size-ratio <X>` | length-ratio floor for the sequence merge (lower to merge more differing lengths) | `0` (off) |
| `--min-haplotypes <N>` / `--min-alt-af <X>` | drop records below N carriers or below ALT carrier frequency; `--min-alt-af` is not minor-allele frequency | `1` / `0` |
| `--cn` | enable the inferred `MODULE_BP` and `PEAK` copy-number routes; exact `REP` self-loop counting is always active. With `--gtf`, separable paralogs also receive per-gene CN | off |
| `--allele-vcf` / `--allele-vcf-max-bp <N>` | write one explicit multiallelic record per bubble; optionally skip bubbles whose alleles exceed N bp | off / `0` (unlimited) |
| `--classify-ins` | refine `INS` subtype `NOVEL`/`DUP` via minimap2 | off |
| `--multiallelic-loci` | collapse a bounded locus into one multiallelic record; `--multiallelic-max-bp` (5000) bounds it | off |
| `--gtf <path>` | gene annotation (needs a PanSN (Pangenome Sequence Naming) `--reference-path`); see [the algorithm](../algorithms/call.md#gene-annotation) | — |
| `--bubble-id <N>` / `--no-per-bubble-vcf` / `--no-variant-nodes` / `-q, --quiet` | scope and output toggles | — |


## Outputs

| file | contents |
|------|----------|
| `<prefix>.region.vcf` | all records, coordinate-sorted, unique IDs (`bgzip`+`tabix`-able). This is the interpreted output: records are merged so a reader can see what varies |
| `<prefix>.alleles.vcf` | with `--allele-vcf`: one record per emitted bubble whose ALT column spells every distinct allele at that site, each haplotype's `GT` indexing its own. Nothing is merged or summarized, so emitted alleles reproduce their graph walks exactly |
| `<prefix>.bubble_<id>.vcf` | per-bubble VCF (unless `--no-per-bubble-vcf`) |
| `<prefix>.variant_nodes.tsv` | per-variant node set — the `describe --variant-nodes` handoff and the `benchmark` input (unless `--no-variant-nodes`) |
| `<prefix>.node_genes.tsv`, `<prefix>.dup_gene_cn.tsv` | with `--gtf` (see below) |

VCF 4.2; samples are haplotypes. Core `FORMAT` fields are `GT` (`1` carrier / `0` reference-like / `.` does not traverse), `CN` (per-sample copy number on `DUP`) and `CNBP` (`DUP` only: the actual signed base-pair change across the module relative to the reference). `REP` records also carry `CNRESID`; `MODULE_BP` records carry the unrounded dosage and rounding margin in `CNR_RAW` and `CNR_MARGIN`. Key `INFO` fields are:

| field | meaning |
|-------|---------|
| `SVTYPE`, `SVLEN` | event type and length difference (ALT − REF) |
| `END`, `BUBBLE_ID`, `START_NODE`, `END_NODE`, `EVENT_NODES` | locus span + graph anchors + the variant's node set |
| `IMPRECISE` | the affected interval is wider than the event bases because within-haplotype coalescing joined separated pieces |
| `AN`/`AC`/`AF` | haplotypes traversing the bubble/carriers/carrier freq `AC/AN` (ALT freq, not folded) |
| `NMERGED`, `SVLEN_RANGE`, `MERGE_JACCARD`/`MERGE_SEQID`/`MERGE_SIZE_RATIO`, `MERGE_DIAMETER` | merge count, size range and evidence; the diameter exposes long single-linkage chains |
| `EVENTID` | links a co-located `DEL`+`INS` substitution |
| `INS_SUBTYPE` | `NOVEL`/`DUP` (`--classify-ins`) |
| `REF_CN` | `DUP` only: reference copy number |
| `CN_METHOD` / `CN_SCOPE` | `DUP` only: which of the three CN routes measured this — `REP` (exact traversal count of a panphorte self-loop), `MODULE_BP` (folded-node bp over a reference-calibrated unit), `PEAK` (highest interior-node multiplicity) — and whether a copy is a copy of a `REPEAT_UNIT` or of a `COLLAPSED_MODULE` |
| `RU_LEN` | `CN_METHOD=REP` only: the literal repeat unit, where `(CN − REF_CN) × RU_LEN` is the haplotype's size |
| `CN_UNIT_BP` | `CN_METHOD=MODULE_BP` only: the calibrated unit CN was divided by. This is the shared per-copy content, not a whole copy, so it understates a carrier's true gain or loss |
| `CN_SHARED_BP` / `CN_REF_FOLD` / `CN_MODULE_REF_BP` | the calibration inputs (`CN_UNIT_BP = CN_SHARED_BP / CN_REF_FOLD`) and the reference's total bp across the module |
| `CN_DOSAGE_MODEL`, `CN_ROUND_AMBIGUOUS_FRAC`, `CN_STEP_*` | `MODULE_BP` model and diagnostics; the VCF header defines each field and records whether panel spacing was used |
| `GENES` | `--gtf`: gene(s) overlapped (whole folded module for a `DUP`) |
| `NALLELES` | `--multiallelic-loci`: number of alleles |
| `INSSEQ`/`DELSEQ`/`INVSEQ` | event sequence (omitted when very long) |


## Limitations

- A merged record carries one representative sequence for every carrier. `SVLEN_RANGE` and `MERGE_DIAMETER` show how far that reaches, but the per-carrier allele is only in the allele VCF.
- `CN_METHOD=PEAK` takes copy number from the highest interior-node multiplicity, which is a heuristic; those records carry `CN_CONFIDENCE=HEURISTIC`.
- Under `CN_METHOD=MODULE_BP`, `CN_UNIT_BP` is a calibration constant rather than a repeat-unit length, so `(CN − REF_CN) × CN_UNIT_BP` understates a carrier. `FORMAT:CNBP` is the per-haplotype size on every route.
- Where a module boundary is visited more than once, the span used is the widest one, from the first source occurrence to the last sink. It can enclose sequence lying between separate visits, and those records carry `CN_SPAN_AMBIGUOUS`.
- The event sequence is omitted from very large records, so a consumer of `INSSEQ`/`DELSEQ`/`INVSEQ` alone cannot reconstruct those; the allele VCF carries them.

With `--gtf` and a PanSN (Pangenome Sequence Naming) reference path, genes overlapping each record are named in `INFO/GENES`, `<prefix>.node_genes.tsv` maps every node to the genes covering it, and for a duplication whose paralogs are separable `<prefix>.dup_gene_cn.tsv` splits the module's copy number per gene.

`scripts/plot_vcf_map.R` renders the region VCF as an oncoprint-style map, and `scripts/build_variant_node_data.R` with `scripts/variant_node_heatmap_app.R` serve an interactive per-node view of the same calls. Each documents its own options under `--help`.

## Example

See the [walkthrough](../walkthrough.md) for this module in a full end-to-end run.
