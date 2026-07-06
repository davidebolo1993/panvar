# Module `call`

CLI: `panvar call`

## What it does

Types structural variants on the pangenome graph into a multi-sample VCF. For each bubble it:
- diffs every haplotype's `source → sink` walk against a reference walk and types each difference (DEL/INS/INV/DUP)
- reduces fragmentation by coalescing events within a haplotype and clustering equivalent events across haplotypes
- recovers sub-threshold carriers by force-calling them against the merged record
- reads copy number (CN) directly off the walk (`--cn`)

Algorithm and worked trace: [algorithms/call.md](../algorithms/call.md).

Call on the `panphorte` graph (the `bubble` → `panphorte` → `call` path). `panphorte` folds genuine tandem repeats into `REP` self-loops — which a tandem needs for an exact count — and leaves paralog clusters and rare duplications untouched. Calling on the `bubble` graph also works for a non-tandem paralog cluster and should give identical results. 

#### Event types

| type | meaning |
|------|---------|
| DEL | reference-only nodes (deleted from the haplotype) |
| INS | haplotype-only nodes (inserted vs reference); `--classify-ins` adds `INS_SUBTYPE=NOVEL\|DUP` |
| INV | a haplotype run that is the reverse-complement node-walk of a reference run |
| DUP | a copy-number gain/loss (`--cn`) |

## Required inputs

- `-i, --gfa <graph.gfa>` — the call substrate from `panphorte` (or `bubble`)
- one of `-b, --bubble-prefix-in <prefix>` or `-c, --bubbles-csv <path>`.
- `-r, --reference-path <name>` — the diff baseline (full name or unique case-insensitive substring).
- `-o, --out-prefix <prefix>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--min-sv-bp <N>` | minimum size of a reported (merged) event | `50` |
| `--rescue-min-bp <N>` | floor for sub-threshold events kept for merge/[rescue](../algorithms/call.md#worked-trace) | `min-sv-bp/2` |
| `--merge-distance-bp <N>` | coalesce nearby same-type events (reference or haplotype space) | `100` |
| `--merge-jaccard <X>` | node-set [Jaccard](../algorithms/call.md#merge-keys--jaccard-vs-sequence-identity) to merge events | `0.80` |
| `--merge-seq-identity <X>` | event-sequence identity to merge | `0.80` |
| `--merge-size-ratio <X>` | length-ratio floor for the sequence merge (lower to merge more differing lengths) | `0` (off) |
| `--min-haplotypes <N>` / `--min-maf <X>` | drop records below N carriers / carrier frequency `AF=AC/AN` | `1` / `0` |
| `--cn` | emit the total copy number of the folded module; with `--gtf` it is reassigned to per-gene CN where the paralogs are separable. See [algorithms/call.md](../algorithms/call.md#copy-number) for the routes and [calling without `--cn`](../algorithms/call.md#copy-number) for what surfaces without it | off |
| `--classify-ins` | refine INS subtype NOVEL/DUP via minimap2 | off |
| `--multiallelic-loci` | collapse a bounded locus into one multiallelic record; `--multiallelic-max-bp` (5000) bounds it | off |
| `--gtf <path>` | gene annotation (needs PanSN `--reference-path`); see [below](#gene-annotation---gtf) | — |
| `--bubble-id <N>` / `--no-per-bubble-vcf` / `--no-variant-paths` / `-q, --quiet` | scope and output toggles | — |


## Outputs

| file | contents |
|------|----------|
| `<prefix>.region.vcf` | all records, coordinate-sorted, unique IDs (`bgzip`+`tabix`-able) |
| `<prefix>.bubble_<id>.vcf` | per-bubble VCF (unless `--no-per-bubble-vcf`) |
| `<prefix>.variant_paths.tsv` | per (variant, carrier) sub-walk provenance (unless `--no-variant-paths`) |
| `<prefix>.variant_nodes.tsv` | per-variant node set — the `describe --variant-nodes` handoff |
| `<prefix>.node_genes.tsv`, `<prefix>.dup_gene_cn.tsv` | with `--gtf` (see below) |

VCF 4.2; samples = haplotypes. `FORMAT`: `GT` (`1` carrier / `0` ref-like / `.` doesn't traverse), `CN` (per-sample copy number on DUP). Key `INFO`:

| field | meaning |
|-------|---------|
| `SVTYPE`, `SVLEN` | event type and length difference (ALT − REF) |
| `END`, `BUBBLE_ID`, `START_NODE`, `END_NODE`, `EVENT_NODES` | locus span + graph anchors + the variant's node set |
| `AN`/`AC`/`AF` | haplotypes traversing the bubble/carriers/carrier freq `AC/AN` (ALT freq, not folded) |
| `NMERGED`, `SVLEN_RANGE`, `MERGE_JACCARD`/`MERGE_SEQID`/`MERGE_SIZE_RATIO` | merge count, size span, merge evidence (≥2-event records only) |
| `EVENTID` | links a co-located DEL+INS substitution |
| `INS_SUBTYPE` | `NOVEL`/`DUP` (`--classify-ins`) |
| `REF_CN` / `RU_LEN` | DUP only: reference copy number / one-copy repeat-unit length |
| `GENES` | `--gtf`: gene(s) overlapped (whole folded module for a DUP) |
| `NALLELES` | `--multiallelic-loci`: number of alleles |
| `INSSEQ`/`DELSEQ`/`INVSEQ` | event sequence (omitted when very long) |

## Plotting

`scripts/plot_vcf_map.R` draws the headline oncoprint-style map (rows = haplotypes, columns = variants grouped by bubble; DEL red, INS-NOVEL green, INS-DUP purple, INV orange, multiallelic yellow-amber, DUP shaded blue by `FORMAT:CN`). It needs only the VCF. Script flags (need `Rscript` + `ggplot2`):

- `--vcf <file>` — the `call` region VCF (required).
- `--out <prefix>` — output prefix (required).
- `--clusters` / `--cluster-by <tsv>` — plot only cluster representatives / group and order rows by cluster (mirrors the inspect heatmaps).
- `--max-paths <N>` — cap the number of haplotype rows.
- `--flip` — transpose (variants as rows).
- `--scale` / `--scale-transform <raw|sqrt|log1p>` — size each cell by `|SVLEN|` (`RU_LEN` for a DUP), with an optional transform.
- `--reference-path <name>` — order columns by reference coordinate.
- `--title <text>`, `--width` / `--height` / `--dpi` — labelling and output resolution.


## Gene annotation (`--gtf`)

A reference-coordinate GTF (Ensembl/GENCODE) projected onto the graph via a PanSN reference (`sample#hap#chrom:start-end`); else skipped with a warning. It adds: `INFO=GENES` per record; `<prefix>.node_genes.tsv` (`node_id → gene(s)`, the bridge for `associate --node-genes`); and `<prefix>.dup_gene_cn.tsv` — the per-haplotype per-gene copy number split of a folded paralog cluster, resolved **from the GTF alone by k-mer dosage** (no per-haplotype realignment). Each gene's coding sequence (its gene span if it has no CDS) supplies a set of k-mers private to it vs its paralogs; a haplotype's per-copy count of those k-mers is the gene's copy number. Near-identical paralog pairs — where a plain count blurs under gene conversion — are split instead by per-site consensus over their divergent coding columns. Each row carries the evidence (`hits`, `priv_kmers`, `dosage`) so the call is auditable, plus a `reliable` flag (`0` only for a gene with no private k-mers, reported as the module total). See [algorithms/call.md](../algorithms/call.md#gene-annotation---gtf).

## Example

```bash
./build/panvar call \
  -i results/real_data/lpa/panphorte/panphorte.normalized.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/panphorte/panphorte \
  --reference-path "grch38#1" \
  -o results/real_data/lpa/call/call \
  --cn \
  --gtf tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz
```
