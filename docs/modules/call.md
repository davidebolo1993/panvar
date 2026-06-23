# Call Module (Module 3 — graph-native SV calling)

CLI: `panvar call`

## What it does

Types structural variants on the pangenome graph into a multi-sample VCF. For each bubble it diffs every
haplotype's `source → sink` walk against a reference walk and types each difference (DEL/INS/INV/DUP). It
then reduces fragmentation by coalescing events within a haplotype and clustering equivalent events
across haplotypes, and recovers sub-threshold carriers. Copy number is read directly off the walk.

Algorithm, copy-number mechanics, merge keys, and worked traces:
[algorithms/call.md](../algorithms/call.md).

Which graph to call on depends on locus topology — tandem repeats on the `panphorte` graph, PGGB-folded
paralog clusters on the `bubble` graph. See the
[CN-topology table](../algorithms/call.md#copy-number-one-method-per-locus-topology).

### Event types

| type | meaning |
|------|---------|
| DEL | reference-only nodes (deleted from the haplotype) |
| INS | haplotype-only nodes (inserted vs reference); `--classify-ins` adds `INS_SUBTYPE=NOVEL|DUP` |
| INV | a haplotype run that is the reverse-complement node-walk of a reference run |
| DUP | a copy-number gain/loss; three sources by precedence (self-loop `REP`, `--cn-from-coverage`, `--cn-from-multiplicity`) — see [algorithms/call.md](../algorithms/call.md#copy-number--three-ways) |

## Required inputs

- `-i, --gfa <graph.gfa>` — the call substrate (`panphorte` `.normalized.sorted.gfa` for tandems, else the
  `bubble` `.sorted.gfa`).
- one of `-b, --bubble-prefix-in <prefix>` or `-c, --bubbles-csv <path>`.
- `-r, --reference-path <name>` — the diff baseline (full name or unique case-insensitive substring).
- `-o, --out-prefix <prefix>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--min-sv-bp <N>` | minimum size of a reported (merged) event | `50` |
| `--rescue-min-bp <N>` | floor for sub-threshold events kept for merge/[rescue](../algorithms/call.md#worked-trace--diff-coalesce-cluster-rescue) | `min-sv-bp/2` |
| `--merge-distance-bp <N>` | coalesce nearby same-type events (reference or haplotype space) | `100` |
| `--merge-jaccard <X>` | node-set [Jaccard](../algorithms/call.md#merge-keys--jaccard-vs-sequence-identity) to merge events | `0.80` |
| `--merge-seq-identity <X>` | event-sequence identity to merge | `0.80` |
| `--merge-size-ratio <X>` | length-ratio floor for the sequence merge (lower to merge differing STR lengths) | `0` (off) |
| `--min-haplotypes <N>` / `--min-maf <X>` | drop records below N carriers / carrier frequency `AF=AC/AN` | `1` / `0` |
| `--cn-from-coverage` | total-module CN on folded paralog clusters the reference traverses ≥2× | off |
| `--cn-from-multiplicity` | `DUP` from peak node multiplicity for folded bubbles with no self-loop | off |
| `--classify-ins` | refine INS subtype NOVEL/DUP via minimap2 | off |
| `--multiallelic-loci` | collapse a bounded locus into one multiallelic record ([mechanics](../algorithms/call.md#multiallelic-mechanics---multiallelic-loci)); `--multiallelic-max-bp` (5000) bounds it | off |
| `--gtf <path>` | gene annotation (needs PanSN `--reference-path`); see [below](#gene-annotation---gtf) | — |
| `--bubble-id <N>` / `--no-per-bubble-vcf` / `--no-variant-paths` / `-q, --quiet` | scope & output toggles | — |

Both CN flags compose safely (disjoint topologies); pass both for widest recall.

## Outputs

| file | contents |
|------|----------|
| `<prefix>.region.vcf` | all records, coordinate-sorted, unique IDs (`bgzip`+`tabix`-able) |
| `<prefix>.bubble_<id>.vcf` | per-bubble VCF (unless `--no-per-bubble-vcf`) |
| `<prefix>.variant_paths.tsv` | per (variant, carrier) sub-walk provenance (unless `--no-variant-paths`) |
| `<prefix>.variant_nodes.tsv` | per-variant node set — the `describe --variant-nodes` handoff |
| `<prefix>.node_genes.tsv`, `<prefix>.dup_gene_cn.tsv` | with `--gtf` (see below) |

VCF 4.2; samples = haplotypes. `FORMAT`: `GT` (`1` carrier / `0` ref-like / `.` doesn't traverse),
`CN` (per-sample copy number on DUP). Key `INFO`:

| field | meaning |
|-------|---------|
| `SVTYPE`, `SVLEN` | event type and length difference (ALT − REF) |
| `END`, `BUBBLE_ID`, `START_NODE`, `END_NODE`, `EVENT_NODES` | locus span + graph anchors + the variant's node set |
| `AN`/`AC`/`AF` | haplotypes traversing the bubble / carriers / carrier freq `AC/AN` (ALT freq, not folded) |
| `NMERGED`, `SVLEN_RANGE`, `MERGE_JACCARD`/`MERGE_SEQID`/`MERGE_SIZE_RATIO` | merge count, size span, merge evidence (≥2-event records only) |
| `EVENTID` | links a co-located DEL+INS substitution |
| `INS_SUBTYPE` | `NOVEL`/`DUP` (`--classify-ins`) |
| `REF_CN` / `RU_LEN` | DUP only: reference copy number / one-copy repeat-unit length |
| `GENES` | `--gtf`: gene(s) overlapped (whole folded module for a DUP) |
| `NALLELES` | `--multiallelic-loci`: number of alleles |
| `INSSEQ`/`DELSEQ`/`INVSEQ` | event sequence (omitted when very long) |

## Plotting

`scripts/plot_vcf_map.R` draws the headline oncoprint-style map (rows = haplotypes, columns = variants
grouped by bubble; DEL red, INS-NOVEL green, INS-DUP purple, INV orange, multiallelic yellow-amber, DUP
shaded blue by `FORMAT:CN`). Needs only the VCF. Example in the per-gene scripts (`scripts/genes/`).

Script flags (needs `Rscript` + `ggplot2`):

- `--vcf <file>` — the `call` region VCF (required).
- `--out <prefix>` — output prefix (required).
- `--clusters` / `--cluster-by <tsv>` — plot only cluster representatives / group & order rows by cluster (mirrors the inspect heatmaps).
- `--max-paths <N>` — cap the number of haplotype rows.
- `--flip` — transpose (variants as rows).
- `--scale` / `--scale-transform <raw|sqrt|log1p>` — size each cell by `|SVLEN|` (`RU_LEN` for a DUP), with an optional transform.
- `--reference-path <name>` — order columns by reference coordinate.
- `--title <text>`, `--width`/`--height`/`--dpi` — labelling and output resolution.

Copy-number concordance against ground truth (used by `scripts/regen_results.sh`) is plotted by
`scripts/plot_cn_correlation.R` (`--table <cn_table.tsv>`, `--out <prefix>`, optional `--width`/`--height`/`--dpi`),
which writes `<prefix>.loci.png` (the four loci as total counts) and `<prefix>.genes.png` (the resolved
CYP2D6/CYP2D7 split).

## Gene annotation (`--gtf`)

A reference-coordinate GTF (Ensembl/GENCODE) projected onto the graph via a PanSN reference
(`sample#hap#chrom:start-end`); else skipped with a warning. Adds: `INFO=GENES` per record;
`<prefix>.node_genes.tsv` (`node_id → gene(s)`, the bridge for `associate --node-genes`);
`<prefix>.dup_gene_cn.tsv` (per-haplotype per-gene CN, `reliable` flag — paralogs that competitive
realignment can separate vs near-identical ones reported as a collapsed module total). Mechanics + trace:
[algorithms/call.md](../algorithms/call.md#gene-annotation-trace---gtf).

## Example

Matches `scripts/genes/lpa.sh` (LPA = tandem → panphorte graph + `--cn-from-multiplicity`; other genes use
`--cn-from-coverage` on the `bubble` graph — see `scripts/genes/`):

```bash
./build/panvar call \
  -i results/real_data/lpa/panphorte/panphorte.normalized.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/panphorte/panphorte \
  --reference-path GRCh38 -o results/real_data/lpa/call/call \
  --cn-from-multiplicity --gtf tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz
```

Algorithm & worked examples: [algorithms/call.md](../algorithms/call.md). References:
[references.md](../references.md#call).
