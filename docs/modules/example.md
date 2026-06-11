# Worked End-to-End Example

This page walks the full `panvar` pipeline on one locus: bubbles → tandem normalization →
graph-native variant calling.

## 0. Inputs

- a GFA with P/W paths (`graph.gfa`), one of which is a usable reference (e.g.
  `GRCh38#0#chr6:31891045-32123783`)
- precomputed snarls: `vg snarls -R -j ... > graph.snarls.jsonl`

## 1. Call bubbles (Module 1)

```bash
panvar bubble -i graph.gfa -o out/bubble --snarls-in graph.snarls.jsonl
```

Produces `out/bubble.bubbles.csv` (one row per site: source, sink, inside nodes).

## 2. (Optional) Inspect a site

```bash
panvar inspect -i graph.gfa --bubble-prefix-in out/bubble --bubble-id 1 -o out/inspect/bubble_1
```

Writes a path FASTA plus node- and edge-traversal matrices; a tandem self-loop shows as an edge cell
`> 1`. Heatmaps: `scripts/plot_node_coverage_heatmap.R`, `scripts/plot_edge_coverage_heatmap.R`.

## 3. Normalize tandem repeats (panphorte)

```bash
panvar panphorte -i graph.gfa --bubble-prefix-in out/bubble -o out/panphorte --min-similarity 0.90
# re-snarl + re-bubble the normalized graph so the caller sees the REP self-loops
vg snarls -R -j out/panphorte.normalized.gfa > out/panphorte.normalized.snarls.jsonl
panvar bubble -i out/panphorte.normalized.gfa -o out/bubble2 \
  --snarls-in out/panphorte.normalized.snarls.jsonl
```

A tandem array collapses to one repeat-unit (`REP`) node looped N times, so copy number is explicit.

## 4. Call structural variants (Module 2)

```bash
panvar call \
  -i out/panphorte.normalized.gfa \
  --bubble-prefix-in out/bubble2 \
  --reference-path "GRCh38#0#chr6:31891045-32123783" \
  -o out/call \
  --merge-distance-bp 100 --merge-jaccard 0.80 --classify-ins
```

Each haplotype's bubble walk is diffed against the reference walk into DEL/INS/INV/DUP events;
fragmented same-type events are coalesced within a bubble (`--merge-distance-bp`) and equivalent events
are merged across haplotypes by node-set Jaccard (`--merge-jaccard`).

## 5. Outputs

- `out/call.bubble_<id>.vcf` — one multi-sample VCF per bubble
- `out/call.region.vcf` — concatenated

Samples are the haplotypes (haploid `GT`: `1` carrier / `0` reference-like / `.` not traversing). DUP
records carry `REF_CN` and per-sample `CN`. `INFO` reports `START_NODE`/`END_NODE` and the
`EVENT_NODES` set. Parse/sort with `bcftools`:

```bash
bcftools sort out/call.region.vcf -Oz -o out/call.region.sorted.vcf.gz
bcftools view out/call.region.sorted.vcf.gz | head
```
