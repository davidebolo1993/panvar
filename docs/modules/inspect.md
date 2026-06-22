# Inspect Utility

Date: 2026-06-07

CLI entrypoint:

- `panvar inspect`

## What it does

`inspect` is a small utility for checking one called bubble, or all called bubbles, before moving downstream.

Given a GFA and a module-1 bubble CSV, it writes:

1. a compressed multi-FASTA with one record per path crossing that bubble
2. a TSV matrix describing how each path traverses the bubble-internal **nodes**
3. a TSV matrix describing how each path traverses the bubble-internal **edges** (adjacencies)

This is meant for sanity-checking questions like:

- which paths actually cross this bubble?
- how long is each path through this bubble?
- which internal nodes and edges are reused, and in which orientation?

## Required inputs

- `--gfa <graph.gfa>`
- one of:
  - `--bubble-prefix-in <module1-prefix>` (auto uses `<module1-prefix>.bubbles.csv`)
  - `--bubbles-csv <module1.bubbles.csv>`

Optional:

- `--bubble-id <N>`: restrict to one bubble

If `--bubble-id` is omitted, `inspect` writes outputs for every bubble in the input CSV.

## Key options

```text
panvar inspect -i <graph.gfa> (-b <prefix> | -c <bubbles.csv>) [--bubble-id <N>] [options]
```

- `--cluster` — group paths by source→sink walk and write `<prefix>.bubble_<N>.clusters.tsv`
- `--cluster-similarity <f>` — walk-similarity threshold for `--cluster` (default `0.90`)

## Outputs

Default output paths are derived from `--out-prefix` and the bubble ID:

- `<out-prefix>.bubble_<N>.paths.fa.gz`
- `<out-prefix>.bubble_<N>.node_counts.tsv`
- `<out-prefix>.bubble_<N>.edge_counts.tsv`
- `<out-prefix>.bubble_<N>.node_lengths.tsv` — per-node bp lengths
- `<out-prefix>.bubble_<N>.clusters.tsv` — only when `--cluster` is set

You can override the first three explicitly with:

- `--fasta-out <path>`
- `--table-out <path>` (node counts)
- `--edge-table-out <path>` (edge counts)

Explicit output paths are only accepted with `--bubble-id`, because all-bubble mode needs one output set per bubble. When inspecting all bubbles, a progress bar is printed to stderr; the machine-readable summary stays on stdout. With a single `--bubble-id`, the per-bubble detail is printed instead.

## Node-length sidecar

`<out-prefix>.bubble_<N>.node_lengths.tsv` lists, in the **same order as the `node.*` columns** of the node-count table, each internal node and its bp length:

```text
node_id	length_bp
```

It exists so the node coverage heatmap can scale its x-axis by node length (many bubble nodes are SNP-sized) without having to reorder columns. It is always written.

## FASTA output

Each FASTA record is the spelled source-to-sink bubble interval for one path.

Headers include:

- path name
- bubble id
- source/sink node ids
- interval length in bp
- whether the original best interval was source-to-sink or sink-to-source
- path-step interval coordinates

A bubble has a fixed `source` and `sink`, but an individual path can walk *through* it in either direction — some paths enter at the `source` and leave at the `sink`, others run the opposite way (enter at the `sink`, leave at the `source`). `inspect` searches both directions and records which one matched in the `source_to_sink=` header field. When a path crossed **sink-to-source**, `inspect` reorients that path's *node walk* — it reverses the step order and flips each node's strand — so the emitted allele always reads `source → sink`. Because the spelled DNA follows the node strands, flipping them yields the reverse-complement of that path's raw sequence. This is driven purely by traversal direction, not by sequence content.

## Node-count table

Columns:

- `path_name`
- `path_length_bp`
- one `node.<id>` column per bubble-internal node from `inside_nodes`

Each node cell is:

```text
total:forward:reverse
```

For example, `3:2:1` means that, inside this bubble interval only, the path traverses that node 3 times total: 2 times in forward orientation and 1 time in reverse orientation. Counts are **interval-local**. If the same path reuses a node elsewhere outside the selected source/sink interval, that outside traversal is not counted for this bubble. Orientation is counted after the canonical `source → sink` normalization above. A `reverse` count therefore means the node is traversed against the canonical frame in this bubble allele — typically a real local inversion — rather than just a path that happened to run the other way. If the graph has no
`-`/`<` path steps in the selected intervals, reverse counts can legitimately be zero.

## Edge-count table

Columns:

- `path_name`
- `path_length_bp`
- one `edge.<from><sign>><to><sign>` column per bubble-internal adjacency, e.g. `edge.1886+>1887+`

Each edge cell is an integer: how many times the path traverses that **directed, orientation-aware adjacency** inside the canonical source-to-sink bubble interval. The edge key encodes the orientation of both endpoints (`+`/`-`), so a forward and a reverse traversal of the same node pair are distinct columns. The edge matrix is **adjacency-aware**, complementing the node matrix: for instance, a tandem repeat that loops over the same unit repeats the *same edge*, so a self-loop / back-edge shows up as a cell value `> 1`. Edge columns are the union of adjacencies observed across all paths, sorted for a stable, deterministic column order. As with nodes, counts are interval-local.

## Path clustering (`--cluster`)

`--cluster` groups the paths crossing a bubble by their `source → sink` walk, so structurally identical haplotypes collapse to one representative. It is an inspection aid only — nothing downstream depends on it — but it is handy for plotting (see the heatmap `--clusters` option). The shingles, MinHash sketch and
connected-component grouping are traced on a tiny worked dataset in [algorithm_example.md](../algorithm_example.md). Identical walks collapse first. The distinct walks are then clustered by **connected components**: each distinct walk gets a **MinHash sketch** over oriented node-step shingles, and two walks' similarity is the sketch-estimated identity. A threshold graph connects every pair at least `--cluster-similarity` similar (default `0.90`), and clusters are its **connected components** — so the grouping is **transitive and order-independent** (if A~B and B~C, then A, B, C land in one cluster). For very short walks that cannot be shingled it falls back to an exact **bp-weighted Jaccard** (each step contributes its node's bp length to an oriented `<node_id><strand>` token multiset, compared as
`sum(min)/sum(max)`). The representative is the member minimizing max-then-mean intra-cluster distance (ties → most-supported walk, then signature).
The sketch is **multiplicity-aware**: a shingle seen *k* times contributes *k* distinct sketch elements, so the sketch Jaccard tracks the shingle **multiset**, not just the set: copies of the same unit share the same shingle *set* regardless of copy number, so a set-based sketch would merge all copy numbers into one cluster — the multiset sketch instead separates them by copy number, matching the bp-weighted fallback. 

Output `<out-prefix>.bubble_<N>.clusters.tsv`:

```text
cluster_id	n_paths	representative_path	members
```

`members` is a `;`-separated list of path names. Lower `--cluster-similarity` to merge more aggressively,
raise it to split.

## Example

Use the **sorted GFA from `bubble`** (its node ids match the bubble CSV).

```bash
./build/panvar inspect \
  -i results/real_data/lpa/bubble/bubble.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/bubble/bubble \
  --bubble-id 7 \
  -o results/real_data/lpa/inspect/inspect \
  --cluster \
  --cluster-similarity 0.97
```

This writes:

- `results/real_data/lpa/inspect/inspect.bubble_7.paths.fa.gz`
- `results/real_data/lpa/inspect/inspect.bubble_7.node_counts.tsv`
- `results/real_data/lpa/inspect/inspect.bubble_7.edge_counts.tsv`
- `results/real_data/lpa/inspect/inspect.bubble_7.clusters.tsv`

To inspect all bubbles, drop `--bubble-id`; `inspect` writes one FASTA/table set per bubble (`inspect.bubble_1.*`, `inspect.bubble_2.*`, …):

```bash
./build/panvar inspect \
  -i results/real_data/lpa/bubble/bubble.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/bubble/bubble \
  -o results/real_data/lpa/inspect/inspect
```

## Node Coverage Heatmap

Use `Rscript scripts/plot_node_coverage_heatmap.R` to visualize an `inspect` node-count table.

Dependency:

- `Rscript` with `ggplot2` 

Example:

```bash
Rscript scripts/plot_node_coverage_heatmap.R \
  --table results/real_data/lpa/inspect/inspect.bubble_7.node_counts.tsv \
  --out results/real_data/lpa/inspect/inspect.bubble_7.node_coverage
```

This writes:

- `results/real_data/lpa/inspect/inspect.bubble_7.node_coverage.png`
- `results/real_data/lpa/inspect/inspect.bubble_7.node_coverage.pdf`


To scale the x-axis by node length (thin tiles for SNP-sized nodes, wide tiles for long nodes) while keeping the column order, pass the node-length sidecar:

```bash
Rscript scripts/plot_node_coverage_heatmap.R \
  --table results/real_data/lpa/inspect/inspect.bubble_7.node_counts.tsv \
  --out results/real_data/lpa/inspect/inspect.bubble_7.node_coverage \
  --node-lengths results/real_data/lpa/inspect/inspect.bubble_7.node_lengths.tsv
```

The `--cluster` output (`clusters.tsv`) can drive the heatmap in two ways. To plot **only the cluster representatives** (one row per cluster), use `--clusters`:

```bash
Rscript scripts/plot_node_coverage_heatmap.R \
  --table results/real_data/lpa/inspect/inspect.bubble_7.node_counts.tsv \
  --out results/real_data/lpa/inspect/inspect.bubble_7.node_coverage
  --node-lengths results/real_data/lpa/inspect/inspect.bubble_7.node_lengths.tsv \
  --clusters results/real_data/lpa/inspect/inspect.bubble_7.clusters.tsv 
```

To **group/order the rows by cluster** (representative first, with a thin separator line between clusters), use `--cluster-by` instead:

```bash
Rscript scripts/plot_node_coverage_heatmap.R \
  --table results/real_data/lpa/inspect/inspect.bubble_7.node_counts.tsv \
  --out results/real_data/lpa/inspect/inspect.bubble_7.node_coverage \
  --node-lengths results/real_data/lpa/inspect/inspect.bubble_7.node_lengths.tsv \
  --clusters results/real_data/lpa/inspect/inspect.bubble_7.clusters.tsv \
  --cluster-by results/real_data/lpa/inspect/inspect.bubble_7.clusters.tsv
```

`--cluster-by` overrides the coverage-profile `--cluster-rows` ordering (rows follow the inspect clusters instead). Paths absent from the cluster file sort last.

Useful plotting options:

- `--value total|forward|reverse`
- `--transform raw|log1p`
- `--node-lengths <path>` (length-scale x; tile widths via `--length-transform raw|sqrt|log1p`, default `sqrt`)
- `--clusters <path>` (keep only `--cluster` representative paths)
- `--cluster-by <path>` (keep all or a subset of paths, group/order rows by `--cluster` assignment)
- `--cluster-rows`
- `--cluster-cols`
- `--max-paths <N>`
- `--max-nodes <N>`

## Edge Coverage Heatmap

Use `Rscript scripts/plot_edge_coverage_heatmap.R` to visualize an inspect edge-count table — the same dependency and the same interface as the node heatmap, minus `--value` (edge cells are plain traversal counts) and `--node-lengths`

```bash
Rscript scripts/plot_edge_coverage_heatmap.R \
  --table results/real_data/lpa/inspect/inspect.bubble_7.edge_counts.tsv \
  --out results/real_data/lpa/inspect/inspect.bubble_7.edge_coverage \
  --clusters results/real_data/lpa/inspect/inspect.bubble_7.clusters.tsv \
  --cluster-by results/real_data/lpa/inspect/inspect.bubble_7.clusters.tsv
```