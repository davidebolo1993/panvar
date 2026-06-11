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
- which internal nodes are reused, and in which orientation?
- which adjacencies (edges) are reused — e.g. a tandem self-loop shows up as a high edge count?

## Required inputs

- `--gfa <graph.gfa>`
- one of:
  - `--bubble-prefix-in <module1-prefix>` (auto uses `<module1-prefix>.bubbles.csv`)
  - `--bubbles-csv <module1.bubbles.csv>`

Optional:

- `--bubble-id <N>`: restrict to one bubble

If `--bubble-id` is omitted, `inspect` writes outputs for every bubble in the input CSV.

## Outputs

Default output paths are derived from `--out-prefix` and the bubble ID:

- `<out-prefix>.bubble_<N>.paths.fa.gz`
- `<out-prefix>.bubble_<N>.node_counts.tsv`
- `<out-prefix>.bubble_<N>.edge_counts.tsv`
- `<out-prefix>.bubble_<N>.node_lengths.tsv` — per-node bp lengths (see below)
- `<out-prefix>.bubble_<N>.clusters.tsv` — only when `--cluster` is set (see below)

You can override the first three explicitly with:

- `--fasta-out <path>`
- `--table-out <path>` (node counts)
- `--edge-table-out <path>` (edge counts)

Explicit output paths are only accepted with `--bubble-id`, because all-bubble mode needs one output set per bubble.

Output directories are auto-created when missing.

When inspecting all bubbles, a progress bar is printed to stderr; the machine-readable summary stays on
stdout. With a single `--bubble-id`, the per-bubble detail is printed instead.

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

### Why some records are reverse-complemented

A bubble has a fixed `source` and `sink`, but an individual path can walk *through* it in either
direction — some paths enter at the `source` and leave at the `sink`, others run the opposite way
(enter at the `sink`, leave at the `source`). `inspect` searches both directions and records which one
matched in the `source_to_sink=` header field.

When a path crossed **sink-to-source**, `inspect` reorients that path's *node walk* — it reverses the
step order and flips each node's strand — so the emitted allele always reads `source → sink`. Because
the spelled DNA follows the node strands, flipping them yields the reverse-complement of that path's
raw sequence. This is driven purely by traversal direction, not by sequence content.

The point is comparability: every record (and every node/edge orientation count below) is reported in
one common `source → sink` frame, so the FASTA alleles line up and can be aligned/clustered directly.

## Node-count table

Columns:

- `path_name`
- `path_length_bp`
- one `node.<id>` column per bubble-internal node from `inside_nodes`

Each node cell is:

```text
total:forward:reverse
```

For example, `3:2:1` means that, inside this bubble interval only, the path traverses that node 3 times total: 2 times in forward orientation and 1 time in reverse orientation.

Important: counts are interval-local. If the same path reuses a node elsewhere outside the selected source/sink interval, that outside traversal is not counted for this bubble.

Orientation is counted *after* the canonical `source → sink` normalization above. A `reverse` count
therefore means the node is traversed against the canonical frame in this bubble allele — typically a
real local inversion — rather than just a path that happened to run the other way. If the graph has no
`-`/`<` path steps in the selected intervals, reverse counts can legitimately be zero.

## Edge-count table

Columns:

- `path_name`
- `path_length_bp`
- one `edge.<from><sign>><to><sign>` column per bubble-internal adjacency, e.g. `edge.1886+>1887+`

Each edge cell is an integer: how many times the path traverses that **directed, orientation-aware
adjacency** inside the canonical source-to-sink bubble interval. The edge key encodes the orientation
of both endpoints (`+`/`-`), so a forward and a reverse traversal of the same node pair are distinct
columns.

The edge matrix is **adjacency-aware**, complementing the node matrix: a tandem repeat that loops over
the same unit repeats the *same edge*, so a self-loop / back-edge shows up as a cell value `> 1` (the
copy number). Edge columns are the union of adjacencies observed across all paths, sorted for a stable,
deterministic column order. As with nodes, counts are interval-local.

## Path clustering (`--cluster`)

`--cluster` groups the paths crossing a bubble by their `source → sink` walk, so structurally identical
haplotypes collapse to one representative. It is an inspection aid only — nothing downstream depends on
it — but it is handy for plotting (see the heatmap `--clusters` option).

Each path's canonical walk is summarized as an **oriented, bp-weighted token multiset**: every step
contributes its node's bp length to the token `<node_id><strand>`. Two walks are compared with weighted
Jaccard — `sum(min) / sum(max)` over the shared tokens — which captures inversions (strand is part of
the token) and copy number (repeated nodes add weight) while ignoring exact step order. Walks are
clustered greedily: identical walks collapse first, then each distinct walk joins the most similar
existing cluster whose representative is at least `--cluster-similarity` similar (default `0.90`), or
starts a new cluster. The representative is the cluster **medoid** (the walk with the highest mean
similarity to the rest; ties go to the most-supported walk).

Output `<out-prefix>.bubble_<N>.clusters.tsv`:

```text
cluster_id	n_paths	representative_path	members
```

`members` is a `;`-separated list of path names. Lower `--cluster-similarity` to merge more aggressively,
raise it to split.

## Example

```bash
./build/panvar inspect \
  -i tests/real_data/c4.gfa \
  --bubble-prefix-in tests/results/c4/bubble \
  --bubble-id 1 \
  -o tests/results/c4/inspect/bubble_1
```

This writes:

- `tests/results/c4/inspect/bubble_1.bubble_1.paths.fa.gz`
- `tests/results/c4/inspect/bubble_1.bubble_1.node_counts.tsv`
- `tests/results/c4/inspect/bubble_1.bubble_1.edge_counts.tsv`

To inspect all bubbles:

```bash
./build/panvar inspect \
  -i tests/real_data/c4.gfa \
  --bubble-prefix-in tests/results/c4/bubble \
  -o tests/results/c4/inspect/all
```

This writes one FASTA/table pair per bubble, for example:

- `tests/results/c4/inspect/all.bubble_1.paths.fa.gz`
- `tests/results/c4/inspect/all.bubble_1.node_counts.tsv`
- `tests/results/c4/inspect/all.bubble_1.edge_counts.tsv`
- `tests/results/c4/inspect/all.bubble_2.paths.fa.gz`
- `tests/results/c4/inspect/all.bubble_2.node_counts.tsv`
- `tests/results/c4/inspect/all.bubble_2.edge_counts.tsv`


## Node Coverage Heatmap

Use `scripts/plot_node_coverage_heatmap.R` to visualize an inspect node-count table.

Dependency:

- `Rscript` with `ggplot2` 

Example:

```bash
scripts/plot_node_coverage_heatmap.R \
  --table tests/results/c4/inspect/bubble_1.bubble_1.node_counts.tsv \
  --out tests/results/c4/inspect/bubble_1.node_coverage
```

This writes:

- `tests/results/c4/inspect/bubble_1.node_coverage.png`
- `tests/results/c4/inspect/bubble_1.node_coverage.pdf`

By default, the heatmap uses total node coverage per path. To check only reverse-orientation traversals:

```bash
scripts/plot_node_coverage_heatmap.R \
  --table tests/results/c4/inspect/bubble_1.bubble_1.node_counts.tsv \
  --out tests/results/c4/inspect/bubble_1.reverse_coverage \
  --value reverse
```

To scale the x-axis by node length (thin tiles for SNP-sized nodes, wide tiles for long nodes) while
keeping the column order, pass the node-length sidecar:

```bash
scripts/plot_node_coverage_heatmap.R \
  --table tests/results/c4/inspect/all.bubble_4.node_counts.tsv \
  --node-lengths tests/results/c4/inspect/all.bubble_4.node_lengths.tsv \
  --out tests/results/c4/inspect/bubble_4.node_coverage
```

To plot only the cluster representatives, pass the `--cluster` output:

```bash
scripts/plot_node_coverage_heatmap.R \
  --table tests/results/c4/inspect/all.bubble_4.node_counts.tsv \
  --clusters tests/results/c4/inspect/all.bubble_4.clusters.tsv \
  --out tests/results/c4/inspect/bubble_4.representatives
```

Useful plotting options:

- `--value total|forward|reverse`
- `--transform raw|log1p`
- `--node-lengths <path>` (length-scale x; tile widths via `--length-transform raw|sqrt|log1p`, default `sqrt`)
- `--clusters <path>` (keep only `--cluster` representative paths)
- `--cluster-rows`
- `--cluster-cols`
- `--max-paths <N>`
- `--max-nodes <N>`

## Edge Coverage Heatmap

Use `scripts/plot_edge_coverage_heatmap.R` to visualize an inspect edge-count table — the same
dependency (`Rscript` + `ggplot2`) and the same interface as the node heatmap, minus `--value` (edge
cells are plain traversal counts):

```bash
scripts/plot_edge_coverage_heatmap.R \
  --table tests/results/c4/inspect/bubble_4.edge_counts.tsv \
  --out tests/results/c4/inspect/bubble_4.edge_coverage \
  --transform log1p --cluster-rows --cluster-cols
```

This writes `<out>.png` and `<out>.pdf`. Options: `--transform raw|log1p`, `--clusters <path>`
(keep only `--cluster` representatives), `--cluster-rows`, `--cluster-cols`, `--max-paths <N>`,
`--max-edges <N>`, `--width`, `--height`.
