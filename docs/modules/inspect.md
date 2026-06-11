# Inspect Utility

Date: 2026-06-07

CLI entrypoint:

- `panvar inspect`

## What it does

`inspect` is a small utility for checking one called bubble, or all called bubbles, before moving to allele clustering.

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

You can override them explicitly with:

- `--fasta-out <path>`
- `--table-out <path>` (node counts)
- `--edge-table-out <path>` (edge counts)

Explicit output paths are only accepted with `--bubble-id`, because all-bubble mode needs one output set per bubble.

Output directories are auto-created when missing.

## FASTA output

Each FASTA record is the spelled source-to-sink bubble interval for one path.

Headers include:

- path name
- bubble id
- source/sink node ids
- interval length in bp
- whether the original best interval was source-to-sink or sink-to-source
- path-step interval coordinates

If the best interval is found as sink-to-source, `inspect` reverse-complements it so the FASTA is still reported in canonical source-to-sink orientation.

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

The orientation is counted after canonical source-to-sink normalization, matching allele clustering. A `reverse` count therefore means the node is traversed in reverse orientation in the canonical bubble allele. If the graph has no `-`/`<` path steps in the selected intervals, reverse counts can legitimately be zero.

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

## Node Coverage Heatmap

Use `scripts/plot_node_coverage_heatmap.R` to visualize an inspect node-count table.

Dependency:

- `Rscript` with `ggplot2` (`conda install -y -c conda-forge r-base r-ggplot2`)

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

Useful plotting options:

- `--value total|forward|reverse`
- `--transform raw|log1p`
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

This writes `<out>.png` and `<out>.pdf`. Options: `--transform raw|log1p`, `--cluster-rows`,
`--cluster-cols`, `--max-paths <N>`, `--max-edges <N>`, `--width`, `--height`.

## Interval selection

For each path, `inspect` uses the same source/sink interval logic as module 2:

1. find source-to-sink and sink-to-source path intervals
2. prefer the interval containing the most bubble-internal node traversals
3. tie-break by shorter span
4. tie-break again by earlier path interval

Paths without a valid source/sink interval crossing at least one internal node are skipped.

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
