# Inspect Utility

CLI: `panvar inspect`

## What it does

A sanity-check utility for one called bubble (or all of them) before going downstream. Given a GFA and a
`bubble`/`panphorte` bubble CSV, it writes, per bubble: a multi-FASTA of each crossing path's `source → sink`
allele, a **node-count** matrix (how each path traverses the internal nodes, with orientation), an
**edge-count** matrix (adjacencies — the tandem signal), and node lengths. With `--cluster` it also groups
structurally identical haplotypes (used for representative-only plots).

Clustering mechanism + worked trace: **[algorithms/inspect.md](../algorithms/inspect.md)**.

## Required inputs

- `-i, --gfa <graph.gfa>` — use the **sorted** GFA from `bubble`/`panphorte` (node ids match the CSV).
- one of `-b, --bubble-prefix-in <prefix>` (auto-uses `<prefix>.bubbles.csv`) or `-c, --bubbles-csv <path>`.
- `--bubble-id <N>` optional; omit to inspect every bubble (one output set each).

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-prefix <p>` | output prefix | — |
| `--bubble-id <N>` | restrict to one bubble (enables the explicit `--*-out` overrides) | all bubbles |
| `--cluster` | group crossing paths by `source → sink` [walk](../algorithms/inspect.md#terms) → `<prefix>.bubble_<N>.clusters.tsv` | off |
| `--cluster-similarity <f>` | walk-similarity threshold for `--cluster` (higher = finer copy-number bands) | `0.90` |
| `--fasta-out` / `--table-out` / `--edge-table-out <path>` | override the FASTA / node-count / edge-count paths (single `--bubble-id` only) | derived |

## Outputs (per bubble)

| file | contents |
|------|----------|
| `<prefix>.bubble_<N>.paths.fa.gz` | one record per crossing path = its spelled `source → sink` allele (canonicalized; `source_to_sink=` header notes if the path ran reversed) |
| `<prefix>.bubble_<N>.node_counts.tsv` | `path_name, path_length_bp`, one `node.<id>` col = `total:forward:reverse` (interval-local) |
| `<prefix>.bubble_<N>.edge_counts.tsv` | `path_name, path_length_bp`, one `edge.<from±>to±>` col = traversal count (self-loops/back-edges show as > 1) |
| `<prefix>.bubble_<N>.node_lengths.tsv` | `node_id, length_bp` in node-column order (lets the heatmap scale x by length) |
| `<prefix>.bubble_<N>.clusters.tsv` | (with `--cluster`) `cluster_id, n_paths, representative_path, members` |

## Plotting

Two R helpers (need `Rscript` + `ggplot2`) visualize the count tables; both take `--clusters` (plot only
cluster representatives) and `--cluster-by` (group/order rows by cluster):

```bash
Rscript scripts/plot_node_coverage_heatmap.R --table <…>.node_counts.tsv \
  --node-lengths <…>.node_lengths.tsv --cluster-by <…>.clusters.tsv --out <…>.node_coverage
Rscript scripts/plot_edge_coverage_heatmap.R --table <…>.edge_counts.tsv \
  --cluster-by <…>.clusters.tsv --out <…>.edge_coverage
```

Node-heatmap extras: `--value total|forward|reverse`, `--transform raw|log1p`, `--length-transform
raw|sqrt|log1p`, `--max-paths/--max-nodes`. (The per-gene scripts in `scripts/genes/` list these commented,
ready to uncomment.)

## Example

Matches `scripts/genes/lpa.sh` (LPA uses `--cluster-similarity 0.95`; the other genes use `0.97`):

```bash
./build/panvar inspect \
  -i results/real_data/lpa/panphorte/panphorte.normalized.sorted.gfa \
  --bubble-prefix-in results/real_data/lpa/panphorte/panphorte \
  --bubble-id 7 -o results/real_data/lpa/inspect/inspect \
  --cluster --cluster-similarity 0.95
```

Drop `--bubble-id` to inspect all bubbles. Clustering algorithm & worked example:
[algorithms/inspect.md](../algorithms/inspect.md). References: [references.md](../references.md#inspect).
