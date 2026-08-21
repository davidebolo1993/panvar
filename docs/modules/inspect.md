# Module `inspect`

CLI: `panvar inspect`

## What it does

A sanity-check utility for one called bubble (or all of them) before going downstream. Given a `bubble`/`panphorte` GFA (Graphical Fragment Assembly) and its bubble CSV, it writes, per bubble:
- a multi-FASTA of each crossing path's `source → sink` allele
- a node-count matrix (how each path traverses the internal nodes, with orientation)
- an edge-count matrix (adjacencies)
- node lengths

With `--cluster` it also groups structurally identical haplotypes (used for representative-only plots).

Algorithm and worked trace: [algorithms/inspect.md](../algorithms/inspect.md).

## Required inputs

- `-i, --gfa <graph.gfa>` — the sorted GFA from `bubble`/`panphorte`/`refine` (node ids should match the CSV to `-b`).
- one of `-b, --bubble-prefix-in <prefix>` (auto-uses `<prefix>.bubbles.csv`) or `-c, --bubbles-csv <path>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-prefix <p>` | output prefix | — |
| `--bubble-id <N>` | restrict to one bubble | all bubbles |
| `--cluster` | group crossing paths by how similarly they traverse the bubble (`source → sink` walk; `<prefix>.bubble_<N>.clusters.tsv`) | off |
| `--cluster-similarity <f>` | walk-similarity threshold for `--cluster` | `0.90` |
| `--fasta-out` / `--table-out` / `--edge-table-out <path>` | override the FASTA/node-count/edge-count paths (single `--bubble-id` only) | derived |

## Outputs (per bubble)

| file | contents |
|------|----------|
| `<prefix>.bubble_<N>.paths.fa.gz` | one record per crossing path = its spelled `source → sink` allele (canonicalized; `source_to_sink=` header notes if the path ran reversed) |
| `<prefix>.bubble_<N>.node_counts.tsv` | `path_name, path_length_bp`, one `node.<id>` col = `total:forward:reverse` (interval-local) |
| `<prefix>.bubble_<N>.edge_counts.tsv` | `path_name, path_length_bp`, one `edge.<from±>to±>` col = traversal count (self-loops/back-edges show as > 1) |
| `<prefix>.bubble_<N>.node_lengths.tsv` | `node_id, length_bp` in node-column order (lets the heatmap scale x by length) |
| `<prefix>.bubble_<N>.clusters.tsv` | (with `--cluster`) one row per cluster: `cluster_id` (group id), `n_paths` (paths in it), `representative_path` (the exemplar plotted), `members` (`;`-separated path names) |

## Limitations

- Clustering compares every pair of distinct walks over a dense matrix, so it is region-scale rather than cohort-scale. It warns above two thousand distinct walks and refuses above twenty-five thousand.
- Walk similarity is estimated from a sketch of the walk's node steps when a walk is long enough to exceed the sketch, and computed exactly when it is not. The estimate carries sampling error near a threshold.
- Clusters are connected components at the similarity threshold, so membership is transitive: two walks below the threshold can land in one cluster through a chain of intermediates.
- Sequences are spelled by concatenating whole nodes, so a graph whose links carry a non-zero overlap is refused rather than mis-measured.

## Plotting

`scripts/plot_node_coverage_heatmap.R` and `scripts/plot_edge_coverage_heatmap.R` render the count tables as heatmaps, optionally grouped by the cluster assignment. Each documents its own options under `--help`.

## Example

See the [walkthrough](../walkthrough.md) for this module in a full end-to-end run.
