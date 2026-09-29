# Module `inspect`

CLI: `panvar inspect`

## What it does

A sanity-check utility for one called bubble (or all of them) before going downstream. Given a `bubble`/`panphorte`/`refine` GFA (Graphical Fragment Assembly) and its bubble CSV, it writes, per bubble:
- a multi-FASTA of each crossing path's `source → sink` allele, and the same records as node walks (`1+,7-,8+`)
- a node-count matrix (how each path traverses the internal nodes, with orientation)
- an edge-count matrix (adjacencies)
- node lengths
- with `--genotype-table`, an allele table (one row per haplotype, or per sample with a [cosigt](https://github.com/davidebolo1993/cosigt) table, one column per distinct allele), the allele FASTA whose ids are its columns, its node walks, and a node dictionary (sequence and reference position of every node on a crossing walk)

With `--cluster` it also groups sufficiently similar walks and selects a representative for each group.

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
| `--genotype-table` | write the per-bubble allele table and allele FASTA (see Outputs) | off |
| `--samples <cosigt.tsv>` | cosigt `sample <tab> hap1 <tab> hap2` table, as for [`describe`](describe.md); table rows become samples. Implies `--genotype-table` | — |
| `-r, --reference-path <name>` | orient `--genotype-table` alleles and their walks to this path's strand, and position the node dictionary on this path (exact name, or a unique case-insensitive substring) | bubble `source → sink` |
| `--fasta-out` / `--table-out` / `--edge-table-out <path>` | override the FASTA/node-count/edge-count paths (single `--bubble-id` only) | derived |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.bubble_<N>.paths.fa.gz` | one record per crossing path = its spelled `source → sink` allele (canonicalized; `source_to_sink=` header notes if the path ran reversed) |
| `<prefix>.bubble_<N>.walks.fa.gz` | the paths FASTA's records as walks: the same header, then the oriented nodes the path takes, `1+,7-,8+`, on one line. Source and sink are included and the walk runs `source → sink`, so concatenating its nodes (reverse-complementing `-` steps) spells exactly the record's sequence. Beside an explicit `--fasta-out X.fa.gz` it is `X.walks.fa.gz` |
| `<prefix>.bubble_<N>.node_counts.tsv` | `path_name, path_length_bp`, one `node.<id>` col = `total:forward:reverse` (interval-local) |
| `<prefix>.bubble_<N>.edge_counts.tsv` | `path_name, path_length_bp`, one `edge.<from±>to±>` col = traversal count (self-loops/back-edges show as > 1) |
| `<prefix>.bubble_<N>.node_lengths.tsv` | `node_id, length_bp` in node-column order (lets the heatmap scale x by length) |
| `<prefix>.bubble_<N>.alleles.fa.gz` | (with `--genotype-table`) one record per distinct `source → sink` sequence, `>alleleK hap1,hap2,…`: the id is the table's column, the description lists the haplotypes spelling it. Alleles are numbered by descending haplotype count, ties by first path in the GFA. With `-r` every allele is reverse-complemented when the reference crosses the bubble `sink → source` |
| `<prefix>.bubble_<N>.alleles.walks.fa.gz` | (with `--genotype-table`) the allele FASTA's records as walks, under the same `>alleleK hap1,hap2,…` header and oriented as the allele is, so each walk spells its allele. Different walks can spell the same bases: such an allele gets one record per walk, `>alleleK <haplotypes taking it> walk=i/n` |
| `<prefix>.bubble_<N>.nodes.json.gz` | (with `--genotype-table`) `{"bubble", "source", "sink", "reference_path", "reference_crosses_bubble", "reference_strand", "haplotypes_crossing", "nodes": {<id>: {...}}}`, one entry per node on any crossing walk: `seq` (forward strand), `length`, `haplotypes` (crossing haplotypes visiting it), `order` (rank along the reference; first appearance when there is none) and `on_reference`. A node the `-r` path visits has `ref_offset`, the 0-based bp offset of its first base from the first base of the reference's walk through the bubble (source included), along that path, plus 1-based `ref_start`/`ref_end` when the path name carries an interval (`chr:start-end`); a node visited more than once also has `ref_occurrences`, and its position is its first visit. A node the reference does not visit has no coordinate of its own: `between` names the reference nodes flanking it on the first walk visiting it, with `between_ref_offsets` (and `between_ref_positions`) giving the last base of the one before and the first base of the one after. Without `-r` every position is `null` |
| `<prefix>.bubble_<N>.genotypes.tsv.gz` | (with `--genotype-table`) `haplotype` (or `sample`) then one `alleleK` column per allele. A haplotype row is 0/1; a sample row counts how many of its haplotypes spell each allele, so a diploid row sums to 2 and a sample assigned the same haplotype twice reads 2. `NA` fills a row whose haplotype, or any of a sample's haplotypes, does not cross the bubble |
| `<prefix>.bubble_<N>.clusters.tsv` | (with `--cluster`) one row per cluster: `cluster_id` (group id), `n_paths` (paths in it), `representative_path` (the exemplar plotted), `members` (`;`-separated path names) |

`scripts/plot_node_coverage_heatmap.R` and `scripts/plot_edge_coverage_heatmap.R` render the count tables as heatmaps, optionally grouped by the cluster assignment. Each documents its own options under `--help`.

## Limitations

- Clustering compares every pair of distinct walks over a dense matrix, so it is region-scale rather than cohort-scale. It warns above two thousand distinct walks and refuses above twenty-five thousand.
- Walk similarity is estimated from a sketch of the walk's node steps when a walk is long enough to exceed the sketch, and computed exactly when it is not. The estimate carries sampling error near a threshold.
- Clusters are connected components at the similarity threshold, so membership is transitive: two walks below the threshold can land in one cluster through a chain of intermediates.
- `--genotype-table` alleles are distinct by exact sequence, so a single base makes a new allele; a hypervariable site comes out with close to one allele per haplotype (459 alleles over 466 haplotypes at bubble 7 of the LPA bubble stage).
- Orientation under `-r` is decided per bubble by how the reference crosses it. A bubble the reference does not cross stays `source → sink`, with a warning.
- Under `--genotype-table`, a path name containing whitespace or a comma is refused, since the allele FASTA lists haplotypes comma-separated after the id.
- Sequences are spelled by concatenating whole nodes, so a graph whose links carry a non-zero overlap is refused rather than mis-measured.

## Example

See the [walkthrough](../walkthrough.md) for this module in a full end-to-end run.
