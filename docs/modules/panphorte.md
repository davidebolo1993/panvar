# Module `panphorte`

CLI: `panvar panphorte`

## What it does

Normalizes contiguous, tandem repeat (TR) bubbles into a compact, copy-number-explicit form. It:
- collapses each detected tandem array to a single repeat-unit (`REP`) node carrying a self-loop, so a haplotype's copy number (`CN`) is recorded as the number of self-loop traversals rather than as a spurious insertion
- writes a new GFA (Graphical Fragment Assembly) — plus a report and `CN`-provenance table — for `inspect` and `call`

Collapse is exact by default (byte-identical copies) and can be made approximate to fold divergent copies. Only genuine population TR are folded, not rare duplications.

Detection normally measures the repeat period in node steps, which requires the graph to split nodes at repeat-unit boundaries. When a bubble yields no step period — as happens on graphs whose node boundaries were set by alignment rather than by the repeat — the unit is seeded from the spelled sequence instead, so an array split across arbitrary node boundaries is still folded.

Algorithm and worked trace: [algorithms/panphorte.md](../algorithms/panphorte.md).

## Required inputs

- `-i, --gfa <graph.gfa>` — the sorted GFA from `bubble`.
- one of `-b, --bubble-prefix-in <bubble-prefix>` (auto-uses `<prefix>.bubbles.csv`) or `-c, --bubbles-csv <path>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-prefix <p>` | output prefix | — |
| `-r, --reference-path <name>` | sort and re-snarl the normalized graph along this reference | — |
| `--min-similarity <f>` | identity to treat a block as a copy of the unit; `1.0` = exact (sequence-preserving), `< 1.0` = approximate/lossy collapse of divergent copies | `1.0` |
| `--min-unit-bp <N>` | minimum repeat-unit span to normalize | `50` |
| `--min-copies <N>` | tandem copies needed (in some haplotype) to treat a bubble as an array; once an array, every haplotype with ≥1 copy folds | `2` |
| `--min-array-prevalence <f>` | min fraction of bubble-traversing haplotypes that must carry a ≥`min-copies` array for the bubble to fold; separates a true population TR from a rare/private duplication | `0.5` |
| `--max-interruption-frac <f>` | max fraction of an array's bp that may be interruptions | `0.25` |
| `--threads <N>` | workers for the approximate seed scan/`CN` detection (`0` = auto) | `0` |
| `--gtf <path>` | after re-sorting, project genes ( `<prefix>.bandage_genes.csv`, needs a PanSN (Pangenome Sequence Naming) `--reference-path`); separate from `bubble --gtf` because collapse, if applicable, renumbers nodes | — |
| `-q, --quiet` | disable progress/logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.normalized.sorted.gfa` | sorted, call-ready normalized graph (with `--reference-path`) |
| `<prefix>.normalized.gfa` | unsorted normalized graph (without `--reference-path`; must be re-sorted before `call`) |
| `<prefix>.bubbles.csv` | re-snarled bubbles (with `--reference-path`) |
| `<prefix>.panphorte.report.tsv` | one row per bubble (columns below) |
| `<prefix>.panphorte.copies.tsv` | (approximate mode) one row per (haplotype, array) — the `CN` provenance for `call` (columns below) |
| `<prefix>.bandage_nodes.csv` | Bandage node colors |
| `<prefix>.bandage_genes.csv` | Bandage gene track (with `--gtf` and PanSN `--reference-path`) |

`<prefix>.panphorte.report.tsv` columns:

| column | meaning |
|--------|---------|
| `bubble_id` | the bubble considered |
| `normalized` | whether it was folded into a `REP` self-loop (`1`/`0`) |
| `unit_bp` | length of the detected repeat unit (bp) |
| `paths_normalized` | how many paths were rewritten |
| `min_copies`, `max_copies` | smallest / largest copy count seen across haplotypes |
| `nodes_collapsed` | number of original nodes folded into the unit |

`<prefix>.panphorte.copies.tsv` columns:

| column | meaning |
|--------|---------|
| `path_name`, `sample` | the haplotype path and its sample |
| `bubble_id` | the array's bubble |
| `copies` | detected tandem copy count for this haplotype |
| `unit_bp`, `region_bp` | repeat-unit length and total array span (bp) |
| `orientations` | strand(s) the copies were found on |
| `mean_identity` | mean identity of the copies to the unit (approximate collapse) |
| `from_node`, `to_node` | graph span of the array |

## Example

See the [LPA walkthrough](../walkthrough.md) for this module in a full end-to-end run.

## Re-snarling under `--reference-path`

With `--reference-path`, the normalized graph is sorted and re-snarled so the output is call-ready. That re-snarl applies its own interior-span filter, `--resnarl-min-variant-bp` (default 50). It is a separate decision from the one that produced the input bubbles: a CSV built with a different threshold can otherwise lose bubbles that normalization never touched. Bubble ids are reassigned by the re-snarl, so an id in the output does not correspond to the same id in the input; the run reports the count change rather than pretending the ids match.

## Report columns

`<prefix>.panphorte.report.tsv` carries, per bubble: `normalized`, `unit_bp`, `paths_normalized`, `min_copies`, `max_copies`, `interruptions_bp`, `nodes_collapsed`, and the diagnostics `n_traversing`, `n_motif_carriers`, `prevalence`, `n_motifs`, `copies_declined_partial_boundary`, `paths_with_partial_boundary` and `status` (`normalized`, `below_prevalence`, `no_tandem_detected`, `no_seed`, `partial_boundary`). `normalized` is `yes`/`no`.

`<prefix>.panphorte.rep_provenance.tsv` maps each REP node to the site and motif it stands for (`rep_node`, `input_bubble_id`, `canonical_motif`, `phase_unit`, `unit_bp`, `copy_quantum`). Two phase-rotated units at one site cannot share an unsplit REP node while exact spelling is preserved, so they become separate nodes; without this table a consumer counting REP occurrences sees two independent DUPs instead of one site's copy number. The ids are the ones panphorte created — under `--reference-path` the sort renumbers them, so a consumer must map through the sorted graph. `call` does not yet aggregate by site. `<prefix>.panphorte.copies.tsv` names its graph columns `input_bubble_id`, `input_from_node` and `input_to_node`: they describe the graph handed in, and sorting renumbers nodes and reassigns bubble ids.
