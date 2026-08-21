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

The CSV must describe the graph handed in and its sites must be disjoint: a node named by no segment, a duplicate bubble id, or two bubbles claiming the same interior node is refused before any work starts. Overlapping sites are not a formality — they describe one piece of sequence twice, so folding both rewrites one span inside the other and the same array is folded at two scales. `bubble` no longer emits such a pair: it runs a deterministic conflict pass before assigning ids, so the sites it emits are pairwise disjoint. ANKRD36C used to fail here, with a site enclosing all ten others folding the same 5616 bp unit as the site inside it; `bubble` now keeps the ten smaller sites, which reconstruct the haplotypes far better than the one enclosing site did. The check remains because the CSV can come from anywhere, and it applies to the *selected* sites, so `--bubble-id` is still the way to work with a hand-made overlapping set.

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
| `--bubble-id <N>` | restrict to these sites (repeatable); overlap and disjointness are then judged among the selected ones only | all |
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
| `<prefix>.panphorte.rep_provenance.tsv` | one row per REP node: the site and motif it stands for (columns below) |
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
| `nodes_collapsed` | distinct original nodes the rewrite replaced at this bubble. A node still walked at another site survives the run, so this counts what was replaced here, not what was deleted from the graph |

`<prefix>.panphorte.copies.tsv` columns:

| column | meaning |
|--------|---------|
| `path_name`, `sample` | the haplotype path and its sample |
| `input_bubble_id` | the array's bubble, numbered as in the CSV handed in |
| `copies` | detected tandem copy count for this haplotype |
| `unit_bp`, `region_bp` | repeat-unit length and total array span (bp) |
| `orientations` | strand(s) the copies were found on |
| `mean_identity` | mean identity of the copies to the unit (approximate collapse) |
| `input_from_node`, `input_to_node` | graph span of the array, in the input graph's numbering |

The three `input_` columns describe the graph handed in. Under `--reference-path` the output is sorted, which renumbers nodes and reassigns bubble ids, so the same numbers still exist afterwards and refer to different sequence.

`<prefix>.panphorte.rep_provenance.tsv` columns:

| column | meaning |
|--------|---------|
| `created_rep_node` | the id panphorte created the REP under, in the input graph's numbering |
| `output_rep_node` | the id it is delivered under — the join key against the GFA beside this file |
| `input_bubble_id` | the site it stands for, numbered as in the CSV handed in |
| `canonical_motif` | the primitive motif, least rotation, min over both strands — the same whichever way the node ends up oriented, so it is the key for grouping phase-rotated units at one site |
| `created_phase_unit` | the sequence panphorte built the node with |
| `output_phase_unit` | the sequence the delivered node actually spells; sorting may flip a REP, so this is not always the same string |
| `unit_bp`, `copy_quantum` | unit length, and how many copies one REP step stands for (always 1) |

## Example

See the [LPA walkthrough](../walkthrough.md) for this module in a full end-to-end run.

## Re-snarling under `--reference-path`

With `--reference-path`, the normalized graph is sorted and re-snarled so the output is call-ready. That re-snarl applies its own interior-span filter, `--resnarl-min-variant-bp` (default 50). It is a separate decision from the one that produced the input bubbles: a CSV built with a different threshold can otherwise lose bubbles that normalization never touched. Bubble ids are reassigned by the re-snarl, so an id in the output does not correspond to the same id in the input; the run reports the count change rather than pretending the ids match.

## Report columns

`<prefix>.panphorte.report.tsv` carries, per bubble: `normalized`, `unit_bp`, `paths_normalized`, `min_copies`, `max_copies`, `interruptions_bp`, `nodes_collapsed`, and the diagnostics `n_traversing`, `n_motif_carriers`, `prevalence`, `n_motifs`, `copies_declined_partial_boundary`, `paths_with_partial_boundary` and `status` (`normalized`, `below_prevalence`, `no_tandem_detected`, `no_seed`, `partial_boundary`). `normalized` is `yes`/`no`.

Two phase-rotated units at one site cannot share an unsplit REP node while exact spelling is preserved, so they become separate nodes; without the provenance table a consumer counting REP occurrences sees two independent DUPs instead of one site's copy number. Group them on `canonical_motif` and join to the graph on `output_rep_node`. `call` does not yet aggregate by site.
