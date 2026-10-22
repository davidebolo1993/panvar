# Module `panphorte`

CLI: `panvar panphorte`

## What it does

Normalizes contiguous tandem-repeat bubbles into a compact, copy-number-explicit form. It:
- detects, per bubble, whether the haplotypes carry a tandem array and what its repeat unit is;
- collapses each detected array to a single repeat-unit node carrying a self-loop, so a haplotype's copy number is the number of self-loop traversals rather than a spurious insertion;
- optionally re-sorts and re-snarls the result along the reference, so the output is ready for `call`;
- writes the normalized GFA (Graphical Fragment Assembly), a per-bubble report, and the tables that let a consumer trace a copy number back to the site it came from.

Collapse is exact by default, folding only byte-identical copies, and can be made approximate to fold divergent ones. Only arrays carried by a substantial share of the panel are folded, so a rare private duplication is left for `call` to type as an ordinary event.

Detection normally measures the repeat period in node steps, which needs the graph to split nodes at unit boundaries. Where a bubble yields no step period, because node boundaries were set by alignment rather than by the repeat, the unit is seeded from the spelled sequence instead.

Algorithm and worked trace: [algorithms/panphorte.md](../algorithms/panphorte.md).

## Required inputs

- `-i, --gfa <graph.gfa>` — the sorted GFA from `bubble`
- one of `-b, --bubble-prefix-in <prefix>` (auto-uses `<prefix>.bubbles.csv`) or `-c, --bubbles-csv <path>`

The CSV must describe the graph handed in, and the sites it selects must be pairwise disjoint. A node named by no segment, a duplicate bubble id, or two sites claiming the same interior node is refused before any work starts: overlapping sites describe one piece of sequence twice, so folding both would rewrite one span inside the other and fold the same array at two scales. Disjointness is judged among the selected sites, so `--bubble-id` is the way to work with a set that overlaps.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-prefix <p>` | output prefix | — |
| `-r, --reference-path <name>` | sort and re-snarl the normalized graph along this reference, producing call-ready output | — |
| `--min-similarity <f>` | identity at which a block counts as a copy of the unit; `1.0` is exact and sequence-preserving, below that the collapse of divergent copies is lossy | `1.0` |
| `--min-unit-bp <N>` | minimum repeat-unit span to normalize | `50` |
| `--min-copies <N>` | tandem copies some haplotype must carry for a bubble to count as an array; once it does, every haplotype with at least one copy folds | `2` |
| `--min-array-prevalence <f>` | fraction of the bubble's haplotypes that must carry an array before it is folded, which is what separates a population repeat from a private duplication | `0.5` |
| `--max-interruption-frac <f>` | fraction of an array's bases that may be interruptions between copies | `0.25` |
| `--resnarl-min-variant-bp <N>` | interior-span filter applied when re-snarling under `--reference-path` | `50` |
| `--bubble-id <N>` | restrict to these sites (repeatable) | all |
| `--allow-partial-boundary` | override the safety refusal if a copy cannot be bounded by a step range; normally unnecessary because fragment nodes preserve partial-node flanks | off |
| `--no-flip` | with `--reference-path`, skip reorienting to the reference strand | off |
| `--threads <N>` | workers for approximate detection (`0` = auto) | `0` |
| `--refine` | after folding, realign the flanks around each folded array (see below); needs `--reference-path` | off |
| `--refine-max-poa-bp <N>` | longest flank `--refine` aligns | `5000` |
| `--refine-max-walks <N>` | most distinct flanks per side `--refine` aligns | `500` |
| `--refine-max-poa-work <N>` | longest × total flank bp per site (`0` = no limit) | `0` |
| `--gtf <path>` | after re-sorting, project genes onto reference nodes; separate from `bubble --gtf` because collapsing renumbers nodes | — |
| `-q, --quiet` | disable progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.normalized.sorted.gfa` | sorted, call-ready normalized graph (with `--reference-path`) |
| `<prefix>.normalized.gfa` | unsorted normalized graph (without `--reference-path`; re-sort before `call`) |
| `<prefix>.bubbles.csv` | bubbles re-snarled on the normalized graph (with `--reference-path`) |
| `<prefix>.panphorte.report.tsv` | one row per bubble |
| `<prefix>.panphorte.copies.tsv` | one row per haplotype and array, in approximate mode |
| `<prefix>.panphorte.rep_provenance.tsv` | one row per repeat-unit node created |
| `<prefix>.bandage_nodes.csv` | Bandage node colors |
| `<prefix>.panphorte.refine.tsv` | one row per input site, with `--refine` |

`panphorte.report.tsv` columns:

| column | meaning |
|--------|---------|
| `bubble_id` | the site, as numbered in the input CSV |
| `normalized` | `yes` if the site was folded |
| `unit_bp` | length of the repeat unit found |
| `paths_normalized` | how many haplotypes were rewritten onto the folded unit |
| `min_copies`, `max_copies` | fewest and most copies any haplotype carries |
| `interruptions_bp` | bases lying between copies rather than in them |
| `nodes_collapsed` | how many nodes the fold removed |
| `n_traversing` | haplotypes crossing the site at all |
| `n_motif_carriers` | haplotypes carrying the detected unit |
| `prevalence` | `n_motif_carriers` over `n_traversing`, the quantity `--min-array-prevalence` gates on |
| `n_motifs` | distinct units detected at the site |
| `copies_declined_partial_boundary`, `paths_with_partial_boundary` | copies, and haplotypes, refused because a copy could not be bounded |
| `status` | `normalized`, or why not: `below_prevalence`, `no_tandem_detected`, `no_seed`, `partial_boundary` |

`panphorte.copies.tsv` gives, per haplotype and array, its `copies`, `unit_bp`, `mean_identity`, the `orientations` of its copies, the `region_bp` they span and the input site they came from. `panphorte.rep_provenance.tsv` maps each created repeat-unit node to the site and `canonical_motif` it stands for, and gives the `copy_quantum` each traversal represents.

`panphorte.refine.tsv` columns (node ids are the folded graph's, before the sort renumbers them):

| column | meaning |
|--------|---------|
| `input_bubble_id`, `input_source`, `input_sink` | the site, as in the input CSV |
| `refine_status` | `refined`; `kept_original:not_simpler`; `not_folded`; or `skipped:<reason>` with `unfolded_dup`, `partial_paths`, `ambiguous_anchor_traversal`, `reference_not_traversing`, `not_a_clean_chain`, `no_carrier`, `segment_too_long`, `poa_work`, `too_diverse` or `poa_not_lossless` |
| `detail` | the node, path or sizes behind a skip |
| `paths`, `carriers` | haplotypes crossing the site, and those walking its repeat-unit node |
| `distinct_flanks_before`, `distinct_flanks_after` | distinct flank sequences aligned on each side |
| `nodes_*`, `steps_*`, `single_path_nodes_*` | before and after: nodes, steps (summed over haplotypes) and nodes used by one haplotype only, outside the array blocks |

### `--refine`

Folding leaves the bases around an array as the input graph spelled them, so one flank can be walked through different node chains by different haplotypes. `--refine` rebuilds them per folded site:

- each carrier's **array block** (its first to last repeat-unit step, with whatever lies between copies) is kept verbatim, so copy counts never change;
- the flanks before the block are aligned together (abPOA), and so are the flanks after it; they are never aligned across the block;
- haplotypes that walk no repeat-unit node at the site are left as folded;
- the rebuild is kept only if it has no more nodes and no more steps than before, and fewer of one; otherwise the site stays as folded (`kept_original:not_simpler`).

A site is skipped, untouched, when a flank walks a node twice (an unfolded duplication that realigning would flatten), when a path touches it without crossing it once from anchor to anchor, or when a path's traversal is ambiguous. Before anything is written, every haplotype must spell what it spelled after folding, walk the same repeat-unit nodes the same number of times, step only across existing links, and no rewritten haplotype's folded route may remain walkable; identical carriers must walk identical routes. Without `--refine` the output is unchanged.

## Limitations

- Below `--min-similarity 1.0` the collapse is lossy: divergent copies are folded onto one consensus unit, so a haplotype no longer spells its own sequence at that site.
- Bubble ids are reassigned when the graph is re-snarled under `--reference-path`, so an id in the output does not refer to the same site as that id in the input. The run reports the count change rather than implying the ids correspond.
- Exact folding seeds from runs of repeated node steps, so two byte-identical copies split at different node boundaries can be missed even though folding them would be lossless.
- Approximate detection looks for copies sharing an exact short seed, so a copy above `--min-similarity` whose substitutions are spread evenly can go undetected.
- Where one site carries two phase-rotated units, each becomes its own repeat-unit node. A consumer counting those nodes sees two independent duplications rather than one site's copy number unless it groups them on `canonical_motif` using the provenance table.
- `--refine` does not realign haplotypes that carry no folded copy. Where those haplotypes still carry the array as literal sequence (copies too divergent to fold), aligning them with the carriers' flanks fragments the site instead of simplifying it, so they keep their folded representation.
- `--allow-partial-boundary` is a safeguard override for an unobserved edge case, not a normal operating mode; using it can create a mixed folded/literal representation that downstream CN cannot interpret safely.

## Example

See the [walkthrough](../walkthrough.md) for this module in a full end-to-end run.
