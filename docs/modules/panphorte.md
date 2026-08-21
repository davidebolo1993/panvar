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

## Limitations

- Below `--min-similarity 1.0` the collapse is lossy: divergent copies are folded onto one consensus unit, so a haplotype no longer spells its own sequence at that site.
- Bubble ids are reassigned when the graph is re-snarled under `--reference-path`, so an id in the output does not refer to the same site as that id in the input. The run reports the count change rather than implying the ids correspond.
- Exact folding seeds from runs of repeated node steps, so two byte-identical copies split at different node boundaries can be missed even though folding them would be lossless.
- Approximate detection looks for copies sharing an exact short seed, so a copy above `--min-similarity` whose substitutions are spread evenly can go undetected.
- Where one site carries two phase-rotated units, each becomes its own repeat-unit node. A consumer counting those nodes sees two independent duplications rather than one site's copy number unless it groups them on `canonical_motif` using the provenance table.
- `--allow-partial-boundary` is a safeguard override for an unobserved edge case, not a normal operating mode; using it can create a mixed folded/literal representation that downstream CN cannot interpret safely.

## Example

See the [walkthrough](../walkthrough.md) for this module in a full end-to-end run.
