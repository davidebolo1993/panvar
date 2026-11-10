# Module `panphorte`

CLI: `panvar panphorte`

## What it does

Folds tandem-repeat sites onto their repeat unit, so a haplotype's copy number is the number of times its walk passes through the unit, instead of a pile of insertions of different sizes.

- Exact mode (default) folds only identical copies; approximate mode (`--min-similarity` below 1) folds near-identical copies, keeping each copy's differences. With `--lossy-units`, approximate mode keeps only the reference's differences, for a simpler graph.
- Only arrays most haplotypes carry are folded; a rare duplication is left for `call` to report as an event.
- With `-r`, the folded graph is sorted and its sites found again.

## How it works

The idea of folding a tandem array onto its unit comes from [Panphorte](https://github.com/GenoGra/Panphorte); copies are aligned with [edlib](https://github.com/Martinsos/edlib), and `--refine` realigns repeat flanks with [abPOA](https://github.com/yangao07/abPOA). As in `rebuild`'s fold, the copies of a unit of `--min-unit-bp` or more are made to walk one shared unit (a unit graph).

1. The unit is found per site. In exact mode it is the shortest block of nodes a haplotype walks twice in a row; a node holding whole copies (two, or a copy and a soft-masked copy) counts each. In approximate mode each haplotype proposes the units of its arrays and the most proposed wins; it is then rotated to start where the copies start. Copies every haplotype walks just outside the site, as shared backbone, are also part of its array, and the site is widened over them.
2. Copies are found along each haplotype: identical node blocks in exact mode, alignments at `--min-similarity` or more in approximate mode. Bases between copies stay in the walk.
3. The site is folded only if at least `--min-array-prevalence` of its haplotypes carry `--min-copies` copies or more.
4. Every haplotype is rewritten onto the unit. In exact mode the unit is one node, stored as its copies read along the reference, so copies on the reference's strand walk it forward even where the reference has none; in approximate mode it is a unit graph: the unit cut wherever a copy differs, each copy walking the shared pieces and its own differences, so every haplotype spells exactly what it did. Copies are counted at an anchor, a piece no copy changes. With `--lossy-units`, the unit graph holds only the reference's differences: the reference still spells exactly, so no coordinate moves, and every other haplotype walks the plain unit once per copy, keeping its copy number but not its own differences between copies. Without `-r`, no haplotype keeps its sequence.
5. With `--refine`, the flanks around each folded array are realigned (abPOA), and kept only if the site gets simpler.
6. With `-r`, the graph is sorted and its sites found again, as `bubble` does but without merging nearby sites, so a folded site keeps its own records. Site ids change.

## Required inputs

- `-i, --gfa <graph.gfa>`: the sorted GFA from `bubble` (a regular file, not a pipe)
- `-b, --bubble-prefix <prefix>` (reads `<prefix>.bubbles.csv`) or `-c, --bubbles-csv <path>`: the sites, made from this graph; a table from another graph is refused
- `-o, --out-prefix <prefix>`: output prefix

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-r, --reference-path <name>` | sort the folded graph and find its sites again, ready for `call`. Without it the output is for viewing only | — |
| `--bubble-id <N>` | only this site; repeatable | all |
| `--min-similarity <f>` | identity for a block to count as a copy: `1` identical only; below 1, approximate mode; above 1 a percentage | `1.0` |
| `--lossy-units` | approximate mode: keep only the reference's differences between copies; other haplotypes keep their copy number but not their own differences | off |
| `--min-unit-bp <N>` | shortest unit folded, at least 1 (approximate mode: at least 16) | `50` |
| `--min-copies <N>` | copies that make an array | `2` |
| `--min-array-prevalence <f>` | share of the site's haplotypes that must carry an array | `0.5` |
| `--max-interruption-frac <f>` | share of an array's bases that may lie between copies | `0.25` |
| `--refine` | realign the flanks of each folded array (needs `-r`); the `refine` module does this and more | off |
| `--refine-max-poa-bp`, `--refine-max-walks`, `--refine-max-poa-work <N>` | limits for `--refine` | `5000`, `500`, `0` |
| `--resnarl-min-variant-bp <N>`, `--merge-nearby-bp <N>`, `--no-fuse-ambiguous`, `--no-flip` | with `-r`: `bubble`'s options for the new sites | `50`, `0`, off, off |
| `--gtf <path>` | with `-r`: genes per node of the folded graph | — |
| `-t, --threads <N>` | worker threads (`0` = all cores) | `0` |
| `-q, --quiet` | no progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.normalized.sorted.gfa`, `<prefix>.bubbles.csv` | with `-r`: the folded graph, sorted, and its sites, for `call` |
| `<prefix>.normalized.gfa` | without `-r`: the folded graph, unsorted |
| `<prefix>.panphorte.report.tsv` | one row per input site: folded or why not |
| `<prefix>.panphorte.unit_graph.tsv` | one row per fold: its unit, nodes and anchor; `call` counts copies with it |
| `<prefix>.panphorte.copies.tsv` | approximate mode: copies per crossing of a site |
| `<prefix>.panphorte.refine.tsv` | with `--refine`: one row per input site |
| `<prefix>.bandage_nodes.csv`, `<prefix>.bandage_genes.csv` | with `-r`: Bandage colours, and genes with `--gtf` |

`panphorte.report.tsv`, one row per input site:

| column | meaning |
|--------|---------|
| `bubble_id` | the site, as numbered in the input table |
| `normalized` | `yes` if the site was folded |
| `unit_bp` | the unit's length; for a site not folded, the unit tried (0 if none; for `short_motif`, the unit made of the short motif) |
| `paths_normalized` | haplotypes rewritten onto the unit |
| `min_copies`, `max_copies` | fewest and most copies a haplotype carries, among those with one |
| `interruptions_bp` | bases between folded copies (0 if not folded) |
| `nodes_collapsed` | nodes the fold removed |
| `n_traversing` | haplotypes crossing the site |
| `n_motif_carriers` | haplotypes carrying an array of the unit |
| `prevalence` | `n_motif_carriers` / `n_traversing`, what `--min-array-prevalence` gates |
| `n_motifs` | exact mode: distinct motifs found; approximate mode: 1 when a unit was chosen |
| `status` | `normalized`, or why not (below) |

Statuses for a site not folded (`n_traversing`, `n_motif_carriers` and `prevalence` are `.` when nothing was counted: the site was declined, when `n_motifs` is `.` too, or had no seed):

- `below_prevalence`: too few haplotypes carry an array;
- `short_motif`: the unit is copies of a motif shorter than `--min-unit-bp`;
- `no_tandem_detected` (exact mode), `no_seed` (approximate mode): no unit found;
- `unit_under_16bp`: approximate mode cannot find copies of a unit this short;
- `declined_shared_route`: after folding, other walks still spelled the old route;
- `no_copies_rewritten`: no walk could be rewritten onto the unit.

`panphorte.unit_graph.tsv`, one row per fold, read by `call` to count copies:

| column | meaning |
|--------|---------|
| `input_bubble_id` | the site, as numbered in the input table |
| `unit_bp`, `unit` | the unit as folded |
| `canonical_motif`, `copy_quantum` | the unit's shortest repeating motif, and how many of it one copy holds |
| `anchor_pos`, `anchor_nodes` | a unit position no copy deletes, and the nodes covering it: each copy passes exactly one of them once, so copies are counted there |
| `n_nodes`, `nodes` | the unit graph's nodes |

`panphorte.copies.tsv` (approximate mode), one row per crossing, so a haplotype crossing a site twice has two rows: `path_name`, `sample`, `input_bubble_id`, `copies`, `unit_bp`, `orientations` (of its copies), `mean_identity` (to the unit), `region_bp` (the bases its copies span), `input_from_node`, `input_to_node`.

### Visualization

With `-r`, the folded graph opens in [BandageNG](https://github.com/asl/BandageNG) coloured by site, as `bubble`'s: File › Load CSV data with `<prefix>.bandage_nodes.csv` (and `<prefix>.bandage_genes.csv` for the genes). To see what a fold changed, compare the graph before and after with the scripts in [rebuild's Visualization](rebuild.md#visualization), and draw the repeat site before and after with `inspect` ([inspect's Visualization](inspect.md#visualization)):

```
python3 scripts/compare_graphs.py -i <input.gfa> -i <prefix>.normalized.sorted.gfa -l input,folded -o <compare prefix>
Rscript scripts/plot_compare_graphs.R --prefix <compare prefix>
```

## Limitations

- `call` counts copies on a unit graph only through `.panphorte.unit_graph.tsv` (from its `-b` prefix or `--panphorte-prefix`); from the graph alone it sees only one-node folds.
- Exact mode needs identical copies cut into the same nodes; real copies usually differ, so approximate mode is the one for real data.
- Units under 16 bp (short tandem repeats) are not folded in approximate mode, whose copies come from 16-base seeds. Exact mode folds them with `--min-unit-bp` lowered only where each copy is its own node, not where a node holds a whole allele, as most graphs store one.
- A copy within `--min-similarity` whose differences are spread evenly can be missed, since copies are found from short exact seeds.

## Example

`scripts/panphorte_examples.sh` folds the graphs in `tests/examples/bubble` (the tandem-repeat locus also with `--lossy-units`), inspects its repeat site before and after each fold, and compares the graphs, in `tests/examples/panphorte`.
