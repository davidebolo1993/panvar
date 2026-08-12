# Module `bubble`

CLI: `panvar bubble`

## What it does

Turns a pangenome graph (GFA — Graphical Fragment Assembly — from [pggb](https://github.com/pangenome/pggb)) into bubble sites for the downstream modules. It:
- sorts/flips the graph along the reference;
- finds snarls internally — a snarl being a boundary-node pair whose removal isolates an interior subgraph (a vendored cactus/3-edge-connected decomposition matching [vg](https://github.com/vgteam/vg) `snarls`);
- infers each bubble's internal nodes from path intervals, and scores path support and internal span;
- applies size/support filters;
- writes the sorted GFA, a bubble CSV and [Bandage](https://github.com/asl/BandageNG)-ready visualization files.

A snarl can contain cycles and inversions, which is where much pangenome variation sits, so snarls are the default; `--superbubbles` restricts to the acyclic subset. Algorithm and worked trace: [algorithms/bubble.md](../algorithms/bubble.md).

## Required inputs

- `-i, --gfa <graph.gfa>` — `pggb` GFA
- `-r, --reference-path <name>` — reference path name (or unique, case-insensitive substring). Orders the internal sort/flip and snarl finder. Not needed with `--snarls-in`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-prefix <p>` | output prefix | `bubble_calls` |
| `-s, --superbubbles` | emit only acyclic superbubbles (single-source/single-sink, no cycle or inversion) instead of all snarls | off (all snarls) |
| `--min-variant-bp <N>` | keep a bubble only if some path's internal span ≥ N (`0` = off) | `50` |
| `--max-variant-bp <N>` | largest variant to keep: drop a bubble if any path's internal span > N (`0` = off). | `0` (off) |
| `--min-path-support <N>` | keep only bubbles crossed by ≥ N paths | `0` (off) |
| `--merge-nearby-bp <N>` | merge consecutive bubbles ≤ N bp apart (after filters) | `0` (off) |
| `--no-flip` | sort but don't reorient nodes to the reference strand | off |
| `--snarls-in <path>` | use an external `vg snarls` JSONL (JSON Lines; graph used as-is, no sort) | — |
| `--emit-snarls-jsonl <path>` | also write the internal snarls as `vg`-style JSONL | — |
| `--snarl-debug-tsv <path>` | write per-snarl-candidate diagnostics (`candidate_id, source, sink, inside_node_count, n_paths, min_inside_bp, long_path_support, inversion_signal, accepted`; merged bubbles get an extra `accepted=1` row) | — |
| `--gtf <path>` | project the genes of a reference-coordinate GTF (Gene Transfer Format) onto reference nodes  (`<prefix>.bandage_genes.csv`, needs a [PanSN](https://github.com/pangenome/PanSN-spec) (Pangenome Sequence Naming) `--reference-path`) | — |
| `-q, --quiet` | disable the progress bar | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.bubbles.csv` | the bubble table consumed by `inspect`/`panphorte`/`call` |
| `<prefix>.sorted.gfa` | reference-sorted/flipped graph — the input for downstream `panphorte`/`call` |
| `<prefix>.bandage_nodes.csv` | Bandage node colors (blue = candidate context, red = retained bubbles) |
| `<prefix>.bandage_genes.csv` | (with `--gtf`) `Name,Colour,Gene` per bubble |

`bubbles.csv` columns:

| column | meaning |
|--------|---------|
| `bubble_id` | unique id for the bubble site |
| `source`, `sink` | the two boundary node ids that delimit the bubble |
| `inside_node_count` | number of interior nodes (strictly between the boundaries) |
| `total_node_count` | interior nodes plus the two boundary nodes |
| `path_support` | how many paths cross the bubble (visit both boundaries) |
| `min_inside_bp`, `max_inside_bp` | smallest/largest interior span (bp) across the supporting paths |
| `inside_nodes` | the interior node ids (`;`-separated) |

## Example

See the [LPA walkthrough](../walkthrough.md) for this module in a full end-to-end run.

## Path support counts pure-deletion alleles

A bubble's supporting paths are found with the shared `bubble_steps()` walk extractor — the same one `call`, `describe` and `genotype` use. `bubbles.cpp` previously carried a second, structurally identical path index with its own interval search, and the two had drifted: the local one required at least one **declared interior node** between the boundaries, so a direct `source → sink` allele — a pure deletion, and usually the most interesting allele in the bubble — was invisible to bubble scoring while every other module handled it correctly.

On the minimal case (paths `1,2,3` and `1,3`) that reported `path_support=1, min_inside_bp=4`; it now reports `path_support=2, min_inside_bp=0`, and the bubble survives `--min-path-support 2` instead of being dropped.

**A consequence, and what is reported instead.** With deletions counted, on a fully-typed panel nearly every haplotype supports nearly every bubble simply by crossing it. On C4 and LPA `path_support` is uniformly 131 and 466 — the panel size. It is *traversal* support and says nothing about any particular allele, so the CSV also reports what those traversals contain:

| column | meaning |
|--------|---------|
| `path_support` | paths that cross the bubble at all — traversal support |
| `distinct_alleles` | distinct `source → sink` walks |
| `ref_allele_support` | paths taking the reference's walk |
| `alt_allele_support_max` | the best-supported non-reference walk |
| `alt_allele_support_min` | the least-supported non-reference walk |

`--min-alt-support <N>` filters on the fourth of these, which is what a support filter is usually wanted for. The difference is not subtle. On C4, where `path_support` is 131 for every bubble:

```
--min-alt-support   0 -> 5     --min-path-support 100 -> 5
--min-alt-support   5 -> 4     --min-path-support 131 -> 5
--min-alt-support  20 -> 2     --min-path-support 200 -> 0
--min-alt-support  90 -> 0
```

One filter is a gradient; the other is a step at the panel size. C4 has a bubble with **127 distinct alleles across 131 paths** (`alt_allele_support_max = 2`, so no alternate is replicated) sitting beside one with `alt_allele_support_max = 80` — `--min-path-support` cannot tell them apart.

## What the bp filters measure

`--min-variant-bp` and `--max-variant-bp` measure the **interior span** between the boundaries — the summed length of the interior nodes a path traverses — not the size of the difference between alleles. A 1 bp substitution inside a 1 kb allele has a 1 kb span and passes `--min-variant-bp 50`. `--min-interior-bp` and `--max-interior-bp` are accepted as aliases that say so; the CSV columns `min_inside_bp` / `max_inside_bp` were always named correctly. Computing true allele-to-reference divergence is not implemented.

## Boundaries carry reference order

A cactus snarl is an **unordered** pair of boundaries — nothing in the decomposition says which is left. But every consumer reads them as an interval in reference order: `call` anchors coordinates on the source, and merging joins one bubble's sink to the next bubble's source. `bubble` therefore orients each pair so `source` is the reference-left boundary, before ids are assigned, so `bubble_id` also increases along the reference. A bubble the reference does not traverse has no order to take and is left alone.

This was not a rare edge case. Measured on the reference loci, the fraction of bubbles stored **reversed** against the reference was:

| locus | reversed |
|-------|----------|
| C4 | 5/5 |
| GSTM1 | 4/4 |
| ACOT | 9/9 |
| LPA | 0/9 |

The bubble *set* is unchanged by orienting — same endpoints, same interiors — but what `call` makes of it is not. On C4 the same five bubbles yield 19 records instead of 15, and the reconstruction benchmark says the extra resolution is real rather than noise: the region VCF closes **9.2%** of the baseline-to-graph gap where it closed **2.5%** before, with residual falling from 1,376,483 to 1,282,382 bases. The events are re-partitioned, not merely multiplied — total allele count falls from 347 to 318 as a few coarse events become several finer ones.

## `--superbubbles` tests the graph, not the paths

A superbubble is a snarl whose interior is a directed acyclic graph. Cyclicity used to be inferred entirely from what the stored paths happen to do — an interior self-loop, a path revisiting an interior node, a node used in both orientations. Those are real signals, but none of them looks at the graph: **an interior cycle that no stored path walks was reported as acyclic and survived `--superbubbles`**, which exists precisely to exclude it.

The interior is now searched directly for a directed cycle, over oriented **handles** rather than node names — leaving a node forward departs its end, leaving it reversed departs its start, so collapsing to node names would call `A+ → B+ → A-` a cycle when it is an ordinary traversal of two distinct handles. The path-derived signals are kept as additional evidence.

On C4 and LPA this changes nothing (3 and 7 superbubbles before and after): the path signals already caught everything cyclic there. It closes a hole those graphs happen not to exercise.

## An empty result is a result

A graph with no snarl — a single node, a linear path — emits an empty bubble table and exits 0. It previously raised `no snarls — provide ...`, because an empty result from the internal finder was indistinguishable from no finder having been supplied.
