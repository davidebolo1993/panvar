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

**A consequence to be aware of:** with deletions counted, on a fully-typed panel nearly every haplotype supports nearly every bubble simply by crossing it. On C4 and LPA `path_support` is now uniformly 131 and 466 — the panel size. It is *traversal* support, not evidence that an alternate allele has multiple observations, and `--min-path-support` should be read that way. Distinct-allele and alternate-allele support are not yet reported separately.
