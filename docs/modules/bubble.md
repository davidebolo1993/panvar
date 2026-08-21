# Module `bubble`

CLI: `panvar bubble`

## What it does

Turns a pangenome graph (GFA — Graphical Fragment Assembly — from [pggb](https://github.com/pangenome/pggb) or [minigraph](https://github.com/lh3/minigraph)) into bubble sites for the downstream modules. It:
- sorts/flips the graph along the reference;
- finds snarls internally — a snarl being a boundary-node pair whose removal isolates an interior subgraph (a vendored cactus/3-edge-connected decomposition matching [vg](https://github.com/vgteam/vg) `snarls`);
- determines each bubble's interior, and scores how many paths cross it, what those crossings contain, and how long the interior is;
- orders each boundary pair along the reference, then applies size and support filters and optional merging;
- writes the sorted GFA, a bubble CSV and [Bandage](https://github.com/asl/BandageNG)-ready visualization files.

A snarl can contain cycles and inversions, which is where much pangenome variation sits, so snarls are the default; `--superbubbles` restricts to the acyclic subset. Algorithm and worked trace: [algorithms/bubble.md](../algorithms/bubble.md).

## Required inputs

- `-i, --gfa <graph.gfa>` — the input graph
- `-r, --reference-path <name>` — reference path name (or unique, case-insensitive substring). Orders the internal sort/flip and snarl finder. Not needed with `--snarls-in`, but see Limitations.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `-o, --out-prefix <p>` | output prefix | `bubble_calls` |
| `-s, --superbubbles` | emit only acyclic superbubbles (single-source/single-sink, no cycle or inversion) instead of all snarls | off (all snarls) |
| `--min-variant-bp <N>` | keep a bubble only if some path's interior span ≥ N (`0` = off). Alias: `--min-interior-bp` | `50` |
| `--max-variant-bp <N>` | drop a bubble if any path's interior span > N (`0` = off); excludes hypervariable tangles. Alias: `--max-interior-bp` | `0` (off) |
| `--min-path-support <N>` | keep only bubbles crossed by ≥ N paths | `0` (off) |
| `--min-alt-support <N>` | keep only bubbles whose best-supported non-reference allele has ≥ N paths | `0` (off) |
| `--merge-nearby-bp <N>` | merge consecutive bubbles whose facing boundaries are ≤ N bp apart along the reference; filters are re-applied afterwards | `0` (off) |
| `--no-flip` | sort but don't reorient nodes to the reference strand | off |
| `--snarls-in <path>` | use an external `vg snarls` JSONL (JSON Lines; graph used as-is, no sort) | — |
| `--emit-snarls-jsonl <path>` | also write the internal snarls as `vg`-style JSONL | — |
| `--sorted-gfa-out <path>` | explicit path for the sorted GFA | `<prefix>.sorted.gfa` |
| `--bubbles-csv <path>` / `--bandage-csv <path>` | explicit output paths | derived from `-o` |
| `--snarl-debug-tsv <path>` | write per-snarl-candidate diagnostics (`candidate_id, source, sink, inside_node_count, n_paths, min_inside_bp, long_path_support, inversion_signal, accepted`; merged bubbles get an extra `accepted=1` row) | — |
| `--gtf <path>` | project the genes of a reference-coordinate GTF (Gene Transfer Format) onto reference nodes (`<prefix>.bandage_genes.csv`, needs a [PanSN](https://github.com/pangenome/PanSN-spec) (Pangenome Sequence Naming) `--reference-path`) | — |
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
| `bubble_id` | unique id for the bubble site; increases along the reference |
| `source`, `sink` | the two boundary node ids, ordered so `source` is the one the reference reaches first |
| `source_orient`, `sink_orient` | the strand on which the reference reads each boundary (`+` forward, `-` reverse) |
| `inside_node_count` | number of interior nodes (strictly between the boundaries) |
| `total_node_count` | interior nodes plus the two boundary nodes |
| `path_support` | how many paths cross the bubble, that is visit both boundaries |
| `distinct_alleles` | how many different `source → sink` walks those paths take |
| `ref_allele_support` | how many take the same walk as the reference |
| `alt_allele_support_max` | paths carrying the best-supported walk that is not the reference's |
| `alt_allele_support_min` | paths carrying the least-supported such walk |
| `min_inside_bp`, `max_inside_bp` | smallest and largest interior span (bp) across the crossing paths |
| `inside_nodes` | the interior node ids (`;`-separated) |

## Limitations

- `--min-variant-bp` and `--max-variant-bp` measure the span between the boundaries, not how much an allele differs from the reference. A small edit inside a long allele is filtered by the length of the allele carrying it. `--min-interior-bp` and `--max-interior-bp` are aliases that name this directly.
- `path_support` counts any path that crosses the site, including one carrying the reference allele, so on a densely typed panel it approaches the panel size. `--min-alt-support` filters on the alternate alleles instead.
- Only top-level snarls are emitted. A large tangle is reported as one site rather than decomposed, so variants nested inside it do not appear separately.
- `--snarls-in` without `--reference-path` is diagnostic-only: imported boundaries carry no reference order, coordinate merging is skipped and no sorted graph is written. Do not pass that output to `call`, which anchors coordinates on the source boundary.
- Interior discovery from the graph is bounded. On a very large or poorly separated site it uses only the interior the paths reveal, and says so on stderr.

## Example

See the [walkthrough](../walkthrough.md) for this module in a full end-to-end run.
