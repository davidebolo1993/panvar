# Module `refine` - algorithm

Mechanism for the `refine` module. For usage/flags see [modules/refine.md](../modules/refine.md); References in [references.md](../references.md#refine).

`refine` POA-realigns the interior of bubbles on the panphorte normalized graph to remove graph-builder alignment artifacts (a spurious insertion-plus-deletion pair where one clean indel belongs), without touching copy number. It works on a region: one or more bubbles fused into a single `source→sink` span and refined as a unit. A `REP` node — panphorte's folded tandem unit, an interior node carrying an L self-edge — is held fixed throughout, and only the residual segments (the maximal runs of non-`REP` interior steps between and around the `REP` blocks) are re-aligned.

## How it works

### 1. Form regions

Working on the panphorte graph, the set of `REP` nodes (those with an L self-edge) is precomputed, then bubbles that share a boundary node are fused into regions. `panphorte` emits its bubbles without merging nearby ones, but directly-adjacent snarls in a chain already share a boundary node — the sink of one is the source of the next — and `refine` fuses exactly those. For a fused set, the two endpoints that appear once across the pooled `source`/`sink` are the outer anchors `a`, `b`; every other node becomes interior (the shared boundaries included, so POA spans across the original bubble boundaries). A set that does not reduce to exactly two anchors is not a clean linear chain and is skipped. This endpoint-based fusion is independent of `bubble --merge-nearby-bp`, which merges by graph distance.

### 2. Classify by copy-number content

For each region, the oriented steps between `a` and `b` are taken for every path (the reference must traverse the region), and the region is handled by its copy-number content:

- No DUP — the interior has no `REP` node and no path revisits a non-`REP` interior node. The whole interior between `a` and `b` is re-aligned.
- Folded DUP — the interior has a `REP` node, held fixed: each path's interior is split at every `REP` block into residual segments and verbatim `REP×n` runs; only the residual segments are re-aligned, and the `REP` runs are copied through unchanged (per-haplotype copy count and orientation preserved). All traversing paths must share the same ordered `REP` skeleton, or the region is skipped.
- Unfolded DUP — the reference or any haplotype revisits a non-`REP` interior node two or more times. The whole region is skipped: POA would linearize the copies and destroy the copy-number signal `call` reconstructs from the revisits. This is what panphorte's prevalence gate leaves unfolded (a private duplication, or a paralog cluster collapsed onto shared nodes).

### 3. Re-align and splice

For the no-DUP and folded-DUP cases, each residual-segment position has its per-haplotype sequences collapsed to their exact distinct set and run through abPOA (one sequence, or none, is trivial). The region is skipped if a segment's median length exceeds `--max-poa-bp` or its distinct-walk count exceeds `--max-walks`. The resulting column MSA is compacted into a minimal sub-graph — a run of columns that all and only the same haplotypes traverse becomes one node, with fresh non-colliding ids — mapping each distinct residual sequence to a rebuilt node path.

Each haplotype's interior is then reassembled as `residual-path, REP-block, residual-path, …` (the `REP` blocks verbatim) and spliced between the anchors: the old interior nodes are dropped, the `REP` nodes and self-loops kept, the rebuilt nodes added, each traversing path's `a…b` sub-walk rewritten, and any now-missing path adjacency added as a `0M` link. Because every rebuilt node spells exactly the bases of its MSA columns and `REP` runs are copied verbatim, each haplotype's spelled sequence over the region is unchanged: the rebuild is sequence-lossless, only the node structure changes.

### 4. Re-sort and re-snarl

Finally the graph is re-sorted and flipped along the reference and re-snarled with the cactus finder (the same step panphorte runs), then written as `<prefix>.normalized.sorted.gfa`, `<prefix>.bubbles.csv`, `<prefix>.bandage_nodes.csv` (and `<prefix>.bandage_genes.csv` with `--gtf`). Running `call` on this graph re-derives the now-clean records; a folded DUP is byte-identical because its `REP` node, self-loop and per-haplotype multiplicity never changed.


## Worked trace

A region between anchors `S` and `E`. The reference interior is `X, R, Y` where `R` is a `REP` node (self-loop, 100 bp unit) and `X`, `Y` are residual flanks. Two haplotypes:

| haplotype | interior walk | copies of R |
|-----------|---------------|-------------|
| reference | `X, R, Y` | 1 |
| H1 | `X, R, R, Y′` (one extra copy; `Y′` is `Y` with a spurious builder-split INS+DEL) | 2 |

Split at `R`: the skeleton is `(R)` for both, so it is consistent. Residual segment 0 is `[X]` for both; segment 1 is `[Y]` (reference) and `[Y′]` (H1). The `REP` block is `[R]` (reference) and `[R, R]` (H1), copied verbatim.

POA the residuals: segment 0 is `X` versus `X`, one node, no change. Segment 1 aligns `Y` against `Y′`; the builder's split insertion and deletion collapse to a single clean indel, or, if the net length is below `--min-sv-bp`, to sub-threshold variation.

Reassemble: the reference interior becomes `x-nodes, R, y-nodes`; H1 becomes `x-nodes, R, R, y′-nodes`. The `R, R` run is untouched.

Result: `call` on the refined graph emits one DUP at this region for H1 (`REF_CN=1`, `CN=2`), byte-identical to before because `R`, its self-loop and H1's two traversals are unchanged, and one clean indel (or nothing) for the `Y′` flank instead of the spurious INS+DEL pair. Benchmark identity for H1 is unchanged, since the same sequence is reconstructed from fewer, cleaner records.
