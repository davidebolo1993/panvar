# Module `refine` - algorithm

Mechanism for the `refine` module. For usage/flags see [modules/refine.md](../modules/refine.md); References in [references.md](../references.md#refine).

## Terms

- **region** — one or more bubbles fused into a single `source…sink` span and refined as a unit. Bubbles fuse when they share a boundary node; a region is refined only if the fused set reduces to exactly two outer anchor nodes (a clean linear chain).
- **residual segment** — a maximal run of non-`REP` interior steps between or around the `REP` blocks; the only part POA re-aligns.
- **`REP` node** — panphorte's folded tandem unit, an interior node carrying an L self-edge (see [panphorte](panphorte.md)). Held fixed during refine.

## How it works

Operate on the panphorte normalized graph. Precompute the set of `REP` nodes (nodes with an L self-edge).

Form regions by fusing bubbles that share a boundary node. `panphorte` emits the bubble decomposition without merging nearby bubbles, but directly-adjacent snarls in a chain already share their boundary node — the sink of one is the source of the next. `refine` fuses exactly those. For a fused set, the two endpoints that appear once across the pooled `source`/`sink` are the outer anchors `a`, `b`; every other node is interior (the shared boundary nodes become interior, so POA spans across the original bubble boundaries). A set that does not reduce to exactly two anchors is not a clean linear chain and is skipped. This endpoint-based fusion is independent of `bubble --merge-nearby-bp`, which merges by graph distance.

For each region, take the oriented steps between `a` and `b` for every path (the reference must traverse the region), and act by copy-number content:

- No DUP — the interior has no `REP` node and no path revisits a non-`REP` interior node. Re-align the whole interior between `a` and `b`.
- Folded DUP — the interior has a `REP` node. Hold it fixed: split each path's interior at every `REP` block into residual segments and verbatim `REP×n` runs, re-align only the residual segments, and copy the `REP` runs through unchanged (per-haplotype copy count and orientation preserved). All traversing paths must share the same ordered `REP` skeleton, otherwise the region is skipped.
- Unfolded DUP — the reference or any haplotype revisits a non-`REP` interior node two or more times. Skip the whole region: POA would linearize the copies and destroy the copy-number signal `call` reconstructs from the revisits. This is the case panphorte's prevalence gate leaves unfolded (a private duplication, or a paralog cluster collapsed onto shared nodes).

Re-align and splice (the no-DUP and folded-DUP cases). For each residual-segment position, collapse the per-haplotype sequences to their exact distinct set and run abPOA over them (one sequence, or none, is handled trivially). Skip the region if a segment's median length exceeds `--max-poa-bp` or its distinct-walk count exceeds `--max-walks`. Compact the resulting column MSA into a minimal sub-graph: a run of columns that all and only the same haplotypes traverse becomes one node, with fresh non-colliding ids. This maps each distinct residual sequence to a rebuilt node path.

Reassemble each haplotype's interior as `residual-path, REP-block, residual-path, …` (the `REP` blocks verbatim) and splice it between the anchors: drop the old interior nodes, keep the `REP` nodes and their self-loops, add the rebuilt nodes, rewrite each traversing path's `a…b` sub-walk, and add any path adjacency the graph now lacks as a `0M` link. Because every rebuilt node spells exactly the bases of its MSA columns and `REP` runs are copied verbatim, each haplotype's spelled sequence over the region is unchanged: the rebuild is sequence-lossless, only the node structure changes.

Finally, re-sort and flip along the reference and re-snarl with the cactus finder (the same step panphorte runs), then write `<prefix>.normalized.sorted.gfa`, `<prefix>.bubbles.csv`, `<prefix>.bandage_nodes.csv` (and `<prefix>.bandage_genes.csv` with `--gtf`). Running `call` on this graph re-derives the now-clean records; a folded DUP is byte-identical because its `REP` node, self-loop and per-haplotype multiplicity never changed.

The rebuild is coordinate-free — it works on step tokens and node-set membership — so a reverse-oriented region needs no special treatment: each haplotype's residual is spelled and reassembled in its own path direction, so its sequence is preserved regardless.

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
