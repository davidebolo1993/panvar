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

For the no-DUP and folded-DUP cases, each residual-segment position has its per-haplotype sequences collapsed to their exact distinct set and run through the aligner; one sequence, or none, is trivial.

Three guards bound the cost, and all three are measured over that distinct set rather than over every carrier, so adding an identical haplotype cannot change a decision without changing what the aligner is actually handed. `--max-poa-bp` skips a segment whose longest sequence exceeds it, `--max-walks` skips one carrying more distinct sequences than it allows, and `--max-poa-work` bounds the estimated alignment cost, the longest sequence against the total bases. A guard firing skips the region, and the reason is recorded rather than silently dropped.

The resulting column alignment is compacted into a minimal sub-graph: a run of columns that all and only the same haplotypes traverse becomes one node, with fresh non-colliding ids, mapping each distinct residual sequence to a rebuilt node path.

Each haplotype's interior is then reassembled as `residual-path, REP-block, residual-path, …` (the `REP` blocks verbatim) and spliced between the anchors: the old interior nodes are dropped, the `REP` nodes and self-loops kept, the rebuilt nodes added, each traversing path's `a…b` sub-walk rewritten, and any now-missing path adjacency added as a `0M` link. Because every rebuilt node spells exactly the bases of its MSA columns and `REP` runs are copied verbatim, each haplotype's spelled sequence over the region is unchanged: the rebuild is sequence-lossless, only the node structure changes.

### 4. Check losslessness

Losslessness is the module's central claim, so it is verified rather than argued. Every haplotype's spelled sequence is compared against what it spelled on the way in, and every consecutive step pair must be joined by a link that exists in the orientation walked. Both run before anything is written, so a violation cannot reach disk, and the output family is staged and committed only once they pass.

The link check is not redundant with the sequence one. A rewrite can leave every haplotype spelling the same bases while walking an adjacency the graph no longer contains, which no sequence comparison can see.

### 5. Re-sort, re-snarl and emit

Finally the graph is re-sorted and flipped along the reference and re-snarled with the cactus finder, applying `--resnarl-min-variant-bp` as its own interior-span filter rather than inheriting whatever produced the input sites. Bubble ids are reassigned by that decomposition. The outputs are written as `<prefix>.normalized.sorted.gfa`, `<prefix>.bubbles.csv`, `<prefix>.bandage_nodes.csv` (and `<prefix>.bandage_genes.csv` with `--gtf`). Running `call` on this graph re-derives the now-clean records; a folded DUP is byte-identical because its `REP` node, self-loop and per-haplotype multiplicity never changed.


## Worked trace

The steps below follow the five above, one for one. A region between anchors `S` and `E`. The reference interior is `X, R, Y`, where `R` is a folded repeat-unit node carrying a self-loop and `X`, `Y` are residual flanks. Two haplotypes:

| haplotype | interior walk | copies of `R` |
|-----------|---------------|---------------|
| reference | `X, R, Y` | 1 |
| H1 | `X, R, R, Y'` — one extra copy, and `Y'` is `Y` carrying a builder-split insertion and deletion | 2 |

1. Form regions. The bubbles between `S` and `E` share boundary nodes, so they fuse into one region with `S` and `E` as its outer anchors.

2. Classify by copy-number content. The interior holds a repeat-unit node with a self-loop, so this is the folded case: `R` is held fixed and only the residual flanks around it are re-aligned. Both haplotypes split at `R` into the same ordered skeleton, which is what the folded case requires. Residual segment 0 is `X` for both; segment 1 is `Y` for the reference and `Y'` for H1. The repeat block is `R` for the reference and `R, R` for H1.

3. Re-align and splice. Segment 0 aligns `X` against `X`: one sequence after collapsing to the distinct set, so nothing changes. Segment 1 aligns `Y` against `Y'`, and the builder's split insertion and deletion collapse into a single clean indel, or into sub-threshold variation if the net length is small. Neither segment approaches the guards. The interiors are reassembled as `x-nodes, R, y-nodes` and `x-nodes, R, R, y'-nodes`, with the repeat runs copied through untouched.

4. Check losslessness. Both haplotypes spell exactly what they spelled before — the rebuilt flank nodes carry the same bases, and the repeat run was never rewritten — and every consecutive step pair is joined by an existing link. The rewrite is accepted.

5. Re-sort, re-snarl and emit. The graph is sorted along the reference and decomposed again, and the family is committed. Calling on it now yields one duplication for H1, unchanged because `R`, its self-loop and H1's two traversals were never touched, and one clean indel in place of the spurious insertion-and-deletion pair at the flank.
