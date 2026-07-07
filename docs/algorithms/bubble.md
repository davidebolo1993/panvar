# Module `bubble` - algorithm

Mechanism for the `bubble` module. For usage/flags see [modules/bubble.md](../modules/bubble.md); References in [references.md](../references.md#bubble).

## Terms

- **snarl** — a pair of boundary nodes whose removal separates an internal subgraph from the rest, defined on the bidirected graph; it can contain cycles and inversions.
- **superbubble** — the acyclic special case: a single-source/single-sink subgraph that is directed and acyclic (no cycle, no inversion).

`bubble` consumes snarls because pangenome variation often lives in the cyclic sites a superbubble omits; `--superbubbles` restricts to the acyclic subset. Internally it mirrors `vg`'s cactus/3-edge-connected decomposition, so the default top-level snarls should match `vg`'s.

## Algorithm overview

For each top-level snarl candidate:

1. get `(source, sink)` boundary from the internal cactus finder (or `--snarls-in` JSONL);
2. find path intervals crossing `source→sink` (or `sink→source`, then canonicalize);
3. collect internal nodes seen between boundaries;
4. measure number of crossing paths and internal bp support across them;
5. apply base filters: `--min-variant-bp` (keep a candidate only if some supporting path has internal span ≥ N bp — the size filter, default 50), `--min-path-support` (keep only candidates with ≥ N supporting paths, default 0);
6. optionally merge nearby surviving bubbles via `--merge-nearby-bp`: consecutive survivors merge when the shortest-path distance from the previous `sink` to the next `source` is ≤ N bp; assign bubble IDs.


## Worked trace 

Input graph (one bubble; `+` = forward strand):

```text
S  1  T                       #flank
S  2  GG                      #source boundary
S  3  AC                      #interior allele "ref"  (2 bp)
S  4  ATG                     #interior allele "alt"  (3 bp)
S  5  CC                      #sink boundary
S  6  A                       #flank
L  1  +  2  +  0M
L  2  +  3  +  0M
L  3  +  5  +  0M
L  2  +  4  +  0M
L  4  +  5  +  0M
L  5  +  6  +  0M
L  1  +  6  +  0M             #bypass edge (hC's deletion)
P  ref  1+,2+,3+,5+,6+  *
P  hA   1+,2+,4+,5+,6+  *
P  hB   1+,2+,3+,5+,6+  *
P  hC   1+,6+           *     #deletion that bypasses the locus
```

1. The cactus finder returns boundary (`source`, `sink`) = (2, 5) — removing 2 and 5 separates `{3,4}`.
2. Per path, find the interval crossing `2 → 5` and collect interiors: `ref → {3}`, `hA → {4}`, `hB → {3}`; `hC` goes `1 → 6`, never visiting both boundaries (it does not support this snarl).
3. Aggregate: `inside_nodes = {3,4}`, `path_support = 3` (ref, hA, hB; `< 4` total because hC bypasses — expected), `min_inside_bp = 2` (`AC`), `max_inside_bp = 3` (`ATG`).

Resulting `bubbles.csv` row:

```text
bubble_id  source  sink  inside_node_count  path_support  min_inside_bp  max_inside_bp  inside_nodes
1          2       5     2                  3             2              3              3,4
```
