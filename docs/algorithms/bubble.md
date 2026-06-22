# bubble — algorithm & worked example

Mechanism and a hand-traced example for **Module 1**. For usage/flags see [modules/bubble.md](../modules/bubble.md);
for citations see [references.md](../references.md#bubble).

## Terms

- **snarl** — a pair of boundary nodes whose removal separates an internal subgraph from the rest, defined
  on the **bidirected** graph; it can contain cycles and inversions. What `vg snarls` emits.
- **superbubble** — the acyclic special case: a single-source/single-sink subgraph that is directed and
  acyclic (no cycle, no inversion). `vg`'s "ultrabubble" ≈ superbubble.

`bubble` consumes **snarls** by default because pangenome variation (tandem expansions, inversions) often
lives in the cyclic sites a superbubble omits; `--superbubbles` restricts to the acyclic subset. Internally
it mirrors `vg`'s cactus / 3-edge-connected decomposition, so the default top-level snarl set matches `vg`.

## Algorithm overview

For each top-level snarl candidate:

1. get `(source, sink)` boundary from the internal cactus finder (or `--snarls-in` JSONL);
2. find path intervals crossing source→sink (or sink→source, then canonicalize);
3. collect internal nodes seen between boundaries;
4. measure number of crossing paths and internal bp support across them;
5. apply base filters (`--min-path-support`, `--min-variant-bp`);
6. optionally merge nearby surviving bubbles (`--merge-nearby-bp`); assign bubble IDs.

**Filters.** `--min-variant-bp N` (default 50): keep a candidate only if some supporting path has internal
span ≥ N (`0` disables). `--min-path-support N` (default 0): keep only candidates with ≥ N crossing paths.
Two derived metrics drive/override the size filter: `long_path_support` (paths with span ≥ `--min-variant-bp`)
and `inversion_signal` (some internal node seen in both orientations → kept even if no path reaches the size).

**Nearby merge.** `--merge-nearby-bp N` (default 0, off), applied after the filters: consecutive survivors
merge when the shortest-path distance from the previous `sink` to the next `source` is ≤ N bp.

## Worked trace — snarl boundary + path support

Input graph (one bubble; `+` = forward strand):

```text
S 3  AC         interior allele "ref"   (2 bp)
S 4  ATG        interior allele "alt"   (3 bp)
L 1 2 , 2 3 , 3 5 , 2 4 , 4 5 , 5 6     (edges; nodes 1..6)
```

Paths: `ref: 1 2 3 5 6`, `hA: 1 2 4 5 6`, `hB: 1 2 3 5 6`, `hC: 1 6` (a deletion that bypasses the locus).

1. The cactus finder returns boundary **(source, sink) = (2, 5)** — removing 2 and 5 separates `{3,4}`.
2. Per path, find the interval crossing `2 → 5` and collect interiors: `ref → {3}`, `hA → {4}`, `hB → {3}`;
   `hC` goes `1 → 6`, never visiting **both** boundaries → it does **not** support this snarl.
3. Aggregate: `inside_nodes = {3,4}`, `path_support = 3` (ref, hA, hB; `< 4` total because hC bypasses —
   expected), `min_inside_bp = 2` (`AC`), `max_inside_bp = 3` (`ATG`).

Resulting `bubbles.csv` row:

```text
bubble_id  source  sink  inside_node_count  path_support  min_inside_bp  max_inside_bp  inside_nodes
1          2       5     2                  3             2              3              3,4
```

## Debug TSV (`--snarl-debug-tsv`)

One row per snarl candidate: `candidate_id, source, sink, inside_node_count, n_paths, min_inside_bp,
long_path_support, inversion_signal, accepted`. Bubbles created by `--merge-nearby-bp` have no candidate
row, so one extra `accepted=1` row is appended per merged bubble.
