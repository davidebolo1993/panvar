# Module `bubble` - algorithm

Mechanism for the `bubble` module. For usage/flags see [modules/bubble.md](../modules/bubble.md); References in [references.md](../references.md#bubble).

`bubble` turns a graph into the variant sites the rest of the pipeline works on. A site is a snarl: a pair of boundary nodes whose removal separates an interior subgraph from the rest, defined on the bidirected graph so it can contain cycles and inversions. Pangenome variation often lives in exactly those cyclic sites, so snarls are the default unit rather than superbubbles — the acyclic special case, a single-source, single-sink subgraph with no cycle or inversion. `--superbubbles` restricts the output to that acyclic subset. The snarls come from a vendored cactus / 3-edge-connected-component decomposition that mirrors [vg](https://github.com/vgteam/vg) `snarls`, so the top-level sites should match vg's.

## How it works

### 1. Sort and flip along the reference

The graph is reordered so numeric node id follows the reference path's order, and nodes the reference traverses in reverse are reverse-complemented to read forward (unless `--no-flip`). This canonical reference-forward frame is the invariant the rest of panvar relies on; it is written as `<prefix>.sorted.gfa` and is what downstream `panphorte`/`call` consume. With `--snarls-in` the graph is taken as-is and this step is skipped.

### 2. Find snarls

The internal cactus / 3-edge-connected finder returns each top-level snarl as a `(source, sink)` boundary pair — the two nodes whose removal isolates the interior. `--snarls-in` substitutes an external `vg snarls` JSONL (JSON Lines) for this step.

### 3. Assign interior and support

For each snarl, the paths crossing `source → sink` (or `sink → source`, then canonicalized) are found, and the nodes each visits strictly between the boundaries are collected as the interior. Two quantities are measured across those crossing paths: path support, the number of paths visiting both boundaries, and the interior span in bp, kept as its smallest and largest value. A path that does not visit both boundaries — a deletion that bypasses the locus, say — does not support the snarl.

### 4. Filter and merge

Base filters drop the uninteresting sites. `--min-variant-bp` keeps a snarl only if some supporting path has an interior span at least that large (the size filter, default 50); `--min-path-support` keeps only snarls crossed by at least that many paths. `--merge-nearby-bp` then optionally coalesces consecutive survivors when the shortest-path distance from one `sink` to the next `source` is within N bp. Surviving sites are assigned bubble ids.

### 5. Emit

The surviving bubbles are written to `<prefix>.bubbles.csv` (consumed by `inspect`/`panphorte`/`call`), alongside the sorted GFA and the [Bandage](https://github.com/asl/BandageNG) visualization files.

## Worked trace

Input graph, one bubble (`+` = forward strand):

```text
S  1  T                       flank
S  2  GG                      source boundary
S  3  AC                      interior allele "ref"  (2 bp)
S  4  ATG                     interior allele "alt"  (3 bp)
S  5  CC                      sink boundary
S  6  A                       flank
L  1  +  2  +  0M
L  2  +  3  +  0M
L  3  +  5  +  0M
L  2  +  4  +  0M
L  4  +  5  +  0M
L  5  +  6  +  0M
L  1  +  6  +  0M             bypass edge (hC's deletion)
P  ref  1+,2+,3+,5+,6+  *
P  hA   1+,2+,4+,5+,6+  *
P  hB   1+,2+,3+,5+,6+  *
P  hC   1+,6+           *     deletion that bypasses the locus
```

1. Find snarls. The cactus finder returns boundary `(source, sink) = (2, 5)`: removing nodes 2 and 5 isolates `{3, 4}`.

2. Assign interior and support. Per path, take the interval crossing `2 → 5` and collect its interior — `ref → {3}`, `hA → {4}`, `hB → {3}`. `hC` goes `1 → 6` and never visits both boundaries, so it does not support this snarl. Aggregating: `inside_nodes = {3, 4}`, `path_support = 3` (ref, hA, hB — fewer than the 4 total paths because hC bypasses, as expected), `min_inside_bp = 2` (`AC`), `max_inside_bp = 3` (`ATG`).

3. Filter and emit. The site is kept (here with a `--min-variant-bp` low enough to admit a 3 bp interior) and one row is written.

Resulting `bubbles.csv` row:

```text
bubble_id  source  sink  inside_node_count  path_support  min_inside_bp  max_inside_bp  inside_nodes
1          2       5     2                  3             2              3              3,4
```
