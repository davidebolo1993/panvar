# Module `bubble` - algorithm

Mechanism for the `bubble` module. For usage/flags see [modules/bubble.md](../modules/bubble.md); references in [references.md](../references.md#bubble).

`bubble` turns a graph into the variant sites the rest of the pipeline works on. A site is a snarl: a pair of boundary nodes whose removal separates an interior subgraph from the rest, defined on the bidirected graph so it can contain cycles and inversions. Pangenome variation often lives in exactly those cyclic sites, so snarls are the default unit rather than superbubbles — the acyclic special case, a single-source, single-sink subgraph with no cycle or inversion. `--superbubbles` restricts the output to that acyclic subset. The sites come from a vendored cactus / 3-edge-connected-component decomposition compatible with the representation used by [vg](https://github.com/vgteam/vg) `snarls`.

## How it works

### 1. Sort and flip along the reference

The graph is reordered so numeric node id follows the reference path's order, and nodes the reference traverses in reverse are reverse-complemented to read forward (unless `--no-flip`). This canonical reference-forward frame is the invariant the rest of panvar relies on; it is written as `<prefix>.sorted.gfa` and is what downstream `panphorte`/`call` consume. With `--snarls-in` the graph is taken as-is and this step is skipped.

### 2. Find snarls

The internal cactus / 3-edge-connected finder returns each top-level snarl as a boundary pair — the two nodes whose removal isolates the interior. `--snarls-in` substitutes an external `vg snarls` JSONL (JSON Lines) for this step.

Snarls nest: a large site can contain smaller ones. `bubble` emits the top level only, so each region of the graph is described once and the emitted sites are non-overlapping, which is what the normalization and calling steps downstream require.

### 3. Order the boundaries along the reference

A cactus snarl is an unordered pair — the decomposition does not say which boundary comes first. Everything downstream reads them as an interval in reference order: coordinates are anchored on `source`, and merging joins one site's `sink` to the next site's `source`. Each pair is therefore ordered so `source` is the boundary the reference reaches first, before ids are assigned, so `bubble_id` also increases along the reference. A site the reference does not traverse has no order to take and keeps the decomposition's pair as given.

A boundary is a node together with the side the site lies on, since the same node can bound a site on either side. `source_orient` and `sink_orient` record the strand on which the reference reads each boundary, which is what tells merging which sides face one another.

### 4. Determine the interior and score it

For each snarl, the paths crossing between the boundaries are found and canonicalized to run `source → sink`. The interior is taken from two sources: the nodes those paths visit strictly between the boundaries, and the nodes the graph itself carries between them, found by walking from the inner side of one boundary to the other over oriented node sides. Using both means a site is described by the structure the graph contains, not only by the subset of it this particular panel happens to walk — which matters for a graph whose paths have been subset, or one built by a tool that records alleles no stored path uses. The graph walk is bounded; past that bound the path-derived interior is used alone and the run reports it.

A path crossing directly from one boundary to the other, carrying nothing between them, is a deletion of the whole interior. It is an allele of the site like any other and counts toward `path_support`, with an interior span of 0. A path visiting neither boundary does not cross the site and does not count.

Scoring then records both how many paths cross and what those crossings contain. `path_support` is the count of crossing paths; because a path carrying the reference allele also crosses, this number tends toward the panel size on a densely typed panel. `distinct_alleles`, `ref_allele_support`, `alt_allele_support_max` and `alt_allele_support_min` describe the walks themselves: how many different ones there are, how many paths follow the reference's, and how well the best- and least-supported of the others are carried. Interior length is kept as `min_inside_bp` and `max_inside_bp`, the shortest and longest interior across the crossing paths.

### 5. Filter and merge

`--min-variant-bp` keeps a site only if some crossing path has an interior span at least that large; `--max-variant-bp` drops one whose longest interior span exceeds a cap, which is how hypervariable tangles are excluded. Both measure the span between the boundaries rather than the difference between alleles, so a small edit inside a long allele is filtered by the length of the allele carrying it.

`--min-path-support` filters on the number of crossing paths and `--min-alt-support` on the best-supported non-reference walk. The second is usually the intended one: since every crossing path counts toward the first, including those carrying the reference allele, it varies little across a typed panel and behaves as a threshold at the panel size, while alternate support varies site by site.

`--merge-nearby-bp` then optionally coalesces survivors that are consecutive along the reference, when the sequence strictly between the facing boundaries is within N bp. Distance is measured between the boundaries themselves, so two sites that abut have a gap of 0 however long the shared boundary node is. A fused site's interior is recomputed from the graph between its new outer boundaries, and the fusion is refused if an interior node with a reference coordinate would fall outside the new span. Filters are re-applied afterwards, since a fusion can exceed a bound that both of its parts satisfied.

Surviving sites are made pairwise disjoint, so no interior node belongs to two of them, and are assigned ids in reference order.

### 6. Emit

The surviving bubbles are written to `<prefix>.bubbles.csv` (consumed by `inspect`/`panphorte`/`call`), alongside the sorted GFA and the [Bandage](https://github.com/asl/BandageNG) visualization files. A graph containing no snarl emits an empty table and exits 0.

## Worked trace

The steps below follow the six above, one for one. Input graph, one bubble (`+` = forward strand):

```text
S  1  T                       flank
S  2  GG                      first boundary
S  3  ATG                     interior allele, 3 bp, carried by hA
S  4  AC                      interior allele, 2 bp, carried by the reference and hB
S  5  CC                      second boundary
S  6  A                       flank
L  1  +  2  +  0M
L  2  +  3  +  0M
L  3  +  5  +  0M
L  2  +  4  +  0M
L  4  +  5  +  0M
L  5  +  6  +  0M
L  1  +  6  +  0M             bypass edge: hC skips the locus entirely
L  2  +  5  +  0M             deletion of the interior: hD
P  ref  1+,2+,4+,5+,6+  *
P  hA   1+,2+,3+,5+,6+  *
P  hB   1+,2+,4+,5+,6+  *
P  hC   1+,6+           *     reaches neither boundary
P  hD   1+,2+,5+,6+     *     crosses the site carrying nothing
```

1. Sort and flip. Node ids already follow the reference and every node reads forward, so the sorted graph is identical to the input and `<prefix>.sorted.gfa` is a copy of it. On a graph that is not already in this frame the ids are reassigned here, and the ids in every output below are the new ones.

2. Find snarls. The cactus finder returns the boundary pair `{2, 5}`: removing nodes 2 and 5 isolates `{3, 4}`. There is no smaller snarl nested inside it.

3. Order the boundaries along the reference. The reference reaches 2 before 5, so the pair is stored `source = 2`, `sink = 5`. The reference reads both forward, so `source_orient` and `sink_orient` are `+`.

4. Determine the interior and score it. The interior from the graph is `{3, 4}`, everything lying between the two boundaries; the crossing paths add nothing beyond it. Per path, what lies between the boundaries is `ref → {4}`, `hA → {3}`, `hB → {4}`, `hD → {}`; `hC` reaches neither boundary and does not cross. So `path_support` is 4, and the walks taken are `{4}` by the reference and hB, `{3}` by hA, and the empty one by hD: `distinct_alleles` 3, `ref_allele_support` 2, and both alternate figures 1. Interior span runs from `min_inside_bp` 0, hD's deletion, to `max_inside_bp` 3, the `ATG` allele.

5. Filter and merge. With the size floor low enough to admit a 3 bp interior the site is kept. Raising `--min-variant-bp` above 3 drops it, since no crossing carries more than that between the boundaries; raising `--min-alt-support` above 1 also drops it, since neither non-reference walk is carried by more than one path. There is only one site, so merging has nothing to join. The interior contains no cycle, so `--superbubbles` keeps it too.

6. Emit. One row is written to `<prefix>.bubbles.csv`, alongside the sorted GFA and the Bandage files:

```text
bubble_id  source  sink  inside_node_count  path_support  distinct_alleles  ref_allele_support  alt_allele_support_max  min_inside_bp  max_inside_bp  inside_nodes
1          2       5     2                  4             3                 2                   1                       0              3              3;4
```
