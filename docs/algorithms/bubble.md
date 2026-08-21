# Module `bubble` - algorithm

Mechanism for the `bubble` module. For usage/flags see [modules/bubble.md](../modules/bubble.md); references in [references.md](../references.md#bubble).

`bubble` turns a graph into the variant sites the rest of the pipeline works on. A site is a snarl: a pair of boundary nodes whose removal separates an interior subgraph from the rest, defined on the bidirected graph so it can contain cycles and inversions. Pangenome variation often lives in exactly those cyclic sites, so snarls are the default unit rather than superbubbles — the acyclic special case, a single-source, single-sink subgraph with no cycle or inversion. `--superbubbles` restricts the output to that acyclic subset. The snarls come from a vendored cactus / 3-edge-connected-component decomposition that mirrors [vg](https://github.com/vgteam/vg) `snarls`, so the top-level sites should match vg's.

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
L  1  +  6  +  0M             bypass edge: hC skips the locus entirely
L  2  +  5  +  0M             deletion of the interior: hD
P  ref  1+,2+,3+,5+,6+  *
P  hA   1+,2+,4+,5+,6+  *
P  hB   1+,2+,3+,5+,6+  *
P  hC   1+,6+           *     reaches neither boundary
P  hD   1+,2+,5+,6+     *     crosses the site carrying nothing
```

1. Find snarls. The cactus finder returns the boundary pair `{2, 5}`: removing nodes 2 and 5 isolates `{3, 4}`.

2. Order the boundaries. The reference reaches 2 before 5, so the pair is stored `source = 2`, `sink = 5`, both read forward.

3. Determine the interior and score it. Per path, take the crossing and collect what lies between the boundaries: `ref → {3}`, `hA → {4}`, `hB → {3}`, `hD → {}`. `hC` reaches neither boundary and does not cross. So `path_support = 4`, and the three walks taken are `{3}` (by ref and hB), `{4}` (hA) and the empty one (hD), giving `distinct_alleles = 3`, `ref_allele_support = 2`, and both alternate-support figures 1. Interior span runs from 0, hD's deletion, to 3, the `ATG` allele.

4. Filter and emit. With the size floor low enough to admit a 3 bp interior the site is kept and one row is written. A `--min-variant-bp` above 3 drops it, since no crossing carries more than that between the boundaries. A `--min-alt-support` above 1 also drops it, since neither non-reference walk here is carried by more than one path.

Resulting `bubbles.csv` row:

```text
bubble_id  source  sink  inside_node_count  path_support  distinct_alleles  ref_allele_support  alt_allele_support_max  min_inside_bp  max_inside_bp  inside_nodes
1          2       5     2                  4             3                 2                   1                       0              3              3;4
```
