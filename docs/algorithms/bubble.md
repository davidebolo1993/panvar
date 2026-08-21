# Module `bubble` - algorithm

Mechanism for the `bubble` module. For usage/flags see [modules/bubble.md](../modules/bubble.md); references in [references.md](../references.md#bubble).

`bubble` turns a graph into the variant sites the rest of the pipeline works on. A site is a snarl: a pair of boundary nodes whose removal separates an interior subgraph from the rest, defined on the bidirected graph so it can contain cycles and inversions. Pangenome variation often lives in exactly those cyclic sites, so snarls are the default unit rather than superbubbles — the acyclic special case, a single-source, single-sink subgraph with no cycle or inversion. `--superbubbles` restricts the output to that acyclic subset. The snarls come from a vendored cactus / 3-edge-connected-component decomposition that mirrors [vg](https://github.com/vgteam/vg) `snarls`, so the top-level sites should match vg's.

## How it works

### 1. Sort and flip along the reference

The graph is reordered so numeric node id follows the reference path's order, and nodes the reference traverses in reverse are reverse-complemented to read forward (unless `--no-flip`). This canonical reference-forward frame is the invariant the rest of panvar relies on; it is written as `<prefix>.sorted.gfa` and is what downstream `panphorte`/`call` consume. With `--snarls-in` the graph is taken as-is and this step is skipped.

### 2. Find snarls

The internal cactus / 3-edge-connected finder returns each top-level snarl as a `(source, sink)` boundary pair — the two nodes whose removal isolates the interior. `--snarls-in` substitutes an external `vg snarls` JSONL (JSON Lines) for this step.

Only top-level sites are kept. Descending into the children of an oversized tangle was measured and rejected: the children are overwhelmingly below the size floor, so almost nothing survives to justify the cost.

### 3. Orient the boundaries

A cactus snarl is an **unordered** pair — the decomposition does not say which boundary is left. Every consumer reads them as an interval in reference order: coordinates are anchored on the source, and merging joins one site's sink to the next site's source. Each pair is therefore oriented so `source` is the reference-left boundary, before ids are assigned, so `bubble_id` also increases along the reference. A site the reference does not traverse has no order to take and is left alone.

`source_orient` / `sink_orient` record the strand on which the reference reads each boundary. A boundary is a *handle*, not a node — the same node can bound a site on either side — and merging needs to know which sides face one another.

### 4. Assign interior and support

For each snarl, the paths crossing `source → sink` (or `sink → source`, then canonicalized) are found. The interior is the union of two sets: the nodes those paths visit strictly between the boundaries, and the nodes the **graph** carries between them — everything reachable from the near boundary's inner handle that also reaches the far boundary, over oriented handles. A branch no panel haplotype walks is still part of the site; omitting it also hid interior cycles from `--superbubbles`, which exists to exclude them. The graph traversal is bounded, and on a site too large or too poorly separated to search it falls back to the panel-derived interior and reports that.

A path crossing directly from source to sink — a deletion of the whole interior — is a real allele of the site and counts as support, with an interior span of 0. A path that visits neither boundary does not support the site at all.

Support is then reported at two granularities, because they answer different questions:

- `path_support` counts paths that cross the site *at all*. On a densely typed panel nearly every haplotype crosses nearly every site, so this approaches the panel size and discriminates little.
- `distinct_alleles`, `ref_allele_support` and `alt_allele_support_max`/`min` describe what those traversals *contain*: how many distinct walks exist, how many paths take the reference's walk, and how well the best and worst non-reference walks are supported.

Interior span in bp is kept as its smallest and largest value across the supporting paths.

### 5. Filter and merge

`--min-variant-bp` keeps a site only if some supporting path has an interior span at least that large; `--max-variant-bp` drops one whose longest interior span exceeds a cap, which is how hypervariable tangles are excluded. Both measure the **span between the boundaries**, not the divergence between alleles, so a small edit inside a long allele is filtered by the length of the allele that carries it.

`--min-path-support` filters on traversal support and `--min-alt-support` on the best-supported alternate. The second is usually the intended one: traversal support is close to constant across a typed panel, so it behaves as a step at the panel size, while alternate support varies site by site and behaves as a gradient.

`--merge-nearby-bp` then optionally coalesces survivors that are consecutive **in reference coordinates**, when the sequence strictly between the facing boundaries is within N bp — so two sites sharing a boundary have a gap of 0 regardless of how long that boundary node is. A fused site's interior is recomputed from the graph between its new outer boundaries rather than assembled from its parts, and the fusion is refused if any interior node with a reference coordinate falls outside the new span. Filters are then re-applied, because a fusion can exceed a bound its parts satisfied.

Surviving sites are made pairwise disjoint — two sites claiming the same interior node describe the same sequence twice, which downstream normalization cannot process — and assigned bubble ids in reference order.

### 6. Emit

The surviving bubbles are written to `<prefix>.bubbles.csv` (consumed by `inspect`/`panphorte`/`call`), alongside the sorted GFA and the [Bandage](https://github.com/asl/BandageNG) visualization files. A graph with no snarl at all emits an empty table and exits 0: an empty result is a result, and is distinct from no snarl source having been supplied.

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
L  1  +  6  +  0M             bypass edge (hC skips the locus entirely)
L  2  +  5  +  0M             deletion of the interior (hD)
P  ref  1+,2+,3+,5+,6+  *
P  hA   1+,2+,4+,5+,6+  *
P  hB   1+,2+,3+,5+,6+  *
P  hC   1+,6+           *     never reaches either boundary
P  hD   1+,2+,5+,6+     *     crosses the site, carrying nothing
```

1. **Find snarls.** The cactus finder returns the boundary pair `{2, 5}`: removing nodes 2 and 5 isolates `{3, 4}`.

2. **Orient.** The reference walks 2 before 5, so the pair is stored `source = 2`, `sink = 5`, both read forward.

3. **Assign interior and support.** Per path, take the interval crossing the boundaries and collect what lies between: `ref → {3}`, `hA → {4}`, `hB → {3}`, `hD → {}`. `hC` visits neither boundary and does not support the site. So `path_support = 4`, and the three distinct walks are `{3}` (taken by ref and hB), `{4}` (hA) and the empty one (hD) — giving `distinct_alleles = 3`, `ref_allele_support = 2`, and both alternate-support figures 1. Interior span runs from `0` (hD's deletion) to `3` (`ATG`).

4. **Filter and emit.** With the size floor low enough to admit a 3 bp interior the site is kept and one row is written. Raising `--min-variant-bp` above 3 drops it; raising `--min-alt-support` above 1 also drops it, since no alternate walk here is carried by more than one path.

Resulting `bubbles.csv` row:

```text
bubble_id  source  sink  inside_node_count  path_support  distinct_alleles  ref_allele_support  alt_allele_support_max  min_inside_bp  max_inside_bp  inside_nodes
1          2       5     2                  4             3                 2                   1                       0              3              3;4
```
