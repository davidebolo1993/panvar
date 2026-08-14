# Module `inspect` - algorithm

Mechanism for the `inspect` utility's `--cluster` mode. For usage/flags see [modules/inspect.md](../modules/inspect.md); References in [references.md](../references.md#inspect).

`--cluster` groups the haplotypes crossing a bubble by how similarly they traverse it, so structurally identical alleles collapse to one representative for plotting. A walk is a path's `source → sink` traversal of the bubble, written as a sequence of `<node_id><strand>` steps; two walks are compared by how much their windowed content overlaps.

## How it works

### 1. Collapse identical walks

Walks that spell the same step sequence fold to one representative up front, so only distinct walks are sketched and compared.

### 2. Sketch each walk

Each distinct walk is broken into shingles — windows of k consecutive steps — and summarized by a MinHash sketch, the bottom-k smallest hashes of those shingles. The sketch is multiplicity-aware: the k-th occurrence of a shingle is salted to a distinct element, so it tracks the multiset (copy number), not just the set of shingles.

### 3. Score similarity

Two walks' similarity is the sketch-estimated identity `2J / (1 + J)`, where `J` is the shingle Jaccard. Very short walks that cannot be shingled fall back to an exact bp-weighted Jaccard over `<node_id><strand>` tokens (`sum(min) / sum(max)`).

### 4. Cluster

A threshold graph connects every pair at or above `--cluster-similarity`; clusters are its connected components, so membership is transitive and order-independent.

### 5. Pick a representative

Each cluster's representative is the member minimizing max-then-mean intra-cluster distance, with the most-supported walk breaking ties.


## Worked trace

Four haplotypes cross one bubble; `a`/`b` are the flanking boundaries and `(u, v)` is the repeat unit:

```text
h1 : a u v u v b              unit (u,v) ×2
h2 : a u v u v b              identical to h1
h3 : a u v u v u v b          unit (u,v) ×3  (copy-number expansion)
h4 : a x y b                  a different allele
```

1. Collapse identical walks. `h1` and `h2` spell the same step sequence, so they fold into one distinct walk `W1` (support 2); `h3` becomes `W3` and `h4` becomes `W4`. Only these distinct walks are sketched and compared.
2. Shingle each walk into windows of length k (here `k = 3`) consecutive steps — the multiset of windows is the walk's fingerprint:
   ```text
   W1 : (a,u,v) (u,v,u) (v,u,v) (u,v,b)
   W3 : (a,u,v) (u,v,u) (v,u,v) (u,v,u) (v,u,v) (u,v,b)
   W4 : (a,x,y) (x,y,b)
   ```
3. Compare shingle multisets, so copy number is visible. As plain sets `W1` and `W3` hold exactly the same shingles, so a set Jaccard would be 1.0 and merge the 2-copy and 3-copy alleles; keeping each shingle's count separates them, because the two middle shingles occur once in `W1` but twice in `W3`:
   ```text
   multiset(W1): auv1 uvu1 vuv1 uvb1
   multiset(W3): auv1 uvu2 vuv2 uvb1
   intersection (min counts) = 4 ;  union (max counts) = 6
   J = 4/6 = 0.667
   identity = 2·0.667/1.667 = 0.80
   ```
4. Estimate that `J` cheaply with a bottom-k MinHash instead of enumerating every shingle. Hash each occurrence-salted shingle to one value — say `auv₁→12, uvu₁→27, vuv₁→41, uvb₁→55`, plus `W3`'s extra copies `uvu₂→08, vuv₂→33`. Each sketch keeps the k smallest of its own walk; with k above 6 nothing is dropped here, so `sketch(W1) = {12, 27, 41, 55}` and `sketch(W3) = {08, 12, 27, 33, 41, 55}`.

   The estimator takes the **k smallest values of the two sketches' union**, with `k = min(|sketch(W1)|, |sketch(W3)|) = 4`, and asks how many of those appear in both. The union sorts to `08, 12, 27, 33, 41, 55`; its bottom 4 are `08, 12, 27, 33`, of which `12` and `27` are in both — giving `2/4 = 0.5`. (This example is small enough that the estimator's sampling error is visible: the multiset Jaccard from step 3 is 0.667. At sketch sizes that matter, hundreds of values, the estimate converges on it.)

   Comparing the two sketches directly as `|shared| / |union|` would be **biased**, and that is why the union is taken first. Each sketch holds the smallest hashes of *its own* shingle set, so when a walk has more than k shingles and is truncated, a shingle both walks carry can sit inside one sketch and below the other's cutoff. The two sketches then sample different regions of the hash space. The error grows with the length ratio, which is exactly the tandem-array case: on a 3000-shingle walk against a 500-shingle one with a true `J` of 0.167, the direct ratio reads **0.085** while the union form reads **0.158**.
5. Cluster by connected components at `--cluster-similarity 0.90`, using the identity derived from the estimated `J`. No pair clears the threshold, so the three walks stay in their own clusters — the 2-copy and 3-copy alleles are kept apart rather than chained into one band.

Resulting `clusters.tsv`:

```text
cluster_id  n_paths  representative_path  members
cluster0    2        h1                   h1;h2
cluster1    1        h3                   h3
cluster2    1        h4                   h4
```
