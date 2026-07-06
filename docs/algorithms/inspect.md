# Module `inspect` - algorithm

Mechanism and a hand-traced example for the `inspect` utility's `--cluster` mode. For usage/flags see [modules/inspect.md](../modules/inspect.md); References in [references.md](../references.md#inspect).

## Terms

- **walk** — a path's `source → sink` traversal of a bubble, as a sequence of `<node_id><strand>` steps.
- **shingle** — a window of `k` consecutive steps (here `k = 3`); the multiset of a walk's shingles is its fingerprint.
- **MinHash sketch** — a small, fixed-size summary (bottom-`k` smallest hashes of the shingles) that estimates the Jaccard similarity of two shingle multisets cheaply, without enumerating them.

## Algorithm Overview

1. Identical walks collapse to one representative (with a support count).
2. Each distinct walk gets a multiplicity-aware MinHash sketch over its shingles: the k-th occurrence of a shingle is salted to a distinct element, so the sketch tracks the multiset (copy number), not just the set.
3. Two walks' similarity is the sketch-estimated identity `= 2J/(1+J)` from shingle Jaccard `J`. A threshold graph connects every pair ≥ `--cluster-similarity` (default 0.90); clusters are its connected components (transitive, order-independent).
4. Very short walks that can't be shingled fall back to an exact bp-weighted Jaccard over `<node_id><strand>` tokens (`sum(min)/sum(max)`).
5. Representative = the member minimizing max-then-mean intra-cluster distance (most-supported walk breaks ties).


## Worked trace

Four haplotypes cross one bubble; `a`/`b` are the flanking boundaries and `(u, v)` is the repeat unit:

```text
h1 : a u v u v b              unit (u,v) ×2
h2 : a u v u v b              identical to h1
h3 : a u v u v u v b          unit (u,v) ×3  (copy-number expansion)
h4 : a x y b                  a different allele
```

1. Collapse identical walks. `h1` and `h2` spell the same step sequence, so they fold into one distinct walk `W1` (support 2); `h3` becomes `W3` and `h4` becomes `W4`. Only these distinct walks are sketched and compared.
2. Shingle each walk into windows of `k = 3` consecutive steps — the multiset of windows is the walk's fingerprint:
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
   J = 4/6 = 0.667  →  identity = 2·0.667/1.667 = 0.80   (< 0.90 → SEPARATE)
   ```
4. Estimate that `J` cheaply with a bottom-k MinHash instead of enumerating every shingle: with illustrative hashes `sketch₄(W1) = {12, 27, 41, 55}` and `sketch₄(W3) = {08, 12, 27, 33}`, the shared bottom-4 hashes give `J ≈ 0.67`, matching the exact multiset Jaccard.
5. Cluster by connected components at `--cluster-similarity 0.90`. No pair clears the threshold (W1–W3 = 0.80 ✗, W1–W4 ✗, W3–W4 ✗), so the three walks stay in their own clusters — the 2-copy and 3-copy alleles are kept apart rather than chained into one band.

Resulting `clusters.tsv`:

```text
cluster_id  n_paths  representative_path  members
cluster0    2        h1                   h1;h2
cluster1    1        h3                   h3
cluster2    1        h4                   h4
```
