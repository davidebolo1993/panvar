# inspect — clustering algorithm & worked example

Mechanism and a hand-traced example for the **inspect** utility's `--cluster` mode. For usage/flags see
[modules/inspect.md](../modules/inspect.md); citations: [references.md](../references.md#inspect).

## Terms

- **walk** — a path's `source → sink` traversal of a bubble, as a sequence of `<node_id><strand>` steps
  (clustering runs on these step tokens, not raw DNA).
- **shingle** — a window of `k` consecutive steps (here `k = 3`); the multiset of a walk's shingles is its
  fingerprint.
- **MinHash sketch** — a small, fixed-size summary (bottom-`k` smallest hashes of the shingles) that
  estimates the **Jaccard** similarity of two shingle multisets cheaply, without enumerating them.

## How `--cluster` works

1. Identical walks collapse to one representative (with a support count).
2. Each distinct walk gets a **multiplicity-aware** MinHash sketch over its shingles: the *k*-th occurrence
   of a shingle is salted to a distinct element, so the sketch tracks the **multiset** (copy number), not
   just the set.
3. Two walks' similarity is the sketch-estimated identity `= 2J/(1+J)` from shingle Jaccard `J`. A threshold
   graph connects every pair ≥ `--cluster-similarity` (default 0.90); clusters are its **connected
   components** (transitive, order-independent).
4. Very short walks that can't be shingled fall back to an exact bp-weighted Jaccard over
   `<node_id><strand>` tokens (`sum(min)/sum(max)`).
5. Representative = the member minimizing max-then-mean intra-cluster distance (ties → most-supported walk).

Multiplicity-awareness is the crux: copies of one unit share the same shingle *set* regardless of copy
number, so a set sketch would merge all copy numbers into one cluster; the multiset sketch separates them.

## Worked trace — four haplotypes through one bubble

```text
h1 : a u v u v b              unit (u,v) ×2
h2 : a u v u v b              identical to h1
h3 : a u v u v u v b          unit (u,v) ×3  (copy-number expansion)
h4 : a x y b                  a different allele
```

**Collapse identical:** h1==h2 → distinct walks `W1`(=h1/h2, support 2), `W3`(=h3), `W4`(=h4).

**Shingles (k=3):**
```text
W1 : (a,u,v) (u,v,u) (v,u,v) (u,v,b)
W3 : (a,u,v) (u,v,u) (v,u,v) (u,v,u) (v,u,v) (u,v,b)
W4 : (a,x,y) (x,y,b)
```

**The copy-number trap (set vs multiset).** As *sets*, `W1` and `W3` are identical → set Jaccard 1.0 →
would merge. Multiset-aware:
```text
multiset(W1): auv1 uvu1 vuv1 uvb1
multiset(W3): auv1 uvu2 vuv2 uvb1
intersection (min counts) = 4 ;  union (max counts) = 6
J = 4/6 = 0.667  →  identity = 2(0.667)/1.667 = 0.80   (< 0.90 → SEPARATE)
```

**Bottom-k MinHash** estimates that `J` without enumerating: with illustrative hashes,
`sketch₄(W1)={12,27,41,55}`, `sketch₄(W3)={08,12,27,33}`, shared bottom-4 ⇒ `J ≈ 0.67`, matching the
multiset Jaccard.

**Connected components @ 0.90:** edges W1–W3 (0.80) ✗, W1–W4 ✗, W3–W4 ✗ → **3 clusters**:
`cluster0 = {h1,h2}` (rep h1, support 2), `cluster1 = {h3}`, `cluster2 = {h4}`.

> Transitivity matters when identities clear the threshold: if W1–W3 and W3–W4 passed but W1–W4 didn't, all
> three would still land in one component (single-linkage). On a copy-number *continuum* (LPA KIV-2) the
> default 0.90 chains many haplotypes into one cluster — raise `--cluster-similarity` to cut finer
> copy-number bands (the per-gene scripts use 0.95 for LPA, 0.97 for the others).
