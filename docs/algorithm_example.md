# Algorithmic Examples (worked traces)

This page complements the per-module docs. Instead of *explaining* each algorithm, it **traces a tiny
concrete input** through it, step by step, so one can see exactly what the module computes. For the
rationale and full option list see the module pages
([bubble](modules/bubble.md), [inspect](modules/inspect.md), [panphorte](modules/panphorte.md),
[call](modules/call.md), [describe](modules/describe.md)).

> Inputs here are deliberately minimal. Where the real code uses 64-bit hashes, the traces use small
> made-up hash values (clearly marked) to keep the arithmetic readable.

---

## bubble — snarl boundary + path support

**Input graph** (one bubble; `+` = forward strand):

```text
S 1  ...        nodes 1..6
S 2  ...
S 3  AC         interior allele "ref"   (2 bp)
S 4  ATG        interior allele "alt"   (3 bp)
S 5  ...
S 6  ...
L 1 2 ,  2 3 ,  3 5 ,   2 4 ,  4 5 ,  5 6     (edges)
```

**Paths (haplotypes):**

```text
ref : 1+ 2+ 3+ 5+ 6+
hA  : 1+ 2+ 4+ 5+ 6+
hB  : 1+ 2+ 3+ 5+ 6+
hC  : 1+ 6+                 (a deletion that bypasses the locus)
```

**Trace:**

1. The internal cactus finder returns the top-level boundary pair **(source, sink) = (2, 5)** — removing
   nodes 2 and 5 separates `{3,4}` from the rest.
2. For each path, find the interval that crosses `2 → 5` and collect the interior nodes between them:
   - `ref` crosses `2 → 3 → 5` → interior `{3}`
   - `hA` crosses `2 → 4 → 5` → interior `{4}`
   - `hB` crosses `2 → 3 → 5` → interior `{3}`
   - `hC` goes `1 → 6`, never visiting **both** boundaries → **does not support** this snarl.
3. Aggregate:
   - `inside_nodes` = union of interiors = `{3, 4}`
   - `path_support` = **3** (ref, hA, hB) — note this is `< 4` total paths, because `hC` bypasses the
     boundary (expected, not an error).
   - `min_inside_bp` = 2 (the `AC` allele), `max_inside_bp` = 3 (the `ATG` allele).

**Resulting `bubbles.csv` row:**

```text
bubble_id  source  sink  inside_node_count  path_support  min_inside_bp  max_inside_bp  inside_nodes
1          2       5     2                  3             2              3              3,4
```

---

## inspect — walk clustering (shingles, multiset MinHash, connected components)

Clustering runs over **node-step walks** (each step is a token `<node_id><strand>`), not raw DNA. Shingle
size is **3** (consecutive step triples). Two distinct walks are connected when their estimated identity
`= 2J / (1 + J)` (from shingle Jaccard `J`) is `≥ --cluster-similarity` (default `0.90`); clusters are the
**connected components** of that graph.

**Input walks** (4 haplotypes through one bubble):

```text
h1 : a+ u+ v+ u+ v+ b+              (unit (u,v) x2)
h2 : a+ u+ v+ u+ v+ b+              (identical to h1)
h3 : a+ u+ v+ u+ v+ u+ v+ b+       (unit (u,v) x3  — a copy-number expansion)
h4 : a+ x+ y+ b+                    (a different allele)
```

### Step 1 — collapse identical walks

`h1` and `h2` are identical → they collapse to one **distinct walk** with support 2. We now cluster three
distinct walks: `W1` (=h1/h2), `W3` (=h3), `W4` (=h4).

### Step 2 — shingles (size 3)

```text
W1 : (a,u,v) (u,v,u) (v,u,v) (u,v,b)
W3 : (a,u,v) (u,v,u) (v,u,v) (u,v,u) (v,u,v) (u,v,b)
W4 : (a,x,y) (x,y,b)        <- only 2 tokens-of-shingle; very short
```

### Step 3 — the copy-number trap (set vs multiset)

Look at `W1` vs `W3`. As **sets** of distinct shingles they are *identical*:

```text
set(W1) = { auv, uvu, vuv, uvb }
set(W3) = { auv, uvu, vuv, uvb }     -> set Jaccard = 4/4 = 1.0  -> identity 1.0  (would MERGE!)
```

A plain set sketch would lump every copy number into one cluster. panvar's sketch is **multiset-aware**:
the *k*-th occurrence of a shingle is salted to a distinct sketch element, so multiplicity is kept:

```text
multiset(W1): auv1  uvu1  vuv1  uvb1
multiset(W3): auv1  uvu2  vuv2  uvb1
intersection (min counts) = auv1+uvu1+vuv1+uvb1 = 4
union        (max counts) = auv1+uvu2+vuv2+uvb1 = 6
multiset Jaccard J = 4/6 = 0.667  ->  identity = 2(0.667)/(1.667) = 0.80   (< 0.90 -> SEPARATE)
```

So `W1` and `W3` end up in **different** clusters, separated by copy number — exactly what we want for a
tandem repeat.

### Step 4 — bottom-k MinHash (how J is estimated cheaply)

The real code does not enumerate full sets; it keeps a bottom-`k` MinHash sketch. With **illustrative**
hash values (real code uses 64-bit hashes; the occurrence index is folded in):

```text
element        hash (made up)
auv#1            12
uvu#1            41
vuv#1            27
uvb#1            55
uvu#2            08      <- only in W3 (2nd copy)
vuv#2            33      <- only in W3 (2nd copy)

sketch_k=4(W1) = {08? no}  bottom-4 of {12,41,27,55}      = {12,27,41,55}
sketch_k=4(W3) = bottom-4 of {12,41,27,55,08,33}          = {08,12,27,33}
shared bottom-4 (after merge of the two sorted sketches)  -> estimates J ~ 0.67
```

The sketch Jaccard tracks the multiset Jaccard above (≈ 0.67), which maps to identity ≈ 0.80.

### Step 5 — connected components

At the default `0.90`: edges `W1–W3` (0.80) ✗, `W1–W4` ✗, `W3–W4` ✗. Result — **3 clusters**:

```text
cluster 0 : h1, h2        (representative = h1; support 2)
cluster 1 : h3
cluster 2 : h4
```

> Transitivity matters when identities clear the threshold. If `W1–W3` and `W3–W4` both passed but
> `W1–W4` did not, all three would still land in **one** component (single-linkage). This is why, on a
> real copy-number *continuum* (e.g. LPA KIV-2), the default threshold chains many haplotypes into one
> cluster — raise `--cluster-similarity` to cut it into finer copy-number bands.

Very short walks (`W4` has fewer than 3 shingles' worth of tokens) cannot be shingled; for those the code
falls back to an exact **bp-weighted Jaccard** over `<node_id><strand>` tokens.

---

## panphorte — tandem-repeat collapse

Detection works on **spelled-sequence tokens**, one per path step: each step's token is the bp sequence
of the node it traverses. Comparing by *sequence* (not node id) means identical copies match even when
they are distinct nodes in the graph. Let `U` denote a 60 bp unit sequence and `x` an 8 bp interruption.

**Input** — a path interior with four copies of `U`, but with a small non-unit node `x` slipped between
the 2nd and 3rd copy (a real array is rarely perfectly clean):

```text
step index :  0    1    2    3    4    5    6
token (seq):  L    U    U    x    U    U    R         (L=flank, U=unit, x=interruption, R=flank)
bp         :  ..   60   60   8    60   60   ..
```

### Step 1 — establish the unit period

Scan left→right. At anchor `i = 1` (`token = U`), look at the other positions where `U` occurs
(`{1,2,4,5}`) and take the **smallest period** `p = j − i` for which the `p`-step block at `i` equals the
block at `i+p` (an *adjacent identical pair*). Here `j = 2 → p = 1`, and `block_equal(1, 2, 1)` is true
(`U == U`). So **period `p = 1`**, **`unit_bp = 60`**. (Up to 8 candidate `j` are tried; if none gives a
clean adjacent pair, advance `i`.)

### Step 2 — extend the run, tolerating interruptions

Start `copy_starts = {1}` and walk right by `p`:

```text
last=1 → adj=2 : block_equal(1,2,1)=U==U  ✓   copy_starts={1,2}
last=2 → adj=3 : token[3]=x ≠ U           ✗   → look for a gap…
          gap g=1: candidate=4, gap_bp=bp[3]=8 ≤ unit_bp=60 ✓ and block_equal(1,4,1)=U==U ✓
                   copy_starts={1,2,4}     (the x is recorded as an interruption)
last=4 → adj=5 : U==U  ✓                       copy_starts={1,2,4,5}
last=5 → adj=6 : token[6]=R ≠ U, no further U → stop
```

Then extend left from `i=1` the same way (here `token[0]=L ≠ U` and no room for a gap → no extension).
A gap is only bridged while it stays short: at most `kMaxGapSteps` steps **and** cumulative `gap_bp ≤
unit_bp`.

### Step 3 — accept test

```text
copies         = 4          ≥ --min-copies (2)              ✓
unit_bp        = 60         ≥ --min-unit-bp (50)            ✓
interruption   = 8 bp;  array_bp = 4·60 + 8 = 248
                8 ≤ --max-interruption-frac (0.25) · 248 = 62  ✓
```

All pass → collapse.

### Step 4 — collapse

One **REP** node carries `U`, with a self-loop; the path is rerouted through REP `copies` times, and the
interruption `x` is kept as literal steps where it sat. Duplicate unit nodes are dropped.

```text
after path :  L+  REP+ REP+  x+  REP+ REP+  R+        (REP self-loop traversed 4× total)
```

**Report row:** `unit_bp=60  copies=4  interruption_bp=8`. Downstream, `call --cn-from-multiplicity`
reads copy number straight off the REP self-loop multiplicity.

### Exact vs approximate (`--min-similarity`)

The trace above is **exact** mode (`--min-similarity 1.0`): `block_equal` requires byte-identical copies.
Real KIV-2 copies differ by SNVs, so the exact pair test fails and nothing is seeded. With
`--min-similarity 0.90`, panphorte switches to a **single-representative, lossy** mode: it picks one unit
`U`, then finds its near-copies by multi-seed **banded alignment** of `U` (and its reverse complement)
against the path's spelled sequence, accepting a copy when `identity ≥ 0.90` and at least half of `U` is
consumed:

```text
copy 1: identity 0.98 ✓   copy 2: 0.95 ✓   copy 3: 0.93 ✓   copy 4: 0.97 ✓   → 4 copies → one REP=U
```

Because alignment (not exact match) bridges a copy, a large *internal indel* between divergent copies is
also bridged when the threshold is low enough — this is what lets a C4 short module (missing the ~6.4 kb
HERV-K) fold onto the long one at `--min-similarity 0.70`. The collapse is lossy: the per-copy SNVs are
discarded in favour of the single representative `U`.

---

## call — diff vs reference, intra-hap coalescing, cross-hap clustering

`call` turns each haplotype's source→sink walk into events vs the reference walk, then merges twice:
**within** a haplotype (coalesce fragments) and **across** haplotypes (cluster the same site into one
record). Copy-number nodes (REP self-loops from panphorte) are dropped from this token walk — they are
typed separately as DUP/CN — so extra repeat copies are never mistyped as insertions.

### Step 1 — anchor, then align between anchors

**Reference walk** `R = A B C D E F`. Take a haplotype `hX = A B D E F` (C missing) plus an inserted run.
Tokens that are **unique in both** `R` and the haplotype are anchor candidates; a monotonic (LIS)
subsequence of them is chosen so anchors are consistent in both walks. Only the **segments between
anchors** are aligned (Needleman–Wunsch), which keeps it fast and localises each edit.

```text
anchors (unique, increasing in both):  A … B … D … E … F
segment between B and D:   ref = [C]      hap = []        → gap block
```

### Step 2 — read events off the gap blocks (`diff_segment`)

A maximal gap block becomes:

```text
ref-only block              → DEL   (nodes = ref nodes,  anchor = first ref node)
hap-only block              → INS   (nodes = hap nodes,  anchored AFTER the last matched ref node)
hap block == revcomp(ref)   → INV   (one event over the reference nodes)
both ref & hap present       → substitution: emit DEL + INS, linked by a shared EVENTID
```

For `hX` the `[C]` ref-only block → **DEL `{C}`**. A haplotype that instead inserted `P Q` between `B`
and `D` would give **INS `{P,Q}`** anchored after `B`.

### Step 3 — coalesce fragments within a haplotype (`coalesce_events`)

Consecutive **same-type** events merge when their gap is `≤ --merge-distance-bp` (default 100) in
**either** of two coordinate spaces:

- **reference bp** — `cur_lo − prev_hi`, using each node's reference start;
- **haplotype bp** — the inserted nodes' own positions in *this* haplotype's walk.

The hap-space clause catches same-type events that are far apart on the reference (a deletion sits
between them) yet contiguous in the sample's own sequence. Example — `hZ = A D E F` deletes both `B` and
`C`:

```text
raw:   DEL{B}  (ref bp [100,160))   DEL{C}  (ref bp [160,210))
gap = 160 − 160 = 0 ≤ 100, same type   →   merge → DEL{B,C}
```

### Step 4 — cluster across haplotypes (transitive single-linkage)

Every surviving event (size `≥` a rescue floor) becomes a node; two are joined when `events_match`:

```text
events_match(a,b) =  same SVTYPE
                  AND |a.ref_pos − b.ref_pos| ≤ merge_distance_bp + min(a.size, b.size)   (position window)
                  AND ( weighted_jaccard(a.nodes, b.nodes) ≥ --merge-jaccard (0.80)        (node-set overlap)
                        OR  seq_identity(a.seq, b.seq) ≥ --merge-seq-identity (0.80),       (sequence)
                            gated by a length ratio so wildly different sizes aren't compared )
```

Clusters are the **connected components** of that graph (union-find), so merging is transitive — *not*
greedy first-fit. Two illustrations:

**(a) identical deletions** — `hX:DEL{C}`, `hW:DEL{C}`: Jaccard `|{C}|/|{C}| = 1.0 ≥ 0.80` → one edge →
one record, carriers `{hX,hW}`, `AC=2`. A third haplotype `hZ:DEL{B,C}` has Jaccard `|{C}|/|{B,C}| = 0.5
< 0.80` and no other link → its **own** record.

**(b) transitive merge of STR alleles** — three haplotypes insert the same motif at the same locus but at
different lengths; the node-set Jaccard is low (different nodes) but `seq_identity` is high:

```text
e1 (30 bp)  ~ e2 (33 bp) : identity 0.97, size ratio 30/33=0.91  → edge
e2 (33 bp)  ~ e3 (36 bp) : identity 0.96, size ratio 33/36=0.92  → edge
e1 (30 bp)  ~ e3 (36 bp) : size ratio 30/36=0.83 → not directly compared/edged
components: {e1, e2, e3}  →  ONE locus  (e1–e2–e3 chain), carriers AC=3
```

Even though `e1` and `e3` are never directly matched, single-linkage unites all three via `e2` — so the
same biological STR does not fragment into separate records across haplotypes.

### Step 5 — copy number

DUP/CN is computed off the dropped REP/multiplicity nodes (e.g. `--cn-from-multiplicity` reads the REP
loop count; `--cn-from-coverage` uses spelled-bp / unit on folded paralog clusters), emitted with
`REF_CN` and per-sample `CN`.

**Resulting records (sketch):**

```text
SVTYPE  EVENT_NODES   carriers       AC   note
DEL     C             hX, hW         2    Jaccard merge (a)
DEL     B,C           hZ             1    coalesced in hZ; distinct site
INS     <STR motif>   h1, h2, h3     3    transitive seq-identity merge (b)
DUP     <REP>         hY (CN=2)      1    from multiplicity
```

---

## describe — k-mer / closed-syncmer sketching

Encoding: `A=0, C=1, G=2, T=3`; each k-mer is **canonicalised** to `min(forward, reverse-complement)` of
its 2-bit code, so a k-mer and its reverse complement count as one feature.

**Canonicalisation example** (k = 4): the k-mer `GACT`

```text
forward  GACT = 2 0 1 3  (2-bit) = 135
revcomp  AGTC = 0 2 3 1  (2-bit) =  45     -> canonical = min(135, 45) = AGTC
```

**Closed-syncmer selection** (k = 4 → `s = max(1, min(11, (4+2)/3)) = 2`). A canonical k-mer is **kept**
when its smallest internal s-mer sits at **either end** (offset `0` or offset `k−s = 2`); dropped if the
minimum is in the middle. Trace on the sequence `GACACG`:

```text
k-mer (canonical)   its 2-mers (offset 0,1,2)     min 2-mer at      keep?
GACA                GA, AC, CA                     AC @ offset 1     DROP (middle)
ACAC                AC, CA, AC                     AC @ offset 0     KEEP (end)
CACG                CA, AC, CG                     AC @ offset 1     DROP (middle)
```

So `GACACG` contributes **one** syncmer feature: `ACAC`.

**Other sampling modes.** `--feature-mode all` keeps every canonical k-mer (3 here: `GACA, ACAC, CACG`).
`--feature-mode minimizer` slides a window of `--minimizer-window` consecutive k-mers and keeps the
**smallest canonical** k-mer in each window (a position chosen by overlapping windows is counted once).
With window `W = 2` over the canonical codes `GACA(132), ACAC(17), CACG(70)`:

```text
window [GACA,ACAC] → min = ACAC (17)
window [ACAC,CACG] → min = ACAC (17)   (same position, counted once)
kept minimizers: { ACAC }
```

### Counting + the discriminative filter

Each retained feature is counted **per haplotype** (counts = multiplicity, so copy-number expansions stay
faithful). A feature is then kept by a **two-part rule** (`--min-paths N`, default `1`):

1. **Copy-number features are always kept** — if the count *varies among the paths that carry it*
   (`min_nonzero_count ≠ max_count`), it tracks copy number and is informative even if present everywhere.
2. **Otherwise a symmetric minor-presence (MAF) cut** — drop when `min(present_paths, absent_paths) ≤ N`.
   This removes features with no contrast: present-in-all-with-one-identical-count, singletons (one path),
   and all-but-one.

Worked example over 5 haplotypes (`--min-paths 1`):

```text
feature   counts per hap        present / absent   min_nz≠max?   verdict
K1=ACAC   1, 1, 3, 1, 0         4 / 1              3≠1 → yes     KEEP  (copy-number: 3× in hap3)
K2        1, 1, 1, 1, 1         5 / 0              1=1 → no      DROP  (no contrast: all equal)
K3        2, 0, 0, 0, 0         1 / 4              2=2 → no      DROP  (singleton: min(1,4)=1 ≤ 1)
K4        1, 1, 1, 1, 0         4 / 1              1=1 → no      DROP  (all-but-one: min(4,1)=1 ≤ 1)
```

`K1` survives on rule 1 (its count varies — it tracks the 3× expansion in `hap3`); the rest carry no
contrast. `--min-paths 0` disables rule 2 (keeps every feature that is non-constant).

### Node provenance

Because k-mers are spelled from the bubble's node sequences, each kept k-mer records the set of
**graph nodes** its occurrences touch (written in the feature map). So a significant marker downstream
maps straight back to a node/edge in the bubble — and to the variant `call` makes there. The pooled
survivors across bubbles are written to the pyseer-ready `fsm_kmers.txt.gz`.
