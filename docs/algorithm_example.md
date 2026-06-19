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
Real repeats (*e.g.*, KIV-2 in LPA) ofter differ by SNVs or small indels, so the exact pair test fails and nothing is seeded. With
`--min-similarity 0.90`, panphorte switches to a **single-representative, lossy** mode: it picks one unit
`U`, then finds its near-copies by multi-seed **banded alignment** of `U` (and its reverse complement)
against the path's spelled sequence, accepting a copy when `identity ≥ 0.90` and at least half of `U` is
consumed:

```text
copy 1: identity 0.98 ✓   copy 2: 0.95 ✓   copy 3: 0.93 ✓   copy 4: 0.97 ✓   → 4 copies → one REP=U
```

Because alignment (not exact match) bridges a copy, a large *internal indel* between divergent copies is
also bridged when the threshold is low enough — this is what lets a C4 short module (missing the ~6.4 kb
HERV-K) fold onto the long one at `--min-similarity 0.70`. The collapse is lossy: the per-copy small events are
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

### Step 2 — read each event off a gap block

Each maximal gap block between two anchors becomes **one typed event**. Reference walk `R = A B C D E F`
in every row; only the haplotype changes:

```text
haplotype          gap block (ref vs hap)         event emitted
A B   D E F        ref=[C]   hap=[]         →      DEL {C}                       (deletes C)
A B P Q C D E F    ref=[]    hap=[P,Q]      →      INS {P,Q}     anchored after B
A B c̅ D E F        hap=[c̅] == revcomp(C)   →      INV over {C}                  (orientation flip)
A B X D E F        ref=[C]   hap=[X]        →      substitution: DEL{C}+INS{X}, one shared EVENTID
```

- **DEL** = reference-only nodes the haplotype skips; **INS** = haplotype-only nodes, pinned to the last
  matched reference node.
- **INV** fires only when the haplotype block is the exact reverse-complement node-walk of the reference
  block — an orientation flip, not merely different sequence.
- **substitution** = a block with content on *both* sides; emitted as a linked DEL + INS so each side
  keeps its own size and sequence, joined by one `EVENTID`.

### Step 3 — coalesce fragments within a haplotype

Consecutive **same-type** events merge when their gap is `≤ --merge-distance-bp` (default 100) in
**either** of two coordinate spaces:

- **reference bp** — `cur_lo − prev_hi`, using each node's reference start;
- **haplotype bp** — the inserted nodes' own positions in *this* haplotype's walk.

**(a) reference space.** `hZ = A D E F` deletes both `B` and `C`:

```text
raw:   DEL{B}  (ref bp [100,160))   DEL{C}  (ref bp [160,210))
gap = 160 − 160 = 0 ≤ 100, same type   →   merge → DEL{B,C}
```

**(b) haplotype space — two insertions "far apart on the reference because a deletion sits between
them."** An insertion has no width on the reference, so each is pinned to a single reference
coordinate (the end of the node it is anchored to). Picture a haplotype that, in one stretch, both
**deletes** a large reference block and makes **two** insertions:

```text
reference:    1[0–10]  2[10–210]  3[210–220]  4 (sink)        node 2 = 200 bp, node 3 = 10 bp
hap walk:     1 — INS5(60bp) — 3 — INS6(60bp) — 4             node 2 deleted; 5,6 inserted
events:       DEL{2}            INS{5}            INS{6}       (deleted block flushes first)
```

`INS5` is anchored after node 1 → reference point **10**; `INS6` is anchored after node 3 → reference
point **220**. The 200 bp deletion of node 2 sits between them and inflates the reference gap:

```text
reference gap : 220 − 10  = 210  > 100   → reference clause FAILS
haplotype gap : node5 ends at hap pos 70, node6 starts at 80 (only node3=10bp between)
                80 − 70   = 10   ≤ 100   → haplotype clause FIRES  → merge → INS{5,6} (120 bp)
```

The two inserted blocks are one near-contiguous insertion in the sample; only the intervening deletion
made them look distant on the reference. (`DEL{2}` stays its own event — the hap-space metric is
defined only for INS, whose nodes actually lie on the haplotype walk; DEL/INV nodes are reference-only,
so they fall back to the reference gap.) On this worked input nodes 5 and 6 emerge as **one**
`INS SVLEN=120 EVENT_NODES=5,6`, with the 200 bp DEL reported separately.

### Step 4 — cluster across haplotypes (transitive single-linkage)

Every surviving event (size `≥` a rescue floor) becomes a node; two are joined when they **match**:

```text
match(a, b) =  same SVTYPE
            AND |a.pos − b.pos| ≤ --merge-distance-bp + min(a.size, b.size)   (position window)
            AND ( node-set overlap (length-weighted Jaccard) ≥ --merge-jaccard (0.80)
                  OR  sequence identity ≥ --merge-seq-identity (0.80),
                      but only when the size ratio clears --merge-size-ratio )
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

### Step 5 — force-call (sub-threshold rescue)

Clustering (Step 4) only seeds from events `≥` a **rescue floor** (`--rescue-min-bp`, default
`min_sv_bp / 2`), and a record is finally kept only if its **representative** — the largest member —
reaches `--min-sv-bp`. That leaves a gap: a haplotype may carry a genuine but *small* version of a
variant that, on its own, never clears the threshold. The force-call pass closes it. For every kept
record, it re-reads **every non-member haplotype's own events** (the full per-haplotype event list is
retained, not just the `≥ floor` subset) and adds the haplotype as a carrier when its walk-diff
**matches** the record's fixed representative.

Worked example — two haplotypes insert the same motif at one anchor, `hap_big` 60 bp and `hap_small`
48 bp (same sequence). With the rescue floor raised to 50, the 48 bp event is **excluded** from the
clustering candidates, so only force-call can recover it:

```text
record seeded from hap_big        : INS 60 bp   (representative, ≥ min_sv 50)
hap_small INS 48 bp               : < floor → not a clustering candidate
force-call: does the 60bp seed match the 48bp event?
   sequence identity high, length ratio 48/60 = 0.80 ≥ gate  → MATCH
   → add hap_small as carrier
result: ONE record, SVLEN=60, SVLEN_RANGE=48,60, NMERGED=2, AC=2  (hap_small GT=1)
```

Two properties make this safe:

- **It demands the haplotype's *own* diff to register the event** — not mere node-set containment. A
  *dissimilar* 48 bp insertion at the same spot does **not** match, so it is **not** rescued (`AC` stays
  `1`, `hap_small` GT=`0`). This is deliberate: on folded paralog loci the inserted nodes are shared
  across nearly all walks, so a pure containment test would force-call essentially every haplotype.
- **It is monotone and single-pass** — the representative is fixed before this pass, so adding carriers
  never changes what any other haplotype is tested against; one sweep reaches the fixpoint with no
  re-alignment.

### Step 6 — copy number

Copy number is read straight off the walk (no re-alignment), one of three ways depending on how the
locus is folded in the graph. `REF_CN` and per-sample `CN` are reported on a `DUP` record.

**(a) clean tandem (self-loop `REP`).** `panphorte` already collapsed an adjacent repeat into a `REP`
node with a self-loop; copy number is just the loop count.

```text
reference walks REP 3×   → REF_CN = 3
sample hY walks REP 2×   → CN = 2   (one copy lost)
```

**(b) folded paralogs, reference has several copies (`--cn-from-coverage`).** Two near-identical genes
(e.g. CYP2D6 / CYP2D7) are folded onto the **same shared nodes**, so a path carrying both genes walks
those nodes more than once — including the reference. We can't spot a *loss* from any single node (a
deleted copy needn't touch the busiest one), so we measure **total sequence** and divide by one copy.
The bp are summed over the **full walk** — the widest source→sink span with every repeated traversal
counted (a node visited twice contributes its length twice). This is essential: the ordinary
minimal-span walk visits each distinct node once and so collapses the folded copies, flattening every
haplotype to the same bp. The size of one copy is calibrated from the reference itself:

```text
reference spells 10,000 bp over its FULL walk, folding ref_fold = 2× over the shared nodes
  unit (one copy) = ref_spelled_bp / ref_fold = 10,000 / 2 = 5,000 bp
  REF_CN          = ref_fold                                = 2

per haplotype:  copies ≈ (bp it spells over the FULL walk) / unit
  hapA spells 15,000 bp  → 15,000 / 5,000 = 3   → CN = 3   (gain)
  hapB spells  5,000 bp  →  5,000 / 5,000 = 1   → CN = 1   (loss)
  hapC spells 10,000 bp  → 10,000 / 5,000 = 2   → CN = 2   (reference-like)
```

Because it uses *all* the traversed sequence, it recovers losses and gains alike, and the reported `CN` is
the haplotype's **absolute** module count (the reference only sets the unit-bp denominator). It reports the
**total module** count (2D6 + 2D7 together), not which paralog is which. When it fires it is the authority
for that bubble — the self-loop and walk-diff paths are skipped.

**(c) folded single extra copy, reference does not fold (`--cn-from-multiplicity`).** Here the reference
is single-copy — every node visited once (`ref_peak = 1`). A haplotype carrying an **extra** copy folds
it back onto the shared nodes, so its busiest node is visited twice. No division: just compare the
**peak** node-visit count.

```text
reference visits every node 1×          → ref_peak = 1, REF_CN = 1
hapD visits its busiest node 2×         → ac_peak = 2 > ref_peak  → DUP, CN = 2
   SVLEN = duplicated content = Σ node_len × excess visits
hapE visits every node 1×               → ac_peak = 1 = ref_peak  → no DUP
```

The **peak** (single busiest node), not per-node excess, is the dosage signal: in a messy cluster many
nodes pick up extra visits just from paralog presence/absence, but the global peak reflects a real
gained copy.

**Resulting records (sketch):**

```text
SVTYPE  EVENT_NODES   carriers       AC   note
DEL     C             hX, hW         2    Jaccard merge (a)
DEL     B,C           hZ             1    coalesced in hZ; distinct site
INS     <STR motif>   h1, h2, h3     3    transitive seq-identity merge (b)
DUP     <REP>         hY (CN=2)      1    self-loop loop count (6a)
DUP     <module>      hA,hB (CN=3,1) 2    coverage total-module CN (6b)
DUP     <peak node>   hD (CN=2)      1    peak multiplicity (6c)
```

### Step 7 — when the merge parameters change the call

Each knob flips one specific decision. When a record actually merges ≥2 events, the evidence is written
back as `MERGE_JACCARD` / `MERGE_SEQID` / `MERGE_SIZE_RATIO`, so you can see *why* it merged.

```text
--merge-distance-bp   two same-type DELs 90 bp apart on the reference (within one haplotype):
                        =100 (default) → 90 ≤ 100 → coalesced into ONE DEL  (Step 3)
                        =50            → 90 > 50  → TWO separate DELs

--merge-jaccard       DEL{2} (60bp) and DEL{2,3} (90bp) share node 2 → length-weighted Jaccard 60/90 = 0.67:
                        =0.80 (default) → 0.67 < 0.80 → two records
                        =0.60          → 0.67 ≥ 0.60 → ONE record   (MERGE_JACCARD=0.67)

--merge-seq-identity  two 60bp insertions, different nodes (Jaccard 0), sequences 92% identical:
                        =0.95          → 0.92 < 0.95 → two records
                        =0.80 (default)→ 0.92 ≥ 0.80 → ONE record   (MERGE_SEQID=0.92)

--merge-size-ratio    same motif inserted at 60 bp and 108 bp (Jaccard 0, identity high):
                        default (floor = --merge-seq-identity = 0.80): 60/108 = 0.56 < 0.80
                            → sizes too different to even compare → two records (size classes kept apart)
                        =0.5           → 0.56 ≥ 0.5 → compared → ONE record spanning both
                                         (SVLEN_RANGE=60,108, MERGE_SIZE_RATIO=0.56)
```

The `MERGE_*` fields appear only on records that merged at least two events; a singleton record has none.

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

**The other mode.** `--feature-mode all` keeps every canonical k-mer (3 here: `GACA, ACAC, CACG`) — the
exhaustive set, versus the one sampled syncmer above.

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
survivors across bubbles are written to the fsm-lite `fsm_kmers.txt.gz`.
