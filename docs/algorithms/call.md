# call — algorithm, copy-number mechanics & worked example

Mechanism and hand-traced examples for **Module 3**. For usage/flags see [modules/call.md](../modules/call.md);
citations: [references.md](../references.md#call).

## Terms

- **anchor** — a node token unique in *both* the reference walk and a haplotype walk; a monotonic (LIS)
  subsequence of anchors aligns the two walks, and only the segments *between* anchors are aligned (fast,
  localised).
- **gap block** — a maximal stretch between two anchors where reference and haplotype differ; each becomes
  one typed event (DEL/INS/INV/substitution).
- **full walk** — the widest source→sink span with every repeated traversal counted (a node visited twice
  contributes its length twice); contrast the *minimal-span* walk that visits each distinct node once.
- **peak multiplicity** — the visit count of a haplotype's single busiest node in the bubble.

## Algorithm

For each bubble:

1. **Alleles.** Group identical source→sink walks into distinct alleles; call once per allele (genotypes
   expand by membership).
2. **Copy-number nodes.** Fold a node's consecutive self-repeats (a `REP` self-loop) into one alignment
   token, so extra copies surface as copy number, not a spurious INS/DEL.
3. **Diff vs reference.** Split both walks at shared anchors and align only between them; each gap block →
   DEL (reference-only nodes), INS (haplotype-only), INV (reverse-complement node-walk of the reference
   run), or a substitution (co-located DEL+INS sharing one `EVENTID`).
4. **Coalesce within a haplotype.** Merge consecutive same-type events whose gap is ≤ `--merge-distance-bp`,
   measured in **either** reference space or the haplotype's own sequence space.
5. **Merge across haplotypes.** Transitive single-linkage (connected components), seeding from events down
   to `--rescue-min-bp`; two events link on same type + a position window + (Jaccard ≥ `--merge-jaccard`
   **or** sequence identity ≥ `--merge-seq-identity`, gated by `--merge-size-ratio`). Largest member
   represents; `MERGE_*` records the evidence. Copy-number records merge separately on shared `REP` identity.
6. **Force-call, then filter.** Re-test every non-carrier against its own diff, adding it as a carrier when
   it matches (sub-threshold events get `GT=1` instead of `0`); one pass (representative is fixed). Keep
   records whose representative reaches `--min-sv-bp` and carrier count reaches `--min-haplotypes`.

## Worked trace — diff, coalesce, cluster, rescue

Reference walk `R = A B C D E F`. Copy-number nodes (REP self-loops) are dropped from this token walk —
typed separately as DUP — so extra repeat copies are never mistyped as insertions.

**Read events off gap blocks** (R fixed; only the haplotype changes):

```text
haplotype          gap block (ref vs hap)      event
A B   D E F        ref=[C]  hap=[]        →     DEL {C}
A B P Q C D E F    ref=[]   hap=[P,Q]     →     INS {P,Q}   (anchored after B)
A B c̅ D E F        hap=[c̅]=revcomp(C)     →     INV over {C}
A B X D E F        ref=[C]  hap=[X]       →     DEL{C}+INS{X}, one EVENTID (substitution)
```

**Coalesce within a haplotype** (`--merge-distance-bp 100`):

- *reference space:* `hZ = A D E F` → `DEL{B}` (ref bp [100,160)) + `DEL{C}` ([160,210)); gap 0 ≤ 100 →
  `DEL{B,C}`.
- *haplotype space:* a haplotype deletes node 2 (200 bp) and inserts `INS5`, `INS6` around it. On the
  reference they are 210 bp apart (the deletion inflates the gap → reference clause FAILS), but on the
  haplotype walk node5 ends at pos 70 and node6 starts at 80 → gap 10 ≤ 100 → **merge → `INS{5,6}` (120 bp)**.
  (The 200 bp `DEL{2}` stays its own event — the hap-space metric is defined only for INS.)

**Cluster across haplotypes (transitive single-linkage):**

- *identical deletions:* `hX:DEL{C}`, `hW:DEL{C}` → Jaccard 1.0 → one record, `AC=2`. A third `hZ:DEL{B,C}`
  has Jaccard `|{C}|/|{B,C}|=0.5 < 0.80` → its own record.
- *STR allele chain:* `e1(30bp)~e2(33bp)` identity 0.97, ratio 0.91 → edge; `e2~e3(36bp)` 0.96, 0.92 →
  edge; `e1~e3` ratio 0.83 not compared. Single-linkage unites `{e1,e2,e3}` via e2 → one locus `AC=3`.

**Force-call (sub-threshold rescue).** Clustering seeds only from events ≥ `--rescue-min-bp` and a record is
kept only if its representative reaches `--min-sv-bp`. With the floor at 50: `hap_big` INS 60 bp seeds the
record; `hap_small` INS 48 bp (same sequence) is < floor, excluded from clustering. Force-call re-reads
`hap_small`'s own diff: identity high, length ratio 48/60 = 0.80 ≥ gate → MATCH → add carrier. Result: one
record `SVLEN=60, SVLEN_RANGE=48,60, NMERGED=2, AC=2`. It demands the haplotype's *own* diff to register the
event (not mere node-set containment), so a *dissimilar* 48 bp insertion at the same spot is **not** rescued
— deliberate on folded paralog loci where inserted nodes are shared across nearly all walks.

**When the merge parameters change the call:**

```text
--merge-distance-bp   two same-type DELs 90 bp apart:   =100 → one DEL ; =50 → two DELs
--merge-jaccard       DEL{2}(60bp) vs DEL{2,3}(90bp), J=0.67:  =0.80 → two ; =0.60 → one (MERGE_JACCARD=0.67)
--merge-seq-identity  two 60bp INS, different nodes, 92% identical:  =0.95 → two ; =0.80 → one (MERGE_SEQID=0.92)
--merge-size-ratio    same motif at 60 vs 108 bp:  default floor 0.80 → 0.56 too different, two ;
                                                   =0.5 → compared → one (SVLEN_RANGE=60,108)
```

## Copy number: one method per locus topology

How a locus is represented in the pangenome decides how `call` reads its copy number — and which graph it
reads. This table is the canonical reference (the module pages link here):

| Region | Topology | Call substrate | CN method | Concordance vs truth |
|--------|----------|----------------|-----------|----------------------|
| **LPA** (KIV-2) | tandem repeat | `panphorte` graph | self-loop `REP` (always on) | **465/465 = 100%** |
| **C4** (RCCX) | PGGB-collapsed paralog | `bubble` graph | `--cn-from-coverage` (full-walk) | **131/131 = 100%** |
| **GSTM1** | deletion/CNV (segdup) | `bubble` graph | `--cn-from-multiplicity` (peak) | **159/159 = 100%** |
| **CYP2D6** | PGGB-collapsed paralog | `bubble` graph | `--cn-from-coverage` (full-walk) | concordant vs D6+D7; residual = unannotated CYP2D8P/hybrid |

The principle: PGGB collapses **identical** copies onto shared nodes (copy number = node multiplicity);
`panphorte` collapses a **variable tandem** into one REP node. So tandem loci are called on the `panphorte`
graph, PGGB-folded paralog clusters on the unfolded `bubble` graph; the mechanics of each method are in
[Copy number — three ways](#copy-number--three-ways) below.

## Copy number — three ways

Copy number is read straight off the walk (no re-alignment), tried per bubble in a fixed precedence so it is
never double-counted. The choice depends on locus topology — see the
[table above](#copy-number-one-method-per-locus-topology).

**(a) Clean tandem → self-loop `REP` (always on).** `panphorte` collapsed an adjacent tandem into a `REP`
self-loop; CN is the loop count. `REF_CN` = reference loops; a sample below it is a loss, above it a gain.
panphorte folds single copies too, so a one-copy haplotype reads `CN=1`; only a haplotype too divergent to
fold reads `CN=0`.
```text
reference REP ×3 → REF_CN=3 ;  sample hY REP ×2 → CN=2 (one copy lost)
```

**(b) PGGB-collapsed paralogs, reference folds ≥2× (`--cn-from-coverage`).** Near-identical paralogs
(CYP2D6/2D7) fold onto shared nodes, so we can't spot a loss from any single node; we measure **total
sequence over the full walk** and divide by one copy (calibrated from the reference):
```text
one-copy bp = ref_spelled_bp / ref_fold = 10,000/2 = 5,000 ;  REF_CN = ref_fold = 2
hapA 15,000 bp → 3 (gain) ;  hapB 5,000 → 1 (loss) ;  hapC 10,000 → 2 (ref-like)
```
Reports the **absolute** total-module count (2D6+2D7 together), recovers losses and gains alike; when it
fires it is the authority for that bubble.

**(c) Single folded extra copy, reference does not fold (`--cn-from-multiplicity`).** Reference visits every
node once (`ref_peak=1`); a haplotype with an extra copy folds it back, so its busiest node is visited
twice. No division — compare the **peak**:
```text
ref_peak=1, REF_CN=1 ;  hapD peak 2 > 1 → DUP, CN=2, SVLEN = Σ node_len × excess ;  hapE peak 1 → no DUP
```
The peak (not per-node excess) isolates real dosage from cluster background.

**Precedence/composition.** Per bubble: coverage CN (if `--cn-from-coverage` and the reference folds ≥2×) is
the authority; else self-loop `DUP`; else peak-multiplicity `DUP` (`--cn-from-multiplicity`). Passing both
flags is safe — they cover disjoint topologies. With no CN flag, an extra copy surfaces as an `INS`
(`INS_SUBTYPE=DUP` under `--classify-ins`). The two folded-cluster detectors report **absolute** per-haplotype
CN (reference only sets the unit-bp denominator), unlike the reference-relative DEL/INS/INV diff.

## Merge keys — Jaccard vs sequence identity

Two events merge across haplotypes on the position window + *either* node-set Jaccard (`--merge-jaccard`) or
sequence identity (`--merge-seq-identity`; tried only if Jaccard misses). They fail in orthogonal ways:

- **same content, different nodes** → Jaccard low, sequence high (one allele threaded through different
  graph nodes — a microsatellite tangle). The sequence gate rescues it (its main reason to exist).
- **same nodes, poorly-aligning sequence** → Jaccard high, sequence low (shared backbone with a big
  internal indel / low-complexity content). Jaccard rescues it (length-weighted, order/orientation-blind).
- **different sizes** → only the sequence path is gated by `--merge-size-ratio`; Jaccard has no size gate.

So lower `--merge-jaccard` to merge events sharing a backbone; lower `--merge-seq-identity`/`--merge-size-ratio`
to merge similar content threading different nodes.

## Multiallelic mechanics (`--multiallelic-loci`)

Collapses a small, bounded bubble that varies mainly by *which sequence* a haplotype carries (e.g. a short
STR with several length alleles) into **one** record: `REF` + `ALT1,ALT2,…` (one per distinct interior
spelling), per-sample `GT` indexing the allele, `NALLELES` counting them. Pure DEL/INS/INV bubbles only — a
bubble that yields a copy-number record (self-loop/coverage/multiplicity `DUP`) is left typed, so DUPs keep
`REF_CN`/`CN`. Fires only when the locus is ≤ `--multiallelic-max-bp` (5000) and the largest allele differs
from `REF` by ≥ `--min-sv-bp`; otherwise it falls back to per-event records so large SVs stay typed.

## Gene annotation trace (`--gtf`)

The bridge is the reference path's PanSN name (`grch38#1#chr6:31891045-32123783`), giving chromosome + start.

1. **Project genes onto reference nodes.** Walk the reference accumulating node lengths; node *k* spans
   `[start+prefix[k], start+prefix[k+1])`; intersect with gene intervals → `node_id → genes`. A node the
   reference visits at two coordinates (a folded paralog) is tagged with **both** genes — which is why
   multiplicity alone can't separate them.
2. **Per-gene CN by competitive realignment (`call`).** Only for a DUP overlapping **≥2 genes** — where
   graph multiplicity reports the module total but not the per-gene split. Each gene's reference sequence is
   realigned (minimap2, all hits) to each haplotype; same-gene hits that are *target-adjacent but
   query-disjoint* chain into one copy (a copy split around a missing HERV counts once, via a gap-compressed
   identity floor); each copy competes across genes by block identity and the winner claims the locus. A DUP
   over a single gene skips realignment entirely: the module total *is* that gene's copy number (so e.g. LPA,
   one DUP over `LPA` across all haplotypes, costs no per-haplotype alignment).
3. **Reliability (collapse groups).** Cluster genes where one aligns to another at > 98% block identity over
   most of its length:
   ```text
   CYP2D6 vs CYP2D7 : 0.95 < 0.98 → separable        C4A vs C4B : 0.9992 > 0.98 → collapse {C4A,C4B}
   ```
   Singleton group ⇒ `reliable=1` (report the per-gene split; CYP2D6 96.9% / CYP2D7 98.4% exact vs truth).
   Group > 1 ⇒ `reliable=0`: emit one row with the collapsing genes (`genes=C4A;C4B`) and the **module
   total** (`FORMAT:CN`) — validated 131/131; only the A-vs-B label is untrustworthy, not the total.
