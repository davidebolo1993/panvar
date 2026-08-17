# Module `call` - algorithm

Mechanism for the `call` module. For usage/flags see [modules/call.md](../modules/call.md); References in [references.md](../references.md#call).

`call` types each haplotype's differences from the reference within a bubble into a multi-sample VCF, and with `--cn` reads copy number off the walk as a `DUP`. It pins the two walks together at their anchors — node tokens occurring exactly once in both the reference and the haplotype walk, so each must mark the same position, chained in the order they appear — and aligns only the gap blocks between consecutive anchors, the maximal stretches where the walks differ. Each gap block becomes one typed event: a `DEL`, `INS`, `INV`, or a substitution (a co-located `DEL` and `INS` sharing one `EVENTID`).

## How it works

For each bubble:

### 1. Group alleles

Identical `source → sink` walks are grouped into distinct alleles and called once per allele; genotypes expand by membership.

### 2. Fold copy-number loops

A node's consecutive self-repeats (a `REP` self-loop) fold to one alignment token, so extra copies surface as copy number rather than a spurious `INS`/`DEL`.

### 3. Diff at anchors into typed events

Both walks are split at their shared anchors and only the gap block between consecutive anchors is aligned. Each becomes one typed event: a `DEL` where only the reference has nodes, an `INS` where only the haplotype does, an `INV` where the haplotype's run is the reverse-complement node-walk of the reference's, or a substitution — a co-located `DEL` and `INS` that share one `EVENTID`.

### 4. Coalesce within a haplotype

Consecutive same-type events whose gap is within `--merge-distance-bp` are joined into one, measured in either reference space or the haplotype's own sequence space.

### 5. Merge across haplotypes

Events merge into one record per site by transitive single-linkage (connected components), seeding from the events down to `--rescue-min-bp`. Two events link when they are the same type, sit within a position window, and either their node sets overlap or their sequences align. The window is `--merge-distance-bp` widened by the smaller event's size, so one large indel whose breakpoint wanders by kilobases across haplotypes still reads as one site; the node/sequence gate below keeps that fuzziness from over-merging genuinely distinct events. The overlap test is a length-weighted node-set Jaccard ≥ `--merge-jaccard`; if that misses, the fallback is sequence identity ≥ `--merge-seq-identity` (length-gated by `--merge-size-ratio`). The largest member represents the record and `MERGE_*` records the evidence. Copy-number (`DUP`) records do not use this path — they merge separately on a shared `REP` self-loop node.

### 6. Force-call and filter

Every non-carrier is re-tested against its own diff and added as a carrier when it matches (sub-threshold events get `GT=1` instead of `0`), in one pass since the representative is fixed. A record is kept only if its representative reaches `--min-sv-bp` and its carrier count reaches `--min-haplotypes`.

#### Merge keys — Jaccard vs sequence identity

The two link tests in step 5 fail on opposite kinds of event, which is why both exist:

- same content, different nodes — node Jaccard low, sequence identity high: one allele threaded through different graph nodes. The sequence gate rescues it — its main reason to exist.
- same nodes, poorly-aligning sequence — Jaccard high, sequence identity low: a shared backbone with a large internal indel or low-complexity content. The Jaccard gate rescues it (length-weighted, order- and orientation-blind).
- different sizes — only the sequence path is size-gated (`--merge-size-ratio`); Jaccard has no size gate.

So lower `--merge-jaccard` to merge events that share a backbone, and lower `--merge-seq-identity` / `--merge-size-ratio` to merge similar content threaded through different nodes.

#### What single-linkage can hide, and how to see it

Merging is transitive: A and C join whenever some B links to both, even when A and C would never have merged directly. That is deliberate — greedy first-fit fragments a real event across haplotypes — but it means a record can span members far more different than any threshold allows. `MERGE_JACCARD` reports the **strongest** edge and cannot show this. `MERGE_DIAMETER` reports the **weakest** pairwise node Jaccard in the component, which can:

| record | `NMERGED` | `MERGE_JACCARD` | `MERGE_DIAMETER` |
|---|---|---|---|
| lpa bubble8_INS | 132 | 1.0000 | **0.0000** |
| gstm1 bubble7_INS | 363 | 1.0000 | 0.4937 |
| ankrd36c bubble11_INS | 11 | 1.0000 | 0.3333 |
| cyp2d6 bubble5_INS | 59 | 1.0000 | 0.9016 |

The first row is the case to know about: 132 events merged into one record whose strongest edge is perfect and two of whose members **share no nodes at all**. The second is the record that hands 363 carriers an 80 bp representative when their true divergence averages 57 bp — the diameter says those members are genuinely heterogeneous, not a tight cluster.

All-pairs comparison is quadratic, so above 128 members the value is measured against the representative only and `MERGE_DIAMETER_EXACT` is absent, marking it an upper bound.

### 7. Place the record

`POS` is the base **before** the event, not its first affected base — the symbolic convention, so `REF` carries one real reference base and `ALT` is the symbol. The event then occupies `POS+1 … END`, and `END` is the last reference base it spans, inclusive:

| record | `POS` | `END` | `END − POS` |
|---|---|---|---|
| `DEL` / `INV` | base before the event | last deleted/inverted base | `\|SVLEN\|` |
| `INS` | base before the insertion | `POS` | 0 — an insertion spans no reference |
| `DUP`, `CN_SCOPE=REPEAT_UNIT` | last base of the upstream flank | last base of the first copy | `RU_LEN` |
| `DUP`, `CN_SCOPE=COLLAPSED_MODULE` | last base of the near boundary | base before the far boundary | `CN_MODULE_REF_BP` |

The two `DUP` rows differ because the routes count different things. A `REP` unit is one node, so its span is that node's length; a collapsed module's span is the bubble **interior**, which is what `FORMAT:CNBP` sums over — both boundary nodes are excluded from that sum, so `END` stops where the far boundary begins. That identity (`END − POS = CN_MODULE_REF_BP`) holds on every module record across the reference loci and is the cheapest check that a module record is placed correctly.

A `DEL` or `INV` whose first affected base is the region's first base has no preceding base to anchor on and keeps its original anchor.

## Worked trace (non-`DUP` events)

Reference interior `R = A B C D E F` — single-letter node tokens (a `REP` self-loop, if present, is dropped from this walk and typed separately as a `DUP`, so extra copies are never mistyped as an insertion). Six haplotypes cross the bubble:

```text
ref  : A B C D E F
hX   : A B   D E F          C deleted
hW   : A B   D E F          C deleted (same walk as hX)
hZ   : A     D E F          B and C deleted
hI   : A B P Q C D E F      P, Q inserted after B
hS   : A B C X E F          D replaced by X   (substitution)
hV   : A B d c E F          C D inverted      (lowercase = reverse-complement run, order reversed)
```

1. Diff each haplotype against the reference. The tokens that appear once in both walks are the anchors (`A B … E F`, in order); only the gap between two consecutive anchors is aligned, and each gap becomes one typed event:
   - `hX`/`hW`: the gap between anchors `B` and `D` holds only the reference's `C`: `DEL{C}`.
   - `hZ`: the gap between `A` and `D` holds the reference's `B C`: one `DEL{B,C}` (both fall in a single gap block).
   - `hI`: the gap after `B` holds the haplotype-only `P Q`: `INS{P,Q}`.
   - `hS`: the gap between `C` and `E` has reference-only `D` and haplotype-only `X` at once: it is a substitution, emitted as a co-located `DEL{D}` and `INS{X}` sharing one `EVENTID`.
   - `hV`: the gap between `B` and `E` is the reference run `C D` traversed as its reverse-complement node-walk: `INV{C,D}`.
2. Coalesce within a haplotype. Same-type events that sit in different gap blocks but within `--merge-distance-bp` are joined into one, measured either in reference space or in the haplotype's own sequence space (the latter catches insertions that a large intervening deletion would otherwise push apart). None of these trigger it.
3. Cluster across haplotypes by transitive single-linkage. `hX:DEL{C}` and `hW:DEL{C}` have identical node sets (Jaccard 1.0), so they become one record with `AC = 2`; `hZ:DEL{B,C}` shares only `C` (Jaccard `1/2 < --merge-jaccard`), so it stays its own record; `hS`'s substitution pair and `hV`'s `INV{C,D}` match nothing else, so each stays a singleton. Events also link when their sequences match (not just their nodes), which rescues one allele threaded through different nodes — see [Merge keys](#merge-keys--jaccard-vs-sequence-identity).
4. Force-call, then filter. Every non-carrier is re-tested against its own diff and added if it matches — so a sub-threshold copy of an event gets `GT = 1` rather than `0` — and a record is kept only if its representative reaches `--min-sv-bp` and its carriers reach `--min-haplotypes`.

Resulting records:

```text
SVTYPE=DEL  EVENT_NODES=C     AC=2   carriers hX, hW
SVTYPE=DEL  EVENT_NODES=B,C   AC=1   carrier  hZ
SVTYPE=INS  EVENT_NODES=P,Q   AC=1   carrier  hI
SVTYPE=DEL  EVENT_NODES=D     AC=1   carrier  hS   EVENTID=e1
SVTYPE=INS  EVENT_NODES=X     AC=1   carrier  hS   EVENTID=e1   (paired with the DEL above)
SVTYPE=INV  EVENT_NODES=C,D   AC=1   carrier  hV
```

Without `--cn`: the `DEL` / `INS` / `INV` records above are emitted either way — copy number is the only thing that changes. The always-on self-loop route still calls a genuine tandem a `DUP` (a `REP` self-loop needs no flag), but the coverage and peak-multiplicity routes are off: a folded extra copy then surfaces as an ordinary `INS` (flagged `INS_SUBTYPE=DUP` under `--classify-ins`) and a lost copy as a `DEL`, both without a per-sample `CN`. Add `--cn` to turn those into `DUP` records carrying `CN` / `REF_CN`.

## Worked trace (`DUP` events)

`call --cn` reads copy number straight off the walk and emits it as a `DUP` record carrying a per-sample `CN`. How it reads the number depends on how the locus is folded in the graph. 

#### A tandem repeat (self-loop `REP`)

`panphorte` has folded the tandem array into one `REP` node with a self-loop, so a walk's copy number is simply how many times it loops that node (`L`/`R` are flanks):

```text
reference : L REP REP REP R          loops REP ×3
hapX      : L REP REP R              loops REP ×2
hapY      : L REP REP REP REP R      loops REP ×4
```

1. The unit is already a single `REP` self-loop node (≥ `--min-sv-bp`), so its integer loop count is the copy number, exactly. The reference loops it three times, so `REF_CN = 3`.
2. Count each haplotype's loops: `hapX` traverses `REP` twice (`CN = 2`, one copy lost) and `hapY` four times (`CN = 4`, one gained); a haplotype at three loops is reference-like.
3. Emit one `DUP` at the bubble with `REF_CN = 3` and a per-sample `CN` from the loop count. `GT = 1` marks only the carriers whose `CN` differs from `REF_CN`, so `AC`/`AF` stay meaningful.

```text
SVTYPE=DUP  REF_CN=3  CN_METHOD=REP  RU_LEN=<unit bp>   hapX CN=2  hapY CN=4  (reference-like CN=3)
```

### Why a collapsed module carries no SVLEN

The `MODULE_BP` and `PEAK` routes count copies of a module that may hold several distinct paralogs, so a "copy" is not one uniform sequence. Two consequences, both visible in the record:

`RU_LEN` is not emitted there. The unit those routes calibrate against is the **shared** per-copy content — the sequence the reference revisits — and not the paralog-private sequence that travels with each copy. It is reported as `CN_UNIT_BP` so it cannot be read as a per-copy size. At GSTM1 `CN_SHARED_BP=22830` over `CN_REF_FOLD=3` gives 7610, while `CN_MODULE_REF_BP=58249` puts the total per-copy content near 19.4 kb: `(CN − REF_CN) × 7610` understates every carrier by about 2.4×. The ratio is not a constant to correct for — across the reference loci it runs from ×1.01 to ×36.75.

#### What the panel-spacing model rests on

`--cn-unit-spacing` takes the per-copy step from the gaps between the panel's own copy-state clusters instead of from `ref_bp/ref_fold`. That is the right shape — the reference-derived unit is affine, not proportional, at about 1.45 units per real copy — but a gap is one copy only if the two clusters are **adjacent** copy states, and nothing inside the estimator can confirm that. `CN_STEP_SUPPORT=clusters,gaps,dropped` and `CN_STEP_MAX_MULTIPLE` expose how thin the evidence is. Measured on the reference loci:

| locus | step / unit | clusters, gaps, dropped | max gap multiple |
|---|---|---|---|
| acot bubble5 | 1.000 | 3, 2, 0 | 1.00 |
| acot bubble7 | 0.999 | 4, 3, 0 | 1.00 |
| gstm1 bubble6 | 1.454 | 4, 2, 1 | 1.00 |
| cyp2d6 bubble5 | 1.453 | 4, 2, 1 | 1.00 |
| ankrd36c bubble7 | 1.104 | 4, **1**, 2 | 1.00 |
| ankrd36c bubble8 | 1.000 | 3, **1**, 1 | 1.00 |
| **c4 bubble6** | **0.197** | 7, 6, 0 | **4.07** |
| lpa bubble5 | — | estimator declines | — |

Three things to read off it. At acot the step equals the reference unit, so the dosage there really is proportional and both models agree. The affine 1.45 appears only at the two paralog modules. And **c4 is where the assumption fails**: one gap is four times the chosen step, so copy states are missing between observed clusters, and the resulting step is a fifth of the reference unit. C4's copy number is exact on all 131 haplotypes under the reference-ratio model, which is why the spacing model is opt-in rather than the default.

Two records rest on a single gap (`gaps=1`), where the median is one difference with nothing to cross-check it. Singleton clusters are dropped — a lone walk length is as likely to be a mis-folded outlier as a real copy state — and `dropped` counts the haplotypes lost that way.

Where the panel supports no estimate at all, `--cn-unit-spacing` **fails** rather than quietly computing the reference-ratio model under a flag set to avoid it.

#### One module span, and when it is a choice

Everything that measures the module — the folded-set construction, reference and haplotype module bp, `FORMAT:CNBP`, `CN_MODULE_REF_BP`, and the record's `POS`/`END` — takes its interval from one resolver. It deliberately selects the **widest** oriented span (first source to last sink) rather than the tight interval used for ordinary allele diffing, because a module's copies are exactly what lies between the outermost boundaries.

That span is well defined only while each boundary is visited once. When a boundary recurs, first-source-to-last-sink is a *choice*: it can enclose two separate visits and everything between them rather than one module, and every quantity measured over it inherits that. `CN_SPAN_AMBIGUOUS` counts the paths where this applies.

It does not fire on any of the six reference loci, and the reason is structural rather than luck: a node the reference revisits cannot be a cut vertex, so the cactus decomposition promotes it into the enclosing snarl's interior and it never becomes a boundary. The counter is reachable through `--snarls-in` or a hand-written bubbles CSV, where boundaries carry no such guarantee — which is exactly the door that needs the warning.

#### Is the step a property of the locus or of the cohort?

The obvious objection to taking the step from the panel is that it makes the answer depend on who is in the panel. `scripts/validate_spacing_stability.py` tests that from outside the caller, by rewriting the panel rather than instrumenting the estimator: it resamples 70% of the non-reference haplotypes ten times, and separately re-runs with a different path as the reference.

| locus | step across 10 subsamples | CN changes | reference switch (GRCh38 → CHM13) |
|---|---|---|---|
| acot bubble5 | 475..475 (0.0%) | 0 / 1308 | — |
| acot bubble7 | 14477..14481 (0.0%) | 0 / 1308 | step 14476 → 14480, CN 0/467 differ |
| gstm1 bubble6 | 11067..11068 (0.0%) | 0 / 1956 | step 11068 → **711**, CN **1/466** differ |

Cohort dependence is not the problem it looked like: the step does not move under resampling, and no copy number changes. Four to six of every ten replicates **refuse** spacing mode outright, which is the refusal conditions above doing their job on a thinner panel rather than the estimator guessing.

The reference switch is the interesting one. At GSTM1 the calibration constants change completely — `CN_UNIT_BP` 7610 → 710, the step 11068 → 711 — because the two references revisit different node sets (`CN_REF_FOLD` 3 over 22830 bp against 2 over 1420 bp), so they are calibrating against different modules. The copy number nonetheless agrees on 465 of 466 haplotypes. **The calibration constants are reference-dependent; the integer copy number is very nearly not.** That is the reassuring version of the result, and it is worth stating in that order rather than reporting only the second half.

#### The integer fit is not a correctness signal

`CN_ROUND_AMBIGUOUS_FRAC` is the share of haplotypes whose dosage sat more than 0.4 units from a whole number, so the integer came from near a coin flip. It replaces the mean residual as what `--max-cn-model-residual` gates on, because the per-sample residuals are bimodal and the mean describes neither mode — at c4 the mean is 0.13 while **no** sample is ambiguous, and at GSTM1 the mean is 0.30 while **two thirds** are.

The gate is off by default, and the reason is GSTM1: it has the worst fit of any reference locus and its copy number is exact against pangene truth on all 466 haplotypes. The unit is a calibration constant for a heterogeneous paralog module, not one real copy, so sitting near a half-integer is what it does when it is working. Declining on a bad fit would discard correct calls.

`SVLEN` is not emitted there either. At GSTM1 the carriers span −32030, −18396 and +18504 bp around `REF_CN=3`; no single record-level size describes them, so the field is absent rather than misleading and **`FORMAT:CNBP` is the per-sample size**. `END` is the module's reference span, the interval CNBP is measured over. Confirmed externally: svim-asm on the assemblies agrees with CNBP to a median of 6 bp at GSTM1 and 2 bp at LPA KIV-2.

This is the only reading that fires without `--cn` too: the self-loop route is always on, so a `REP` tandem is a `DUP` regardless of the flag.

#### A collapsed paralog cluster (coverage)

Here `pggb` collapsed near-identical paralog copies onto shared nodes, and the reference itself traverses them more than once. No single node can show a loss, so we measure the total sequence over the folded nodes and calibrate one copy from the reference:

```text
folded nodes (reference visits each ≥2×):  reference total = 10,000 bp,  ref folds 2×
hapA : 15,000 bp over the folded nodes
hapB :  5,000 bp
hapC : 10,000 bp
```

1. Calibrate the one-copy unit from the reference: `one-copy bp = reference folded-node bp / ref_fold = 10,000 / 2 = 5,000`, and `REF_CN = ref_fold = 2`.
2. Divide each haplotype's folded-node bp by that unit: `hapA` gives 3 (a gain), `hapB` gives 1 (a loss), `hapC` gives 2 (reference-like). One measure recovers both losses and gains.
3. Emit one `DUP` with `REF_CN = 2` and per-sample `CN` — the absolute total of all the collapsed paralogs together. Measuring over the folded nodes only is what keeps an edit to unique content in the bubble (interstitial sequence, a single-visit paralog) from being misread as copy number.

```text
SVTYPE=DUP  REF_CN=2      hapA CN=3   hapB CN=1   hapC CN=2
```

#### A single-copy-reference duplication (peak multiplicity)

The reference here visits each inside node once — it does not fold — so there is no reference fold for coverage to calibrate against. This reading covers a haplotype that folds an extra copy back onto those nodes (a duplication where the reference is the single/short allele):

```text
reference : each inside node ×1            ref_peak = 1
hapD      : busiest inside node ×2         one extra copy folded back
hapE      : each inside node ×1            no extra copy
```

1. `REF_CN = ref_peak = 1`. There is no unit to divide by, so we compare peak multiplicities directly.
2. Take each haplotype's busiest inside node: `hapD` peaks at 2 (`> 1`) and is interpreted as a gain; `hapE` stays at 1. Using the global peak, not a per-node excess, keeps scattered cluster background from inflating the call.
3. Emit a `DUP` for `hapD` with `REF_CN = 1`, `CN = 2`, and `SVLEN = Σ node_len × excess` (the duplicated bp). This reading is gains only — a loss leaves no extra fold to count.

```text
SVTYPE=DUP  REF_CN=1      hapD CN=2   (hapE is not a carrier)
```


## Gene annotation trace (`--gtf`)

With `--gtf`, `call` records which genes each variant touches and splits a `CN` call that spans several genes into per-gene `CN` — from the GTF (Gene Transfer Format) alone, by k-mer dosage. The reference path's PanSN (Pangenome Sequence Naming) name (`sample#hap#chrom:start-end`) supplies the chromosome and start coordinate the projection needs.

First the genes are projected onto the reference nodes: walking the reference and accumulating node lengths gives each node a coordinate span, intersected with the GTF gene intervals to tag every node with the gene (or genes) it falls in. Where the reference revisits the same nodes at two coordinates — a folded paralog — that node is tagged with both genes, which is precisely why graph multiplicity alone cannot tell the paralogs apart. When a `DUP` overlaps two or more genes, each gene's discriminative sketch is its merged coding sequence (CDS; the k-mers private to it versus its paralogs), and per-haplotype copy number is the dosage of those private k-mers over the module sub-walk.


## Worked trace (assign genes in `DUP` events)

#### Divergent paralogs (private-k-mer dosage)

A folded cluster carries two coding paralogs `A` and `B` plus a pseudogene `P`. Coding sequence is where paralogs differ, so each gene's sketch is its merged CDS — except `P`, a gene with no CDS, which falls back to its gene span:

```text
gene   marker   private k-mers (unique vs the other paralogs)
A      CDS      1200
B      CDS       800
P      span     4800     (no CDS -> gene span)
```

1. Build each gene's canonical k-mer set from its marker and keep the k-mers private to it (absent from every paralog in the cluster). These sets discriminate from the reference alone.
2. Scan a haplotype's module sub-walk and count each gene's private k-mer **occurrences**. A gene present in `c` copies carries each private k-mer about `c` times, so `dosage = hits / occurrences-per-reference-copy ≈ c`, rounded to the copy number. The denominator counts occurrences rather than distinct k-mers, because a k-mer appearing twice inside one copy contributes two hits from that copy and a distinct count would read it as extra copies.
3. Emit one row per gene with its evidence. Keeping a CDS-less gene `P` in the set matters even when its own copy number is not separately validated: its private k-mers absorb its locus, so one of its copies is not mis-counted as a coding paralog.

`reliable` is decided per haplotype, not per gene pair. A pair that is separable somewhere in the panel is not necessarily separable in a given haplotype, and reporting `CN=0` from `hits=0` over an empty private set as a confident call was wrong; such a haplotype falls back to the module total instead.

**Measured and rejected: screening private k-mers against the whole locus.** Privacy is judged against the sibling genes only, so a k-mer occurring elsewhere in the locus still counts as discriminative. Widening it — dropping any k-mer whose locus occurrences are not all inside the gene's own marker — is the obvious correction and it is worse against pangene truth:

| variant | GSTM1/2/4/5 exact | CYP2D6 exact |
|---|---|---|
| sibling-only (shipped) | 466/466 each | 126/127 |
| + locus screen | 466/466 each | 123/127 |
| + locus screen, keeping CDS-junction k-mers | 464/466 | 122/127 |

It removes about a third of the markers, and the dosage ratio loses more to the smaller denominator than it gains in specificity. The second row exists because a marker is a gene's *merged CDS*, so a k-mer spanning an exon-exon junction occurs nowhere in the genomic reference; keeping those rather than discarding them as absent made the result worse still, which rules out the obvious repair. Not retained in any form.

```text
one haplotype:  A  hits=990/1200  dosage=0.83  -> CN=1
                B  hits=560/800   dosage=0.70  -> CN=1
                P  hits=5500/4800 dosage=1.15  -> CN=1
```


#### Near-identical pair (per-site consensus)

Pooled dosage suffices for divergent paralogs but blurs on near-identical pairs: two paralogs `A` and `B` that differ at only a handful of coding sites, where gene conversion makes individual copies mosaics — so two converted `A` copies count like one `A` + one `B`. When the pairwise k-mer Jaccard flags a pair as near-identical, `call` splits it by per-site consensus. Take a haplotype whose module total (from coverage) is 2 copies:

```text
per-site A-allele vs B-allele k-mer fraction, one haplotype (module total = 2):
  site 1  A=1.0 B=0.0
  site 2  A=1.0 B=0.0
  site 3  A=0.5 B=0.5   <- a copy gene-converted at this site
  ...     (more sites at A=1.0 B=0.0)
```

1. Align the two coding sequences once (reference-side) to enumerate the divergent columns; at each column collect the k-mers carrying the A-allele and the B-allele.
2. Per haplotype, count each site's A-allele vs B-allele k-mers. Because every copy carries some allele at every site, `A + B` sums to the module total, so each site independently reports the split as the fraction `A/(A+B)`.
3. Take the median of that fraction across all sites: the converted site (`0.5`) is an outlier out-voted by the clean sites (`1.0`), so the median is `1.0`. Split the reliable module total by it → `A = round(2 × 1.0) = 2`, `B = 0`. A pooled count would have read the mosaic as one copy each.

```text
SVTYPE=DUP ...  gene A CN=2   gene B CN=0     (per-site recovers the mosaic; pooled dosage would miss)
```

Every `dup_gene_cn.tsv` row carries `hits`, `priv_kmers`, and `dosage`, so the split is auditable and the copy number carries through the whole spectrum — a deletion reads `0`, a single copy `1`, a duplication `2`. A gene with no CDS is sketched from its gene span instead of its coding sequence; a gene with no private k-mers at all (indistinguishable from a paralog) is reported as the module total with `reliable=0`.

