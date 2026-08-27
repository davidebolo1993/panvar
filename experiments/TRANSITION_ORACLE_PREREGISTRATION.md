# Pre-registration: does boundary-transition evidence resolve what block-internal evidence cannot?

Written 2026-08-25, **before the oracle was run**. This is an identifiability question, answered with
no emission and no reads, so it can only be wrong about the representation — never about the model.

## Why this and not the alternatives

Three interventions have been closed at the marker-poor blocks, each on measurement:

| intervention | result |
|---|---|
| denser markers (`--all-kmers`) | falsified; same filters, same place |
| keeping ambiguous markers (`--no-region-unique`) | falsified; breaks leave-zero-out |
| smaller k | rejected; creates candidates, all filtered away |
| 2-syncmer adjacency, retained set | adds nothing at cyp2d6 block 5 (10/105 classes, unchanged) |

All four share a defect: the evidence is owned by ONE block, so a marker touching two blocks is
discarded by the confinement rule. Boundary-transition evidence changes the ownership, not just the
sequence: a k-mer spanning a junction belongs to the TRANSITION (left allele -> current allele), which
is a state the chain already carries, rather than to either block's emission.

This is why it does not inherit the earlier falsifications, and why gstm1's lack of an invariant
neighbouring block is not fatal — the variable neighbour is part of the transition state itself.

## Construction, frozen

For block `B` with left neighbour `L`:

1. For every ordered pair (allele `a` of `L`, allele `b` of `B`) that **co-occurs on at least one
   panel haplotype** (from `--dump-haplotype-alleles`; pairs no haplotype carries are not panel
   evidence and are excluded), form `suffix(seq(a), k-1) + seq(b) + prefix(seq(next), k-1)`.
2. Keep only k-mers **crossing a boundary** — those with at least one base on each side. Purely
   internal k-mers are the evidence already shown not to work and must not be counted here.
3. A haplotype's transition vector is the multiset of junction k-mers its own (a, b) pair yields.
4. A diploid transition pair's predicted vector is the SUM of its two haplotype-transition vectors.

## Metrics

Identical instrument to section 12, so results are comparable:

- distinct pair-vectors out of all diploid pairs;
- pairs in non-singleton (ambiguous) classes;
- largest equivalence class;
- whether the truth pair's class contains the pair production actually called.

Reported for: block-internal retained markers (the production baseline), junction markers alone, and
the two combined.

## Declared prediction

Junction markers separate pairs that block-internal retained markers cannot, at gstm1 blocks 9 and 11
(alleles under k, where internal evidence is empty or fully filtered) and at cyp2d6 block 5 (alleles
long enough, where internal evidence exists but is filtered to 10 of 36).

## Gate — with the floor the last pre-registration lacked

The k-size gate passed while **1 of 30 observations changed**. A mean that improves on one moved
observation is not evidence. Therefore:

1. **Minimum informative count: at least 20 (sample, block) observations must CHANGE**, across at
   least two blocks and at least two loci. Below that the arm is reported as inconclusive, whatever
   the means do.
2. Ambiguity must fall at the target blocks: fewer pairs in non-singleton classes, and a smaller
   largest class, at gstm1 9/11 and cyp2d6 5.
3. Control blocks (marker-rich, e.g. cyp2d6 block 7 at 1,253 retained markers) must not get worse.
4. The truth pair's class must shrink, and where production called a wrong pair inside the truth's
   class, that pair must leave it.

Only if 1-4 hold does this become a pairwise HMM transition emission. An oracle win is a licence to
implement, not a result about accuracy.

## Declared failure modes

- **Junction markers fail the same filters.** They contain flank sequence shared with neighbours, so
  under the CURRENT rules they would be dropped as multi-block. The oracle deliberately bypasses the
  filters to ask whether the INFORMATION is there; if it is, the filter rule must change alongside,
  and that change is then part of the implementation rather than a surprise.
- **Combinatorial growth.** Transition states are pairs of alleles, so a block with 21 alleles beside
  one with 48 has up to 1,008 states. Restricting to pairs observed on a panel haplotype is what
  keeps this finite, and the observed count is reported as part of the result.
- **Short alleles may still yield nothing.** gstm1 block 11 has an empty bypass allele; its junction
  k-mer is the direct left->right join. That case is in the construction by design, and if it is the
  only evidence recovered, that is a small result and will be reported as one.

---

# RESULT (2026-08-25, corrected after review)

## Correction first: an allele-omission bug voided the first run's denominators

The oracle built its allele list from the `--dump-block` FASTA. `--dump-block` skipped alleles with
an empty sequence, so the **empty bypass allele — the one a deletion takes, and the case the
construction was written to reach — was absent from every count.** Blocks were scored as 20/16/13
alleles instead of 21/16/14.

The symptom was visible in the first run (13 alleles against section 12's 14) and was recorded as a
caveat instead of chased. A count that does not reconcile is a bug until shown otherwise.

Fixed at the source: `--dump-block` now emits empty alleles, and the oracle derives its allele
universe from the haplotype-allele matrix and asserts it against the dump.

## The first run also measured the wrong factorisation

It combined the left and right junctions of a block into one vector and attached the result to the
centre allele as if it were a property of that allele. That is a **three-block factor**, not the
first-order pairwise transition being proposed, and 12 of 21 alleles at gstm1 block 11 occur beside
more than one neighbour, so "the allele's context" is a choice rather than a fact.

Rebuilt as the actual proposal: states are (left, centre, right); each boundary contributes
`suffix(left) + prefix(right)` k-mers separately; a diploid pair sums two state vectors.

## What is actually being called is the CENTRE, not the transition label

Several transitions share a centre allele, so a vector class that cannot name the transition may
still name the allele. Reporting label ambiguity alone would have read as a negative result and been
wrong. Both are reported; **CENTRE ok** is the quantity a caller needs.

| block | centre alleles | state set | evidence | states | pairs | CENTRE ok | label ok |
|---|---:|---|---|---:|---:|---:|---:|
| gstm1 11 | 21 (0-21 bp) | panel triples | boundary only | 94 | 4,465 | **100.0%** | 2.7% |
| | | panel triples | boundary + internal | 94 | 4,465 | **100.0%** | 2.7% |
| | | sampled all | boundary + internal | 400 | 80,200 | 99.7%* | 6.3% |
| gstm1 9 | 16 (14-42 bp) | panel triples | boundary only | 65 | 2,145 | 76.7% | 22.6% |
| | | panel triples | boundary + internal | 65 | 2,145 | **88.2%** | 24.1% |
| | | sampled all | boundary + internal | 400 | 80,200 | 96.6%* | 42.5% |
| cyp2d6 5 | 14 (261-851 bp) | panel triples | boundary only | 17 | 153 | **9.8%** | 3.9% |
| | | panel triples | boundary + internal | 17 | 153 | 98.0% | 76.5% |

\* Sampling STATES removes collision opportunities, so these are **upper bounds**, not estimates.
They are here only to show that allowing novel recombinations does not collapse the result — which is
required, since restricting to panel-observed transitions permanently would make mosaics
unrepresentable by construction. Note block 9's sampled figure EXCEEDS its panel-triple figure purely
from this bias; the panel-triple number is the trustworthy one.

Block-internal evidence alone, for scale: gstm1 11 gives **1 of 231** distinct pair-vectors (nothing),
gstm1 9 gives 21 of 136, cyp2d6 5 gives **103 of 105**.

## Verdict, by block

- **gstm1 block 11 — build it.** From no identifiability whatever (1/231) to the centre genotype
  determined in 100% of 4,465 diploid pairs, on boundary evidence alone. Internal evidence adds
  nothing because there is none.
- **gstm1 block 9 — build it, but boundary evidence is not sufficient alone.** 76.7% on boundaries,
  88.2% adding internal. The two evidence classes are complementary, and 11.8% of pairs remain
  ambiguous with both.
- **cyp2d6 block 5 — do NOT add junction evidence.** Boundary-only pins the centre in 9.8% of pairs
  against internal k-mers' 103/105. With alleles of 261-851 bp the boundaries are a negligible
  fraction of the sequence. Its 98.0% combined figure is internal evidence carrying the result.
  cyp2d6 needs the FILTER changed, not new evidence.

## The two-locus gate CANNOT be satisfied on this panel set — stated, not worked around

The pre-registration required changed observations across at least two loci. Auditing the other four
loci for the same failure mode:

| locus | blocks with >= 8 alleles | fewest informative markers |
|---|---:|---|
| gstm1 | — | **0** (block 11), **7** (block 9) |
| lpa | 21 | 15 |
| ankrd36c | 18 | 18 |
| acot | 15 | 22 |
| c4 | 10 | 47 |

**No second locus in this cohort exhibits marker starvation.** The positive evidence rests on one
locus, and the gate is not met. Two honest readings, both worth holding:

- against: a mechanism seen at one locus may be a property of that locus's block partition rather
  than of short alleles generally, and both gstm1 blocks are BACKBONE blocks carrying 16 and 21
  alleles, which is itself irregular and may be the real defect;
- for: the target is extremely concentrated. gstm1 blocks 9 and 11 are **36 of 2,296 block
  observations (1.6%) and carry 11.1% of the entire cohort's gap mass** — 34% of gstm1's.

Recommendation: implement, scoped to blocks whose alleles are shorter than k, and treat the
single-locus evidence as the main risk. Adding a locus with short-allele blocks to the panel set
would retire that risk and is cheaper than the implementation.

## Binding constraint for the implementation

**Hundreds of overlapping junction k-mers are not hundreds of independent observations.** They come
from the same reads. gstm1 block 11 shows 155 markers over a 0-21 bp allele — every one of them
inside a single fragment. The existing `marker_clumps` ESS discount exists for exactly this and must
be applied to junction evidence before any accuracy claim; without it the factor will be
overconfident by roughly the clump size. Nothing measured here has read support, depth at a junction,
or sequencing-error sensitivity in it.
