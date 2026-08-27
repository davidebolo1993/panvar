# Plan: the blocks that carry the error

Written 2026-08-25 from the six-locus cohort (`results/real_data/<locus>/genotype/`).

## The observation that makes a plan possible

Leave-zero-out is **almost, but not entirely, exact on simulated reads**: 7 of 2,296 block
observations (0.30%) carry non-zero error, and all 7 also called the wrong allele pair. So the model
nearly always finds the answer when the answer is in the panel — but not always, and the exceptions
matter more than their weight suggests (see "The seven clean failures" below). The remaining error is
off-panel generalisation, and it is not spread out:

| | |
|---|---|
| blocks carrying measurable leave-one-out error | **22 of 115** |
| top 5 blocks | **42%** of all error |
| top 10 | **67%** |
| top 15 | **80%** |
| blocks contributing under 1% each | **93 of 115** |

So this is not "the genotyper is 45% accurate". It is "twelve blocks are broken and the other
hundred are fine". gstm1 is the clearest case: blocks 0-7 and 13-19 are flat at zero in both
regimes and both read sources, and blocks 8-12 carry the locus.

## The seven clean failures — the most diagnostic observations in the cohort

| locus | block | sample | error |
|---|---|---|---:|
| cyp2d6 | 3 (bubble) | NA19238 | 0.0745 |
| acot | 5 (bubble) | HG00350 | 0.0227 |
| cyp2d6 | 4 (backbone) | NA19238 | 0.0088 |
| lpa | 8 (backbone) | HG00253 | 0.0002 |
| acot | 6 (backbone) | HG00128 | 0.00005 |
| acot | 4 (backbone) | HG00350 | 0.00001 |
| cyp2d6 | 10 (backbone) | NA19238 | 0.000007 |

These are leave-zero-out on simulated reads: **the truth pair is in the panel and the reads are
noiseless copies of it.** Read acquisition, mapping, depth, sequencing error and panel coverage are
all excluded by construction. Whatever makes the model choose wrongly here is the emission or the
chain, isolated.

Three of the seven are one donor (NA19238) and two are another (HG00350), which suggests a per-sample
property rather than a per-block one — a hypothesis, not a finding, since n=7.

**These are the cheapest failures in the whole cohort to debug and should be looked at first**, before
the top-15 classification: they are few, fully controlled, and reproducible in seconds.

**The panel floor at these blocks is ~0.** At cyp2d6 block 5 it is 0.0001 against 0.138 achieved; at
lpa 15 and lpa 9 it is exactly 0.0000. The panel contains a near-perfect pair and we do not pick it.
These are model failures, not representation failures.

## Three candidate mechanisms — and the classifier is NOT marker count

| | mechanism | what fixes it |
|---|---|---|
| **A** | filter discards markers that would have been sufficient | change the filter |
| **B** | alleles shorter than k, so no candidate exists at all | boundary-transition evidence |
| **C** | retained markers ARE sufficient, and the call is still wrong | change the emission |

The tempting classifier is "what fraction of informative markers survived". It is wrong. cyp2d6
block 5 retains 10 of 36 — a 28% survival rate that looks healthy — but those 10 markers separate
**10 of 105** diploid pairs while all 36 separate **103 of 105**. It is squarely category A.

The correct classifier is **identifiability of the retained set**: do the markers the model actually
scores distinguish the truth pair from the others? That is measured, not inferred, and the
instrument already exists (`--ledger-block` plus the pair-vector class analysis).

Measured so far, on three blocks only:

| block | alleles | informative | retained | retained separates | all separate | class |
|---|---:|---:|---:|---|---|---|
| cyp2d6 5 | 14 | 36 | 10 | **10/105** | 103/105 | **A** |
| gstm1 11 | 21 | 0 | 0 | — (no evidence) | — | **B** |
| gstm1 9 | 16 | 7 | 0 | — (all filtered) | 21/136 | **A + B** |

Everything else in the top 12 is **unclassified**. Assigning them by marker count would repeat the
mistake this project has already made three times: reading a proxy as the quantity.

## Phase 0 — classify the top 15 blocks. Do this before anything else.

Run the ledger and the pair-vector class analysis on each of the 15 blocks carrying 80% of the
error, and record, per block:

1. informative candidates, retained, and the fate split (multi-block / over-expected / both);
2. diploid pair-vector classes under **retained** markers and under **all informative** markers;
3. whether the truth pair sits in a class containing the pair production actually called;
4. allele length distribution against k.

Then A / B / C follows from the measurement:

- retained separates poorly **and** all-informative separates well → **A**;
- no informative candidates, alleles under k → **B**;
- retained separates well and the call is still wrong → **C**.

Cost: minutes per block, no reads needed. **Nothing downstream should start until this is done**,
because A, B and C need three different fixes and the current split of the top 12 is a guess.

## Phase 1 — one workstream per class, gated

### A. The filter discards sufficient information

Both region rules drop markers *hard*: `blocks_with > 1` (confinement) or `actual > expected`
(over-expected). At the poor blocks nearly every dropped marker fails **both**, so relaxing either
alone rescues nothing — measured at cyp2d6 block 5, where relaxing confinement alone rescues
**zero** markers.

Globally keeping them is already falsified (`--no-region-unique` breaks leave-zero-out, 13/13 samples
worse). The reason is sound: a marker occurring elsewhere accumulates counts from everywhere, so its
observed count stops meaning what the model assumes.

**The fix is to stop choosing between keep and drop.** Model the contamination instead: a marker
appearing in blocks *b1..bn* has expected count `lambda * sum_i m_i(g_i)` over the blocks that carry
it, not `lambda * m(g)`. Retain it with a *background term* for the other blocks' contribution.
Cheapest honest version: treat the other blocks' contribution as a nuisance offset estimated from
their own anchors, so the marker still discriminates within this block.

- Gate: cyp2d6 blocks 3/5 and lpa 9/15 improve; **leave-zero-out stays exact** (this is what killed
  the last attempt); marker-rich blocks do not regress; direct and indexed agree.
- Risk: this is a real modelling change, not a flag. Prototype in the oracle first — does the truth
  pair become separable under retained-plus-background? — before touching the emission.

### B. Alleles shorter than k

Only gstm1 9/10/11 today, ~14% of the top-12 error mass. The boundary-transition oracle already
answered the identifiability question: centre genotype determined in **100%** of 4,465 diploid pairs
at block 11 (against 1/231 for internal k-mers) and 88% at block 9.

Blocked on three things, all recorded in `TRANSITION_ORACLE_PREREGISTRATION.md`:

1. **one locus only** — no other locus in the cohort has marker starvation, so the pre-registered
   two-locus gate cannot be met. Adding a locus with short-allele blocks is cheaper than the
   implementation and retires the main risk;
2. **effective sample size** — 155 junction markers over a 0-21 bp allele all sit in one fragment.
   They must go through `marker_clumps` or the factor will be overconfident by roughly the clump
   size;
3. **it is a transition factor, not an emission** — at block 9 the residual signal is the
   neighbour's identity, so folding it into one block's emission would consume as evidence a
   quantity the model is inferring.

### C. Markers sufficient, call still wrong

If Phase 0 puts blocks here, the evidence is fine and the *likelihood* is wrong — which is the same
shape as the long-standing KIV-2 result, where the certified pair loses under **noiseless** evidence.
That rules out read sampling, depth, and pruning, and leaves the emission's geometry: the presence
veto, where a handful of absent-marker mismatches outvote dozens of multiplicity agreements.

This is the least-understood class and probably the largest. Do not propose a fix before Phase 0
says which blocks are in it; the existing negative-results ledger already contains a dozen emission
tweaks that failed.

## What NOT to do

Retired on measurement, with the evidence in `GENOTYPE_STATUS_2026-08-25.md`:

| | |
|---|---|
| denser markers (`--all-kmers`) | falsified; same filters, same place |
| keep everything (`--no-region-unique`) | falsified; breaks leave-zero-out |
| smaller k | rejected; creates candidates, all filtered away |
| adjacency on the retained set | adds nothing at cyp2d6 5 (10/105 unchanged) |
| junction evidence at long-allele blocks | 9.8% centre-identifiable against internal 103/105 |
| mapped-coverage pipeline | `mass_bp` beats it and needs no aligner |

## Sequencing

0. **The seven clean leave-zero-out failures** — hours, fully controlled, and they isolate the
   emission from everything else.
1. **Phase 0 classification** of the top 15 blocks — days, no new code, decides everything else.
2. **A** in the oracle (does retained-plus-background separate the truth pair?) — no emission change.
3. **B** implementation, once a second short-allele locus exists.
4. **C** only after Phase 0 says which blocks are in it.

Read acquisition is a separate track and is not in this plan: real reads cost 9.5 points of
leave-zero-out exactness (99.7% simulated against 90.2% real) with the answer sitting in the panel,
and at gstm1 they cost error up to 0.456 at a block that is exact on simulated reads. That is a
larger single lever than anything above and it is input plumbing, not modelling.
