# Pre-registration: clump-normalized scoring, offline, before any production change

Supersedes an earlier draft of this file. That draft's S2/S3 were an **aggregate-Huber** construction
-- multiplicities summed within a clump, then a Huber loss on the aggregate residual. It is renamed
here so it cannot be confused with the corrected score, and its negative result is kept in
`PHASE1_AGGREGATE_HUBER_NEGATIVE.md`.

Written before any score below was run, and after the holdout was frozen
(`results/genotype_clump_oracle_holdout.tsv`).

## Why a clump-level score

Phase 0: support mismatch is 96.9% of the gross loss and 108% of the net, concentrated in 1-45 clumps
per donor with the top five carrying 53.6-99.2%. A global `rho` cannot address this -- multiplying
every candidate by one constant rescales confidence and chain influence but leaves the local ranking
untouched. Evidence has to be normalized or bounded **before** clumps are combined.

## The three scores. Exactly these.

Marker evidence for candidate pair `q`, relative to the background-only model:

```
e_j(q) = log NB( y_j | lambda * m_j(q) + mu )  -  log NB( y_j | mu )
```

```
S1(q) = the current production score                       (control, not a candidate)
S2(q) = sum over clumps of  mean_{j in clump} e_j(q)
S3(q) = sum over clumps of  clip( mean_{j in clump} e_j(q), -tau, +tau )
```

**Mean, not median, and that is fixed here.** Choosing between them after seeing results would be an
undeclared tuning dimension. The mean keeps every marker's direction while stopping 20 syncmers in one
clump from carrying 20 times the evidence of a clump holding one.

**The bound is symmetric** because Phase 0 found both failure directions: the certified pair is
penalised for lacking truth-supported sequence AND for carrying sequence the truth lacks.

Clump membership is `node_first_pos / fragment_len` -- panel-derived, independent of the candidate and
of the truth, so it is computable at genotyping time. Markers with no recorded position are their own
clump; the count of these is reported, not silently pooled.

**Not tested:** any dosage term, `mass_bp`, `mass_weight`, `mass_window`, `compositional`,
`marker_outlier`, `robust_c`, `--all-kmers`. Phase 0 says the immediate error is support geometry, and
every one of those is in the negative-results ledger.

## tau

Chosen on the 30-donor **selection** set, from a grid written down now and not extended afterwards:

```
tau in { 0.5, 1, 2, 5, 10, 20, inf }      (inf = S2)
```

Objective, fixed: minimum total certified sequence excess over the selection set, with ties broken
toward the larger tau (less intervention).

No "one-fragment" formula is used. Inventing a theoretically-labelled constant whose scale is not
justified after clump averaging would be tuning with extra steps.

## Reported for every donor and score

Excess alone can be gamed by declining to decide, so the tie structure is reported alongside it:

- selected pair; certified-pair rank and tie count;
- exact number of tied top pairs, and whether the certified pair is among them;
- best and worst sequence excess **within** the tied top set;
- total certified sequence excess;
- donor-level wins / losses / unchanged against S1;
- number of contributing clumps and their concentration.

First pass: noiseless truth multiplicities, all candidates, no HMM. That isolates the score.

## Gate, numerically

Written before any selection-set output was inspected. "Materially" and "mainly" are removed: each
condition below is a number, so a result cannot be argued past the gate after the fact.

Let `X(score)` be total certified sequence excess over the set being scored.

1. **Control.** Every on-panel control ranks its own source pair 1, uniquely. 30/30 on selection.
2. **No exact donor breaks.** No donor whose excess is 0 under S1 becomes nonzero. This is checked on
   the five zero-excess DEVELOPMENT donors as well as on the selection set, and S2 and every S3 arm
   must be run on those five -- S1 matching the binary in 15/15 validates the implementation and says
   nothing about this condition.
3. **Size of the win.** `X(candidate) <= 0.90 * X(S1)` -- at least a 10% reduction.
4. **Not one donor.** Removing the single most-improved donor, the improvement is still positive.
5. **Not bought with ties.** The worst-case excess among the tied top pairs is reported, and
   condition 3 must still hold when each donor is scored at that worst case rather than at its best.
6. **Direction.** Donor-level wins strictly exceed donor-level losses.
7. **Final holdout.** Conditions 1-6 hold on the 50-donor set, opened once, after the score and tau
   are frozen.
8. **No collateral damage.** Ordinary non-array blocks unchanged or better.

## Declared outcomes

- Normalization works and bounding adds nothing -> ship normalization alone.
- Bounding is necessary -> report sensitivity across the declared tau grid.
- Neither works under noiseless evidence -> the retained feature representation must change, and no
  further likelihood parameter is justified.

## Order, fixed

1. Phase 0 write-up and terminology committed.
2. Holdout manifest and this file committed. **Done before any run.**
3. S1/S2/S3 implementation validated on the 15 development donors.
4. 30-donor selection set run; one score and one tau frozen.
5. 50-donor final holdout run once.
6. Only on a pass: integrate, then observed reads, GQ recalibration, equivalence-set output, HMM.

Production inference is untouched until step 6.
