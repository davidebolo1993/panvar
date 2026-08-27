# Finishing `genotype`: point-by-point to a clean release

Written 2026-08-27, replacing the genotype sections of the 2026-08-21 plan. Everything below is
either DONE, RUNNING, READY (unblocked, not started) or BLOCKED (waiting on something named).

The central decision this plan turns on: **the module does not have to solve KIV-2 to ship.** It has
to be honest about what it cannot resolve. Sections B and C are shippable whatever Phase 1 returns.

---

## A. Close the scientific line

### A1. Phase 1 selection run -- RUNNING

30 untouched donors, frozen in `results/genotype_clump_oracle_holdout.tsv`. Scores S1 (control), S2
(clump-normalized), S3 (symmetrically bounded) over the declared tau grid. Definitions and gates in
`PHASE1_ORACLE_PREREGISTRATION.md`; nothing there may change now.

### A2. Five zero-excess controls -- QUEUED behind A1

Gate condition 2. S1 matching the binary in 15/15 validates the implementation and says nothing about
whether S2 or any S3 arm keeps an exact donor exact.

### A3. Freeze one score and one tau -- BLOCKED on A1, A2

Selected on the 30-donor set by the pre-registered objective. Written down before A4 is run.

### A4. Final holdout, opened once -- BLOCKED on A3

50 donors. Opened exactly once, after the score and tau are frozen. A failure here ends the line; it
does not license a re-tune.

### A5. Branch on the outcome

- **Pass** -> A6.
- **Fail** -> the retained feature representation cannot support sequence-nearest off-panel
  projection. Stop tuning the likelihood. Section B becomes the whole deliverable, and that is a
  legitimate result, not a defeat: it converts a confidently wrong caller into a correctly
  uncertain one.

### A6. Production prototype -- BLOCKED on A5 passing

Integrate as an emission behind a flag, default off. Then, in order: observed reads end to end,
interaction with the chain, GQ recalibration, and only then a default change.

### A7. Outstanding diagnostic debt

- Map dominant clumps to sequence spans. Partly done: median 26 markers over 268 bp, 19.6% of markers
  recurring elsewhere in the array, and **15 clumps implicated in more than one donor carrying 44.6%
  of the dominant-class loss**. What is missing is the sequence-level count -- how many distinct
  differences a clump holds -- which needs the allele sequences and is READY.
- Whether any of this transfers to cyp2d6 and gstm1, whose errors are a different bucket. READY.

---

## B. The output contract. Shippable regardless of A.

This is the part that makes the module honest, and none of it depends on Phase 1.

### B1. Equivalence sets -- READY, highest value

Measured: 44% of KIV-2 loss is dosage-split -- the diploid total is right and the split between
homologues is wrong. The emission's mean is `lambda * (m_a + m_b) + mu`, a function of the SUM alone,
so two allocations with the same summed vector are **exactly tied by construction**. Where short
reads carry no information separating them, no parameter can recover it.

Emit, per block: the best pair, the set of pairs within a calibrated score margin, posterior mass on
that set, and an explicit unresolved-allocation flag. A correctly reported six-pair class beats a
confident wrong one.

### B2. Separate the three answers -- READY

Total dosage, composition given dosage, and allocation between homologues are different questions
with different evidence. Report them separately with their own confidences rather than collapsing
them into one argmax and one GQ.

Caveat that must be respected: `mass_bp` is NOT the dosage answer. Its pre-registered replacement
rule failed -- B = 0.00 at six of seven array blocks, mean absolute length error 7,981 bp against
`called_bp`'s 412. Any dosage estimator needs its own validation before it is promoted.

### B3. GQ calibration -- READY

At GQ 40-60 the observed error rate is about 1,000x what the number claims, and the curve saturates
near 0.5% beyond GQ 60. Fit a monotone map on held-out donors and emit `GQ_CAL` alongside `GQ`; do
not overwrite the raw value. Gate: observed error within 2x of claimed in every calibrated decile.

### B4. Reliability fields -- READY

Surface `markers/allele` (below ~2 the mean gap is 0.039-0.045; above it, <=0.0012), the truth-class
size where computable, and a reason string so a declined call says why.

---

## C. Release blockers, none of them scientific

### C1. `docs/modules/genotype.md` and `docs/algorithms/genotype.md` -- READY

Owed since the first plan and still absent. Everything else in `docs/` has both pages. This is the
single most visible gap in the tree.

### C2. Non-numeric node name silently becomes 0 -- READY

`genotype_blocks.cpp:115`. A malformed name is read as node 0 rather than refused.

### C3. Numeric range validation -- READY

One validated option in the whole command.

### C4. Quantized depth fallback -- READY

The pooled median of integer anchor counts cannot express a non-half-integer whatever the sample size.

### C5. CTest does not register `genotype_stats` when the module is off -- READY

A green CTest can currently mean the genotype tests never ran. The suite should say so.

### C6. `n_scored_alleles` reads -1 at non-representable blocks -- READY

Same nesting mistake as `local_best_*`, already fixed there.

---

## D. What "clean release" means

1. A and B resolved, with B shipped whether or not A passes.
2. Every C item closed.
3. No default changed without a held-out result: improvement on donors the method never saw, and no
   improvement bought by declining more calls.
4. The negative-results ledger carried forward intact. Twelve interventions failed before this line;
   the record of why is worth more than any one of them.
5. `--certified-oracle`, `--noiseless-counts`, `--probe-pair` and the ledger dumps stay in the shipped
   binary. They are how the next person reproduces any of this.

## E. Explicitly NOT in scope

- Graph-flow synthesis for unrepresentable blocks.
- Long-read or linked-read allocation evidence. Measured: at 350 bp the order evidence favours the
  wrong allele, and the crossover is 2-5 kb against a ~5.5 kb repeat unit. Allocation at this locus
  is very likely not identifiable from these reads, which is why B1 is the answer and not a
  consolation prize.
- Any return of `--all-kmers`, `marker_outlier`, `mass_window`, `compositional` or read-derived
  routing except as a new, pre-registered, clump-level model.
