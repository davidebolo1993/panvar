# M11-E result: all five catastrophic blocks are one-unit copy-number errors that the objective prefers in the winner's own context; no search miss, no numerical inconsistency

**Run.** Branch `m11-e-catastrophic-attribution`; freeze `fc6d5ab`; result after.
- **Results:** `tests/results/m11e_catastrophic_attribution_20261021/`.
- **Budget:** 196.4 of 900 s: 120 s of disclosed development plus 76.4 s of ledger wall.
- **Units:** 15 of 15 passed first time, with no amendment and no watchdog or envelope event; peak 1.57 GiB.
- **Probe:** `genotype_probe` sha256 `d889c7d5…`, the same as M11-D's.

**Not done:**
- no caller rerun, new reads or new inference;
- no fit, cutoff, metric, floor search, validation, clustering, prior change, promotion, merge, push or cleanup.

The 18 ACOT validation families are untouched. M11-D's frozen results are unchanged; one correction is **appended** to its
RESULT.

## 1. Reporting and eligibility corrections

**Artefact totals.** `report.a4` summed each interval's artefact bp over both of its flagged blocks.
- **Corrected,** once per interval: **CYP2D6 1,770 bp** (5 intervals) and **ACOT 16 bp** (5 intervals). C4 has no
  degeneracy interval.
- **The re-count reproduces report.a4's 3,540 / 32 exactly,** which confirms the cause.
- **Record:** `runs/fixes/artefact_correction.tsv` and `artefact_intervals.tsv`. The original report is preserved, and an
  addendum is appended to M11-D's RESULT.

**Calibration-eligibility annotations** (`runs/fixes/eligibility.tsv`). Labels are copied unchanged, nothing is relabelled,
and nothing is refitted.

| eligibility | blocks | non-top | catastrophic | fit-informative |
|---|---|---|---|---|
| BLOCK_OUTCOME | 395 | 14 | 5 | yes (if assessable) |
| ALL_EMPTY_COUNTED | 34 | 0 | 0 | no; counted in coverage |
| REPRESENTATION_DEPENDENT_INTERVAL_EXACT | 11 | 7 | 0 | interval level only. The concatenated predicted interval equals truth exactly (CYP2D6 HG02769, NA19650, NA19983; ACOT HG02083). The block labels are **not** relabelled correct |
| REPRESENTATION_DEPENDENT_UNRESOLVED | 14 | 6 | 0 | no. The interval is not exact (CYP2D6 HG03456: 59 edits; NA19705: 60; ACOT HG01975, HG02178, HG02074, NA18971: 1–6), and its floor was never computed |
| NO_PREDICTION / UNDEFINED_TRUTH | 19 / 1 | | | no |

**QVfrac ≥ 0.8 can accept a 12-kb-excess prediction.** The case is CYP2D6 NA19983 sv:9: excess 12,075 bp, QVfrac 0.864.

| homologue | edits / length | QV predicted | QV best-available |
|---|---|---|---|
| truth1 | 4 / 39,639 | 39.96 | 39.96 (at its best) |
| truth2 | 13,696 / 27,464 | 3.02 | 9.76 |

- **Why it passes:** QVs are summed (42.98 / 49.72), so the accurate homologue supplies 93% of the numerator.
- **Why the denominator is small:** the truth2 best-available allele is itself 1,621 edits from truth, which leaves little
  for the failure to lose.
- **The same pattern** appears in ACOT NA19909 (0.749) and C4 NA19434 (0.683): one homologue at its best, the other
  27 kb wrong.
- **Status:** QVfrac stays descriptive. Excess and normalized excess stay primary.

**Runner.**
- **Fix:** a unit passes only when it exits 0 **and** the watchdog verdict is `ok`. A dependent runs only if its
  dependencies passed, and `next` returns non-zero on a failure or NOT RUN.
- **`gate.runner`: 8 of 8 checks pass:**
  - a failing analysis under verdict `ok` stops its dependent, and so does a timeout;
  - an independent unit still runs;
  - the sealed M11-D ledger shows the old precondition admitted `fit` after `outcomes` exited 1.
- **This batch:** run with the fixed runner; the loop stopped only at "all recorded".

## 2. Attribution (frozen method; `runs/trace/attribution.tsv`)

**Method:** in each reported winner's context, only the affected block is replaced by the evaluator's best-available pair,
in both homologue assignments, and scored by `score-set` under the run's exact objective and settings.

**Checks — all pass in all 8 cases:**
- **reproduction:** the incumbent total is reproduced within 1e-6;
- **sequence:** the incumbent's spelled bytes equal the run's recorded haplotype SHA-256, flanks are byte-identical, the
  block is exactly the declared allele, and the probe's scored bytes match (FNV-64 and length);
- **origin re-summation:** equals the scorer's fragment term in every decomposition;
- **run record:** where the run itself scored a replacement (10 comparisons, including the controls), its recorded total
  equals score-set's within 4.7e-8.

| case | block | call vs truth length (bp) | Δfragment | Δexposure | Δprior | **Δtotal** (better assignment) | best pair: local rank (list = 32) / scored by the run? | dosage audit | classification |
|---|---|---|---|---|---|---|---|---|---|
| K1 C4 NA24385 | sv:7 | call **+26,370** (one extra unit) | −385.3 | +329.6 | −19.3 | **−75.0** | 1 / yes (repair round 3) | CONTENT_REPAIRED: 5u → **4u** (count features); the repaired content scores −26.0 by likelihood | OBJECTIVE_PREFERS_FARTHER_PREDICTION |
| K2 ACOT NA19909 | sv:13 | call **−27,889** | −5,326.7 | −343.5 | −0.5 | **−5,670.6** | 8,808 / no | **FAILED** ("KMAX > 3 is not enumerated") | OBJECTIVE_PREFERS_FARTHER_PREDICTION |
| K3 C4 NA19434 | sv:7 | call **−26,371** (one unit short) | −576.3 | −329.6 | −9.5 | **−915.3** | 1,521 / no | CONTENT_REPAIRED: 5u → **6u**; the repaired content scores **+221.3** by likelihood but is not panel-spellable (prior −inf) | OBJECTIVE_PREFERS_FARTHER_PREDICTION |
| K4 CYP2D6 NA19650 | sv:9 | call **−12,145** | **+127.1** | **−151.9** | −6.7 | **−31.5** | 24 / yes (repair round 5) | **FAILED** (unit library: flanks disagree in sv:9) | OBJECTIVE_PREFERS_FARTHER_PREDICTION |
| K5 CYP2D6 NA19983 | sv:9 | call **+13,699** | −284.5 | +151.7 | −16.5 | **−149.3** | 73 / no | **FAILED** (same) | OBJECTIVE_PREFERS_FARTHER_PREDICTION |

**Every catastrophic error is a one-unit copy-number error** in a `repeated_unit` block: about 26 kb (C4), 28 kb (ACOT)
and 13.7 kb (CYP2D6).

**No search miss.** No best-available replacement scores higher in the winner's context.
- K1 and K4: the run scored the replacement in repair and rejected it correctly by its own objective.
- K2, K3 and K5: the pair was outside the 32-pair local list and never scored in that context, but it loses anyway.

**No numerical or implementation inconsistency** was found.

**What decides each case:**
- **K1, K4, K5 — the unit count is set by fragment placement against exposure.** The prior is small and decisive in no case.
  - **K1 and K5:** the reads favour the haplotype with the extra unit. K1: 1,138 fragments against 212, almost all
    interior to the block. Exposure only partly offsets this, so the net likelihood is −56 / −133 for the truth-length
    pair.
  - **K4, the reverse:** the reads favour the truth-length pair (+127; 528 fragments against 28), but exposure for the
    extra 12 kb (−152) outweighs them.
  - **So the objective's implicit read-depth copy-number balance errs in both directions,** by margins of 25–133 nats.
- **K2 and K3 — the panel's best-available allele is read-inconsistent.**
  - **K2:** 16 fragments place in the incumbent's block and nowhere in the best pair's haplotypes, costing −4,280. Allele
    320 is the nearest by edit distance (619 edits) but lacks read-supported sequence.
  - **K3:** the best pair loses on fragment and exposure. The dosage audit's off-panel content, with one more unit and
    in the truth's direction, **gains +221 nats of likelihood**. It cannot be reported because it has no panel prior.
  - **These are availability and representation limits** on top of the copy-number error. No panel allele has both the
    right unit count and read-consistent content.

**The dosage audit** is evidence, not proof.
- **Where it ran (K1, K3),** its count-based correction moves the same way the truth does.
- **K1 is a disagreement:** the count features say one unit fewer, while the likelihood prefers the extra unit (the
  repaired content scores −26).
- **It FAILED in 3 of 5 catastrophic cases,** and in the CYP2D6 control, so the independent count evidence was missing
  exactly where it was needed.

**Controls.** Each control already calls its best pair, so one assignment is identical (Δ = 0 in every term).
- **The other assignment is the homologue swap at the block:** fragment Δ = 0, and exposure Δ ≤ 5e-12.
- **The swap loses only through the prior:** −30.2, −21.1 and −40.2. Phase across these long blocks is invisible to the
  reads.
- **About their label:** the frozen rule labels them `CONTROL:OBJECTIVE_PREFERS_FARTHER_PREDICTION` mechanically, but the
  swap is not a farther prediction. The controls serve to show that reproduction, byte identity, record consistency and
  decomposition reconciliation all hold.

**Limits:**
- A losing single-block replacement says nothing about the global optimum.
- It is not by itself a likelihood bug.
- Different neighbouring blocks were not explored, and joint multi-block alternatives were not explored.

## 3. Recommendation: the single correction to investigate next (no implementation)

**Investigate the copy-number model for `repeated_unit` blocks.** Two terms set the unit count:
- the fragment term's multi-placement mass, which grows with each added unit;
- the exposure term.

Check whether their balance is calibrated, against an independent count-based dosage estimate.

**Why this, and not search, numerics or the prior:**
- **All five catastrophic errors are one-unit dosage errors.**
- **In K1, K4 and K5 the decomposition isolates the fragment-versus-exposure balance as the deciding term,** erring in
  both directions by only 25–133 nats.
- **In K3 the likelihood favours the truth-direction unit count once read-consistent content exists.**
- **The independent count evidence points the truth's way wherever it ran (2 of 2),** but the audit cannot currently run
  on CYP2D6 sv:9 (unit-library flank disagreement) or ACOT sv:13 (KMAX > 3).

**Evidence against the alternatives:**
- **Search:** no case shows a better-scoring replacement.
- **Numerics:** every reconciliation and record check holds.
- **Prior:** never decisive. In the better assignment |Δprior| ≤ 19.3 (≤ 35.4 in either assignment), while the
  likelihood alone (Δfragment + Δexposure) already favours the incumbent by ≥ 24.8 in every case.

**Starting point:** M6-A's 31-mer count fit, which corrected CYP2D6 LOO dosage but not the C4 3+1 case. It is a natural
reference.

**The first investigation step:** in K1, K4 and K5, compare the per-unit depth the likelihood implies against
`marker_depth_rate_per_copy` and the audit's robust features. This should happen before any change to the objective.

**Second, separate thread** (not recommended first): **K2 and K3 are availability limits.** K3's off-panel, likelihood-better
content is blocked by the panel prior.

## Files

- `runs/fixes/`: corrections, eligibility, the QVfrac breakdown.
- `runs/prepare/`: candidates, sequence checks, spellings.
- `runs/score.*/`: score-set tables, local ranks, decompositions.
- `runs/trace/attribution.tsv`.
- `MANIFEST.tsv`, `ledger.tsv`.

## Closure (2026-10-21): local scope

**M11-E is closed with its local scope explicit:**
- **What was shown:** the tested best-available replacements lose in each winner's own context.
- **What it does not exclude:** better multi-block alternatives, or other correct-copy-number alternatives.
- **K2 and K3:** they show that the **tested** alleles are inconsistent with the reads. They do not show the exhaustive
  absence of every read-consistent, correct-copy-number panel allele.
