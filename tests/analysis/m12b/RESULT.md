# M12-B result: dosage-first selection — INCONCLUSIVE by the frozen rule; the count-based copy state matches the truth in 27 of 29 blocks, but sequence-level benefit is shown in only 1 of 3 scored catastrophic blocks

**Run.** Branch `m12-b-dosage-first`; freeze `85dfdf8`; amendments A1 `83e8812`, A2 `42bc0b2`,
A3 `14998cf`; result after.
- **Results:** `tests/results/m12b_dosage_first_20261021/`.
- **Budget:** 1,423.6 of 1,800 s: 200 s of disclosed development plus 1,223.6 s of ledger wall.
- **Units:** 60 recorded — 52 passed, 2 failed, and 6 NOT RUN on the budget reserve.
- **Failed units (both mechanical, amended, preserved):**
  - `gate.features`: a dispatch-name defect (A1);
  - `evaluate`: the CSV field-size limit (A3).
- **Watchdog:** no event. Peak 1.42 GiB (the ACOT count units).

**Not done:**
- no new reads, simulation, caller rerun, integration, cap increase, repaired FASTA, prior change, fit or clustering;
- no merge, push or cleanup; the original calls and every earlier output are untouched.

**The 18 ACOT validation families are untouched.** M12-A's K1 qualification was appended to its RESULT.

**What the selections are:** experimental point selections, frozen in `select` and the score units **before** truth was
opened in `evaluate`. They are not calls, posteriors or confidence.

## 1. Gates

| gate | result |
|---|---|
| selection logic (synthetic) | 7 of 7 pass: total-state grouping (3+1 and 2+2 share a state), both assignments, exact ties retained, order invariance, missing features or losses → unsupported, multi-feature sums |
| score-set candidate order | totals and fragment-term bits identical with the CYP2D6 preflight list reversed |
| copy states = the caller's own census | the probe's pooled features equal the caller's `dosage_feature_copies.tsv` in **22 of 22** blocks where the audit wrote one (C4 and ACOT; the CYP2D6 audit had failed before writing). A real C4 3+1 and 2+2 pair share one total state, so total dosage alone does not resolve allocation |
| every scored block (12) | incumbent total reproduced; incumbent distance = the evaluator's; every scored haplotype's length = flanks + the declared allele; FNV identity for the incumbent and the first and last candidates (6 per block); every candidate the run itself scored carries the same total (e.g. 63/63, 1/1); no alias inconsistencies (no aliased alleles); coverage = 100% of each selected state; no exact ties |

## 2. Coverage and bounds

**All 29 blocks are SUPPORTED:** every block has pooled features and count losses, so no structural alphabet was needed.

**Preflight:** 0.24 s per unique haplotype (C4), 0.31 (CYP2D6), 0.64 (ACOT). The admission bound was an estimate of
≤ 240 s.

| locus | blocks | scored | not scored |
|---|---|---|---|
| CYP2D6 | 6 | 6 (43–45 s each) | — |
| C4 | 12 | 6 (B07, B08, B09, B11, B15, B18; 7–35 s) | **6 on the A2 budget reserve: B10, B12, B13, B14, B16 and B17 (K1)** |
| ACOT | 11 | 0 | **11: LIMITATION_STATE_TOO_LARGE** (306–672 unique haplotypes, 1,224–55,080 candidates; estimates 300–649 s). Never shortlisted |

**Catastrophic blocks:** three were scored — B04 (K4), B06 (K5) and B18 (K3). Two were not: B17 (K1, budget) and B21 (K2,
bound).

## 3. Results

### 3a. Copy-state selection (truth-free) against the truth's copy state, all 29 blocks

The truth state is assessable in all 29 blocks: each feature's k-mers sit at one integer count in the truth.

| | blocks | the count rule selects the truth's total state |
|---|---|---|
| original state already correct | 24 | **24 / 24 kept**: the rule changed no correct state |
| original state wrong (the five catastrophic blocks) | 5 | **3 / 5 fixed**: B04 K4 (4→5 units), B17 K1 (5,4,0 → 4,4,0), B21 K2 (3,2,2,1,1,1 → 4,2,2,2,2,2) |
| | | 2 / 5 moved to a different wrong state |

**The two misses:**
- **B06 K5:** (5,3,3,2,2) → (4,2,2,2,2); the truth is (4,2,2,1,2). **No panel pair has the truth's state**, and no allele
  matches truth homologue 2's state, so this is an availability limit.
- **B18 K3:** (5,3,0) → (6,2,0); the truth is (6,3,0). **40 panel pairs have the truth's state.** The count loss preferred
  a wrong one, so this is a count-model miss on the second feature.

**Content sensitivity.** Every block except B07, B15 and B17 is **CONTENT_SENSITIVE**: the state margin is smaller than the
median-to-best spread of losses inside the selected state. **So the count preference is not a pure copy-number
measurement.** Sequence content materially shapes it, even where it selects the right state.

### 3b. Sequence-level outcome, 12 scored blocks (exhaustive within the selected state; unchanged fragment + exposure + prior)

| block | original pair → selected | state: original → selected (truth) | excess (bp): original → selected | Δtotal vs incumbent (fragment / exposure / prior) | outcome |
|---|---|---|---|---|---|
| **B04 K4** CYP2D6 NA19650 | 23,73 → **22,80** | (4,2,2,2,2) → **(5,3,3,2,2)** (= truth; allocation matches) | **12,159 → 17** | −3.6 (+155.8 / −152.1 / −7.4) | **FIXED_LARGE** |
| **B06 K5** CYP2D6 NA19983 | 95,97 → 104,97 | (5,3,3,2,2) → (4,2,2,2,2) (truth (4,2,2,1,2), not available) | **12,075 → 24,240** | −37.0 (−170.3 / +152.1 / −18.8) | RETAINED_LARGE, **worse** |
| **B18 K3** C4 NA19434 | 102,60 → 102,27 | (5,3,0) → (6,2,0) (truth (6,3,0)) | **26,210 → 6,263** | −205.9 (+53.6 / −250.1 / −9.5) | RETAINED_LARGE, reduced |
| B01 B02 B03 B05 (CYP2D6), B07 B08 B09 B11 B15 (C4) | unchanged (B05 is the homologue swap 20,8 → 8,20; same genotype) | unchanged (= truth) | 0, 0, 2, 27, 0, 0, 0, 9, 2 — all unchanged | 0 | UNCHANGED |

**Totals over the scored blocks:**
- 1 fixed large error; 2 retained large errors (1 reduced, 1 worse);
- 0 introduced large errors;
- **0 of 8 near-floor controls damaged**, and every scored control unchanged.

**Total copy number, allocation, sequence and phase are not the same thing:**
- **B04:** total state and allocation both match the truth; the sequence is near the floor (excess 17).
- **B18:** the unit total is now right but the second feature is wrong, and the sequence is still 6.3 kb in excess.
- **B06:** the state is wrong and the sequence is worse.
- **Phase** across the block was not assessed.
- **Getting total dosage right does not establish correct haplotypes.**

## 4. Decision (frozen rule): **INCONCLUSIVE**

- **Why not MERITS:** 3 catastrophic blocks were scored, but only **1** was fixed; MERITS needs at least 3.
- **Why not FAILS:** nothing was introduced, so introduced (0) < fixed (1), and fixed is not 0.

**What the evidence supports:**
- **For the copy-state step:** truth-free count selection keeps all 24 correct states and recovers 3 of the 5 wrong ones.
  The other 2 are an availability limit and one count miss.
- **For the full rule:** one clear sequence-level fix (K4), one reduction (K3) and one worsening (K5). K5's worsening
  follows from forcing a state the panel cannot realise.
- **Untested:** K1 and K2, whose copy states the count step fixes, were not scored, because of the budget (K1) and the
  exhaustive bound (K2).

**Honest failures:**
- the rule worsened K5;
- its count loss mis-chose K3's second feature even though the truth state was available;
- most count preferences are content-sensitive;
- exhaustive within-state scoring does not fit ACOT at all under these bounds, and fits only part of C4 inside 1,800 s.

**No post-hoc change** was made to features, loss, copy-state rules, ties or bounds.

**Status of the evidence:** inspected development evidence on 29 blocks, not validation or calibrated confidence.

**What a resolution would need** (not run): the 6 C4 blocks lost to budget, including K1, scored under the same frozen rule
with more time. That would test whether K1's correct state also yields a correct sequence.

## Corrections (M12-C, 2026-10-21): interpretation only; results and decision unchanged

1. **K5 (B06).** "No panel pair has the truth's state … an availability limit" means only that **no available pair matches
   the truth's full pooled-feature vector**.
   - **What it does not mean:** that the correct total copy number is unavailable.
   - The selected pair 104/97 and the evaluator-best pair 95/10 have the **same** state, (4,2,2,2,2,0,0). Both differ from
     the truth only in feature 3, which every panel allele carries at least once.
   - Their feature-0 totals equal the truth's (4).
2. **K5's worsening.** "Follows from forcing a state the panel cannot realise" is incorrect as an explanation.
   - **The evaluator-best pair 95/10 was in the selected state and was scored.**
   - The unchanged objective preferred 104/97 by 112.3 nats, and the count loss preferred it too.
   - The worsening comes from the rule's **within-state** choice.
3. **Safety.** "0 introduced large errors" is not evidence of safety. K5, an already-catastrophic prediction, got worse
   (12,075 → 24,240). M12-C adds K1 (39,155 → 65,462).

Details: `tests/analysis/m12c/RESULT.md`.
