# M11-D result: unified-configuration development pilot. Neither confidence candidate beats a constant baseline out of family; all five catastrophic blocks receive ≥ 0.91 confidence

**Run.** Branch `m11-d-dev-pilot`.
- **Commits:** freeze `19fabee`; amendments A1 `ece6e4e`, A2 `4ec52ed`, A3 `77f185a`, A4 `fde8c84`; result after.
- **Results:** `tests/results/m11d_dev_pilot_20261020/`.

**Budget:** 2,002.0 of 3,600 s. That is 1,747.0 s of ledger wall plus 255 s of disclosed development: 60 pre-freeze, 15 for
A2, 120 for A3, 60 for A4.
- **Units:** 101 recorded; every caller, evaluation and truth unit passed first time.
- **Watchdog:** no timeout or envelope event.
- **Peaks:** C4 2.53 GiB; CYP2D6 1.21 GiB; ACOT 2.87 GiB.
- **Call wall time:** C4 223 s, CYP2D6 138 s, ACOT 786 s; longest single call 129 s (ACOT HG02723), against its 300 s
  watchdog.

**Failed analysis units (kept recorded):**
- `outcomes`: a caller no-call's NA alleles (A1);
- `fit`: ran without an outcomes table, so its output is invalid (A1);
- `outcomes.a1`: the consistency check on NO_PREDICTION rows (A2); `fit.a1` NOT RUN;
- `report`: zero-length blocks' undefined normalized excess (A4).

None of these fixes changed a cell, the configuration, the label, a candidate or the comparison.

**Untouched:** the 18 ACOT validation families. No new simulation. **Every number here is a development estimate.**

## 1. Inputs, completeness and coverage

| locus | cells | families | completed ranking | resource-incomplete | caller status | assessable runs |
|---|---|---|---|---|---|---|
| C4 | 12 (hash rule; HG00096 ineligible as reference donor) | 12 | 12 | 0 | COMPLETE 12 | 12 |
| CYP2D6 | 6 | 6 | 6 | 0 | COMPLETE 5, CONDITIONAL_GENOTYPE_EQUIVALENT 1 (NA19705) | 6 |
| ACOT | 12 | 12 | 12 | 0 | COMPLETE 10, CONDITIONAL_GENOTYPE_EQUIVALENT 1 (HG02083), **RESCORE_INCOMPLETE 1 (NA21106)** | 11 |
| all | 30 | 29 (NA19705 = family 2368 in C4 and CYP2D6) | 30 | **0** | | 29 |

**No search was resource-incomplete.**

**One run is a completed ambiguous ranking:** ACOT NA21106.
- Caller no-call `divergent_tie`: two candidates are tied on total score and disagree at ctx:18 (`SEQUENCE_TIED`).
- **Not new:** M9-C (`eval.LOO.d5`) has the same no-call with the same note.
- Its 19 blocks are counted as NO_PREDICTION.

**The two CONDITIONAL_GENOTYPE_EQUIVALENT runs** resolve the genotype with a sequence set. They are predicted and assessable.

**Block coverage:**

| | count |
|---|---|
| model blocks | 474 |
| truth-defined | 473 (1 undefined: CYP2D6 NA19650) |
| **predicted** | **454 (96.0%)** |
| NO_PREDICTION | 19 |
| fit rows (assessable, defined C1 input) | 454 |

**C1 input source:** in all 454 rows it is `error_mass_upper`; none printed 0.

**GQ correspondence:** −10·log10(error_mass_upper) vs `legacy_representative_GQ` on 458 blocks with lower = upper. Max |diff|
0.0048; all within 0.01.

## 2. Outcomes (primary: excess and normalized excess)

| locus | predicted | excess sum / mean / median / p90 / max (bp) | normalized excess mean / median / max | floor attained | catastrophic ≥ 1 kb | top (QVfrac ≥ 0.8) | QVfrac quintiles 0–.2 / .2–.4 / .4–.6 / .6–.8 / .8+ | zero-length |
|---|---|---|---|---|---|---|---|---|
| C4 | 132 | 65,606 / 497 / 0 / 6 / 39,155 | 0.0043 / 0 / 0.249 | 83 (62.9%) | 2 | 128 (97.0%) | 0 / 1 / 0 / 3 / 128 | 0 |
| CYP2D6 | 113 | 26,226 / 232 / 0 / 118 / 12,159 | 0.079 / 0 / 0.908 | 70 (61.9%) | 2 | 100 (88.5%) | 2 / 2 / 6 / 3 / 100 | 12 |
| ACOT | 209 | 28,976 / 139 / 0 / 4 / 27,270 | 0.010 / 0 / 0.355 | 140 (67.0%) | 1 | 199 (95.2%) | 0 / 0 / 3 / 7 / 199 | 22 |
| all | 454 | 120,808 / 266 / 0 / 6 / 39,155 | 0.025 / 0 / 0.908 | 293 (64.5%) | **5** | 427 (94.1%) | 2 / 3 / 9 / 13 / 427 | 34 |

**Zero-length blocks** (truth and called alleles all empty: CYP2D6 sv:1, sv:11 and sv:15; ACOT sv:3 and sv:11):
- excess 0; normalized excess undefined (0/0, excluded from its summary); QVfrac = 1 by the L = 0 → QV 60 convention;
- they are trivially top and are part of the top rate.

**QVfrac > 1** in 2 blocks: CYP2D6 NA19983 ctx:12 (1.000012) and ACOT HG02178 ctx:18 (1.000001).
- Both are at the floor, and the predicted edits equal the best edits (3 + 1 and 3 + 3).
- The alignment lengths differ by one or two columns: for example, 9,692 + 9,689 predicted vs 9,690 + 9,689 best. A longer L
  gives a marginally higher QV.
- They are not clipped, and the ratio is not a probability.

**The label hides large errors in long blocks.** CYP2D6 NA19983 sv:9 has excess 12,075 bp but QVfrac 0.864, so it is
"top". Excess and the catastrophic flag catch it; the categorical label does not.

## 3. Development out-of-family comparison (frozen; `runs/fit.a2/performance.tsv`)

**Method:** leave-one-family-out, families pooled across loci. Δlog loss is shown with a family-bootstrap 95% interval
(2,000 resamples). **Negative Δ = better than the baseline.**

| locus | blocks / families | top rate | log loss B0 / C1 / C2 | Brier B0 / C1 / C2 | Δ C1 − B0 | Δ C2 − B0 | Δ C2 − C1 |
|---|---|---|---|---|---|---|---|
| all | 454 / 28 | 0.941 | 0.2276 / 0.2320 / 0.2224 | 0.0562 / 0.0564 / 0.0557 | +0.0044 [−0.0014, +0.0148] | −0.0052 [−0.0160, +0.0065] | −0.0096 [−0.0231, +0.0033] |
| C4 | 132 / 12 | 0.970 | 0.1458 / 0.1500 / 0.1391 | 0.0303 / 0.0307 / 0.0303 | +0.0042 [−0.0009, +0.0134] | −0.0067 [−0.0174, +0.0064] | −0.0109 [−0.0258, +0.0039] |
| CYP2D6 | 113 / 6 | 0.885 | 0.3836 / 0.3812 / 0.3644 | 0.1054 / 0.1052 / 0.1025 | −0.0025 [−0.0034, −0.0014] | −0.0192 [−0.0450, +0.0015] | −0.0168 [−0.0384, +0.0033] |
| ACOT | 209 / 11 | 0.952 | 0.1949 / 0.2032 / 0.1982 | 0.0459 / 0.0463 / 0.0466 | +0.0083 [−0.0021, +0.0275] | +0.0033 [−0.0148, +0.0187] | −0.0050 [−0.0308, +0.0169] |

**Reading:**
- **C1 does not beat the constant baseline pooled.** Its point estimate is worse; the interval includes 0.
- **C2 is nominally better than B0 and C1 pooled, in C4 and in CYP2D6, but every interval includes 0.**
- **The only interval excluding 0 is CYP2D6 C1 − B0 (−0.0025), and it is tiny, from 6 families.**
- **The fitted C2 is essentially a READS / PRIOR_ONLY switch.** Out of family it gives READS blocks about 0.98 and
  PRIOR_ONLY blocks about 0.91–0.92; GQ moves it by about 0.01.

**AUCs of the LOFO predictions** (B0 0.239, C1 0.318, C2 0.457 pooled) are **not interpretable as discrimination**. Leaving
a family out shifts the training prevalence against that family's own outcome, so B0, a constant per fold, already has AUC
0.24. C1's small slope does not overcome that shift.

**In-sample direction** (descriptive only, `report.a4/descriptive_association.tsv`):

| locus | GQ AUC | median GQ top / non-top | PRIOR_ONLY share, top / non-top |
|---|---|---|---|
| all | 0.581 | | 228/427 (53%) / 22/27 (81%) |
| CYP2D6 | 0.61 | | |
| ACOT | 0.61 | 57.0 / 28.2 | |
| C4 | **0.40 (reversed)** | 51.3 / **69.9** | |

C4's non-top blocks carry *higher* GQ than its top blocks.

## 4. Catastrophic failures (excess ≥ 1,000 bp), each with its out-of-family confidence

| block | excess (bp) | normalized | QVfrac | GQ | read support | census route | JP tie | p B0 / C1 / C2 |
|---|---|---|---|---|---|---|---|---|
| C4 NA24385 sv:7 | 39,155 | 0.249 | 0.209 | 276.6 | PRIOR_ONLY | repeated_unit | not eligible | 0.941 / **0.965** / 0.919 |
| ACOT NA19909 sv:13 | 27,270 | 0.240 | 0.749 | 442.7 | PRIOR_ONLY | repeated_unit | not eligible | 0.940 / **0.991** / 0.949 |
| C4 NA19434 sv:7 | 26,210 | 0.148 | 0.683 | 85.0 | PRIOR_ONLY | repeated_unit | not eligible | 0.941 / 0.941 / 0.912 |
| CYP2D6 NA19650 sv:9 | 12,159 | 0.181 | 0.509 | 19.0 | PRIOR_ONLY | repeated_unit | not eligible | 0.944 / 0.942 / 0.922 |
| CYP2D6 NA19983 sv:9 | 12,075 | 0.180 | 0.864 (top) | 159.2 | **READS** | repeated_unit | not eligible | 0.944 / 0.948 / **0.979** |

**Every catastrophic block is a `repeated_unit` census route:**
- 5 of the 29 `repeated_unit` blocks;
- 0 of the 368 `not_sensitive` blocks;
- 0 of the 57 `pair` blocks.

**Neither candidate flags them:**
- C1 gives 0.94–0.99; two of them carry GQ 277 and 443;
- C2 gives 0.91–0.98, and the READS catastrophic block gets C2's highest confidence (0.979).

Q2 is a diagnostic only, and this pilot neither fits nor selects it. **This is the main development finding.**

## 5. Synonymous-boundary artefacts (identified separately; not treated as proven sequence failures)

**Measure:** the interval's summed block distances minus its concatenated-interval distance. **Intervals:** CYP2D6 sv:3..sv:5
and ACOT sv:5..ctx:6 (M9-I inventories). C4: this pilot's model inventory found no degeneracy interval in these 12 models.

| locus | flagged intervals (blocks) | artefact bp, per interval | non-top blocks that are artefacts |
|---|---|---|---|
| CYP2D6 | 5 (10: sv:3 + sv:5 in HG02769, HG03456, NA19650, NA19705, NA19983) | 236, 472, 236, 236, 590 = **1,770** | 10 of 13 |
| ACOT | 5 (10: sv:5 + ctx:6 in HG01975, HG02178, HG02074, NA18971, HG02083) | 1, 4, 1, 4, 6 = **16** | 3 of 10 (sv:5 in NA18971, HG01975, HG02083) |
| C4 | 0 (the model inventory finds no coalescence interval in any of the 12 models) | 0 | 0 of 4 |
| all | 10 (20 blocks) | **1,786** | **13 of 27** |

**`report.a4/blocks_summary.tsv` overstates artefact bp.** Its `artefact_bp` column sums the interval value over every
flagged block in the interval, so each interval is counted twice (3,540 / 32 / 3,572). The per-interval values above are
correct.

**CYP2D6:** in 4 of 5 donors the interval's artefact bp equals the pair's summed excess. For example, HG03456:
177 + 295 = 472. The concatenated interval matches as well as the floor allows, so the block-level excess there is a
boundary-placement effect.
- **NA19705 is the exception:** 236 of its 296 bp of summed excess are artefact, leaving a 60 bp residual.

**ACOT HG00639 and HG02723 sv:5** (artefact bp 0) are genuine small sequence differences, not artefacts.

As frozen, artefact blocks stay in the comparison with their label. Almost half the non-top labels are therefore boundary
effects that no per-block confidence input could be expected to predict. No re-fit without them was run, because none was
declared.

## 6. Diagnostics Q2 / Q3 (cross-tabs only; not fitted)

| input | value | blocks | top | floor | catastrophic |
|---|---|---|---|---|---|
| census route | not_sensitive | 368 | 354 | 247 | 0 |
| | pair | 57 | 48 | 36 | 0 |
| | repeated_unit | 29 | 25 | 10 | **5** |
| JP emission tie | NOT_ELIGIBLE | 437 | 415 | 288 | 5 |
| | TIED_ELIGIBLE | 17 | 12 | 5 | 0 |
| read support | PRIOR_ONLY | 250 | 228 | 133 | 4 |
| | READS | 204 | 199 | 160 | 1 |

## 7. What this does and does not establish

**Established (development, 28–29 families, one configuration):**
- per-block `error_mass_upper` is not, by itself, an out-of-family predictor of the QVfrac ≥ 0.8 label beyond its base rate;
- adding PRIOR_ONLY gives a small, non-significant gain;
- all five catastrophic errors receive high confidence from both candidates.

**Not done:**
- no new scientific no-call, no production CONFIDENT cutoff, no default promoted;
- no caller-scoring, prior, clustering, alignment or floor-search change; no validation simulation, merge, push or cleanup.

**Not claimed:**
- calibrated confidence;
- cross-locus reliability, since CYP2D6 has 6 families and pooled LOFO mixes loci;
- any validation result.

**Decisions for review:**
- whether the evaluation label should exclude zero-length blocks and/or identified synonymous-boundary artefacts. Each would
  have to be declared before any further comparison;
- whether the catastrophic-error concentration in `repeated_unit` blocks justifies declaring Q2 as a candidate in a separate,
  frozen comparison. It was **not** selected here;
- whether a continuous primary target (excess or normalized excess), and not the near-saturated top label (94% positive),
  should drive any later confidence work.
