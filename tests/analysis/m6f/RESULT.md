# M6-F result: whole-component copying beats independent sites decisively; switching adds little and inconsistently, so C1 suffices

**Run.** Branch `m6-f-copying-predictive`, freeze `abdeb14`, no amendments. Results in
`tests/results/m6f_copying_predictive_20261019/`.

**Budget.** 29.0 of 600 s, 4 threads. All 4 units exited 0 on the first attempt.

**Scope.** Panel sequences and the IGSR pedigree only. No reads, read likelihood, truth-selected constructions, production
change, flag, reconstruction, merge or push. The M6-E qualifications are appended to its RESULT.md.

**What is measured.** Predictive probability of held-out component bases **conditional on** the supplied role, length class
and declared alignment mask. It is not full-sequence probability and not genotyping accuracy.

## Holdout (`prepare`; all checks pass)

**Pedigree.** IGSR `20130606_g1k_3202_samples_ped_population.txt` (SHA-256 `4e164b71…`; URL and retrieval time in
`inputs/SOURCE.tsv`).
- It confirms SH032 (HG00512, HG00513 → HG00514), PR05 (HG00731, HG00732 → HG00733) and Y117 (NA19239, NA19238 → NA19240).
- **60 families**: 3 with several donors; 4 donors absent from the pedigree are their own families.

**Training-only families.**
- Every development family: 54 donors with sources in `DEVELOPMENT_DONORS.tsv`, including the template donors of all ten
  mosaic plans.
- chm13v2 and grch38.
- HG00096 (coordinate reference).

**Evaluation families.** 14 families were eligible. The first 8 by sha256("M6-F:" + id) were taken before any scoring:
GB40, GB66, BB05, HG00864, NG98, PEL38, NA19036, 1328.

**Weights.** Founder policy: the trio children HG00514, HG00733 and NA19240 have weight 0, everyone else 1. Identical
sequences in unrelated donors keep their full multiplicity.

**Coordinate references** come from training-only families:

| stratum | reference | evaluation copies identical to the reference (kept) |
|---|---|---|
| period/long | period_u13 | 0 |
| period/short | period_u27 | 0 |
| prefix | prefix_u6 | 5 |
| tail/long | tail_u74 | 0 |
| tail/short | tail_u9 | 1 |

**Evaluation copies:** 47 in total, none unassessable.

| stratum | copies | families | coverage (scored / reference bp) | excluded by the training mask (median) | excluded by the held-out alignment (max) | training-variable sites scored (median) | novel substitutions at training-monomorphic positions |
|---|---|---|---|---|---|---|---|
| period/long | 16 | 8 | 0.987 | 414 | 31 | 284 | 6 |
| period/short | 1 | 1 | 0.989 | 289 | 0 | 125 | 1 |
| prefix (78 bp) | 15 | 8 | 0.692 | 24 | 0 | 4 | 0 |
| tail/long | 4 | 3 | 0.988 | 382 | 0 | 186 | 5 |
| tail/short | 11 | 6 | 0.991 | 241 | 0 | 194 | 2 |

## Gates (`gates`: all pass)

| gate | result |
|---|---|
| G1 normalisation (4⁵ sequences, every model and grid point) | max error ≤ 5.1e-15 |
| G2 forward = exhaustive hidden paths; compressed = uncompressed | 7.1e-15; 5.7e-14 |
| G3 ρ = 0 reproduces C1 | 3.4e-13 |
| G4 redraw at every site reproduces C0 | 1.1e-13 |
| G5 held-out isolation (26 family × stratum folds, held-out copies scrambled) | fits identical |
| G6 aliases (fixture; real training-only library) | 1.1e-13; 1.1e-13 |
| G7 same observations (2,413 `score_all` calls) | one digest per call |

## Selected parameters (training-family inner validation only)

- **ε = 1e-4 for every model in every non-prefix fold.** This is the **lower edge of the frozen grid**, so the residual
  rate is not bracketed and calibration is limited by the grid. The comparison is still paired on identical observations
  and identical grids.
- **C2's ρ = 1e-4 per bp** in every long and short stratum: interior to the grid, about 3 redraws per 32-kb unit.
- **Prefix:** C1's ε = 1e-3; C2's ε = 1e-4 with ρ = 1e-3 (upper edge of the ρ grid).
- Parameters are in `runs/evaluate/m6f.fits.tsv`.

## Paired predictive log-score differences (nats; positive = second model better)

**Per stratum:**

| stratum | families | copies | C1 − C0: sum, per copy, families positive | C2 − C0: per copy, families positive | C2 − C1: sum, per copy, per base, families positive, largest share | stratum flags for C2 − C1 |
|---|---|---|---|---|---|---|
| period/long | 8 | 16 | 790.1, **49.4**, 8/8 | 51.6, 8/8 | 36.1, **2.25**, 1.1e-4, **6/8**, 0.46 | meaningful, consistent |
| period/short | 1 | 1 | 49.2, 49.2, 1/1 | 57.9, 1/1 | 8.6, 8.6, 3e-4, 1/1, 1.00 | one family only |
| prefix | 8 | 15 | 1.2, 0.08, 6/8 | 0.12, 6/8 | 0.66, 0.04, 8e-4, 8/8, 0.15 | not meaningful |
| tail/long | 3 | 4 | 118.4, 29.6, 2/3 | 36.0, 3/3 | 25.5, 6.37, 2e-4, 2/3, 0.73 | not consistent (one family dominates) |
| tail/short | 6 | 11 | 474.9, **43.2**, 6/6 | 43.9, 6/6 | 7.9, **0.72**, 3e-5, **2/6**, 0.83 | neither |
| **pooled** | **8** | **47** | **1,433.9, 30.5, 8/8** (largest share 0.17) | 32.2, 8/8 | **78.7, 1.67, 5/8**, 0.41 | **not meaningful, not consistent** |

Mean per-base log scores (C0 / C1 / C2):
- period/long −0.00219 / −0.00066 / −0.00059;
- tail/short −0.00234 / −0.00068 / −0.00066;
- tail/long −0.00217 / −0.00125 / −0.00106.

**Every evaluated family, pooled over strata:**

| family | copies | C1 − C0 | C2 − C0 | C2 − C1 |
|---|---|---|---|---|
| 1328 | 6 | 119.8 | 143.3 | +23.5 |
| BB05 | 6 | 149.4 | 146.5 | −2.9 |
| GB40 | 6 | 212.6 | 218.1 | +5.5 |
| GB66 | 6 | 205.1 | 202.9 | −2.2 |
| HG00864 | 4 | 174.9 | 185.8 | +11.0 |
| NA19036 | 6 | 239.4 | 233.4 | −5.9 |
| NG98 | 6 | 162.9 | 200.0 | **+37.2** |
| PEL38 | 7 | 169.8 | 182.4 | +12.7 |

- C2 − C1 per family in period/long: +5.0, +0.2, +10.1, −5.4, +5.3, −5.8, +21.9, +4.9.
- C2 − C1 per family in tail/short: −2.3, −3.2, −4.8, +3.1, −0.2, +15.2.

The per-copy table is `runs/evaluate/m6f.copies.tsv`. Evaluation took 8.9 s in total, about 0.8 s per period/long fold
including inner validation.

## Decision (frozen rule): **C1_SUFFICES**

1. **Does linked copying beat independent sites? Yes, decisively and consistently.**
   - C1 − C0 is +30.5 nats per copy pooled, positive in **8/8 families**, with the largest family share 0.17.
   - Per stratum: +49 per period/long copy and +43 per tail/short copy, all families positive.
   - Only the 78-bp prefix (4 variable sites) shows no meaningful difference.
2. **Does switching improve on a whole-component mixture? Only slightly, and not consistently.**
   - Pooled C2 − C1 is +1.67 nats per copy (below the 2.0 threshold), with **5/8 families positive**. Three families are
     worse with switching: BB05, GB66 and NA19036.
   - Within period/long alone, C2 − C1 meets the stratum flags: +2.25 per copy, 6/8, largest share 0.46. That is about
     1e-4 nats per base, so the gain is small.
   - In tail/short, switching is worse in 4 of 6 families.
3. **Consistent, or driven by a few families?**
   - The C1 gain is broad.
   - The C2 gain is concentrated: NG98 (+37.2) and 1328 (+23.5) provide 61 of the 78.7 nats.
   - **Post hoc** (not part of the rule): excluding the prefixes raises the pooled mean to about 2.4 per copy, but the same
     three families stay negative. The decision is unchanged.

**Recommendation.**
- Prefer **C1**, a weighted mixture over whole training components with residual substitutions, as the sequence-prior
  ingredient.
- Do not pursue within-unit switching (C2) now. Its benefit is small (≈ 1e-4 nats/base), family-dependent, and positive in
  only one stratum.
- **No next test is recommended by this batch**, by the rule's own terms.
- If C1 is carried forward, its residual rate needs a grid that brackets the optimum. Every fold chose the lowest ε
  (1e-4). This would have to be declared as new work, not as tuning of this result.

**Not claimed.**
- A positive C1 result validates one sequence-prior ingredient, under the stated conditioning (role, length class,
  declared mask; substitutions only; no indels, segments or structure).
- It does not show recovery of C4 NA19650, any genotyping improvement, or that κ, a structural prior or an escape
  distribution would work.
- Eight families are a small evaluation, and the short-period (1 copy) and long-tail (4 copies) strata are especially thin.
