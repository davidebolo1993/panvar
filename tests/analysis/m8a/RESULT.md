# M8-A result: copying-prior calibration feasibility. Direct predictive fitting is not supported; a masked-emission conditional likelihood is the smallest workable alternative

**Run.** Branch `m8-a-copying-calibration`: freeze `77c3207`, plus `masked.py` (declared after `support.full`; no
fitting). Results in `tests/results/m8a_copying_calibration_20261020/`.

**Execution.** 34.2 of 300 s: audit_gates 0.2, support.full 13.2, support.loo 19.0, masked 1.8.

**Scope.** No fitting, genotyping, caller or default change, clustering, residual editing, merge or push.

**M7 status, recorded at Sol's review.** M7 is closed with M7-C accepted. No further Locityper experiments.

## 1. Reuse (not reimplemented)

| component | role |
|---|---|
| `build_copying_transitions` / `CopyingTransitionSpec` | `constant_v1` (r per boundary; default r = 0.01) and `distance_ref_core_v1_experimental` (r_b = 1 − exp(−d_b / λ)) |
| `LocusGeometry` and the M4-A geometry outputs | block ids and reference core midpoints |
| `ls_hmm.hpp` | `log_haploid_mosaic_prior` (fixed allele path) and `log_haploid_sequence_prior` (summed over every allele path that spells a sequence; −inf if none) |
| M4-A gates G1–G10, M4-B ctests | existing correctness checks |
| `genotype_probe block-sequences` | allele carriers and truth projections used in this batch |

## 2. Estimand

The transition is T_b = (1 − r_b) I + r_b U, where U is uniform over all n template rows: the current row included and
duplicate rows counted. The two candidate parameterisations:

| parameter | definition |
|---|---|
| **constant redraw parameter r** | the same donor-redraw probability at every block boundary |
| **copying length λ** (distance mode) | r_b = 1 − exp(−d_b / λ), with d_b the distance between consecutive reference core midpoints; d_b = 0 means no redraw |

Both are **copying-model parameters, not biological crossover probabilities.**

**Donor redraws vs observable allele changes.** A redraw is visible only when the newly drawn template carries a
different allele at a later informative block.
- It is invisible when it lands on the same row (probability r/n), on a duplicate or identical path, or on any template
  sharing the held-out allele at every subsequent block up to the next difference.
- So the observable change rate is below r, and it depends on allele diversity and panel composition.

**What changes the interpretation of a fitted value:**
- **Panel size n:** U includes the current row. Larger panels make redraws more often visible and give the held-out
  haplotype closer donors.
- **Duplicate or identical haplotypes:** these add weight to a lineage and make redraws among them invisible. C4 and
  ACOT each contain one pair of identical complete paths.
- **Ancestry composition:** a well-represented population yields long copying tracts, and therefore smaller fitted r or
  larger λ.
- **Relatedness:** relatives give near-trivial copying and must leave a fold together.
- **Block and allele granularity:** blocks are locus-specific, and alleles are exact byte sequences, so fitted values are
  not portable across loci or graph builds.

## 3. Family-aware leave-out inventory (`support.full`; IGSR pedigree; `DEVELOPMENT_DONORS.tsv`, 63 donors)

| locus (model) | blocks | templates | donors | families (multi-donor) | uninspected | inspected (development) | reference genomes |
|---|---|---|---|---|---|---|---|
| C4 (`c4.gfa`) | 11 | 131 | 67 | 61 (3 trios) | **6 families / 12 haplotypes** | 53 | 2 |
| CYP2D6 (`cyp2d6.gfa`, curated) | 19 | 126 | 67 | 61 (3 trios) | **6 / 12** (the same 6 families as C4) | 53 | 2 |
| ACOT (`acot.gfa`) | 19 | 467 | 235 | 235 (0) | **225 / 450** | 7 | 3 |

- The uninspected C4 and CYP2D6 families are VN047 (HG02018), NG109 (HG03520), ST116 (HG04036), 2368 (NA19705), 2424
  (NA19836) and 2484a (NA20355).
- All three trios are inspected development families. No relative crosses a fold.
- No template has an incomplete or ambiguous block projection.

## 4. Can held-out sequences be assigned probabilities?

**Mostly no.** Support is judged against templates outside the haplotype's own family, by allele label and by allele
bytes; unsupported positions are never omitted silently. Table: `support.tsv` / `support_summary.tsv`.

| locus | status | haplotypes | **fully supported** | unsupported blocks (total; per haplotype) | label-only unsupported | identical full path in training |
|---|---|---|---|---|---|---|
| C4 | uninspected | 12 | **0** | 51 (4.3) | 0 | 0 |
| C4 | inspected | 117 | **0** | 529 (4.5) | 0 | 0 |
| CYP2D6 | uninspected | 12 | **0** | 70 (5.8) | 0 | 0 |
| CYP2D6 | inspected | 112 | **0** | 596 (5.3) | 0 | 0 |
| ACOT | uninspected | 450 | **4 (0.9 %)** | 2,285 (5.1) | 0 | 0 |
| ACOT | inspected | 14 | 0 | 76 (5.4) | 0 | 0 |

**Ambiguity types:**
- **Sequence ambiguity:** none at block level. No block holds two distinct alleles with identical bytes, and no case is
  supported by bytes but not by label.
- **Cross-block spelling ambiguity** (one sequence spelled by allele paths with different block boundaries) is not
  measured here. `log_haploid_sequence_prior` aggregates it by construction.

**Frame-leakage check (`support.loo`).** Existing LOO models cover 5 of 6 uninspected C4 donors and 1 CYP2D6 donor; ACOT
has none. Projecting the true haplotypes onto the LOO block chain:
- every truth block has a frame (0 projection failures);
- **every haplotype still has 3–7 unsupported blocks.**
- The unsupported blocks are mostly large context blocks and the structural block: C4 ctx:6, ctx:10 and sv:7; CYP2D6
  ctx:2, sv:9, ctx:10, ctx:12 and ctx:14.

So lack of support is a property of the panel, not of the graph frame.

**Consequence.** With no mutation or emission model, the copying-only prior gives **probability 0 to 100 % of C4 and
CYP2D6 held-out haplotypes, and to 99.1 % of ACOT's.**
- The held-out log score is −∞, so **direct predictive fitting of r or λ is undefined.**
- That holds whether the probability is computed on allele paths or with `log_haploid_sequence_prior`.

## 5. The smallest alternative: masked-emission conditional likelihood

**Definition:**

  L(θ) = P(alleles at the supported blocks | training templates, θ)

The emission at each unsupported block is set to 1 for every template; it is marginalised as missing data. Transitions
across masked blocks compose exactly, because T_b products remain in the same family.

**What it is not:**
- not a mutation model;
- not a nearest-panel projection;
- not clustering;
- not a new inference engine: it is the existing Li–Stephens forward pass with missing emissions.

Masked counts are always reported. They are about 39–41 % of C4 blocks, 28–31 % of CYP2D6 blocks and 27–29 % of ACOT blocks per
haplotype.

**Identifiability signal (`masked`).** The minimum number of template switches over training templates needed to explain
only the supported blocks:

| locus | status | haplotypes | supported blocks (median) | min switches: count | needing ≥ 1 switch |
|---|---|---|---|---|---|
| C4 | uninspected | 12 | 7 | 0: 9, 1: 2, 2: 1 | 3 |
| C4 | inspected | 117 | 6 | 0: 57, 1: 45, 2: 11, 3: 4 | 60 |
| CYP2D6 | uninspected | 12 | 14 | 1: 2, 2: 6, 3: 3, 4: 1 | **12** |
| CYP2D6 | inspected | 112 | 14 | 0: 7, 1: 27, 2: 33, 3: 29, 4: 12, 5: 3, 6: 1 | 105 |
| ACOT | uninspected | 450 | 14 | 0: 56, 1: 109, 2: 147, 3: 101, 4: 25, 5: 10, 6: 2 | **394** |
| ACOT | inspected | 14 | 14 | 0: 1, 1: 3, 2: 3, 3: 3, 4: 3, 5: 1 | 13 |

- **CYP2D6 and ACOT carry switch information,** so r and λ should be identifiable there under L(θ).
- **C4 is weakly informative:** most supported-block paths need no switch, so a fitted C4 value will likely sit near the
  lower grid boundary.
- **Prior calibration cannot address C4's missing content.** Unsupported blocks such as sv:7 are outside what any copying
  parameter can spell.

## 6. Proposed next experiment (M8-B): one frozen predictive comparison (not run)

**Minimal implementation:**
- a masked-emission overload of the existing haploid forward pass in `ls_hmm`. Missing blocks get all-ones emissions;
  compatibility at observed blocks is by sequence identity. It is gated against a brute-force sum over template histories
  on tiny panels, and with no masked blocks it must equal `log_haploid_mosaic_prior`;
- a diagnostic probe that builds a family-excluded template set from the model's `PanelAlleleMatrix` and evaluates L(θ)
  over a grid, using `build_copying_transitions` for both modes.

No caller or default change.

**Models:** (1) legacy `constant_v1`, r = 0.01; (2) fitted constant r; (3) fitted distance λ.

**Grids (frozen):**
- r ∈ {0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3};
- λ ∈ λ_HM × {1/30, 1/10, 1/3, 1, 3, 10, 30}, with λ_HM, the hazard-matched control length (C4 1.87 Mb, CYP2D6 1.13 Mb,
  ACOT 1.22 Mb), as scale only.

**Fitting (per locus).**
- **Objective:** Σ log L(θ) over haplotypes, each scored with templates excluding its whole family.
- **Training/validation:** fitted on **training families only** (inspected development families plus the non-selected
  uninspected families), by inner leave-one-family-out.
- **Ties:** the smaller r, or the larger λ.
- **Boundary check:** a selection at a grid edge is reported as BOUNDARY. The grid is not extended in that batch.

**Common evaluation population (all three models).**
- C4: the 6 uninspected families (12 haplotypes).
- CYP2D6: the same 6 families (12 haplotypes).
- ACOT: 30 of the 225 uninspected families, chosen by sha256("M8-B:" + family) before scoring.
- Every haplotype is scored; masked counts are reported per haplotype.
- The block frame is the full graph with the held-out family's templates removed. Frame leakage is disclosed: §4 shows it
  does not change support.

**Comparisons:** paired Δ log L per haplotype and family, for fitted-constant vs legacy and fitted-distance vs
fitted-constant, separately per locus.
- **Meaningful:** ≥ 1 nat per haplotype.
- **Consistent:** positive in ≥ 75 % of evaluation families, with no family > 50 % of the gain.

**Stopping criteria:**
- no grid extension or re-fitting within the batch;
- if neither fitted mode beats legacy meaningfully and consistently, keep legacy;
- if a locus selects a boundary value, report it and stop that locus.

The CYP2D6 NA19650 failure is a later diagnostic only, never a fitting target.

**Expected cost:** O(blocks × templates) per haplotype and grid point, so seconds per locus. Proposed bound 300 s.

## 7. M7-C audit maintenance (`audit_gates`; all PASS)

`tests/analysis/m7c/audit_gates.py` recomputes the attempt reconstruction from the archived M7-B.1 outputs and
**enforces** these gates in all 11 runs:

| gate | tolerance | observed max |
|---|---|---|
| mean | ≤ 1e-6 | 3.0e-8 |
| lik_sd | ≤ 1e-4 | 1.4e-6 |
| weighted score | ≤ 1e-6 | 1.45e-7 |
| attempts | exactly 1–20 | all exact |

All six perturbation controls are detected in every run:
- +1e-3 on one attempt's `lik`;
- `lik_sd` + 1e-2;
- population instead of sample variance;
- +1e-3 on one attempt's `aln_lik`;
- depth weight 1.84 instead of 1.85;
- one attempt removed.

No genotyping was re-run and no historical package was altered.

## Decision

- **Direct predictive calibration is not supported by the current model:** held-out sequences have probability 0.
- The smallest workable route is the **masked-emission conditional likelihood** (§5), which is informative for CYP2D6 and
  ACOT and weak for C4.
- **Proposed next experiment:** M8-B, the frozen three-model comparison in §6, including its minimal masked-forward
  implementation and gates.
- Stopped before fitting.
