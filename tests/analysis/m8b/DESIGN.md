# M8-B: bounded exploratory copying-prior calibration on family-held-out panel haplotypes (frozen before any unit)

**Ancestry.** Branch `m8-b-copying-calibration` from `8caf884` (M8-A, kept as the historical record).

**Scope.** Exploratory. A favourable masked score would only authorise considering a later genotyping validation. No
parameter or default promotion.
- **Excluded:** caller integration, genotyping, mutation models, clustering, residual editing, merge, push and cleanup.
- Budget: 300 s.

## Code (diagnostic additions only)

- **`ls_hmm`:** `log_haploid_mosaic_prior_masked`, the existing labelled-path forward recursion with all-ones emissions at
  masked blocks. Every transition, including those around masked blocks, is kept with its distance. It is additive, and
  callers are untouched.
- **`genotype_probe copying-heldout`**, three modes:
  - `gates`;
  - `fit`: inner leave-one-family-out training profiles, selection and outer evaluation; `--scramble-eval` for the
    isolation test;
  - `support`: the whole-sequence prior on family-excluded models of the same graph frame.
- **Reused unchanged:** `build_copying_transitions`, `hazard_matched_copying_length`, the M4-B `LocusGeometry`,
  `log_haploid_mosaic_prior`, `log_haploid_sequence_prior`, `spell_panel_template`, `build_locus_input`
  (`excluded_templates`) and `prepare_genotype_model`.

## 1. Corrected inventory (`inventory`, before any parameter scoring)

**Development donors** are the union of:
- authoritative sample-level records under `tests/results`: per-donor directories, LOO model files, mosaic `stage.log`
  donors, simulation donors and campaign signature `donors` rows;
- M8-A's list.

Status expands to whole IGSR pedigree families: any used member makes the family DEVELOPMENT.

**Asserted:** NA19705 and HG03520 are development donors.

**Evaluability (frozen):** a locus is evaluated only with ≥ 4 uninspected families. ACOT evaluates the first 30
uninspected families by sha256("M8-B:" + family). Otherwise the locus has **no outer evaluation (INFEASIBLE)** and only
training profiles are reported.

**Disclosed:** M8-A already examined support and minimum-switch summaries over the then-proposed evaluation population.

## 2. Whole-sequence support (`support`)

**Check.** `log_haploid_sequence_prior` at constant r = 0.01, on a model of the **same graph** with every eval template
and the queried family excluded. A finite value means the whole sequence is supported.
- Bounded at 25 s per locus. Eval families come first; families not reached are reported UNTESTED.
- Whole-sequence absence is never inferred from block absence.

**Negative control (gate G6):** query blocks AB | C against training blocks A | BC.
- The labelled path must be unsupported (−∞).
- The whole sequence must be supported (finite).

## 3. Training isolation

- **Outer eval families** are excluded from fitting queries, the fitting template bank and selection.
- **Inner leave-one-family-out:** each query family's whole family is removed from its bank.
- **Outer evaluation** is against the frozen training bank.
- **Masks** are recomputed against every actual bank.
- **Isolation test:** `fit_scrambled` perturbs every eval row; training profiles and per-family contributions must be
  byte-identical.
- **What stays fixed and is derived from the full panel, evaluation families included:** the graph frame (blocks,
  allele labels, bubbles) and the reference-core geometry used for distances.

## 4. Objective

Σ over query haplotypes of log P(alleles at observed blocks | bank, θ), with masked blocks given all-ones emissions.

- **Labelled-path score, not sequence-aggregated.** Within the full frame no block holds duplicate allele sequences
  (M8-A).
- **It is a score of observed supported blocks, conditional on the declared masking procedure.** It is not a predictive
  probability of novel haplotypes.
- **Selection limitation:** the masks depend on observed content and panel support.
- Haplotypes scoring −∞ (possible only with r_b = 0) are counted and reported.

**Gates (all must pass before fitting):**

| gate | check | tolerance |
|---|---|---|
| G1 | exhaustive tiny-history oracle, 200 fixtures, including duplicate rows | ≤ 1e-10 |
| G2 | no masks = `log_haploid_mosaic_prior` | ≤ 1e-12 |
| G3 | all masked = 0 | |
| G4 | a masked interior block = removing it with its boundaries composed as 1 − (1 − r₁)(1 − r₂) | |
| G5 | donor permutation invariance; duplicating every row leaves the score unchanged (declared uniform row weighting; partial duplicates are covered by G1) | |
| G6 | AB \| C whole-sequence control | |

## 5. Grids, selection, evaluation (unchanged from M8-A)

| model | values |
|---|---|
| legacy constant | r = 0.01 |
| fitted constant | r ∈ {0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3} |
| fitted distance | λ = λ_HM × {1/30, 1/10, 1/3, 1, 3, 10, 30} |

- **Ties** within 1e-9 go to the smaller r / larger λ. A grid edge is reported as BOUNDARY, with no extension.
- **Reported:**
  - complete profiles;
  - flatness;
  - selection change when each training family is removed;
  - every eval family, including all-masked ones, with its mask fraction and paired results;
  - meaningful: ≥ 1 nat per haplotype;
  - consistent: ≥ 75 % of families positive and largest positive share ≤ 0.5.
- **Never used for fitting:** the CYP2D6 NA19650 failure and C4's missing-content problem.
