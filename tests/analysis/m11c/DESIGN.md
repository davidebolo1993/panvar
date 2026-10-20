# M11-C: confidence-input preparation (frozen before execution; no calibration, no genotyping campaign)

- **Branch:** `m11-c-confidence-inputs`, from `m11-a-confidence-design` `e971768`.
- **Results:** `tests/results/m11c_confidence_inputs_20261020/`.
- **Budget:** 900 s, including 150 s of disclosed development (the diagnostic's builds).
- **Not done:** no simulation, calibration, internal whole-haplotype alignment in the caller, clustering, merge or push.

## 1. Q1 semantics, corrected from the implementation

`genotype_support` comes from `block_genotype_read_supported` (`genotype_repair.cpp` pooled path; `genotype_command.cpp`
unpooled path; same rule).
- **Fragment contender.** A **ranked rival** is a fragment contender unless its fragment term (the read likelihood alone,
  without exposure or prior) is certifiably below the winner's:
  `rival.fragment_term + rival.numeric_error_bound < best.fragment_term − best.numeric_error_bound` excludes it.
- **The rule.** A block is `READS` iff **no** fragment contender has a different unordered genotype at that block (from the
  contender group's representative labelling). Otherwise the block is `PRIOR_ONLY`.

**Correction.** PRIOR_ONLY does **not** mean "no fragment informs the block". It means that, **within the explored and ranked
pool**, at least one rival that the reads alone do not place below the winner carries a different genotype there, so the
call's separation from that rival rests on the exposure and prior terms.
- **READS is equally pool-conditional.** It says nothing about unexplored candidates.
- **It depends on representative labellings,** so synonymous labellings can affect it.
- **M11-B's association results stand** (floor attainment lower for PRIOR_ONLY in every stratum). Only their interpretation
  changes.

## 2. Diagnostic outputs (`block_diagnostics off|diagnostic_v1`, default off)

**Writes** `<prefix>.block_diagnostics.tsv`, one row per model block in model order: `block_index`, `block_id`,
`census_route`, `jp_emission_tie`, plus a header comment and the dosage-audit status.

**`census_route`:**
- **Value:** the dosage census route (`repeated_unit` / `pair` / `not_sensitive`).
- **When:** captured through `DosageShadowOptions::census_routes_out`, set immediately after `run_census`, so it **survives
  a later audit-stage failure**.
- **Otherwise:** an explicit `NOT_EVALUATED:{dosage_audit_off, no_genotype_emitted, unsupported_library, no_rank1_call}`, or
  `UNAVAILABLE:census_not_completed`. A route is never inferred.

**`jp_emission_tie`:**
- **Value:** the JP1.1 eligibility condition (every pair emission tied, more pairs than the local-pair limit), as
  `TIED_ELIGIBLE` / `NOT_ELIGIBLE`.
- **Where:** it is evaluated only by repair searches with JP1.1 on, and carried in `JunctionProposalReport::block_eligible`
  (OR across searches).
- **Otherwise:** `NOT_EVALUATED:{junction_proposals_off, no_repair_search_report}`. No search runs just to fill it.

**Output only:**
- nothing reads these fields;
- the census vector is filled only when the key is on;
- the per-block flags are copied from a condition the search already computed;
- the file is run-owned (removed at run start);
- the settings row is written only when the key is declared.

**Not changed:** selection, scores, likelihoods, priors, confidence, filters, predictions and every existing column and its
precision.

## 3. Gates and smoke runs (saved inputs only; no new donors or reads)

| unit | checks |
|---|---|
| build; JP1.1, compound and sequence-report units; e2e `sequence_first` | unchanged behaviour |
| `smoke.c4.{off,on,on.r}` (M10-A unified D configuration, c4_5) | off/on identity; determinism (repeat) |
| `smoke.c4.dosageoff.on` | census `NOT_EVALUATED:dosage_audit_off` everywhere (no route inferred) |
| `smoke.acot.{off,on}` (D, acot_5, JP1.1 on) | off/on identity; the eligible block equals sv:5 and the status count |
| `smoke.cyp.{off,on}` (M1-B CYP2D6 NA19650 5× command; JP off; the dosage audit fails after the census, "unit library") | off/on identity; routes retained despite the audit failure; JP `NOT_EVALUATED:junction_proposals_off` |
| `cmp.*` | s2l comparator v3.2 (operational, profile, checkpoints; declared setting `block_diagnostics`; the M10-A timer allowances) with the diagnostic file skipped, which must exist only in the on arm |
| `diag_check` | one row per model block in model order; legal values only; route and JP expectations above; status-count agreement |

## 4. Offline evaluation quantities (no caller change)

**Primary outcomes (unchanged):** per-block **excess** and **normalized excess** over the existing availability floor
(`scripts/genotype_report/load.py`). No whole-locus floor.

**QVfrac** (Davide; COSIGT interpretation), offline only:
- `QVfrac = (QV₁_pred + QV₂_pred) / (QV₁_best + QV₂_best)`, with QV = −10·log10(max(0.5, d)/L);
- d and L come from `panvar::nw_edit_distance` (the benchmark convention; L = 0 gives 60).

**The pairs:**
- **predicted:** the reported block alleles paired with the truth haplotypes by the evaluator's minimum summed edits (tie:
  p1–t1, p2–t2);
- **best:** per truth haplotype, the **evaluator's recorded closest allele** over the active model's block alleles (its
  allowed space; LZO: the truth itself).

**Checks:**
- the floor equals closest₁ + closest₂ (verified on all 438 truth-defined current-schema development blocks before freezing);
- the recomputed distances must equal the evaluator's.

**Ties:**
- **Best allele:** the evaluator's (increasing length gap, then allele order, strict improvement). QV-maximization among
  equal-edit alleles is **not** evaluated; it would need a new floor search.
- **Predicted pairing:** pairings where summed QV would prefer the other assignment are counted and reported, never
  substituted.

**Never substituted:** QV divided by the perfect-match ceiling.

**Scope:** the current-schema development runs of M11-B (M1-B, M9-C, M9-G JP-on, M9-H). Truth and allele sequences come
from `genotype_probe block-sequences` with each evaluation's own model, curation, truth graph and donor (the fingerprint
must match).
- **Legacy strata** (C4 production, CYP2D6 pilot): what would be needed is stated, not computed.

**Historical interpretation recorded (Davide):** the COSIGT QV categories (`<17 / 17–23 / 23–33 / >33`) and **availability-
relative** QVfrac quintiles (0–0.2, …, 0.8–1.0 of QVfrac). These are **not** QV / QV_max quintiles. Any categorical
acceptance rule remains a **proposal pending approval**.

## 5. Simplified prospective calibration design (replaces the stratified scheme of CALIBRATION_DESIGN.md)

**Two pre-declared candidates:**
- **C1:** the existing per-block confidence quantity alone (`error_mass_upper`);
- **C2:** C1 plus Q1 (`genotype_support`).

Q2 (census route) and Q3 (emission tie) are **diagnostic only** at first. There is no three-way stratification or merge
scheme.

**Outcomes:** excess and normalized excess primary; QVfrac and floor attainment secondary or descriptive.
- **Nothing fitted here.**
- **Not claimed:** calibrated confidence.
- **The small, ACOT-only validation** cannot establish cross-locus reliability or a very low confident-error rate.
