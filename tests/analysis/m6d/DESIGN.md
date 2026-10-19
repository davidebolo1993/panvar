# M6-D: candidate-centred vs common-content nuisance prior, synthetic diagnostic only (frozen before any unit)

**Ancestry.** Branch `m6-d-prior-audit` from `20d0916` (M6-C result), with the M6-C qualifications appended to
`tests/analysis/m6c/RESULT.md`. M6-C's outputs are read only, never rewritten.

**Scope.** One diagnostic executable, `tests/nuisance_prior_audit.cpp`.
- **Copied from M6-C verbatim:** the world, scaffold, read simulator, catalogue and the four fixtures (checked textually
  against `tests/nuisance_structural.cpp`).
- **Not in scope:** real loci, caller or production changes, flags, reconstruction, probability tuning or grids, fitting,
  candidate discovery, additional sites, factorisation, merge, push and cleanup.

## Comparison (frozen)

For each scaffold G, every **absolute** sequence state θ (one allele per (copy, site); 2^(2(nA+nB)) ≤ 1,024 states) is
enumerated and scored **once**:
- with the existing paired-fragment likelihood (`score_diplotypes`: placements summed, exposure, background);
- distinct diplotypes are interned.

Two priors are then applied to the same likelihood table:

| prior | P(θ \| G) | depends on the scaffold's base annotation |
|---|---|---|
| original (M6-C) | Π over (copy, site) of 0.05 if θ ≠ base, 0.95 if θ = base | yes |
| common-content control | Π over (copy, site) of 0.05 if θ = World.alt, 0.95 if θ = World.ref; identical for every scaffold | **no** |

The common-content prior is a **diagnostic control, not a proposed biological prior**.

**Unchanged from M6-C:** the likelihood, reads (M6-C's deterministic seeds), the catalogue of nine allocations and the
uniform allocation weights.
- These weights induce total-copy prior masses of 2/9, 3/9 and 4/9 for totals 3, 4 and 5. The summary reports this, and
  also the per-total evidence (the posterior odds between totals divided by the induced prior odds).

**Attribution identity (exact, reported per scaffold and prior):**

  log m(G) = E_q[log L] − KL(q ‖ P(θ | G)),   q(θ) = P(θ | G) L(θ) / m(G)

The best state is never substituted for it.

**State table.** Each run writes a compact table, one row per (scaffold, θ), sufficient to reproduce both
marginalisations:
- absolute alleles, alt count, flips from base, diplotype slot;
- log L, fragment and exposure terms;
- both log priors.

## Gates (`gates`; M6-C's gate world 5, 200-bp unit, reads seed 99)

| gate | checks |
|---|---|
| HA | both priors normalise within every scaffold: all 9 allocations × both base annotations (≤ 1e-15) |
| HB | common-content: changing **only** the base annotation (2+1 base ref vs base alt; separate tables, different state order) leaves the marginal unchanged (≤ 1e-9). The original prior changes (contrast, ≠ 0) |
| HC | renaming ref/alt (World.ref ↔ World.alt) with probabilities transformed consistently leaves the result unchanged. Common-content: q → 1 − q. Original: base complemented (≤ 1e-9) |
| HD | with a synonymous site, merging labels into unique diplotypes (prior mass summed) and splitting each label in two (half mass each) preserve total prior mass (≤ 1e-15), the marginal (≤ 1e-12) and the KL, for both priors |
| HE | a likelihood constant across 1,024 states returns that constant, with no −log(state count) penalty (≤ 1e-12), for both priors |
| HF | table marginals equal an independent exhaustive oracle, for both priors (≤ 1e-6). The oracle enumerates absolute states, scores each state alone and computes priors in-loop. Scaffolds: 1+1 base ref, 2+1 base alt |
| HG | the KL identity holds (≤ 1e-9), for both priors |

## Runs (frozen: the four existing combinations only)

| unit | fixture | seed | truth | catalogue base content |
|---|---|---|---|---|
| run 1 | F2a missing variants, allocation invisible (600 bp) | 23 | 2+2, all alt | balanced ref, unbalanced alt |
| run 2 | F2b missing variants, allocation informative (200 bp) | 11 | 2+2, all alt | balanced ref, unbalanced alt |
| run 3 | F4 adjacency, true 2+2 (200 bp) | 23 | 2+2, ref | all ref |
| run 4 | F4 adjacency, true 3+1 (200 bp) | 23 | 3+1, ref | all ref |

## Expectations (declared before execution; not pass rules)

**F4: mechanically identical under both priors.**
- Every F4 scaffold has base = ref, so "flip from base" and "alt" coincide state by state: the two priors are identical
  arrays.
- Preservation of the F4 allocation controls is therefore **mechanical**, not a test of the control's effect on
  allocation evidence.
- M6-C recorded both runs as CORRECT (0.997 / 0.9997 for 2+2; 0.997 for 3+1).

**F2: only the base-alt scaffolds' priors change.**
- The balanced 2+2 scaffold (base ref) keeps its prior. The unbalanced (base alt) scaffolds' priors change.
- Under the control, every total-4 allocation pays the same content cost for an all-alt state.
- Expected outcomes:
  - **F2a:** the total-4 allocations become nearly tied, since allocation is invisible at 600 bp.
  - **F2b:** the outcome depends on whether the likelihood evidence for 2+2 over 2+3 exceeds what the control's per-copy
    content costs leave. M6-C's best-state terms suggested about 6 nats; that is a best-state statement, not a marginal
    decomposition.

**What the summary checks** (`summary.py`, frozen):
- C1: recomputation from the state table equals the executable's table (≤ 1e-6).
- C2: the original-prior marginals **reproduce M6-C's recorded marginals** for the same fixture and seed (≤ 1e-6).
- C3: the KL identity holds (≤ 1e-9).

It also reports:
- allocation posteriors and class (as M6-C);
- total-copy posteriors and the dosage class;
- per-total evidence;
- the truth allocation's posterior within its own total;
- for the truth against every other allocation: d log m = d E_q log L − d KL under each prior, and the centring effect
  (d_common − d_original).

Ties are read from the log marginals only. Pooled-content equality is never assumed to mean observation-model equivalence.

## Decision rule (declared)

- **Improvement** in either F2 combination (truth or truth total moves to ≥ 0.95 under the control, or the F2 error
  changes from a wrong structure to retained uncertainty), **with F4 preserved:** this justifies investigating a
  panel-grounded generative prior. It does not justify adopting this control.
- **No improvement:** this argues for parking the nuisance model before any expansion.

## Bounds

- 300 s counting every attempt; 4 threads; existing envelope.
- **Declared timeouts:** gates 60 s; F2a 50 s; F2b and each F4 run 30 s; summary 30 s.
- **Admission rule:** gates and runs are admitted only if spent + timeout + 30 s (summary reserve) ≤ 300.
- Runs execute only after the gates pass. Builds are not counted. Failed units stay recorded.
