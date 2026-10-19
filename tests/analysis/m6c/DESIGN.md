# M6-C: nuisance-marginalised structural likelihood, synthetic feasibility only (frozen before any unit)

**Ancestry.** Branch `m6-c-nuisance-marginal` from `b318c34` (M6-B), plus the M6-B closure.

**Scope.** One diagnostic executable, `tests/nuisance_structural.cpp`.
- **Not in scope:** caller integration, public flags, real-locus fitting, variant discovery, adaptive editing, indels and
  reconstructed output.
- **Count comparator:** diagnostic only and unused here.

## Model (frozen)

  P(reads | G) = Σ_θ P(θ | G) · P(reads | G, θ)

- **G (structural scaffold):**
  - haplotype A = F_A + n_A copies of a unit + G_A;
  - haplotype B = F_B + n_B copies + G_B;
  - unique random 1,200-bp flanks, so allocation is meaningful;
  - each copy carries a base allele at two predeclared biallelic sites, at unit positions ⌊L/3⌋ and ⌊2L/3⌋, with alt =
    the next base in ACGT.
- **θ (nuisance state):** one allele for every (copy, site), **shared by the whole read set**. There are 2^(2(n_A+n_B))
  states, at most 1,024, all enumerated; nothing is truncated.
- **P(θ | G):** an independent flip from the scaffold's base allele at each (copy, site), with **p = 0.05**, declared and
  normalised. It is not tuned.
- **P(reads | G, θ):** the existing paired-fragment likelihood (`score_diplotypes`) of the complete diplotype, with
  placement alternatives summed, exposure, and the outlier/background mixture.
  - Distinct sequences are interned; equivalent states keep their summed prior mass.
- **Marginalisation:** θ is summed by long-double log-sum-exp; the maximum is used only as a negative control.
- **Not applied:** no count score and no copying prior on structures. The structural posterior uses uniform weights over
  the declared catalogue.
- **Library:** insert 350 ± 35 bp (±4 sd), substitution error 0.002, outlier weight 0.001, diploid depth 20
  (λ = 20/(4·100)), 100-bp reads in both orientations.

**Catalogue.** All allocations with total 3–5 and at least one copy per homologue: 1+2, 2+1, 1+3, 2+2, 3+1, 1+4, 2+3,
3+2, 4+1. Each scaffold's base content is declared by the fixture.

## Gates (`gates`)

| gate | checks |
|---|---|
| GA | P(θ \| G) sums to 1 (≤ 1e-15) |
| GB | p = 0 reproduces fixed-sequence scoring exactly |
| GC | the batched, interned implementation equals an independent exhaustive oracle, where each θ is scored alone with no deduplication (≤ 1e-6) |
| GD | negative control: the best-θ value differs from the marginal and from the oracle |
| GE | with a synonymous site, the 16 θ labels collapse to 4 unique haplotypes and still equal the oracle |
| GF | whole-haplotype exchange and site relabelling are invariant (≤ 1e-9) |
| GG | a long unit with no spanning fragment ties 3+1 with 2+2 (≤ 1e-6); a short unit with the same pooled content does **not** tie, so pooled-count equality is not identity |

## Fixtures (frozen; seeds 11, 23, 37)

| fixture | unit | truth | catalogue base content | expectation (declared, not a pass rule) |
|---|---|---|---|---|
| F1 control | 600 bp | 2+2, all ref | all ref | 2+2 under both fixed and marginal |
| F2a missing variants, allocation invisible | 600 bp | 2+2, every copy alt at both sites | balanced = ref, unbalanced = alt | fixed favours an unbalanced alt scaffold. Marginal: allocation unobservable, so the preference is expected to be prior-driven (8 flips cost about 21 nats). |
| F2b missing variants, allocation informative | 200 bp | as F2a | as F2a | marginal recovers 2+2 only if adjacency evidence exceeds the flip cost |
| F3 matched content, true 3+1 / true 2+2 | 1,500 bp (no fragment spans a copy) | 3+1 / 2+2, ref | all ref | ties across same-total allocations: uncertainty retained |
| F4 adjacency, true 3+1 / true 2+2 | 200 bp | 3+1 / 2+2, ref | all ref | correct allocation in both, with no balanced preference |
| F6 dosage | 600 bp | 2+3 | all ref | total 5 preferred; allocations within it unobservable |
| F7 contamination | 200 bp | 2+2 | all ref | plus fragments from a 4%-diverged 3-copy paralog array (5% of fragments) and unrelated sequence (2%) |

**Summary (`summary.py`, frozen rules).** Per run:
- **marginal and fixed classes:** CORRECT (P_truth ≥ 0.95), UNCERTAIN (truth in a tied top set), WRONG, or SPLIT;
- **attribution:** CONTENT_ACCOMMODATION, ALLOCATION_EVIDENCE, PRIOR_DRIVEN (the winner and runner-up have equal
  likelihood terms and differ only through the nuisance prior), UNRESOLVED or MISLEADING;
- the decomposition, runtime and memory.

Failures are reported as they are; there is no tuning, no added states and no fixture repair.

## Bounds

- 600 s counting every attempt, 4 threads.
- **Declared timeouts:** gates 120 s; each fixture run 60 s; summary 30 s.
- The fixtures run only after the gates pass. Builds are not counted.

## Amendment 1 (declared after `run.F1_control_2+2.11`, before any further unit; 30.9 s spent)

**Why.** The first fixture run took 20.9 s. Under the frozen fixture-major order, the 26 remaining runs, plus the summary,
could exhaust the 600-s budget before the last fixtures (F6, F7) ran.

**Change: order only.** The remaining runs follow manifest A1: every fixture at seed 11, then seed 23, then seed 37. The
units, seeds, timeouts, model, fixtures and classification are unchanged, and the completed unit stays recorded. Any
unit the budget refuses is recorded as NOT RUN.

## Amendment 2 (declared after seed 23 completed, before any seed-37 unit; 429.1 s spent)

**Why.** The plain budget rule would admit seed-37 fixture runs that leave less than the summary's declared 30 s.

**Change.** A fixture run is admitted only if spent + its timeout + 30 s ≤ 600. Refused runs are recorded NOT RUN. The
summary rule is unchanged, and so is everything else.

## Amendment 3 (declared after `summary` failed; 533.3 s spent)

**What failed.** `summary` exited 1. A structural posterior underflowed to 0.0, and the tie test took `math.log` of it.

**Fix.** Ties (within 1e-6) are now decided on the log marginals directly. The classification thresholds and attribution
rules are unchanged.

**Dry run.** The fixed script was run once on the recorded outputs (to the scratchpad, not counted as a unit).

**Units.** Manifest A3, `summary.a3` (30 s). `summary` stays recorded as failed.
