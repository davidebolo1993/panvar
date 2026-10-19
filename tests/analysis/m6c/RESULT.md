# M6-C result: nuisance-marginalised structural likelihood, synthetic feasibility

**Run.** Branch `m6-c-nuisance-marginal`: freeze `b5ad696`, plus three amendments:
- A1: seed-major order;
- A2: the summary's time reserved;
- A3: a summary log(0) defect.

Results in `tests/results/m6c_nuisance_marginal_20261019/`.

**Budget.** 533.4 of 600 s. 22 of 27 fixture runs completed. Five seed-37 runs were refused by the budget rule (F3 true
2+2, F4 ×2, F6, F7) and are recorded NOT RUN.

**Scope.** Synthetic only. No caller, flag, real-locus fitting or reconstruction.

## M6-B closure (appended to its record)

The M6-B findings stand, qualified to the tested dictionary, objective and catalogue:
- counts carry useful dosage information in the CYP2D6 example;
- C4's available content biases both tested models;
- matched allocations can be count-indistinguishable;
- paired adjacency can carry information the counts lack.

This does not show that every loss or longer-range model must fail. The count comparator stays diagnostic.

## Gates (all pass)

- **Prior mass:** sums to 1 (≤ 1e-15).
- **Disabled nuisance:** p = 0 reproduces fixed-sequence scoring exactly.
- **Oracle:** the batched, interned marginal equals the independent exhaustive oracle (each θ scored alone): difference 0.
- **Negative control:** maximisation differs from the marginal by 3.15 nats, so it is detected.
- **Equivalent states:** with a synonymous site, 16 labels collapse to 4 unique haplotypes, and the oracle is still
  equal.
- **Invariance:** whole-haplotype exchange and site relabelling are invariant (exactly 0).
- **Ties:**
  - a long unit (no spanning fragment) ties 3+1 with 2+2 exactly;
  - a short unit with the **same pooled content** separates them by 2,873.7 nats, so pooled-count equality is not
    identity.

## Fixture outcomes

Structural posterior over the 9-scaffold catalogue, from the marginal likelihood, with uniform structural weights and no
structural prior. Nuisance states per scaffold: 64 (total 3), 256 (total 4), 1,024 (total 5), enumerated completely.

| fixture (unit) | seeds | marginal outcome | fixed-sequence outcome | reading |
|---|---|---|---|---|
| F1 control 2+2 (600 bp) | 11, 23, 37 | seed 11: totals 4 vs 5 split (P_truth 0.11; total 5 ahead by 0.36 nats); seeds 23, 37: **1+3 / 2+2 / 3+1 tied** (about 0.33 each; total 4 by 4.5 nats) | the same | Allocation is **unobservable** with 600-bp copies; uncertainty retained. Dosage is weakly identified at depth 20: one seed leans to 5 copies. The declared "control" expectation (2+2 unique) was wrong for this unit length. |
| F2a missing variants, allocation invisible (600 bp) | 11, 23, 37 | P_truth ≈ 1e-11. Winner: an alt-content unbalanced scaffold (total 4, or total 5 on seed 11) | the same (fixed 2+2 is 694–739 nats behind) | **Content accommodation works for the likelihood:** marginal 2+2 has exactly the winner's fragment and exposure terms. The ranking is then **prior-driven:** 2+2 pays log(0.05⁸) = −23.97 for 8 corrections; the winners need none. |
| F2b missing variants, allocation informative (200 bp) | 11, 23, 37 | **2+3 / 3+2 wrong** (P ≈ 0.5 each; P_truth about 1e-8 to 2e-7) | 2+3 / 3+2 | Adjacency excludes 3+1 and 1+4 (−2,220 to −4,423). The corrected 2+2 **beats 2+3 on likelihood by about 6.3 nats** (exposure +10.0, fragment −3.6), but the 23.5-nat nuisance-prior difference reverses it. **Prior-driven and misleading.** |
| F3 matched content, true 3+1 (1,500 bp) | 11, 23, 37 | 1+3 / 2+2 / 3+1 **tied** (0.33 each) | the same | Uncertainty retained, as declared. No balanced preference. |
| F3 matched content, true 2+2 (1,500 bp) | 11, 23 | tied (0.32–0.33) | the same | as above |
| F4 adjacency, true 2+2 (200 bp) | 11, 23 | **CORRECT** (0.997, 1.000) | CORRECT | Allocation evidence from flank-to-flank spanning fragments. |
| F4 adjacency, true 3+1 (200 bp) | 11, 23 | seed 11: 3+1 top at **0.82** (SPLIT); seed 23: **CORRECT** (0.997) | the same | Allocation evidence holds for both truths. There is no balanced bias; one seed is under-confident. |
| F6 dosage, true 2+3 (600 bp) | 11, 23 | **total 5 identified** (four allocations at 0.25 each) | the same | Dosage is evidenced; allocation unobservable, uncertainty retained. |
| F7 contamination, 2+2 (200 bp; 5% paralog + 2% off-target) | 11, 23 | **CORRECT** (0.997, 0.999) | CORRECT | Robust to this small contamination. |

**Note on the frozen attribution rule.** `summary.py` compared each winner with its own runner-up, which is the
allocation's mirror-image tie, so it labelled the F2 runs SPLIT. The decompositions above, which compare the winner with
the truth scaffold, show the actual driver: the nuisance prior.

**Prior structure.** With no flips the prior is still 0.95^(2n). Each extra copy costs 0.10 nats, a small declared
preference for fewer copies.

**Cost.** Each catalogue of nine scaffolds took:
| unit | time | states scored |
|---|---|---|
| 200 bp | 8–12 s | up to 1,024 per scaffold, with fresh placements |
| 600 bp | 20–28 s | (as above) |
| 1,500 bp | 47–55 s | (as above) |

Peak RSS was under 0.03 GiB.

## Feasibility decision

**Mathematically feasible; scientifically not yet sufficient.**
- **Improvement from accommodating missing content: real but incomplete.** Marginalising a shared θ absorbs the
  missing-variant likelihood penalty exactly (F2: fixed 2+2 was about 700 nats behind; marginal 2+2 has the
  correct-content likelihood). The structural decision then rests on the nuisance prior.
- **The prior can dominate the structure.** In both missing-variant fixtures the declared prior's cost for the
  corrections (about 24 nats for 8 site-copies) exceeded the available structural evidence (about 6 nats of dosage
  evidence in F2b). The result is a confident **wrong** structure. Content-matched base scaffolds win whenever allocation
  is unobservable. This is the central risk for real C4, where many more sites per copy would need accommodating.
- **Evidence genuinely supporting allocation** comes only from spanning adjacency (F4; F2b's exclusion of 3+1 and 1+4).
  Where none exists, ties are retained honestly (F3, F6, F1). There is no systematic preference for balanced allocation.
- **Dosage** is evidenced (F6), but only weakly at depth 20 with 600-bp copies (F1 seed 11).
- **Unresolved or misleading:** F2a and F2b are misleading (prior-driven); F1, F3 and F6 are unresolved (correctly
  uncertain); F4 true 3+1 seed 11 is under-confident (0.82).

## How this differs from M2

M2 discovered residual edits **from the reads**: data-selected sites, proposal and acceptance, decoy competition, the
best joint hypotheses. That is effectively optimising edits per hypothesis. M6-C instead:
- uses sites **declared before any read**, identical for every structure;
- uses a **normalised** prior;
- **sums** every nuisance state for every structure, so all structures are compared on the same nuisance terms with no
  per-candidate editing.

## What real loci would still require (not attempted)

1. **Candidate-independent nuisance sites:** positions polymorphic among the panel's unit copies, from a unit alignment,
   declared before reads. On C4 this is hundreds to thousands of sites per module, plus indels.
2. **A prior grounded in the panel:** per-site allele frequencies, likely correlated across sites and copies, rather
   than independent flips. Its specification decides structural outcomes (F2), so it would need its own validation.
3. **Factorised computation.** 2^(sites × copies) is intractable. A copy- and site-local decomposition needs the fragment
   likelihood split into site-local terms that preserve exact placement sums.
4. **Copy-permutation equivalence and identifiability analysis,** to report which allocations are unobservable instead of
   letting the prior choose.

**Synthetic success here does not solve real C4.** The F2 failures show that, even in a two-site world, the nuisance
prior can override weak structural evidence.

## Qualifications appended at M6-D (2026-10-19; Sol's review; the record above is not rewritten)

The arithmetic above is accepted. **The model is not approved for caller integration.**

1. **Best-state columns are not a decomposition of the marginal.** `best_state_fragment`, `best_state_exposure` and
   `best_state_log_prior` describe the single nuisance state with the highest posterior. They do not decompose the
   marginal.
   - The F2 readings above ("8 corrections cost −23.97", "2+2 beats 2+3 on likelihood by about 6.3 nats", "the
     decompositions show the actual driver: the nuisance prior") are best-state statements.
   - The exact marginal decomposition is log m = E_q[log L] − KL(q ‖ prior), reported in M6-D.
2. **The no-change state alone does not establish a copy-number preference.** "0.95^(2n), 0.10 nats per extra copy"
   concerns that one state. Each copy's nuisance prior sums to one, so this is not a marginal copy-number preference by
   itself.
3. **Uniform allocation weights induce a prior on total copy number.** Uniform weights over the nine allocations give
   prior masses **2/9, 3/9 and 4/9** to totals 3, 4 and 5. Total 5 carries twice the prior mass of total 3. This bears
   directly on the runs that leaned toward total 5: F1 seed 11 and F2a seed 11.
4. **Exact factorisation is unproved.** "Mathematically feasible" refers to complete enumeration in this tiny state space.
   Exact site-local factorisation of the paired-fragment likelihood remains unproved.

**Unchanged:**
- F2a and F2b remain failures. F2a was prior-driven (truth posterior about 1e-11); F2b chose the wrong dosage (2+3 or 3+2).
- F1, F3 and F6 remain unresolved.
- F4 true 3+1, seed 11, remains under-confident (0.82).
