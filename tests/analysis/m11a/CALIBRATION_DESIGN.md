# Revised block-confidence calibration design (proposal; no genotyping, simulation, fitting or thresholds)

**Requirement.** **One user-facing assessment per reported block genotype.** It may use several inputs, combined inside one
pre-declared calibration. That is different from asking users to combine several independent production filters. The
existing headline filters are not extended.

The primary outcome stays the floor-adjusted per-block excess (bp and diploid-length normalized). The error-mass-only
calibration campaign is **not** frozen.

## 1. Candidate additional truth-free inputs (at most three; from saved development evidence only; validation untouched)

The evidence below comes from the M11-B development observations (`block_observations.tsv`) and the current-schema runs'
`dosage_blocks.tsv`. It is descriptive, with donors as units, and no test is made.

### Q1. Read support of the block genotype

| | |
|---|---|
| **definition** | `genotype_support` ∈ {READS, PRIOR_ONLY} in `<run>.block_genotypes.tsv`. PRIOR_ONLY means no fragment informs the block's genotype (`quality_reason no_read_support`). The same information is `read_support_status` in `dosage_blocks.tsv` |
| **availability** | every run, in both the legacy and the current schema |
| **mechanism** | copying-prior effects (E2, E4); prior-only confident calls |
| **evidence for** | floor attainment is lower for PRIOR_ONLY in **every** stratum: M1-B 0.36 vs 0.75; M9-C 0.63 vs 0.91; M9-G 0.66 vs 0.90; M9-H 0.61 vs 0.81; C4 0.74 vs 0.95; CYP2D6 0.66 vs 0.87. Extreme-confidence errors that are PRIOR_ONLY: NA19909 sv:13 (error mass 1.1e-38, 27 kb); CYP2D6 sv:9 5× (12 kb); 27 of 44 C4 catastrophic blocks; all three ACOT sv:9 ≈ 475 bp errors |
| **evidence against** | C4 NA19650 sv:7 (65 kb) is READS-supported. Most PRIOR_ONLY blocks are correct (many context blocks) |
| **calling-time cost** | none (already emitted) |

### Q2. Block census route

| | |
|---|---|
| **definition** | `census_route` ∈ {repeated_unit, pair, not_sensitive} in `<run>.dosage_blocks.tsv`. It is model-derived from node multiplicity and repeated bases (`run_census`); it is a block property, **not** sample evidence |
| **availability** | current-schema runs **when the dosage audit succeeds**. It is absent when the audit fails (CYP2D6 M1-B 5×: "unit library: alleles disagree on their flanks"; M9-C d1: "KMAX > 3"). It is absent from the legacy C4 production runs |
| **mechanism** | copy-number and allocation errors whose true structure is absent from, or loses within, the pool |
| **evidence for** | both catastrophic development sites with a route are `repeated_unit`: C4 sv:7 (NA19650, and every high-GQ C4 campaign catastrophic block is at sv:7) and ACOT sv:13 (NA19909) |
| **evidence against** | a class label. It can lower the **average** confidence of repeated-unit blocks, but cannot separate good from bad calls within that class. A dosage-audit PASS proves nothing: NA19650 sv:7 is `CONSISTENT` / `NOT_TESTED` (total dosage right, allocation wrong). CYP2D6 sv:9's route is unavailable in the saved runs |
| **calling-time cost** | negligible. The census already runs in the default audit. Making it available when later audit stages fail is an output change, not new computation |

### Q3. Proposal-emission tie at the block (the JP1.1 eligibility condition)

| | |
|---|---|
| **definition** | the block is eligible when **every pair emission ties and more than 32 pairs exist** (`junction_proposals.cpp`): the read-marker proposal cannot rank the block's alleles |
| **availability** | computed in every repair when JP1.1 is on. Emitted only as a **count** (`junction_eligible_blocks`); the block identity is not written |
| **mechanism** | search omissions at marker-uninformative blocks (E1, E3, E6) |
| **evidence for** | in every JP-on ACOT development run (M9-G 5, M9-H 6) exactly one block is eligible, and M9-E identified it as sv:5, where 3 of the 6 documented confident errors were search gaps. JP1.1 reduced, but did not remove, the sv:5 residuals (M9-H) |
| **evidence against** | C4 and CYP2D6 have 0 eligible blocks, so there is no information there. Many sv:5 calls are correct. The evidence comes from one locus |
| **calling-time cost** | none when JP1.1 is on. Emitting a per-block flag is an **output change** of an already-computed condition. Nothing is emitted today |

**Not proposed:** another structural metric, internal whole-haplotype alignment, or completion status as evidence. A
completed bounded search does not prove candidate adequacy.

## 2. Revised calibration design (for approval; to be frozen before any development outcome is examined)

| element | proposal |
|---|---|
| **unit** | one reported block genotype; donors (families) are the independent units for every rate and interval |
| **outcome** | the floor-adjusted excess (bp, normalized). The "acceptable" outcome is Davide's definition: floor attainment, or an approved excess allowance. Catastrophic and floor attainment are reported descriptively |
| **assessability gate (not evidence)** | completed ranking, not provisional, no cap / time / footprint stop, representation guard complete. Otherwise `NOT_ASSESSED` with the reason; the prediction is kept |
| **inputs** | `error_mass_upper` (per block; continuous), Q1 read support, Q2 census route (with an explicit `route_unavailable` level), Q3 emission tie (with an explicit `not_emitted` level until it is emitted). Locus-level support is reported separately, never broadcast |
| **combination** | **one** calibrated estimate per block: the probability of an acceptable outcome. A monotone map of error mass is fitted **within pre-declared strata** of (Q1 × Q2 × Q3), with a pre-declared merge order for strata with too few development donors (Q3, then Q2, then Q1). A stratum below the declared minimum number of donors after merging stays `NOT_ASSESSED (insufficient development evidence)` |
| **development data** | the unified M10-A configuration re-run on exposed families with **existing reads**: CYP2D6 6 pilot donors, C4 12 campaign families by frozen hash, ACOT 12 previously genotyped families. Legacy outputs are descriptive only and never pooled across configurations |
| **predictor selection** | fixed **now**, from development evidence and mechanisms. The validation outcomes are never used to select or drop inputs |
| **validation** | ACOT's 18 never-genotyped families (prior-score-exposed in M8-B, **not** pristine). Truth is sealed until the procedure (strata, merge rule, minimum donors, map form) is frozen. C4 and CYP2D6 have no available validation, which is stated |
| **reported** | per stratum and pooled: calibration (predicted vs observed acceptable rate, by donor-level bootstrap), the confident-error rate at the decision Davide approves, and catastrophic confident errors listed individually. The attribution classes of §Qualifications (search omission, prior, likelihood/model, mixed, unattributed) are kept as columns |

## 3. Limitations

- **The documented failure classes are only partly separable by these inputs.**
  - The C4 sv:7 allocation error is READS-supported, repeated-unit and not eligible for JP1.1. Its inputs resemble correct
    repeated-unit calls, so calibration can only **lower the average confidence** of that class.
  - The unattributed errors (ACOT sv:13 and sv:9; ten C4 sv:7 donors) may share inputs with correct calls.
- **Availability gaps:**
  - Q2 is missing when the dosage audit fails.
  - Q3 needs a per-block output, not present in any saved run.

  Both would be reported as explicit levels, never imputed.
- **Small evidence:** few development donors per locus, correlated blocks within a donor, and one locus (ACOT) for the
  validation and for Q3.
- **Rounding:** printed error masses of 0 and posteriors of 1 lose order. The map must use GQ / log-scale values, or a
  higher-precision output.
- **Synonymous-boundary artefacts** stay visible. How they count is Davide's decision.
- **Configuration drift:** calibration applies only to the configuration it was developed on (the unified M10-A
  configuration).

## 4. Decisions needed (block-level)

1. **The acceptable block outcome.** Floor attainment, or an excess allowance (bp and/or normalized; possibly different for
   context and SV blocks).
2. **The acceptable confident-error rate,** and whether any confident catastrophic error is tolerable (catastrophic
   definition: reuse the 1,000 bp report default?).
3. **The input set:** Q1–Q3 as above, or a subset. Permission for the two **output-only** changes: Q2 available on audit
   failure; Q3 per-block flag.
4. **The minimum development donors per stratum**, and the merge order.
5. **Synonymous-artefact treatment.**
6. **Validation acceptance:** ACOT 18, prior-score-exposed; C4 / CYP2D6 development-only. Read simulation for validation;
   5× first.

---
**Superseded in part (appended by M11-C).**
- **Q1's definition:** see the correction in `RESULT_M11B.md`. PRIOR_ONLY means some fragment contender in the ranked pool
  disagrees at the block.
- **§2's stratified (Q1 × Q2 × Q3) scheme:** replaced by M11-C's two pre-declared candidates, C1 (`error_mass_upper`) and
  C2 (+ Q1). Q2 and Q3 stay diagnostic only, as default-off outputs (`block_diagnostics`). See
  `tests/analysis/m11c/RESULT.md` §4.
