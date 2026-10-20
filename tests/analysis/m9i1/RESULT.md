# M9-I.1 result: contract v2 hardens missing, empty, zero-mass and invalid cases and normalisation; reapplied to the 34 M9-I cases it changes no scientific conclusion

**Run.** Branch `m9-i1-compound-contract`, freeze `42a311a`, with no amendments. Results are in
`tests/results/m9i1_compound_contract_20261020/`.
- **Execution:** 99.3 of 600 s, including 2 s pre-freeze.
- **Read-only:** M9-I's code is reused unchanged. M9-I's tables are untouched, and its corrections are appended to its
  RESULT.

## Corrected contract (v2)

**Target:** the unordered pair of concatenated interval sequences.
- It is exact-sequence, **including within-interval phase**, and it is not structural confidence.
- Over the explored pool, every recorded labelling is projected and each group's mass counted once.

**Inputs:** finite, non-negative masses, otherwise **INVALID**. The original total is recorded, masses are normalised by
it (`fsum`), and the discrepancy is reported.

**Empty pool:** **UNASSESSED**, never zero error.

**Missing projections:** the mass enters both upper bounds and neither lower bound. It is counted, and the status becomes
`ASSESSED_WITH_MISSING_PROJECTIONS`.

**Bounds:** direct sums.

| bound | groups counted |
|---|---|
| support_lower | all projections equal T |
| support_upper | any projection equals T, plus missing groups |
| error_lower | labellings present, none equals T |
| error_upper | any projection differs, plus missing groups |

**Ambiguity:** a group with several projections is reported as ambiguous, and its mass is never split.

**Competitor, decided from projections:**

| status | meaning |
|---|---|
| NO_EVALUATED_ALTERNATIVE | every recorded projection equals T |
| EVALUATED_ALTERNATIVE | a differing projection with positive mass |
| EVALUATED_ALTERNATIVE_ZERO_MASS | differing projections exist, with mass exactly 0 |
| UNKNOWN | missing groups, and no recorded difference |

NO_EVALUATED_ALTERNATIVE is a statement about the explored pool, never evidence of correctness.

**Scope:**
- on saved tables these are approximate masses reconstructed from 12-digit printed values; groups below double range
  print as 0;
- they are not certified bounds on the full-precision posterior;
- they do not cover unsearched alternatives.

## Gates (11/11 PASS)

| gate | checks |
|---|---|
| M1 | missing group: both uppers, neither lower; counted; competitor UNKNOWN |
| M2 | empty or zero-mass pool is UNASSESSED |
| M3 | NaN, infinite and negative masses are INVALID |
| M4 | totals of 1.0000003 and 0.9999997 are normalised, with the discrepancy recorded |
| M5 | a 1e-300 error mass is retained exactly |
| M6 | a zero-mass competitor is EVALUATED_ALTERNATIVE_ZERO_MASS |
| M6b | a wholly supporting pool is NO_EVALUATED_ALTERNATIVE |
| M7 | **physically realisable mixed projection**: "A" + "CG" and "AC" + "G" spell the identical whole haplotype (independent spelling), project differently onto the single-block interval, and give support [0, 0.8] and error [0.2, 1.0] with the mass unsplit |
| M7b | on the enclosing interval the same group is one projection |
| M8 | **phase counterexample**: truth AC/GT against prediction AT/GC; single-position unphased genotypes match, the compound pair does not |

## Reapplication to the 34 saved M9-I cases

| category | finding |
|---|---|
| **Scientific conclusions that change** | **none**. All 34 are ASSESSED, with 0 missing groups and 0 ambiguous groups. Every competitor status is EVALUATED_ALTERNATIVE with positive mass (M9-I's "present"). Compound exactness and distances are unchanged, being prediction facts. Single-block controls still reproduce the saved `sequence_blocks` bounds 32/32 |
| **Numerical changes from normalisation only** | |total − 1| ≤ 5.0e-10 per run; every bound changes by ≤ 5.0e-10 relative |
| **Semantics exercised only by fixtures** | missing projections (M1), UNASSESSED and INVALID pools (M2, M3), zero-mass-only competitors (M6), a realisable mixed projection (M7), and the phase case (M8). Note that every run contains 10⁴–9 × 10⁴ zero-mass groups, but always alongside positive-mass disagreements |

**Interpretation corrections** (appended to M9-I):
- "compound distance never exceeds block sums" is scoped to the observed cases; M8 is a counterexample;
- compound targets include within-interval phase;
- "honest masses" is replaced by descriptive wording;
- exact-sequence confidence stays distinct from structural confidence.

## Recommendation on additive integration

**Ready for additive integration, under contract v2.** Add one annotation row per declared interval beside
`sequence_blocks`:
- **identity:** interval id and source (model inventory and/or representation interval);
- **target:** compound target digests;
- **status:** status and competitor;
- **bounds:** support and error bounds;
- **counts:** ambiguous, missing and zero-mass group counts;
- **accounting:** normalisation discrepancy;
- **scope:** a scope string.

**Conditions:**
1. **Compute from in-memory doubles, not printed tables.** That removes the reconstruction caveat. Underflowed groups
   still contribute exactly 0 and must be reported as zero-mass, not omitted.
2. **Never change predictions, filters or block rows.** Block-level rows stay as provenance and diagnostics.
3. **Declare intervals conservatively** from the model inventory, merged with the run's representation intervals, and
   state that they are not minimal.
4. **Gate the C++ port** against this Python reference on the 34 cases and on the M1–M8 fixtures.
5. **Treat any user-facing confidence decision as a separate step,** needing held-out calibration in the compound unit.
   No threshold follows from these inspected cases.
