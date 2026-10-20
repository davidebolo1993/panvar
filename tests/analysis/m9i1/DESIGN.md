# M9-I.1: hardened compound-annotation contract (read-only)

**M9-I stays the historical record.** Its code is reused unchanged: case list, dumps and projector are loaded by executing
`m9i.py` in a no-op mode. Its tables are untouched, and corrections are appended to its RESULT.

**Not done:** no caller rerun, new family, production change or default change.

**Results:** `tests/results/m9i1_compound_contract_20261020/`.

## Contract v2 (`m9i1.py bounds_v2`)

**Target T:** the unordered pair of concatenated interval sequences.
- It includes **within-interval phase**: for example, truth AC/GT against a prediction AT/GC match at every single
  position unphased, but not as a compound pair.
- It is an **exact-sequence** target. It tolerates no SNP or indel and is not structural confidence.

**Inputs:**
- every group mass must be finite and ≥ 0, otherwise **INVALID**, with no bounds;
- the original total is recorded, masses are **normalised by it** (`math.fsum`), and the discrepancy (total − 1) is
  reported.

**Empty pool** (no group, or total mass 0): **UNASSESSED**, never zero error.

**Missing projections** (a group with mass but no recorded labelling) — the policy chosen: **conservative**.
- Its mass enters both `support_upper` and `error_upper` and **neither** lower bound.
- It is counted (`missing_groups`, `missing_mass`).
- The status becomes ASSESSED_WITH_MISSING_PROJECTIONS.

**Bounds,** each a direct `fsum` of normalised masses and never a subtraction from one:

| bound | groups counted |
|---|---|
| support_lower | every recorded projection equals T |
| support_upper | some projection equals T, plus missing groups |
| error_lower | labellings present, and none equals T |
| error_upper | some projection differs, plus missing groups |

**Ambiguity:** a group with several projections is ambiguous, and its mass is never split or assigned.

**Competitor status, decided from projections, not mass:**

| status | meaning |
|---|---|
| NO_EVALUATED_ALTERNATIVE | every recorded projection equals T and none is missing |
| EVALUATED_ALTERNATIVE | a differing projection with positive mass |
| EVALUATED_ALTERNATIVE_ZERO_MASS | differing projections exist, and their mass is exactly 0 |
| UNKNOWN | no differing projection recorded, but some groups are missing |

NO_EVALUATED_ALTERNATIVE describes the explored pool only. It is **never** evidence of correctness.

**Scope:**
- saved-table masses are printed to 12 significant digits, so these are approximate pool-conditional masses
  reconstructed from printed values;
- they are not certified bounds on the full-precision posterior;
- they cover explored candidates and recorded labellings only.

## Gates (`gates`)

| gate | checks |
|---|---|
| M1 | a missing group adds mass to both upper bounds and neither lower; it is counted; competitor UNKNOWN |
| M2 | an empty pool, or one with zero total mass, is UNASSESSED |
| M3 | NaN, infinite and negative masses are INVALID |
| M4 | totals of 1.0000003 and 0.9999997 are normalised, with the discrepancy recorded |
| M5 | a 1e-300 error mass is retained exactly |
| M6 | a zero-mass competitor is EVALUATED_ALTERNATIVE_ZERO_MASS |
| M6b | a wholly supporting pool is NO_EVALUATED_ALTERNATIVE |
| M7 | a **physically realisable** mixed projection: two labellings spell the same whole haplotype ("A" + "CG" = "AC" + "G", checked by independent spelling), and their compensating boundary gives different projections on the single-block interval; the mass is kept as an interval [0, 0.8] support / [0.2, 1.0] error |
| M7b | the same group on the enclosing interval is one projection |
| M8 | the phase counterexample |

**Disclosed development:** one development run of the gates before the freeze (0.14 s).

## Reapplication (`reapply`)

The 34 M9-I cases are recomputed under v2, and the single-block controls against the saved `sequence_blocks` bounds.

**Reported separately:**
- scientific changes: status, competitor and missing groups;
- numerical changes caused only by normalisation: discrepancy and maximum relative bound change;
- semantics exercised only by fixtures: missing, ambiguous and zero-mass cases that do not occur in the saved data.

## Budget

**600 s,** including 2 s pre-freeze.
