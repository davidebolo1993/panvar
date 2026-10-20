# M11-D amendment A2 (declared before further execution)

**State at declaration:**
- `outcomes.a1` exited 1 at its final evaluator-consistency check, after writing its tables. The check indexed
  `distance_matches_evaluator` on every truth-defined row, including the A1 NO_PREDICTION rows, which have no such column.
- `fit.a1` was NOT RUN (precondition).
- Spent: 1,786.9 of 3,600 s.

**Disclosed A2 development:** 15 s. I checked the structure of `outcomes.a1`'s partial `blocks.tsv`:
- row classes: 454 predicted, 19 NO_PREDICTION (all `acot.NA21106`, the caller no-call), 1 undefined-truth (`cyp.NA19650`);
- consistency columns: 0 distance or floor mismatches among the predicted rows.

No outcome value was used. Nothing in the frozen fitting procedure changes.

**Fix:** the consistency check runs on `predicted` rows only. **Re-run units:** `outcomes.a2` (300 s) and `fit.a2` (120 s).
Admission: 1,786.9 + 15 + 300 + 120 ≤ 3,600. The failed and NOT RUN units stay recorded.
