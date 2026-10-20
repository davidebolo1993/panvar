# M11-D amendment A1 (declared before further execution)

**State at declaration:**
- All 30 cells' call, eval and truth units completed (exit 0, no watchdog or envelope event), plus `inventory`.
- Spent: 1,763.0 of 3,600 s.

**Defects (analysis code only; no caller, configuration, input or comparison change):**
1. **`outcomes` (exit 1).** Cell `acot.NA21106` is a caller **no-call**:
   - `RESCORE_INCOMPLETE` / `divergent_tie`: two candidates are tied on total score and disagree at ctx:18;
   - `sequence_completion = completed_ranking`, `SEQUENCE_TIED`.

   Its evaluator rows carry `called_allele = NA`, and the QVfrac step indexed an `allele_NA` sequence. Such a block is now
   `predicted = False` (NO_PREDICTION): it is counted in coverage and gets no excess, QVfrac or confidence row.
   Artefact intervals require every block in the interval to be predicted.
2. **`fit` (exit 0, output invalid).** It ran after the failed outcomes unit, because its precondition did not require
   outcomes to succeed, and wrote an empty comparison. Now:
   - `fit` reads the outcomes unit with the same suffix and fails on a missing table;
   - it excludes blocks with `predicted ≠ True`;
   - `fit.a1` requires `outcomes.a1` exit 0.
3. **Operator loop.** My shell loop tested the watchdog verdict and not the exit code, so it did not stop at the failed
   unit. Later loops stop on a non-zero exit.

**Re-run units:** `outcomes.a1` (300 s) and `fit.a1` (120 s), declared timeouts unchanged. The failed `outcomes` unit and the
invalid `fit` unit stay recorded in the ledger. Admission: 1,763.0 + 300 + 120 ≤ 3,600.
