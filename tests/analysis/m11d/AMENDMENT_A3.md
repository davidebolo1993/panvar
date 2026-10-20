# M11-D amendment A3 (declared before further execution)

**State at declaration:**
- `outcomes.a2` and `fit.a2` exited 0. Spent: 1,821.6 of 3,600 s.
- The frozen comparison is final: `runs/fit.a2/performance.tsv`.

**Added:** a **reporting-only** unit, `report` (120 s), running `m11d_report.py`. Its tables come from `outcomes.a2` and
`fit.a2`:
- run summary (completeness; ambiguous-completed vs resource-incomplete);
- block summary by locus (excess, normalized excess, floor attainment, catastrophic, top rate, quintiles, QVfrac > 1, artefacts);
- every non-top or catastrophic block with its out-of-family B0/C1/C2 probabilities;
- the pair edits and lengths behind each QVfrac > 1;
- the in-sample **direction** of association of each frozen input;
- Q2/Q3 cross-tabs.

**Why the direction table:** `fit.a2` shows a pooled AUC for C1 below 0.5. The direction table is needed to report that
honestly. It is **descriptive only**: it selects no predictor and changes no model, label, cell or comparison.

**Disclosed A3 development:** 120 s (reading `runs.tsv` and the fit printout; writing the script).

**Admission:** 1,821.6 + 15 (A2, already counted) + 120 + 120 ≤ 3,600.
