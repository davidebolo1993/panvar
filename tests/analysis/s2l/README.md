# S2-L: bounded cross-locus validation and performance-baseline closure

**Tools:**
- comparator v3 → v3.2 (`s2l_compare*.py`, categorised explicit allowances; Amendments 1–2);
- controls (`compare_controls_v3.py`);
- frozen matrix runner (`s2l_run.py`; Amendment 3 for completed ties);
- whole-process summary (`s2l_summary.py`).

**Outcome:**
- Control (ordinary union) vs combined (blocked/fine union + round-1 reuse-only admission) are scientifically identical on every completed pair: C4, CYP2D6 and ACOT at 5×; C4 and CYP2D6 at 30× (CYP2D6 a genuine tie).
- ACOT 30×: the control hits the 500 M origin cap; combined completes (103.8 s, peak 4.54 GB) and repeats identically.

Results: `tests/results/s2l_cross_locus_20261019/README.md`.
