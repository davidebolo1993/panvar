# S2-K: S2-J validation closed; post-refinement resource boundary

- **Comparator v2:** `s2j_compare_v2.py`, with documented allowances for the origin-cache end bytes and the scheduling-dependent union prefilter counters. `compare_controls.py` runs negative and positive controls on relocated copies (Amendment 1).
- **Guard lifecycle:** repair-path outputs are prepared in memory under the armed watchdog, then claimed and published together (`flush_outputs`). Post-refinement footprint and allocator checkpoints are recorded. The test hooks `post_refinement_guard` and `post_refinement_outputs` exercise stops there.
- **6 GiB-envelope ACOT 30× validation:** COMPLETE, peak 4.53 GB, science identical to S2-J's completed run.

Results: `tests/results/s2k_post_refinement_20261019/README.md`.
