# S2-K: close S2-J validation and the post-refinement resource boundary

Branch `s2-k-post-refinement` from `3be807b`. S2-J's admission policy stays experimental and default-off.

## 1. Comparator v2 (`s2j_compare_v2.py`)
v1 (`tests/analysis/s2j/s2j_compare.py`) and its two failed results are unchanged. v2 adds exactly four status keys:
- **`repair_origin_cache_bytes_end`**: changing retention is the intended intervention.
- **`union_cached_prefilters_consumed`, `union_completed_blocks_cached_prefilters_consumed`, `union_prefilter_cache_misses`.** These belong to the union prefilter-decision cache. A decision is a pure function of the (read, window) contents, so scores cannot change. The hit/miss split, however, depends on scheduling:
  - two workers can look up the same key before either stores it; both miss and both run the kernel, and `remember_prefilter` keeps the first;
  - once the table is at `max_held`, which keys were admitted depends on store order.

  Hits count against `max_union_cached_prefilters_consumed`. Near that cap, the same inputs could therefore refuse in one execution and not in another: resource behaviour is not identical, only the science is. The caps are unchanged.

**Validation of the comparator** (`compare_controls.py`):
- S2-J's saved C4 and ACOT 5× outputs are recompared, without rerunning the cells.
- Negative controls, on copies, must fail: a score, a candidate, the trajectory, a filter, the prediction, an undeclared status value, an undeclared new key.
- Positive controls must pass: the two declared allowance kinds.

## 2. Post-refinement checkpoints
These are main-thread `note_footprint` calls, now also recording the allocator's in-use bytes (`unavailable` where unsupported). A checkpoint appears only when its stage runs:
- `post_refinement_start`;
- `pool_decision_and_summary`;
- `representation_guard_done` (includes competitor scoring);
- `library_sensitivity_done` (only when the audit runs);
- `outputs_prepared`, `outputs_published`;
- `sequence_report_published`;
- `status_written:<status>`.

## 3. Guard lifecycle
**Before:** `claim_run()` ran right after refinement and tie branching. That disarmed the watchdog before the rerank, the representation guard, the library-sensitivity audit and output preparation.

**Now:**
- the repair path's final outputs (tie tables, candidates, scores, equivalence set, library-sensitivity blocks, genotypes, block genotypes, FASTA) are prepared in memory (`defer_output`);
- `flush_outputs`, called by `stop_with` on every repair-path exit, claims the run and then writes them all by atomic replacement, followed by the sequence report and the status, as before;
- the claim and a watchdog stop exclude each other under the watch mutex, so nothing is created between a stop's cleanup and its exit;
- a stop before the claim removes every final output and keeps only what the existing contract keeps: the S2-A checkpoint, as PROVISIONAL.

**Honest description:** the footprint limit is a sampled guard (polled every 100 ms while the run is unclaimed), not an allocation-level ceiling. The help text now says so.

**Test hooks:** `test_interrupt_point post_refinement_guard` (a time stop after the guard) and `post_refinement_outputs` (a memory stop with outputs prepared but unpublished).

**Gate** (`test_post_refinement_guard`):
- a normal completion records the stages in order and publishes after preparation;
- both stops end as budget stops with the checkpoint PROVISIONAL, no COMPLETE status and no final output;
- normal completion and publication-failure behaviour are covered by the full e2e (S1-C lifecycle tests).

## 4. One ACOT 30× confirmation (conditional on every gate)
Sol's prospective **6 GiB-envelope validation**:
- `footprint_limit_bytes 6442450944`;
- S2-J's completed arm otherwise unchanged: policy, inputs, scientific settings, cache budgets, every work and consumption cap, the per-run timeout;
- compared scientifically with S2-J's completed run by `compare_6gib.py`.
