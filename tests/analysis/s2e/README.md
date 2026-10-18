# S2-E: exact fragment-blocked union scoring (experimental kernel, measurement only)

**What it is:** `score_diplotypes_indexed_blocked` builds one content index, then processes contiguous blocks of read pairs. Each block gets a fresh, one-shot batch using the existing discovery, pairing and reduction. Its mass columns are folded into every candidate's accumulators in the original fragment order, then released with the batch.

**What stays whole-run** (see `DECOMPOSITION.md`):
- exposure;
- background;
- the prefilter-cache read-length condition;
- the mass-cell cap.

**Caps.** Executed-work and consumption caps are cumulative (`ContentIndexBatch::set_cumulative_baseline`). Memory budgets apply per block.

**Tools:**
- the gate: `tests/genotype_blocked_union.cpp` (ctest `panvar_genotype_blocked_union`);
- the measurement: harness mode `s2c_progressive_harness blocked`, with a 4 GiB guard on bulk and blocked runs;
- the runner: `s2e_run.py`.

Results: `tests/results/s2e_fragment_blocked_20261018/README.md`.
