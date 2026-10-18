# S2-D: bounded emission-cache lifetime and narrow progressive-checkpoint integration (experimental)

## Bounded emission reuse within one batch

`ContentIndexBatch::set_persistent_emission_cache(true)` turns it on. It is off by default; off is S2-C's per-call cache.

**How it works.** After each selective discovery, the emission cache is detached into a batch-owned handoff. The next discovery adopts it.

**Safeguards:**
- **Budget:** the same byte budget applies; it is not enlarged.
- **Identity:** the library, band and recursion backend must match (the existing I2-K check). Keys hold the read and window bytes.
- **Lifetime:** the handoff is destroyed with the batch, so it never crosses samples. It is exclusive with the adaptive handoff.
- **Exceptions:** a guard ends the global cache, so it is never left armed.
- **Refusals:** a refused or failed call poisons the batch. Preconditions are checked first and do not poison.

## Progressive scheduling in the caller

The profile key `union_progressive_schedule experimental_v1` turns it on. It needs `early_incumbent experimental_v1` and the exact indexed route.

**What it does.** The indexed union pass runs as a schedule on one selective batch, keeping every candidate:
1. the S2-A pair;
2. the next 31 pairs by proposal score;
3. the rest.

**The checkpoint.** After the 32-pair prefix, its best fully evaluated pair is checkpointed as `union_prefix_best`. The prefix's certified ties and the pair's projection ambiguity are recorded with it.

**What stays the same.** The pool's scores feed every downstream consumer unchanged, and a refusal goes through the existing union-cap path. S2-A's union check compares against the stored early candidate. The S2-A pair's re-evaluation inside the union index is reported as repeated work (`union_progressive_repeated_early_evaluation_*`).

## Tools

- **Gates:** `genotype_progressive_kernel` (S2-D cases) and e2e `test_progressive_schedule`.
- **Runner:** `s2d_run.py`.
- **Comparison:** `s2d_compare.py`, complete-run output identity.
