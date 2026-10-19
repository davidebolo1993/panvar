# M3-D Amendment 1: cache work counters in the option-off identity check

Declared after `compare` (310.0 of 900 s spent, exit 1) and before any further execution.

## The failure

`compare` failed 3 of 36 identity checks, all on `.mosaic_status.tsv`. In every cell, the only non-volatile lines that
differ between the baseline binary and the M3-D binary with the option off are emission- and prefilter-cache **work
counters**:
- `emission_cache_dp_hits` / `dp_misses` / `executed_prefilters` / `executed_dps`;
- `scorer_emission_cache_arena_bytes`;
- `union_cached_prefilters_consumed`;
- `union_prefilter_cache_misses`.

Every other status line agrees, and so does every other output file: genotypes, block genotypes, scores, candidates,
equivalence set, FASTA and the rounds table's science columns. That is the other 33 checks.

**Why these could vary.** With 4 threads the emission cache is filled concurrently, so hit and miss counts can depend on
scheduling. The caps count executed work and are cache-dependent, as recorded earlier. That is a hypothesis to test, not
a classification to assume.

## Change

1. **`run.mosaic_2x.baseline_repeat`.** The **baseline binary** again, with identical inputs, profile and 4 threads, the
   cheapest cell.
2. **`repeat.check`.** Lists the non-volatile status differences between the two baseline runs.
3. **`compare.a1`.** The frozen comparator with exactly those seven keys treated as volatile
   (`M3D_CACHE_COUNTERS_VOLATILE=1`); nothing else changes.
   - **Interpretation:** if `repeat.check` shows the same keys varying between two baseline runs, they are
     scheduling-dependent and `compare.a1` is the identity verdict. If it does not, the difference is attributed to the
     M3-D binary and reported as a failure.
4. **Evaluation units** (`eval.a1.*`) run only after `compare.a1` freezes the selection.

`compare` stays recorded with exit 1.

**Unit timeouts:**

| unit | declared timeout |
|---|---|
| `run.mosaic_2x.baseline_repeat` | 120 s |
| `repeat.check` | 30 s |
| `compare.a1` | 30 s |
| `eval.a1.<cell>.<arm>`, ×6 | 60 s each |
