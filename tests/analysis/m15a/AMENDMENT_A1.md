# M15-A amendment A1 (declared before further execution)

**State at declaration:**
- `build` passed (23.3 s).
- `ctest` stopped at its 420 s unit watchdog. Spent: 506.6 of 900 s.
  - **Passed:** 58 tests.
  - **Disabled in the suite:** 2 (`panvar_real_regressions`, `panvar_rebuild_stats`).
  - **`panvar_genotype_command_e2e`** reached the declared **120 s per-test cap**. Every sub-check it had printed passed. Its
    recorded cost (8.9 s, CTestCostData) was stale.
  - **`panvar_smoke_real_c4`** started last and was interrupted by the unit watchdog.

**Fix:**
- re-run **exactly those two tests** as `ctest.a1`, with a 280 s per-test timeout under a 300 s unit watchdog;
- `reproduce` runs as frozen.

**Not changed:** the test selection otherwise, and every test.

**The failed `ctest` unit stays recorded.** If either test fails or times out again, it is reported as unverified, not
re-tried.
