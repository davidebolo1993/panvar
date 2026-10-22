# Regression tests

Build with `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`, then
`cmake --build build -j 4` and `ctest --test-dir build --output-on-failure`.
Python 3 is needed for assessment fixtures; reporting tests additionally use
pandas, NumPy and matplotlib. The existing module shell tests and shared graph
fixtures remain in their usual locations. Long real-data tests remain opt-in
with `-DPANVAR_SLOW_TESTS=ON`.

The `genotype_{smoke,cli,index,assessment,dosage,resources,exactness}.sh` scripts
run named CTest groups: e.g. `bash tests/genotype_cli.sh build`.
They fail if the group is missing. All groups also run in the ordinary full suite.
`helpers/` holds C++ assertions and independent Python oracles; `fixtures/`
contains immutable profile parsing examples, not profiles users need to run.

The release suite requires neither campaign outputs nor the general research
probe. Historical experiments, schedules and result analyses were archived
outside this repository before removal, with a source snapshot and coverage map.
Current production code and mathematical assertions were not simplified during
this organization pass. Within-process concurrent callers remain unsupported.
