# Tests

| folder | holds |
|---|---|
| `suites/` | the test scripts, one per module plus smoke and regression runs; CMake registers each one |
| `tools/` | the suites' helpers: fixture writers and the diff and report tools of `call_regress.sh` |
| `synthetic_data/` | small graphs and fixtures the suites run on; `make_synthetic.py` writes `syn.gfa` and `syn_w.gfa` |
| `real_data/` | graphs of human loci and their truth sets (see its README) |
| `examples/` | outputs of `scripts/*_examples.sh` and the `dev/` checks; not committed (see below) |

Run the suites from the build folder:

```bash
ctest --output-on-failure        # every enabled suite
ctest -LE real                   # without the ones on real loci
```

The slow suites on real loci (`rebuild_stats.sh`, `real_regressions.sh`) are off unless CMake is configured with `-DPANVAR_SLOW_TESTS=ON`; then `ctest -L slow` runs them. `call_regress.sh` compares `call` with a pinned earlier build and also needs `-DPANVAR_ORIG_BIN=<that binary>`; every known difference is listed in `suites/call_regress.expected.tsv`.

A suite also runs by hand with the arguments its header lists, for example `bash tests/suites/call_stats.sh build/panvar <output dir>`.

## Examples

The example outputs, `examples/<module>/`, are regenerated from the repository root. Each script reads the one before it, so they run in this order (about 1.5 hours on 2 threads, most of it `rebuild` and the thresholds):

```bash
for s in rebuild bubble inspect panphorte refine call describe benchmark associate; do
  bash scripts/${s}_examples.sh
done
bash scripts/benchmark_thresholds.sh   # after benchmark_examples.sh
bash dev/call_truth.sh                 # copy number against the truth sets, after call_examples.sh
```

The read-recruitment examples run on a synthetic genome set, in seconds, the second after the first:

```bash
bash scripts/background_examples.sh
bash scripts/recruit_examples.sh
```

They need `python3` and `Rscript` (ggplot2, ggrepel, data.table). The `dev/` checks that call svim-asm, truvari or GEMMA, and the read-recruitment benchmarks, are described in `dev/README.md`.
