#!/usr/bin/env python3
"""M15-D.1 runner: budget 1,800 s; one process at a time; 4 threads; 6 GiB. Units: build (records the binary's hash), the
focused reproductions (the e2e sequence-first gates, then binary-level examples), then the full enabled suite on the final
revision. Usage: m15d1_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import runner  # noqa: E402
import m11d  # noqa: E402

WT = m11d.WT
OUT = f"{m11d.R}/m15d1_correctness_20261021"
# disclosed development (before the freeze): three runs of the e2e sequence-first subset (44.7 + 42.2 + 45.5 s), one
# 4 s smoke run of the footprint hook through a profile, and rebuilds
BUDGET, PRE = 1800, 240.0
UNITS = [("build", 300), ("e2e.sequence_first", 150), ("ex.library_preflight", 120), ("ex.term_time", 120),
         ("ex.term_footprint", 120), ("ex.term_checkpoint", 150), ("ex.term_nocall", 120), ("ex.term_locus_failure", 150),
         ("ex.cyp_audit", 150), ("ctest.full", 1200)]


def command(unit, d):
    if unit == "build":
        return ["sh", "-c", f"cmake --build {WT}/build -j 4 && shasum -a 256 {WT}/build/panvar {WT}/build/genotype_command_e2e_test"]
    if unit == "e2e.sequence_first":   # termination, declared-library preflight, the restored API dosage-off provenance (sf21_off)
        return ["env", "PANVAR_E2E_ONLY=sequence_first", f"{WT}/build/genotype_command_e2e_test", d]
    if unit == "ctest.full":
        return ["ctest", "--test-dir", f"{WT}/build", "-j", "1", "--timeout", "900", "--output-on-failure"]
    return [sys.executable, os.path.join(HERE, "examples.py"), unit.split(".", 1)[1], d]


R = runner.Runner(OUT, BUDGET, PRE, UNITS, command, lambda u: [] if u == "build" else ["build"],
                  revision=lambda: runner.git_revision(WT), require_clean=runner.git_clean(WT))

if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        os.makedirs(f"{OUT}/runs", exist_ok=True)
        text = "\n".join([f"# executable\tpanvar\t{m11d.sha(WT + '/build/panvar')}",
                          f"# source\tgit\t{runner.git_revision(WT)}"] +
                         ["unit\tdeclared_timeout_seconds"] + [f"{u}\t{t}" for u, t in UNITS]) + "\n"
        open(f"{OUT}/MANIFEST.tsv", "w").write(text)
        print(text)
    elif sys.argv[1] == "next":
        sys.exit(R.next())
    else:
        sys.exit(R.run_unit(sys.argv[1]))
