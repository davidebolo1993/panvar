#!/usr/bin/env python3
"""M15-G runner: budget 2,400 s; one process at a time; 4 threads; 6 GiB. Units: build (hashes), the assessment
equivalence test (C++ vs specification on every fixture), a read-only replay of every saved development v2 prefix, the
completed / dosage / legacy-tie / provisional / absent / reporting-failure runs, then the full enabled suite on the final
revision. Usage: m15g_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import runner  # noqa: E402
import m11d  # noqa: E402

WT = m11d.WT
OUT = f"{m11d.R}/m15g_assessment_20261022"
# disclosed development (before the freeze): five e2e subset runs (21 + 49 + 21 + 21 + 45 s), one C4 smoke run (20 s), the
# legacy-tie fixture check (<1 s), quick contract tests and spec fixture checks (~6 s), rebuilds (~40 s)
BUDGET, PRE = 2400, 240.0
UNITS = [("build", 300), ("spec_equivalence", 60), ("ex.replay_saved", 400), ("ex.complete_c4", 150), ("ex.complete_cyp2d6", 200),
         ("ex.complete_acot", 300), ("ex.dosage_cyp2d6", 200), ("ex.legacy_tie", 60), ("ex.provisional", 120), ("ex.absent", 120),
         ("ex.reporting_failure", 150), ("ctest.full", 1250)]


def command(unit, d):
    if unit == "build":
        return ["sh", "-c", f"cmake --build {WT}/build -j 4 && shasum -a 256 {WT}/build/panvar {WT}/build/genotype_assessment_test"]
    if unit == "spec_equivalence":
        return ["ctest", "--test-dir", f"{WT}/build", "-R", "panvar_genotype_assessment$", "--output-on-failure", "-V"]
    if unit == "ctest.full":
        return ["ctest", "--test-dir", f"{WT}/build", "-j", "1", "--timeout", "900", "--output-on-failure"]
    return [sys.executable, os.path.join(HERE, "examples.py"), unit.split(".", 1)[1], d]


R = runner.Runner(OUT, BUDGET, PRE, UNITS, command, lambda u: [] if u == "build" else ["build"],
                  revision=lambda: runner.git_revision(WT), require_clean=runner.git_clean(WT))

if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        os.makedirs(f"{OUT}/runs", exist_ok=True)
        text = "\n".join([f"# executable\tpanvar\t{m11d.sha(WT + '/build/panvar')}", f"# source\tgit\t{runner.git_revision(WT)}"] +
                         ["unit\tdeclared_timeout_seconds"] + [f"{u}\t{t}" for u, t in UNITS]) + "\n"
        open(f"{OUT}/MANIFEST.tsv", "w").write(text)
        print(text)
    elif sys.argv[1] == "next":
        sys.exit(R.next())
    else:
        sys.exit(R.run_unit(sys.argv[1]))
