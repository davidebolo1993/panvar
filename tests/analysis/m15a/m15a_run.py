#!/usr/bin/env python3
"""M15-A verification runner: budget 900 s; one process at a time; 4 threads for builds/tests where applicable.
Units: build (every target), ctest (the frozen focused selection below, -j1), reproduce (archived-evidence check, no caller run).
Usage: m15a_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import runner  # noqa: E402
import m11d  # noqa: E402

WT = m11d.WT
OUT = f"{m11d.R}/m15a_release_candidate_20261021"
BUDGET, PRE = 900, 60.0   # disclosed: reading, one reproduction comparison of archived files (no caller run)
# frozen focused selection: every test whose recorded cost is < 30 s (CTestCostData), i.e. all except the three heaviest
# (panvar_genotype_report 105 s, panvar_genotype_branch_lifetime 98 s, panvar_genotype_block_delta 72 s) -- excluded for budget
EXCLUDE = "panvar_genotype_report$|panvar_genotype_branch_lifetime$|panvar_genotype_block_delta$"
UNITS = [("build", 420), ("ctest", 420), ("reproduce", 60)]
NEEDS = {"ctest": ["build"]}


def command(unit, d):
    if unit == "build":
        return ["cmake", "--build", f"{WT}/build", "-j", "4"]
    if unit == "ctest":
        return ["ctest", "--test-dir", f"{WT}/build", "-j", "1", "-E", EXCLUDE, "--output-on-failure", "--timeout", "120"]
    return [sys.executable, os.path.join(HERE, "reproduce.py"), d]


R = runner.Runner(OUT, BUDGET, PRE, UNITS, command, lambda u: NEEDS.get(u, []), revision=lambda: runner.git_revision(WT),
                  require_clean=runner.git_clean(WT))

if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        os.makedirs(f"{OUT}/runs", exist_ok=True)
        text = "\n".join([f"# executable\t{b}\t{m11d.sha(f'{WT}/build/{b}')}" for b in ("panvar", "genotype_probe")] +
                         ["unit\tdeclared_timeout_seconds"] + [f"{u}\t{t}" for u, t in UNITS]) + "\n"
        open(f"{OUT}/MANIFEST.tsv", "w").write(text)
        print(text)
    elif sys.argv[1] == "next":
        sys.exit(R.next())
    else:
        sys.exit(R.run_unit(sys.argv[1]))
