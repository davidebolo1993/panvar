#!/usr/bin/env python3
"""M15-C runner: budget 1,800 s; one process at a time; 4 threads. Units: build (records the help-corrected binary's hash),
the CLI contract test, and the documented examples (examples.py). Usage: m15c_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import runner  # noqa: E402
import m11d  # noqa: E402

WT = m11d.WT
OUT = f"{m11d.R}/m15c_interface_review_20261021"
BUDGET, PRE = 1800, 150.0   # disclosed: two panvar rebuilds for the help text (40 s), one CLI-contract test run, profile composition, reading
EXAMPLES = [("standard_c4", 150), ("standard_cyp2d6", 150), ("standard_acot", 300), ("resource_adjusted", 150), ("time_limit_binding", 150),
            ("invalid", 120)]
UNITS = [("build", 300), ("cli_contract", 60)] + [(f"ex.{n}", t) for n, t in EXAMPLES]


def command(unit, d):
    if unit == "build":
        return ["sh", "-c", f"cmake --build {WT}/build --target panvar genotype_probe genotype_cli_contract_test -j 4 && shasum -a 256 {WT}/build/panvar"]
    if unit == "cli_contract":
        return ["ctest", "--test-dir", f"{WT}/build", "-R", "panvar_genotype_cli_contract$", "--output-on-failure"]
    return [sys.executable, os.path.join(HERE, "examples.py"), unit.split(".", 1)[1], d]


R = runner.Runner(OUT, BUDGET, PRE, UNITS, command, lambda u: ["build"] if u != "build" else [], revision=lambda: runner.git_revision(WT),
                  require_clean=runner.git_clean(WT))

if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        os.makedirs(f"{OUT}/runs", exist_ok=True)
        text = "\n".join([f"# executable\tpanvar\t{m11d.sha(WT + '/build/panvar')}"] +
                         [f"# profile\t{p}\t{m11d.sha(WT + '/profiles/' + p)}" for p in sorted(os.listdir(WT + "/profiles")) if p.startswith("genotype_experimental")] +
                         ["unit\tdeclared_timeout_seconds"] + [f"{u}\t{t}" for u, t in UNITS]) + "\n"
        open(f"{OUT}/MANIFEST.tsv", "w").write(text)
        print(text)
    elif sys.argv[1] == "next":
        sys.exit(R.next())
    else:
        sys.exit(R.run_unit(sys.argv[1]))
