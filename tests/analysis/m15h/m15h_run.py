#!/usr/bin/env python3
"""M15-H runner: budget 1,800 s; one process at a time; 4 threads. Units: build (hashes), the CLI contract (help
organisation, corrected claims, inventory coverage and defaults, the guide's generated inventory), the captured help texts,
preserved outputs (C4, against M15-G), and the full enabled suite on the final binary. Usage: m15h_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import runner  # noqa: E402
import m11d  # noqa: E402

WT = m11d.WT
OUT = f"{m11d.R}/m15h_help_20261022"
# disclosed development (before the freeze): rebuilds of panvar and the contract test, five contract-test runs (<1 s each),
# help previews and the inventory generation
BUDGET, PRE = 1800, 90.0
UNITS = [("build", 300), ("cli_contract", 60), ("ex.help_text", 30), ("ex.preserved_c4", 150), ("ctest.full", 1250)]


def command(unit, d):
    if unit == "build":
        return ["sh", "-c", f"cmake --build {WT}/build -j 4 && shasum -a 256 {WT}/build/panvar"]
    if unit == "cli_contract":
        return ["ctest", "--test-dir", f"{WT}/build", "-R", "panvar_genotype_cli_contract$", "--output-on-failure", "-V"]
    if unit == "ctest.full":
        return ["ctest", "--test-dir", f"{WT}/build", "-j", "1", "--timeout", "900", "--output-on-failure"]
    return [sys.executable, os.path.join(HERE, "examples.py"), unit.split("ex.", 1)[1], d]


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
