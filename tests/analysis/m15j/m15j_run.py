#!/usr/bin/env python3
"""M15-J runner: budget 2,400 s; one process at a time; 4 threads. Units (DESIGN.md section 4): build, the CLI contract, the
ordinary development controls and the dosage check against M15-G, refusals, the panel index against M15-E, the unrepaired view
against the M15-I binary, and the full enabled suite. Usage: m15j_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import runner  # noqa: E402
import m11d  # noqa: E402

WT = m11d.WT
OUT = f"{m11d.R}/m15j_interface_20261022"
# disclosed development (before the harness freeze): incremental rebuilds of panvar and the two contract tests, about ten
# contract-test runs (< 1 s each), two curation-contract runs (3 s), manual refusal checks (< 1 s each), help previews and the
# inventory generation
BUDGET, PRE = 2400, 240.0
UNITS = [("build", 300), ("cli_contract", 60), ("ex.ordinary_c4", 150), ("ex.ordinary_cyp2d6", 200), ("ex.ordinary_acot", 300),
         ("ex.dosage_check_cyp2d6", 250), ("ex.dosage_check_c4", 150), ("ex.refusals", 120), ("ex.index_c4", 200),
         ("ex.unrefined", 300), ("ctest.full", 1250),
         # A1 (declared before the rerun): ex.index_c4 required the index to be byte-identical to M15-E's, but the header records
         # the builder's source revision (40 hex characters) and a 32-byte header digest covering it, so 71 bytes differ by
         # design. 37 of 38 checks had passed (cold = indexed, both equal to M15-E). The check now allows exactly those two
         # fields and requires every other byte to be identical; the unit is rerun (the failure stays recorded)
         ("a1.ex.index_c4", 200)]


def command(unit, d):
    if unit == "build":
        return ["sh", "-c", f"cmake --build {WT}/build -j 4 && shasum -a 256 {WT}/build/panvar"]
    if unit == "cli_contract":
        return ["ctest", "--test-dir", f"{WT}/build", "-R", "panvar_genotype_(cli|curation)_contract$", "--output-on-failure", "-V"]
    if unit == "ctest.full":
        return ["ctest", "--test-dir", f"{WT}/build", "-j", "1", "--timeout", "900", "--output-on-failure"]
    return [sys.executable, os.path.join(HERE, "examples.py"), unit.split("ex.", 1)[1], d]


R = runner.Runner(OUT, BUDGET, PRE, UNITS, command, lambda u: [] if u == "build" else ["build"],
                  revision=lambda: runner.git_revision(WT), require_clean=runner.git_clean(WT))

if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        os.makedirs(f"{OUT}/runs", exist_ok=True)
        text = "\n".join([f"# executable\tpanvar\t{m11d.sha(WT + '/build/panvar')}", f"# source\tgit\t{runner.git_revision(WT)}",
                          f"# baseline\tpanvar.m15i\t{m11d.sha(OUT + '/baseline/panvar.m15i')}"] +
                         ["unit\tdeclared_timeout_seconds"] + [f"{u}\t{t}" for u, t in UNITS]) + "\n"
        open(f"{OUT}/MANIFEST.tsv", "w").write(text)
        print(text)
    elif sys.argv[1] == "next":
        sys.exit(R.next())
    else:
        sys.exit(R.run_unit(sys.argv[1]))
