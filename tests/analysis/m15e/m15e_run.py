#!/usr/bin/env python3
"""M15-E runner: budget 2,400 s; one process at a time; 4 threads; 6 GiB. Units: build (records the hashes), the CLI
contract test and the e2e prepared-locus gates (the panel-index CLI on a tiny locus), then the binary-level examples on
the development cells, then the full enabled suite on the final revision. Usage: m15e_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import runner  # noqa: E402
import m11d  # noqa: E402

WT = m11d.WT
OUT = f"{m11d.R}/m15e_panel_index_20261021"
# disclosed development (before the freeze): five runs of the e2e prepared-locus subset (24 + 25 + 27 + 27 + 27 s), one C4
# smoke build/refusal run (~4 s), two CLI-contract runs and rebuilds
BUDGET, PRE = 2400, 260.0
UNITS = [("build", 300), ("cli_contract", 60), ("e2e.prepared_locus", 150),
         ("ex.build_c4", 120), ("ex.build_cyp2d6", 120), ("ex.build_acot", 180),
         ("ex.index_c4", 240), ("ex.index_cyp2d6", 300), ("ex.index_acot", 400),
         ("ex.dosage_check_c4", 200), ("ex.two_samples_c4", 240), ("ex.refusals_c4", 180),
         ("ctest.full", 1250),
         # A1 (declared before the rerun): ex.refusals_c4's curation manifest was named inside the "curation" refusal prefix,
         # so the nothing-written count included the manifest itself; renamed, the unit is rerun once (original kept)
         ("a1.ex.refusals_c4", 180)]
NEEDS = {"ex.index_c4": ["ex.build_c4"], "ex.index_cyp2d6": ["ex.build_cyp2d6"], "ex.index_acot": ["ex.build_acot"],
         "ex.dosage_check_c4": ["ex.build_c4"], "ex.two_samples_c4": ["ex.build_c4", "ex.index_c4"],
         "ex.refusals_c4": ["ex.build_c4"], "a1.ex.refusals_c4": ["ex.build_c4"]}


def command(unit, d):
    if unit == "build":
        return ["sh", "-c", f"cmake --build {WT}/build -j 4 && shasum -a 256 {WT}/build/panvar {WT}/build/genotype_command_e2e_test"]
    if unit == "cli_contract":
        return ["ctest", "--test-dir", f"{WT}/build", "-R", "panvar_genotype_cli_contract$", "--output-on-failure"]
    if unit == "e2e.prepared_locus":
        return ["env", "PANVAR_E2E_ONLY=prepared_locus", f"{WT}/build/genotype_command_e2e_test", d]
    if unit == "ctest.full":
        return ["ctest", "--test-dir", f"{WT}/build", "-j", "1", "--timeout", "900", "--output-on-failure"]
    return [sys.executable, os.path.join(HERE, "examples.py"), unit.split("ex.", 1)[1], d]


R = runner.Runner(OUT, BUDGET, PRE, UNITS, command, lambda u: [] if u == "build" else ["build"] + NEEDS.get(u, []),
                  revision=lambda: runner.git_revision(WT), require_clean=runner.git_clean(WT))

if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        os.makedirs(f"{OUT}/runs", exist_ok=True)
        text = "\n".join([f"# executable\tpanvar\t{m11d.sha(WT + '/build/panvar')}", f"# source\tgit\t{runner.git_revision(WT)}"] +
                         ["unit\tdeclared_timeout_seconds\tneeds"] +
                         [f"{u}\t{t}\t{','.join(NEEDS.get(u, [])) or '-'}" for u, t in UNITS]) + "\n"
        open(f"{OUT}/MANIFEST.tsv", "w").write(text)
        print(text)
    elif sys.argv[1] == "next":
        sys.exit(R.next())
    else:
        sys.exit(R.run_unit(sys.argv[1]))
