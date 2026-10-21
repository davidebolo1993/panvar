#!/usr/bin/env python3
"""M15-D runner: budget 2,400 s; one process at a time; 4 threads; 6 GiB. Units: build (records the binary's hash), the CLI
contract test, the full enabled suite, then the basic-command examples (examples.py) on development controls only.
Usage: m15d_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import runner  # noqa: E402
import m11d  # noqa: E402

WT = m11d.WT
OUT = f"{m11d.R}/m15d_basic_command_20261021"
# disclosed development (before the freeze): rebuilds of panvar and every target, the focused non-heavy suite (149 s), two
# reruns of the two tests whose old-default expectations were updated, the command e2e prepared-locus subset (45 s)
# Amendments A1/A2 (declared in DESIGN.md section 5 before their units ran) add 360 s of disclosed development: two
# command_e2e development runs (136.8 + 169.8 s), the quick contract tests and two rebuilds.
BUDGET, PRE = 2400, 660.0
EXAMPLES = [("invalid", 120), ("noprofile_c4", 300), ("noprofile_cyp2d6", 400), ("noprofile_acot", 600), ("auto_equivalence_c4", 300),
            ("declared_c4", 150), ("declared_cyp2d6", 200), ("declared_acot", 300), ("dosage_absent_c4", 150),
            ("dosage_absent_cyp2d6", 200), ("dosage_absent_acot", 300), ("unrepaired_view", 150), ("auto_refused", 120),
            ("resource_stop", 150)]
UNITS = [("build", 300), ("cli_contract", 60), ("ctest.full", 1200)] + [(f"ex.{n}", t) for n, t in EXAMPLES]
# A1/A2 reruns: the failed units only (the originals stay recorded)
AMEND = [("a.build", 120), ("a2.ex.auto_equivalence_c4", 60), ("a2.ex.declared_c4", 60), ("a2.ex.declared_cyp2d6", 90),
         ("a2.ex.dosage_absent_c4", 60), ("a2.ex.dosage_absent_cyp2d6", 80), ("a2.ex.resource_stop", 40), ("a1.ctest", 240)]
UNITS += AMEND
NEEDS = {"a2.ex.auto_equivalence_c4": ["ex.noprofile_c4"], "a2.ex.dosage_absent_c4": ["a2.ex.declared_c4"],
         "a2.ex.dosage_absent_cyp2d6": ["a2.ex.declared_cyp2d6"],
         "ex.auto_equivalence_c4": ["ex.noprofile_c4"], "ex.dosage_absent_c4": ["ex.declared_c4"],
         "ex.dosage_absent_cyp2d6": ["ex.declared_cyp2d6"], "ex.dosage_absent_acot": ["ex.declared_acot"]}


def command(unit, d):
    if unit == "a1.ctest":
        return ["ctest", "--test-dir", f"{WT}/build", "-j", "1", "--timeout", "230", "--output-on-failure", "-R",
                "panvar_genotype_command_e2e$|panvar_genotype_cli_contract$|panvar_genotype_sequence_report$|panvar_library_estimate$"]
    if unit == "a.build":
        unit = "build"
    if unit == "build":
        return ["sh", "-c", f"cmake --build {WT}/build -j 4 && shasum -a 256 {WT}/build/panvar"]
    if unit == "cli_contract":
        return ["ctest", "--test-dir", f"{WT}/build", "-R", "panvar_genotype_cli_contract$", "--output-on-failure"]
    if unit == "ctest.full":
        return ["ctest", "--test-dir", f"{WT}/build", "-j", "1", "--timeout", "900", "--output-on-failure"]
    return [sys.executable, os.path.join(HERE, "examples.py"), unit.split("ex.", 1)[1], d]


def needs(u):
    if u in ("build", "a.build"):
        return []
    return (["a.build"] if u.startswith("a") else ["build"]) + NEEDS.get(u, [])


R = runner.Runner(OUT, BUDGET, PRE, UNITS, command, needs, revision=lambda: runner.git_revision(WT), require_clean=runner.git_clean(WT))

if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        os.makedirs(f"{OUT}/runs", exist_ok=True)
        text = "\n".join([f"# executable\tpanvar\t{m11d.sha(WT + '/build/panvar')}"] +
                         ["unit\tdeclared_timeout_seconds\tneeds"] + [f"{u}\t{t}\t{','.join(needs(u)) or '-'}" for u, t in UNITS]) + "\n"
        open(f"{OUT}/MANIFEST.tsv", "w").write(text)
        print(text)
    elif sys.argv[1] == "next":
        sys.exit(R.next())
    else:
        sys.exit(R.run_unit(sys.argv[1]))
