#!/usr/bin/env python3
"""M15-L runner: budget 3,600 s (raised from the design's 2,400 s before this harness was frozen: see DESIGN.md); one process
at a time; 4 threads. Units (DESIGN.md section 5): build, the CLI and curation tests, the layout on C4, CYP2D6 and ACOT, the
dosage check, an indexed run, provisional and failed outcomes, an earlier flat run's reports, and the full enabled suite. Each
example runs the M15-L and the M15-K binaries on the same command. Usage: m15l_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import runner  # noqa: E402
import m11d  # noqa: E402

WT = m11d.WT
OUT = f"{m11d.R}/m15l_diagnostics_dir_20261022"
# disclosed development (DESIGN.md): rebuilds of panvar and the e2e test, one lost background e2e attempt, e2e reruns (some on a
# stale test binary), CLI/curation/report test runs, one C4 trial run and the report-loader check -- about 1,500 s, an upper
# estimate
BUDGET, PRE = 3600, 1500.0
UNITS = [("build", 300), ("tests", 120), ("ex.layout_c4", 200), ("ex.layout_cyp2d6", 250), ("ex.layout_acot", 350),
         ("ex.dosage_check_cyp2d6", 250), ("ex.indexed_c4", 200), ("ex.provisional", 150), ("ex.failed", 300),
         ("ex.stale_flat", 150), ("ctest.full", 1300),
         # A1 (declared before the reruns): ex.provisional compared the primary record byte for byte, but a provisional
         # record carries the checkpoint's timing rows (checkpoint_seconds_*), which differ by run; the record now gets the
         # timing-row allowance every keyed file has. ex.failed's record comparison replaced the primary prefix before the
         # longer diagnostics path it begins with, so the failure message naming the diagnostics locus table differed; the
         # longer path is now replaced first. Every layout and outcome check had passed (12/13, 11/12). Both are rerun; the
         # failures stay recorded.
         ("a1.ex.provisional", 150), ("a1.ex.failed", 300)]


def command(unit, d):
    if unit == "build":
        return ["sh", "-c", f"cmake --build {WT}/build -j 4 && shasum -a 256 {WT}/build/panvar"]
    if unit == "tests":
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
                          f"# baseline\tpanvar.m15k\t{m11d.sha(OUT + '/baseline/panvar.m15k')}", f"# disclosed_development_seconds\t{PRE}"] +
                         ["unit\tdeclared_timeout_seconds"] + [f"{u}\t{t}" for u, t in UNITS]) + "\n"
        open(f"{OUT}/MANIFEST.tsv", "w").write(text)
        print(text)
    elif sys.argv[1] == "next":
        sys.exit(R.next())
    else:
        sys.exit(R.run_unit(sys.argv[1]))
