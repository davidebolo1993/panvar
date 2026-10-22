#!/usr/bin/env python3
"""M15-K runner: budget 2,400 s; one process at a time; 4 threads. Units (DESIGN.md section 4): build, the contract tests (CLI,
curation, source provenance, the e2e unrefined-v2 gate), refined runs against M15-J, unrefined runs against the M15-J binary
(differences classified), the index with producer provenance, a rebuild after a new commit, and the full enabled suite.
Usage: m15k_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import runner  # noqa: E402
import m11d  # noqa: E402

WT = m11d.WT
OUT = f"{m11d.R}/m15k_unrefined_v2_20261022"
# disclosed development (before the harness freeze): incremental rebuilds, one reconfigure (CMakeLists.txt gained the
# provenance target), contract/curation/prepared-locus/provenance test runs (< 5 s each), five focused e2e unrefined-v2 runs
# (about 20 s each, recording and then asserting the outcomes), one unrefined C4 trial run and one unrefined C4 calibration
# against the M15-J binary (about 30 s each), help previews and the inventory generation
BUDGET, PRE = 2400, 300.0
UNITS = [("build", 300), ("contract", 240), ("ex.refined_c4", 150), ("ex.refined_cyp2d6", 200), ("ex.refined_acot", 300),
         ("ex.unrefined_c4", 200), ("ex.unrefined_cyp2d6", 250), ("ex.unrefined_acot", 300), ("ex.index_c4", 200),
         ("ex.provenance_rebuild", 200), ("ctest.full", 1300)]


def command(unit, d):
    if unit == "build":
        return ["sh", "-c", f"cmake --build {WT}/build -j 4 && shasum -a 256 {WT}/build/panvar && cat {WT}/build/generated/panvar_source_provenance.h"]
    if unit == "contract":
        return ["sh", "-c", f"ctest --test-dir {WT}/build -R 'panvar_genotype_(cli|curation)_contract$|panvar_source_provenance$' "
                            f"--output-on-failure -V && PANVAR_E2E_ONLY=unrefined_v2 {WT}/build/genotype_command_e2e_test {d}/e2e"]
    if unit == "ctest.full":
        return ["ctest", "--test-dir", f"{WT}/build", "-j", "1", "--timeout", "900", "--output-on-failure"]
    return [sys.executable, os.path.join(HERE, "examples.py"), unit.split("ex.", 1)[1], d]


R = runner.Runner(OUT, BUDGET, PRE, UNITS, command, lambda u: [] if u == "build" else ["build"],
                  revision=lambda: runner.git_revision(WT), require_clean=runner.git_clean(WT))

if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        os.makedirs(f"{OUT}/runs", exist_ok=True)
        text = "\n".join([f"# executable\tpanvar\t{m11d.sha(WT + '/build/panvar')}", f"# source\tgit\t{runner.git_revision(WT)}",
                          f"# baseline\tpanvar.m15j\t{m11d.sha(OUT + '/baseline/panvar.m15j')}"] +
                         ["unit\tdeclared_timeout_seconds"] + [f"{u}\t{t}" for u, t in UNITS]) + "\n"
        open(f"{OUT}/MANIFEST.tsv", "w").write(text)
        print(text)
    elif sys.argv[1] == "next":
        sys.exit(R.next())
    else:
        sys.exit(R.run_unit(sys.argv[1]))
