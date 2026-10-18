#!/usr/bin/env python3
"""S1-C validation runner (frozen manifest; one process at a time; <= 4 threads; <= 300 s per unit; 1200 s budget,
counted as time spent in this batch plus the declared timeout of the next unit; a unit that does not fit is NOT_RUN).
The implementation under test is this worktree's build; every unit records the revision and requires a clean tree.

Order: the sequence-report unit test (new S1-C cases: orientation, linear interval, publication failure) and the CLI
contract -> Sol's two reproductions (orientation: his unmodified source against this build, and his compiled S1-B binary;
lifecycle: CLI sequences with the S1-B and S1-C binaries) -> fixture gates (permutations, oracle, legacy identity,
persisted-index route) -> the end-to-end suite (S1-B schema test plus the S1-C lifecycle regressions) -> the CYP2D6 30x
comparison closure (saved runs only, with negative controls; no sentinel rerun) -> every other enabled test, one per unit,
heaviest last.

  s1c_run.py manifest | next | <unit>
"""
import csv
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "s1b"))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
import s1b_run  # noqa: E402  (test lists, test-name rule)
import p1_run  # noqa: E402  (run_timed)

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/s1c_reporting_integrity_20261018"
S1B = "/Users/davide.bolognini/panvar/tests/results/s1b_sequence_first_20261018"
S1B_BUILD = "/Users/davide.bolognini/panvar-worktrees/s1-b-sequence-first/build"
SOL = "/private/tmp/panvar-s1b-review.yFD8kk"
P11_CYP = "/Users/davide.bolognini/panvar/tests/results/p11_validation_closure_20261017/runs/o.cyp_30"
F1_READS = "/Users/davide.bolognini/panvar/tests/results/representation_fixtures_20261015/F1/reads"
BUDGET = 1200
tn = s1b_run.test_name
MANIFEST = [("t.genotype_sequence_report", "ctest", "genotype_sequence_report", 30),
            ("t.genotype_cli_contract", "ctest", "genotype_cli_contract", 30),
            ("repro.orientation", "repro_orientation", "", 120),
            ("repro.lifecycle", "repro_lifecycle", "", 120),
            ("fixtures.A", "fixtures", "F1,F2,F3,F4", 240),
            ("fixtures.B", "fixtures", "F5,F6,F7,F7full,F9", 240),
            ("t.genotype_command_e2e", "ctest", "genotype_command_e2e", 300),
            ("cyp30.closure", "cyp_closure", "", 60)]
MANIFEST += [(f"t.{tn(t)}", "ctest", tn(t), 60) for t in s1b_run.LIGHT + s1b_run.REST_SMALL + s1b_run.OTHER]
MANIFEST += [(f"t.{t}", "ctest", t, s) for t, s in s1b_run.HEAVY]
LEDGER = "unit\trevision\tdirty_tree\tdeclared_timeout_seconds\twall_seconds\tpeak_rss_gib\texit\ttimed_out\tnote\n"


def spent():
    p = f"{OUT}/ledger.tsv"
    return sum(float(r["wall_seconds"]) for r in csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else 0.0


def done():
    p = f"{OUT}/ledger.tsv"
    return {r["unit"] for r in csv.DictReader(open(p), delimiter="\t")} if os.path.exists(p) else set()


def state():
    rev = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", WT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
    return rev, len([l for l in dirty.split("\n") if l.strip()])


def write_manifest():
    os.makedirs(OUT, exist_ok=True)
    text = "unit\tkind\ttarget\tdeclared_timeout_seconds\n" + "".join("\t".join(str(x) for x in m) + "\n" for m in MANIFEST)
    p = f"{OUT}/MANIFEST.tsv"
    if os.path.exists(p):
        assert open(p).read() == text, "a different manifest is already frozen"
    else:
        open(p, "w").write(text)
    print(text, end="")


def record(unit, rev, dirty, timeout, wall, rss, code, timed_out, note):
    new = not os.path.exists(f"{OUT}/ledger.tsv")
    with open(f"{OUT}/ledger.tsv", "a") as f:
        if new:
            f.write(LEDGER)
        f.write(f"{unit}\t{rev}\t{dirty}\t{timeout}\t{wall:.2f}\t{rss:.3f}\t{code}\t{int(bool(timed_out))}\t{note}\n")


def orientation(d):
    """Sol's compiled S1-B binary as delivered, then his unmodified source compiled against this build."""
    print("== Sol's binary (compiled against S1-B):", flush=True)
    old = subprocess.run([f"{SOL}/orientation"], capture_output=True, text=True)
    print(old.stdout.strip(), f"\nexit {old.returncode}")
    src = f"{d}/orientation.cpp"
    open(src, "w").write(open(f"{SOL}/orientation.cpp").read())
    cc = subprocess.run(["c++", "-O3", "-DNDEBUG", "-std=c++17", "-arch", "arm64", f"-I{WT}/include", "-c", src, "-o",
                         f"{d}/orientation.o"], capture_output=True, text=True)
    ld = subprocess.run(["c++", "-O3", "-arch", "arm64", f"{d}/orientation.o",
                         f"{BUILD}/CMakeFiles/genotype_sequence_report_test.dir/src/cli_utils.cpp.o", "-o", f"{d}/orientation_s1c",
                         f"{BUILD}/libpanvarlib.a", f"{BUILD}/libedlib.a", f"{WT}/external/abpoa/lib/libabpoa.a",
                         f"{WT}/external/minigraph/libminigraph.a", "-lz"], capture_output=True, text=True)
    if cc.returncode or ld.returncode:
        print("build failed", cc.stderr[-400:], ld.stderr[-400:])
        return 1
    new = subprocess.run([f"{d}/orientation_s1c"], capture_output=True, text=True)
    print("== the same source against S1-C:\n" + new.stdout.strip(), f"\nexit {new.returncode}")
    return 0 if old.returncode == 1 and new.returncode == 0 else 1


def run_unit(unit):
    u, kind, target, timeout = next(m for m in MANIFEST if m[0] == unit)
    assert open(f"{OUT}/MANIFEST.tsv").read() and timeout <= 300 and unit not in done()
    rev, dirty = state()
    assert dirty == 0, "the tree must be committed and clean"
    used = spent()
    if used + timeout > BUDGET:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {used:.1f} s + {timeout} s > {BUDGET} s")
        print(f"{unit}: NOT RUN")
        return 3
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d), f"{d} exists; a run is never repeated into the same place"
    os.makedirs(d)
    me = os.path.abspath(__file__)
    if kind == "ctest":
        cmd = ["nice", "ctest", "--test-dir", BUILD, "-R", f"^panvar_{target}$", "--output-on-failure", "-j1"]
    elif kind == "fixtures":
        cmd = ["nice", sys.executable, os.path.join(HERE, "..", "s1b", "s1b_fixtures.py"), f"{OUT}/fixtures",
               f"{BUILD}/panvar", f"{BUILD}/cyp2d6_block_alleles", f"{BUILD}/genotype_prepared_locus_driver", target]
    elif kind == "repro_lifecycle":
        cmd = ["nice", sys.executable, os.path.join(HERE, "repro_lifecycle.py"), d, f"{S1B_BUILD}/panvar", f"{BUILD}/panvar",
               f"{S1B}/fixtures_a1/F1", F1_READS]
    elif kind == "repro_orientation":
        cmd = [sys.executable, me, "_orientation", d]
    else:
        cmd = [sys.executable, os.path.join(HERE, "cyp_closure.py"), P11_CYP, f"{S1B}/runs/cyp30.v1"]
    cmd = ["/usr/bin/time", "-l", *cmd]
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s) " + " ".join(cmd) + "\n")
        log.flush()
        code, timed_out = p1_run.run_timed(cmd, log, timeout)
    wall = time.time() - t
    rss = 0
    text = open(f"{d}/run.log").read()
    for line in text.splitlines():
        if "maximum resident set size" in line:
            rss = int(line.split()[0])
    note = "-"
    if kind == "ctest":
        m = re.search(r"Test +#\d+: \S+ \.+\s*\**\s*(\w[\w ]*?)\s+([\d.]+) sec", text)
        note = f"ctest: {m.group(1)} {m.group(2)} s" if m else "ctest: no result line"
    record(unit, rev, dirty, timeout, wall, rss / 2**30, code, timed_out, note)
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s, {note} (spent {spent():.1f}/{BUDGET})")
    return 0


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest()
    elif arg == "_orientation":
        sys.exit(orientation(sys.argv[2]))
    elif arg == "next":
        todo = [m[0] for m in MANIFEST if m[0] not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(arg))
