#!/usr/bin/env python3
"""S2-A validation runner (frozen manifest; one process at a time; <= 4 threads; <= 300 s per unit; 1200 s budget,
counted as time spent in this batch -- development runs included -- plus the declared timeout of the next unit; a unit
that does not fit is NOT_RUN).

Order: the sequence-report unit test and the CLI contract -> the S2-A fixtures (s2a_fixtures.py: F1 equivalent
labellings, F5 genuine tie, F7full, F7 provisional; off / on / stopped) -> the end-to-end suite (+ test_early_incumbent:
deterministic interruptions before the first evaluation, after the checkpoint, during union scoring and refinement;
prefix reuse; publication failure; 1 vs 4 threads; independent rescoring of the early candidate) -> real data (P1
development overlay, sequence_first_v2, dosage by its default): C4 and CYP2D6 5x off/on (completed runs compared), ACOT
30x on (unchanged caps: does a checkpoint appear before its resource stop?) -> the summary table.

  s2a_run.py manifest | next | <unit>
"""
import csv
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
import p1_run  # noqa: E402  (cells, profiles, run_timed)

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/s2a_early_incumbent_20261018"
BUDGET = 1200
MANIFEST = [("t.genotype_sequence_report", "ctest", "genotype_sequence_report", 30),
            ("t.genotype_cli_contract", "ctest", "genotype_cli_contract", 30),
            ("fixtures", "fixtures", "F1,F5,F7full,F7", 120),
            ("t.genotype_command_e2e", "ctest", "genotype_command_e2e", 300),
            ("c4.off", "real", "c4_5:off", 90),
            ("c4.on", "real", "c4_5:on", 90),
            ("cyp.off", "real", "cyp_5:off", 90),
            ("cyp.on", "real", "cyp_5:on", 90),
            ("acot30.on", "real", "acot_30:on", 300),
            ("summary", "summary", "-", 30)]
ORDER = [m[0] for m in MANIFEST]
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
    return rev, len([line for line in dirty.split("\n") if line.strip()])


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
    with open(f"{OUT}/ledger.tsv", "a") as f:
        f.write(f"{unit}\t{rev}\t{dirty}\t{timeout}\t{wall:.2f}\t{rss:.3f}\t{code}\t{int(bool(timed_out))}\t{note}\n")


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
    if kind == "ctest":
        cmd = ["nice", "ctest", "--test-dir", BUILD, "-R", f"^panvar_{target}$", "--output-on-failure", "-j1"]
    elif kind == "fixtures":
        cmd = ["nice", sys.executable, os.path.join(HERE, "s2a_fixtures.py"), f"{OUT}/fixtures", f"{BUILD}/panvar"]
    elif kind == "real":
        cell_id, arm = target.split(":")
        cell = p1_run.cells()[cell_id]
        prof = p1_run.profile_for(cell, "opt", d, timeout, cell_id)
        args = p1_run.genotype_args(cell, prof, f"{d}/call")
        i = args.index("--dosage-mode")
        del args[i:i + 2]   # sequence_first_v2 decides dosage: audit by default
        if arm == "on":
            early = f"{d}/profile.early.tsv"
            open(early, "w").write("profile\ts2a-early-incumbent\nearly_incumbent\texperimental_v1\n")
            args[args.index("--profile") + 2:args.index("--profile") + 2] = ["--profile", early]
        cmd = ["nice", f"{BUILD}/panvar", "genotype", *args, "--output-schema", "sequence_first_v2"]
    else:
        cmd = [sys.executable, os.path.join(HERE, "s2a_summary.py"), OUT]
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
    elif arg == "next":
        todo = [u for u in ORDER if u not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(arg))
