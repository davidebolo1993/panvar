#!/usr/bin/env python3
"""S2-B validation runner (frozen manifest; one process at a time; <= 4 threads; <= 300 s per unit; 1200 s budget,
counted as time spent in this batch -- development runs included -- plus the declared timeout of the next unit; a unit
that does not fit is NOT_RUN).

Order: part 1 tests (sequence report, CLI contract, the full e2e with test_early_incumbent and test_early_incumbent_caps
in the normal route) -> the bound's admissibility gate (s2b_preplacement_bound fixtures) -> per 5x cell (C4, CYP2D6,
ACOT; P1 overlay, sequence_first_v2, early_incumbent, diagnostic_union_evidence_dump): the evidence-dump run, then the
bound measurement against it -> ACOT 30x: bound CONSTRUCTION only (the ACOT 5x union's haplotypes, 30x reads, the
incumbent's haplotypes plus the first 40 others; no reference, no exclusion claim).

  s2b_run.py manifest | next | <unit>
"""
import csv
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
import p1_run  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/s2b_preplacement_20261018"
BUDGET = 1200
MANIFEST = [("t.genotype_sequence_report", "ctest", "genotype_sequence_report", 30),
            ("t.genotype_cli_contract", "ctest", "genotype_cli_contract", 30),
            ("t.s2b_preplacement_bound_fixtures", "ctest", "s2b_preplacement_bound_fixtures", 60),
            ("t.genotype_command_e2e", "ctest", "genotype_command_e2e", 300),
            ("dump.c4", "dump", "c4_5", 90),
            ("measure.c4", "measure", "c4_5", 120),
            ("dump.cyp", "dump", "cyp_5", 90),
            ("measure.cyp", "measure", "cyp_5", 120),
            ("dump.acot5", "dump", "acot_5", 150),
            ("measure.acot5", "measure", "acot_5", 300),
            ("construct.acot30", "construct", "acot_30", 300)]
ORDER = [m[0] for m in MANIFEST]
UNIT_OF_DUMP = {"c4_5": "dump.c4", "cyp_5": "dump.cyp", "acot_5": "dump.acot5", "acot_30": "dump.acot5"}


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
    cell = p1_run.cells()[target] if kind != "ctest" else None
    if kind == "ctest":
        cmd = ["nice", "ctest", "--test-dir", BUILD, "-R", f"^panvar_{target}$", "--output-on-failure", "-j1"]
    elif kind == "dump":
        prof = p1_run.profile_for(cell, "opt", d, timeout, target)
        args = p1_run.genotype_args(cell, prof, f"{d}/call")
        i = args.index("--dosage-mode")
        del args[i:i + 2]
        extra = f"{d}/profile.s2b.tsv"
        open(extra, "w").write("profile\ts2b-evidence\nearly_incumbent\texperimental_v1\ndiagnostic_union_evidence_dump\t1\n")
        args[args.index("--profile") + 2:args.index("--profile") + 2] = ["--profile", extra]
        cmd = ["nice", f"{BUILD}/panvar", "genotype", *args, "--output-schema", "sequence_first_v2"]
    else:
        dump = f"{OUT}/runs/{UNIT_OF_DUMP[target]}/call.union_evidence.bin"
        cmd = ["nice", f"{BUILD}/s2b_preplacement_bound", "measure", dump, cell["reads1"], cell["reads2"], "4", f"{d}/measure.tsv"]
        if kind == "construct":
            cmd += ["construction_only", "42"]
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
