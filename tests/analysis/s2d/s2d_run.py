#!/usr/bin/env python3
"""S2-D validation runner (frozen manifest; one process at a time; <= 4 threads; <= 300 s per unit; 1200 s budget,
counted as time spent in this batch -- development runs included -- plus the declared timeout of the next unit; a unit
that does not fit is NOT_RUN).

Order: the kernel gate (S2-D cache lifetime cases) and the full e2e (test_progressive_schedule in the normal route) ->
part 2, the kernel cache-lifetime comparison on S2-C's union inputs (C4, CYP2D6, ACOT 5x; per-call vs persistent, two
repeats each, alternating) -> part 4, end-to-end caller runs (P1 overlay, sequence_first_v2, early_incumbent; control =
the bulk indexed union, progressive = union_progressive_schedule experimental_v1, per-call cache) and their output
comparison -> one bounded ACOT 30x progressive caller run.

  s2d_run.py manifest | next | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/s2d_cache_lifetime_20261018"
S2C = "/Users/davide.bolognini/panvar/tests/results/s2c_demand_driven_20261018/runs"
BUDGET = 1200
CELLS = (("c4_5", "c4", 60, 90), ("cyp_5", "cyp", 60, 90), ("acot_5", "acot5", 120, 150))
MANIFEST = [("t.genotype_progressive_kernel", "ctest", "genotype_progressive_kernel", 60),
            ("t.genotype_command_e2e", "ctest", "genotype_command_e2e", 300)]
for cell, tag, kt, ct in CELLS:
    for rep in ("r1", "r2"):
        for life in ("per_call", "persist"):
            MANIFEST.append((f"kcache.{tag}.{life}.{rep}", "kcache", f"{cell}:{life}", kt))
for cell, tag, kt, ct in CELLS:
    MANIFEST += [(f"caller.{tag}.control", "caller", f"{cell}:control", ct),
                 (f"caller.{tag}.progressive", "caller", f"{cell}:progressive", ct),
                 (f"compare.{tag}", "compare", cell, 30)]
MANIFEST.append(("caller.acot30.progressive", "caller", "acot_30:progressive", 300))
ORDER = [m[0] for m in MANIFEST]
TAG = {c[0]: c[1] for c in CELLS}
TAG["acot_30"] = "acot30"


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
    elif kind == "kcache":
        cell_id, life = target.split(":")
        cell = p1_run.cells()[cell_id]
        inputs = f"{S2C}/dump.{TAG[cell_id]}/call.union_kernel_inputs.bin"
        cmd = ["nice", f"{BUILD}/s2c_progressive_harness", "progressive", inputs, cell["reads1"], cell["reads2"], f"{d}/kernel",
               "31", str(4 << 30), life]
    elif kind == "caller":
        cell_id, arm = target.split(":")
        cell = p1_run.cells()[cell_id]
        prof = p1_run.profile_for(cell, "opt", d, timeout, cell_id)
        args = p1_run.genotype_args(cell, prof, f"{d}/call")
        i = args.index("--dosage-mode")
        del args[i:i + 2]
        extra = f"{d}/profile.s2d.tsv"
        text = "profile\ts2d-caller\nearly_incumbent\texperimental_v1\n"
        if arm == "progressive":
            text += "union_progressive_schedule\texperimental_v1\n"
        open(extra, "w").write(text)
        args[args.index("--profile") + 2:args.index("--profile") + 2] = ["--profile", extra]
        cmd = ["nice", f"{BUILD}/panvar", "genotype", *args, "--output-schema", "sequence_first_v2"]
    else:
        tag = TAG[target]
        cmd = [sys.executable, os.path.join(HERE, "s2d_compare.py"), f"{OUT}/runs/caller.{tag}.control",
               f"{OUT}/runs/caller.{tag}.progressive", f"{d}/compare.tsv"]
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
