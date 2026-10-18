#!/usr/bin/env python3
"""S2-C validation runner (frozen manifest; one process at a time; <= 4 threads; <= 300 s per unit; 1200 s budget,
counted as time spent in this batch -- development runs included -- plus the declared timeout of the next unit; a unit
that does not fit is NOT_RUN).

Order: the progressive kernel's fixture gate, the indexed-route regressions (union placement cache; the full e2e, which
holds the indexed-union tests) -> per 5x cell (C4, CYP2D6, ACOT): the caller run writing union_kernel_inputs.bin (P1
overlay, sequence_first_v2, early_incumbent; the bulk union is the caller's own), then the bulk adapter and the
progressive schedule on the same inputs and revision (separate processes), then the bit-for-bit comparison -> ACOT 30x:
the caller run (writes the inputs before its union pass; its own stop is recorded) and one bounded progressive run under
the 4 GiB footprint guard (completed stages kept; no completion claimed unless it occurs).

  s2c_run.py manifest | next | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/s2c_demand_driven_20261018"
BUDGET = 1200
MANIFEST = [("t.genotype_progressive_kernel", "ctest", "genotype_progressive_kernel", 60),
            ("t.genotype_union_placement_cache", "ctest", "genotype_union_placement_cache", 120),
            ("t.genotype_command_e2e", "ctest", "genotype_command_e2e", 300)]
for cell, tag, dump_t, run_t in (("c4_5", "c4", 90, 60), ("cyp_5", "cyp", 90, 60), ("acot_5", "acot5", 150, 120)):
    MANIFEST += [(f"dump.{tag}", "dump", cell, dump_t), (f"bulk.{tag}", "bulk", cell, run_t),
                 (f"progressive.{tag}", "progressive", cell, run_t), (f"compare.{tag}", "compare", cell, 30)]
MANIFEST += [("dump.acot30", "dump", "acot_30", 300), ("progressive.acot30", "progressive", "acot_30", 300)]
ORDER = [m[0] for m in MANIFEST]
TAG = {"c4_5": "c4", "cyp_5": "cyp", "acot_5": "acot5", "acot_30": "acot30"}


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
    harness = f"{BUILD}/s2c_progressive_harness"
    if kind == "ctest":
        cmd = ["nice", "ctest", "--test-dir", BUILD, "-R", f"^panvar_{target}$", "--output-on-failure", "-j1"]
    else:
        cell = p1_run.cells()[target]
        inputs = f"{OUT}/runs/dump.{TAG[target]}/call.union_kernel_inputs.bin"
        if kind == "dump":
            prof = p1_run.profile_for(cell, "opt", d, timeout, target)
            args = p1_run.genotype_args(cell, prof, f"{d}/call")
            i = args.index("--dosage-mode")
            del args[i:i + 2]
            extra = f"{d}/profile.s2c.tsv"
            open(extra, "w").write("profile\ts2c-kernel\nearly_incumbent\texperimental_v1\ndiagnostic_union_kernel_dump\t1\n")
            args[args.index("--profile") + 2:args.index("--profile") + 2] = ["--profile", extra]
            cmd = ["nice", f"{BUILD}/panvar", "genotype", *args, "--output-schema", "sequence_first_v2"]
        elif kind in ("bulk", "progressive"):
            cmd = ["nice", harness, kind, inputs, cell["reads1"], cell["reads2"], f"{d}/{kind}"]
            if kind == "progressive":
                cmd += ["31", str(4 << 30)]
        else:
            cmd = [harness, "compare", inputs, f"{OUT}/runs/bulk.{TAG[target]}/bulk",
                   f"{OUT}/runs/progressive.{TAG[target]}/progressive", f"{d}/compare.tsv"]
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
