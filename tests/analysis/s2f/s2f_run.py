#!/usr/bin/env python3
"""S2-F validation runner (frozen manifest; one process at a time; <= 4 threads; <= 300 s per unit; 1200 s budget,
counted as time spent in this batch -- development runs included -- plus the declared timeout of the next unit; a unit
that does not fit is NOT_RUN).

Order: the combined gate (genotype_blocked_fine) and the indexed regressions (blocked union, progressive kernel, union
placement cache, full e2e) -> per 5x cell (C4, CYP2D6, ACOT; S2-C's union kernel inputs): ordinary blocked and
blocked+fine (4,096 read pairs, this revision), compared bit for bit with each other and with S2-E's bulk scores ->
ACOT 30x blocked+fine with the ACTUAL 30x pool and unchanged caps (S2-C dump.acot30) under the harness's 4 GiB guard ->
CONDITIONAL: only if that run's union completes, one ordinary-blocked ACOT 30x run with a prospectively declared
1,000,000,000 physical-origin limit -- an extended-resource validation run, not the main arm.

  s2f_run.py manifest | next | <unit> | dev <name> <note> -- <command...>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/s2f_blocked_fine_aggregates_20261018"
S2C = "/Users/davide.bolognini/panvar/tests/results/s2c_demand_driven_20261018/runs"
S2E = "/Users/davide.bolognini/panvar/tests/results/s2e_fragment_blocked_20261018/runs"
BUDGET = 1200
EXTENDED_ORIGIN_CAP = 1_000_000_000   # declared before any S2-F measurement; the conditional control only
MANIFEST = [("t.genotype_blocked_fine", "ctest", "genotype_blocked_fine", 120),
            ("t.genotype_blocked_union", "ctest", "genotype_blocked_union", 60),
            ("t.genotype_progressive_kernel", "ctest", "genotype_progressive_kernel", 60),
            ("t.genotype_union_placement_cache", "ctest", "genotype_union_placement_cache", 120),
            ("t.genotype_command_e2e", "ctest", "genotype_command_e2e", 300)]
for cell, tag, t in (("c4_5", "c4", 60), ("cyp_5", "cyp", 60), ("acot_5", "acot5", 120)):
    MANIFEST += [(f"blocked.{tag}", "blocked", cell, t), (f"blocked_fine.{tag}", "blocked_fine", cell, t),
                 (f"compare.{tag}", "compare", cell, 30), (f"compare_s2e_bulk.{tag}", "compare_s2e_bulk", cell, 30)]
MANIFEST += [("blocked_fine.acot30", "blocked_fine", "acot_30", 300),
             ("control_blocked_1b.acot30", "control_blocked_1b", "acot_30", 300)]
ORDER = [m[0] for m in MANIFEST]
TAG = {"c4_5": "c4", "cyp_5": "cyp", "acot_5": "acot5", "acot_30": "acot30"}


def ledger_rows():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def spent():
    return sum(float(r["wall_seconds"]) for r in ledger_rows())


def done():
    return {r["unit"] for r in ledger_rows()}


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
    os.makedirs(OUT, exist_ok=True)
    p = f"{OUT}/ledger.tsv"
    if not os.path.exists(p):
        open(p, "w").write("unit\trevision\tdirty_tree\tdeclared_timeout_seconds\twall_seconds\tpeak_rss_gib\texit\ttimed_out\tnote\n")
    with open(p, "a") as f:
        f.write(f"{unit}\t{rev}\t{dirty}\t{timeout}\t{wall:.2f}\t{rss:.3f}\t{code}\t{int(bool(timed_out))}\t{note}\n")


def timed(cmd, d, timeout):
    cmd = ["/usr/bin/time", "-l", *cmd]
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s) " + " ".join(cmd) + "\n")
        log.flush()
        code, timed_out = p1_run.run_timed(cmd, log, timeout)
    wall = time.time() - t
    text = open(f"{d}/run.log").read()
    rss = 0
    for line in text.splitlines():
        if "maximum resident set size" in line:
            rss = int(line.split()[0])
    return code, timed_out, wall, rss, text


def tsv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path):
            k, _, v = line.rstrip("\n").partition("\t")
            out[k] = v
    return out


def run_unit(unit):
    u, kind, target, timeout = next(m for m in MANIFEST if m[0] == unit)
    assert open(f"{OUT}/MANIFEST.tsv").read() and timeout <= 300 and unit not in done()
    rev, dirty = state()
    assert dirty == 0, "the tree must be committed and clean"
    if kind == "control_blocked_1b":   # conditional on the combined ACOT 30x union completing
        combined = tsv(f"{OUT}/runs/blocked_fine.acot30/blocked_fine.tsv")
        if combined.get("completed") != "yes":
            record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0,
                   f"NOT RUN (condition): the combined ACOT 30x union did not complete ({combined.get('refusal', 'no result')[:80]})")
            print(f"{unit}: NOT RUN (condition)")
            return 3
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
        tag = TAG[target]
        inputs = f"{S2C}/dump.{tag}/call.union_kernel_inputs.bin"
        if kind in ("blocked", "blocked_fine"):
            cmd = ["nice", harness, kind, inputs, cell["reads1"], cell["reads2"], f"{d}/{kind}", "4096", str(4 << 30)]
        elif kind == "control_blocked_1b":
            cmd = ["nice", harness, "blocked", inputs, cell["reads1"], cell["reads2"], f"{d}/blocked", "4096", str(4 << 30),
                   str(EXTENDED_ORIGIN_CAP)]
        elif kind == "compare":
            cmd = [harness, "compare", inputs, f"{OUT}/runs/blocked.{tag}/blocked", f"{OUT}/runs/blocked_fine.{tag}/blocked_fine",
                   f"{d}/compare.tsv"]
        else:
            cmd = [harness, "compare", inputs, f"{S2E}/bulk.{tag}/bulk", f"{OUT}/runs/blocked_fine.{tag}/blocked_fine", f"{d}/compare.tsv"]
    code, timed_out, wall, rss, text = timed(cmd, d, timeout)
    note = "-"
    if kind == "ctest":
        m = re.search(r"Test +#\d+: \S+ \.+\s*\**\s*(\w[\w ]*?)\s+([\d.]+) sec", text)
        note = f"ctest: {m.group(1)} {m.group(2)} s" if m else "ctest: no result line"
    elif kind == "control_blocked_1b":
        note = f"EXTENDED-RESOURCE VALIDATION RUN: physical-origin limit {EXTENDED_ORIGIN_CAP} (declared before execution)"
    record(unit, rev, dirty, timeout, wall, rss / 2**30, code, timed_out, note)
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s, {note} (spent {spent():.1f}/{BUDGET})")
    return 0


def dev(name, note, cmd):
    """A development run before the tools are committed: counted against the budget, recorded, never a result."""
    used = spent()
    assert used + 300 <= BUDGET, f"development run refused: {used:.1f} s + 300 s > {BUDGET} s"
    d = f"{OUT}/dev/{name}"
    assert not os.path.exists(d), f"{d} exists"
    os.makedirs(d)
    code, timed_out, wall, rss, text = timed(cmd, d, 300)
    record(name, "uncommitted-development", "-", "-", wall, rss / 2**30, code, timed_out,
           f"development run before the tools were committed (counted): {note}")
    print(text[-3000:])
    print(f"{name}: exit {code}, {wall:.1f} s (spent {spent():.1f}/{BUDGET})")
    return code


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest()
    elif arg == "next":
        todo = [u for u in ORDER if u not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    elif arg == "dev":
        sep = sys.argv.index("--")
        sys.exit(dev(sys.argv[2], sys.argv[3], sys.argv[sep + 1:]))
    else:
        sys.exit(run_unit(arg))
