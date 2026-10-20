#!/usr/bin/env python3
"""M8-B runner -- bounded copying-prior calibration on family-held-out panel haplotypes (masked score; no reads) (DESIGN.md; frozen manifests; one process
at a time; 4 threads; 300 s counted as time spent plus the declared timeout of the next unit).

  m8b_run.py manifest P | next P | <unit>
"""
import csv
import hashlib
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
import p1_run  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/m8b_copying_calibration_20261020"
M6F = "/Users/davide.bolognini/panvar/tests/results/m6f_copying_predictive_20261019"
M1D = "/Users/davide.bolognini/panvar/tests/results/m1d_unit_reachability_20261019/runs"
BUDGET = 300
RESERVE = 30
R = "/Users/davide.bolognini/panvar/tests/results"
LOCI = {"C4": (f"{R}/genotype_c4_production_11fe681/graphs/c4.gfa", "HG00096#1#haplotype1-0000024:31848049-32074186", None),
        "CYP2D6": (f"{R}/cyp2d6_shortread_pilot_v1/data/graphs/cyp2d6.gfa", "grch38#1#chr22:42031864-42245566", "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"),
        "ACOT": (f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", "GRCh38#0#chr14:73416151-73668669", None)}
MANIFEST_P = ([("inventory", "inventory", "-", 30), ("gates", "gates", "-", 30)]
              + [(f"{k}.{L}", k, L, 60) for L in LOCI for k in ("fit", "fit_scrambled")]
              + [(f"support.{L}", "support", L, 45) for L in LOCI] + [("summary", "summary", "-", 30)])
# Amendment 1 (declared after `gates` exited 1 at 15.2 s): the G5 fixture gave the doubled rows one shared name, which panel
# validation refuses; the copies now get unique names. gates.a1 replaces gates; the rest of the frozen order is unchanged.
MANIFEST_A1 = [("gates.a1", "gates", "-", 30)] + [m for m in MANIFEST_P if m[0] not in ("inventory", "gates")]
MANIFESTS = {"P": MANIFEST_P, "A1": MANIFEST_A1}
INPUTS = []


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


def write_manifest(part):
    os.makedirs(OUT, exist_ok=True)
    text = "unit\tkind\ttarget\tdeclared_timeout_seconds\n" + "".join("\t".join(str(x) for x in m) + "\n" for m in MANIFESTS[part])
    p = f"{OUT}/MANIFEST_{part}.tsv"
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


def verify():
    for r in INPUTS:
        assert hashlib.sha256(open(r["file"], "rb").read()).hexdigest() == r["sha256"], f"{r['file']} differs from the frozen input"


def command(kind, target, d):
    ok = lambda u: next((r["exit"] for r in ledger_rows() if r["unit"] == u), None) == "0"
    probe = f"{BUILD}/genotype_probe"
    if kind == "inventory":
        return [sys.executable, os.path.join(HERE, "inventory.py"), d], "corrected inventory and frozen splits"
    if kind == "gates":
        return [probe, "copying-heldout", "--mode", "gates", "--out", f"{d}/h"], "gates"
    if kind in ("fit", "fit_scrambled", "support"):
        assert ok("inventory") and (ok("gates") or ok("gates.a1")), "fitting runs only after the inventory and the gates pass"
        gfa, ref, cur = LOCI[target]
        cmd = [probe, "copying-heldout", "--mode", "support" if kind == "support" else "fit", "--out", f"{d}/h", "--reference", ref, "--model-gfa", gfa,
               "--split", f"{OUT}/runs/inventory/split.{target}.tsv"]
        if cur:
            cmd += ["--panel-curation", cur]
        if kind == "fit_scrambled":
            cmd += ["--scramble-eval", "1"]
        if kind == "support":
            cmd += ["--support-seconds", "25"]
        return cmd, kind
    if kind == "summary":
        return [sys.executable, os.path.join(HERE, "summary.py"), f"{OUT}/runs", f"{d}/m8b"], "summary"
    raise ValueError(kind)


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in reversed(list(MANIFESTS.items())) for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and unit not in done()
    rev, dirty = state()
    assert dirty == 0, "the tree must be committed and clean"
    used = spent()
    verify()
    if used + timeout + (RESERVE if kind != "summary" else 0) > BUDGET:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {used:.1f} s + {timeout} s > {BUDGET} s")
        print(f"{unit}: NOT RUN")
        return 3
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d), f"{d} exists; a run is never repeated into the same place"
    os.makedirs(d)
    cmd, note = command(kind, target, d)
    cmd = ["/usr/bin/time", "-l", *cmd]
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s) " + " ".join(cmd) + "\n")
        log.flush()
        code, timed_out = p1_run.run_timed(cmd, log, timeout)
    wall = time.time() - t
    text = open(f"{d}/run.log").read()
    rss = next((int(l.split()[0]) for l in text.splitlines() if "maximum resident set size" in l), 0)
    record(unit, rev, dirty, timeout, wall, rss / 2**30, code, timed_out, note)
    print("\n".join(l for l in text.splitlines() if not l.startswith(" ") and not l.startswith("$"))[-3500:])
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s (spent {spent():.1f}/{BUDGET})")
    return 0


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest(sys.argv[2])
    elif arg == "next":
        todo = [m[0] for m in MANIFESTS[sys.argv[2]] if m[0] not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(arg))
