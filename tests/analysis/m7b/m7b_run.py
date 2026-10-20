#!/usr/bin/env python3
"""M7-B runner -- corrected fixed-pool comparison, Panvar vs Locityper (frozen order; one process at a time) (DESIGN.md; frozen manifests; one process
at a time; 4 threads; 300 s counted as time spent plus the declared timeout of the next unit).

  m7b_run.py manifest P | next P | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m7b_fixed_pool_20261020"
M6F = "/Users/davide.bolognini/panvar/tests/results/m6f_copying_predictive_20261019"
M1D = "/Users/davide.bolognini/panvar/tests/results/m1d_unit_reachability_20261019/runs"
BUDGET = 900
RESERVE = 40   # compare + reporting
CELLS = ["c4_lzo_hg00171_5", "cyp_lzo_hg02769_5", "c4_loo_na19650_10", "cyp_loo_na19650_5"]
MANIFEST_P = ([("provenance", "provenance", "-", 120), ("background", "background", "-", 60),
               ("reference.C4", "reference", "C4", 30), ("reference.CYP2D6", "reference", "CYP2D6", 30)]
              + [(f"preproc.{c}", "preproc", c, 45) for c in CELLS] + [("pools", "pools", "-", 90)]
              + [(f"panvar.D.{c}", "panvar", f"{c}:D", 45) for c in CELLS]
              + [("panvar.X.cyp_loo_na19650_5", "panvar", "cyp_loo_na19650_5:X", 30), ("panvar.O.c4_loo_na19650_10", "panvar", "c4_loo_na19650_10:O", 120)]
              + [(f"loc.{c}.D.s1", "loc", f"{c}:D:20261020", 60) for c in CELLS]
              + [("loc.cyp_loo_na19650_5.DX.s1", "loc", "cyp_loo_na19650_5:DX:20261020", 60), ("loc.c4_loo_na19650_10.DO.s1", "loc", "c4_loo_na19650_10:DO:20261020", 90)]
              + [("loc.c4_lzo_hg00171_5.D.s1r", "loc", "c4_lzo_hg00171_5:D:20261020", 60), ("loc.c4_loo_na19650_10.D.s2", "loc", "c4_loo_na19650_10:D:20261021", 60),
                 ("loc.cyp_loo_na19650_5.D.s2", "loc", "cyp_loo_na19650_5:D:20261021", 60)]
              + [("truthdist", "truthdist", "-", 90), ("compare", "compare", "-", 30)])
# Amendment 1 (declared after `background` exited 1 at 8.2 s): lambda compared at simulation.tsv's printed precision (6 significant
# digits; relative 1e-5) instead of 1e-6. Every other check passed. background.a1 replaces background; downstream paths read background.a1.
MANIFEST_A1 = [("background.a1", "background", "-", 60)] + [m for m in MANIFEST_P if m[0] not in ("provenance", "background")]
# Amendment 2 (declared after `truthdist` exited 1 at 63.4 s): the M6-A summary's second section (block_excess) was not parsed;
# truthdist.a2 and compare.a2 replace truthdist and compare.
MANIFEST_A2 = [("truthdist.a2", "truthdist", "-", 90), ("compare.a2", "compare", "-", 30)]
MANIFESTS = {"P": MANIFEST_P, "A1": MANIFEST_A1, "A2": MANIFEST_A2}
INPUTS = [{"file": "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references/tests/analysis/m7a/CANDIDATES.tsv",
           "sha256": "a8cd9516798a4795" }]   # prefix-checked in verify(); full bindings are re-verified by pools/provenance via M7-A BINDINGS.tsv


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
        assert hashlib.sha256(open(r["file"], "rb").read()).hexdigest().startswith(r["sha256"]), f"{r['file']} differs from the frozen input"
    for r in csv.DictReader(open("/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references/tests/analysis/m7a/BINDINGS.tsv"), delimiter="\t"):
        assert hashlib.sha256(open(r["file"], "rb").read()).hexdigest() == r["sha256"], f"{r['file']} differs from its M7-A binding"


def command(kind, target, d):
    py = lambda script, *a: [sys.executable, os.path.join(HERE, script), *a]
    P = f"{OUT}/runs/pools"
    if kind == "provenance":
        return py("provenance.py", d), "provenance"
    if kind == "background":
        return py("calib.py", "background", d), "background"
    if kind == "reference":
        return py("calib.py", "reference", target, d), "reference"
    if kind == "preproc":
        return py("calib.py", "preproc", target, d), "preproc"
    if kind == "pools":
        return py("pools.py", d), "pools"
    if kind == "panvar":
        cell, tier = target.split(":")
        sys.path.insert(0, HERE)
        from m7b_common import CELLS as CC, LOCI, SCORING_PROBE
        c = CC[cell]
        cmd = [SCORING_PROBE, "score-set", "--out", f"{d}/fs", "--model-gfa", c["model"], "--reference", LOCI[c["locus"]]["ref_path"], "--donor", c["donor"],
               "--run", c["run"], "--reads1", f"{c['reads']}/reads_1.fq", "--reads2", f"{c['reads']}/reads_2.fq", "--library", f"{c['reads']}/library.txt", "--threads", "4"]
        if c["curation"]:
            cmd += ["--panel-curation", c["curation"]]
        cmd += ["--explicit", f"{P}/{cell}.panvar.O.tsv"] if tier == "O" else ["--candidates", f"{P}/{cell}.panvar.{tier}.tsv"]
        return cmd, "panvar score-set"
    if kind == "loc":
        cell, pool, seed = target.split(":")
        return py("loc.py", cell, pool, seed, d), "locityper"
    if kind == "truthdist":
        return py("truthdist.py", d), "truth distance (evaluation only)"
    if kind == "compare":
        return py("compare.py", d), "compare"
    raise ValueError(kind)


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in reversed(list(MANIFESTS.items())) for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and unit not in done()
    rev, dirty = state()
    assert dirty == 0, "the tree must be committed and clean"
    used = spent()
    verify()
    if used + timeout + (RESERVE if kind != "compare" else 0) > BUDGET:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {used:.1f} s + {timeout} s > {BUDGET} s")
        print(f"{unit}: NOT RUN")
        return 3
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d), f"{d} exists; a run is never repeated into the same place"
    os.makedirs(d)
    cmd, note = command(kind, target, d)
    cmd = ["nice", "/usr/bin/time", "-l", *cmd]
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
