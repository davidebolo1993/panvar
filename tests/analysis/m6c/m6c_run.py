#!/usr/bin/env python3
"""M6-C runner -- nuisance-marginalised structural likelihood, synthetic feasibility only (DESIGN.md; frozen manifests; one process
at a time; 4 threads; 600 s counted as time spent plus the declared timeout of the next unit).

  m6c_run.py manifest G|R | next G|R | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m6c_nuisance_marginal_20261019"
BUDGET = 600
FIXTURES = ["F1_control_2+2", "F2a_missing_variants_invisible_allocation", "F2b_missing_variants_informative_allocation",
            "F3_matched_content_true_3+1", "F3_matched_content_true_2+2", "F4_adjacency_true_3+1", "F4_adjacency_true_2+2",
            "F6_dosage_true_2+3", "F7_contamination_2+2"]
SEEDS = [11, 23, 37]   # frozen
MANIFEST_G = [("gates", "gates", "-", 120)]
MANIFEST_R = [(f"run.{f}.{s}", "run", f"{f}:{s}", 60) for f in FIXTURES for s in SEEDS] + [("summary", "summary", "-", 30)]
# Amendment 1 (declared after run.F1_control_2+2.11 at 30.9 s spent): runs take ~21 s, so the frozen fixture-major order could
# exhaust the budget before the last fixtures; the remaining runs go seed-major (every fixture at seed 11, then 23, then 37).
# Same units, timeouts and seeds; only the order changes. The completed unit stays recorded.
MANIFEST_A1 = [(f"run.{f}.{s}", "run", f"{f}:{s}", 60) for s in SEEDS for f in FIXTURES if (f, s) != ("F1_control_2+2", 11)] + [("summary", "summary", "-", 30)]
MANIFESTS = {"G": MANIFEST_G, "R": MANIFEST_R, "A1": MANIFEST_A1}
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
    if kind == "gates":
        return [f"{BUILD}/nuisance_structural", "gates"], "mathematical gates"
    if kind == "run":
        f, seed = target.rsplit(":", 1)
        assert exit_of_gates() == "0", "fixtures run only after the gates pass"
        return [f"{BUILD}/nuisance_structural", "run", f, seed, f"{d}/n"], "synthetic fixture"
    if kind == "summary":
        return [sys.executable, os.path.join(HERE, "summary.py"), f"{OUT}/runs", f"{OUT}/ledger.tsv", f"{d}/m6c"], "summary"
    raise ValueError(kind)


def exit_of_gates():
    return next((r["exit"] for r in ledger_rows() if r["unit"] == "gates"), None)


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in reversed(list(MANIFESTS.items())) for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and unit not in done()
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
