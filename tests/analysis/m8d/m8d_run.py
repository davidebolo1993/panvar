#!/usr/bin/env python3
"""M8-D runner -- continuation: HG01975 LOO r = 0.01 / r = 0.1 with the caller's declared runtime (DESIGN.md; frozen manifests; one process
at a time; 4 threads; 2,100 s counted as time spent plus the declared timeout of the next unit
plus a reserve for every call still pending and for evaluation/analysis).

  m8d_run.py manifest P | next P | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m8d_hg01975_loo_20261020"
M8C = "/Users/davide.bolognini/panvar/tests/results/m8c_copying_challenge_20261020"
BUDGET = 2100
RESERVE = 200   # two evals (60 each) + contrast (60) + summary (20)
WATCHDOG = 930
R = "/Users/davide.bolognini/panvar/tests/results"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
MANIFEST_P = ([("hashcheck", "hashcheck", "-", 30), ("call.r0.01", "call", "r0.01", WATCHDOG), ("call.r0.1", "call", "r0.1", WATCHDOG), ("freeze", "freeze", "-", 10),
               ("eval.r0.01", "eval", "r0.01", 60), ("eval.r0.1", "eval", "r0.1", 60), ("contrast", "contrast", "-", 60), ("summary", "summary", "-", 20)])
MANIFESTS = {"P": MANIFEST_P}
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
    P = f"{M8C}/runs/prepare.a1"
    py = lambda mode: [sys.executable, os.path.join(HERE, "analyse.py"), mode, d]
    if kind == "hashcheck":
        return py("hashcheck"), "frozen-input hash check"
    if kind == "call":
        assert ok("hashcheck")
        reads = f"{P}/reads/HG01975/cov5"
        return (["nice", f"{BUILD}/panvar", "genotype", "--profile", f"{P}/profile.{target}.tsv", "--threads", "4", "--gfa", ACOT_GFA, "--reference-path", ACOT_REF,
                 "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq", "--library", f"{reads}/library.txt", "--out-prefix", f"{d}/call", "--quiet",
                 "--panel-curation", f"{P}/curation.LOO.HG01975.txt", "--dosage-mode", "repair"], "caller (M8-C command, new directory)")
    if kind == "freeze":
        return py("freeze"), "prediction freeze"
    if kind == "eval":
        fz = {l.split("\t")[0]: l.split("\t")[1] for l in open(f"{OUT}/runs/freeze/freeze.tsv").read().splitlines()[1:]}
        assert fz.get(f"HG01975.LOO.{target}") == "COMPLETE", "only completed authoritative outputs are evaluated"
        return ([f"{BUILD}/genotype_probe", "evaluate", "--scheme", "loo", "--donor", "HG01975", "--reference", ACOT_REF,
                 "--run-prefix", f"{OUT}/runs/call.HG01975.LOO.{target}/call", "--model-gfa", ACOT_GFA, "--truth-gfa", ACOT_GFA, "--out", f"{d}/eval",
                 "--panel-curation", f"{P}/curation.LOO.HG01975.txt"], "truth evaluation")
    if kind in ("contrast", "summary"):
        return py(kind), kind
    raise ValueError(kind)


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in reversed(list(MANIFESTS.items())) for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and unit not in done()
    rev, dirty = state()
    assert dirty == 0, "the tree must be committed and clean"
    used = spent()
    verify()
    if kind == "eval":
        fz = {l.split("\t")[0]: l.split("\t")[1] for l in open(f"{OUT}/runs/freeze/freeze.tsv").read().splitlines()[1:]}
        if fz.get(f"HG01975.LOO.{target}") != "COMPLETE":
            record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, "NOT RUN: no completed authoritative output (never evaluated)")
            print(f"{unit}: NOT RUN (incomplete)")
            return 0
    pending_calls = sum(1 for m in MANIFESTS["P"] if m[1] == "call" and m[0] not in done() and m[0] != unit)
    reserve = (RESERVE + WATCHDOG * pending_calls) if kind in ("hashcheck", "call") else (RESERVE if kind in ("freeze", "eval") else 0)
    if used + timeout + reserve > BUDGET:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {used:.1f} s + {timeout} s > {BUDGET} s")
        print(f"{unit}: NOT RUN")
        return 3
    d = f"{OUT}/runs/call.HG01975.LOO.{target}" if kind == "call" else f"{OUT}/runs/{unit}"
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
