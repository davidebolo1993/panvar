#!/usr/bin/env python3
"""M8-C runner -- paired ACOT genotyping challenge r = 0.01 vs r = 0.1 (two donors, LZO/LOO, 5x) (DESIGN.md; frozen manifests; one process
at a time; 4 threads; 300 s counted as time spent plus the declared timeout of the next unit).

  m8c_run.py manifest P | next P | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m8c_copying_challenge_20261020"
M6F = "/Users/davide.bolognini/panvar/tests/results/m6f_copying_predictive_20261019"
M1D = "/Users/davide.bolognini/panvar/tests/results/m1d_unit_reachability_20261019/runs"
BUDGET = 1200
RESERVE = 60   # contrast + summary
R = "/Users/davide.bolognini/panvar/tests/results"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
CALLS = [(i, sch, r) for i in ("d1", "d2") for sch in ("LZO", "LOO") for r in ("r0.01", "r0.1")]
EVALS = [(i, sch, r) for sch in ("LOO", "LZO") for i in ("d1", "d2") for r in ("r0.01", "r0.1")]
MANIFEST_P = ([("prepare", "prepare", "-", 90)] + [(f"call.{i}.{s}.{r}", "call", f"{i}:{s}:{r}", 120) for i, s, r in CALLS]
              + [("checks", "checks", "-", 30)] + [(f"eval.{i}.{s}.{r}", "eval", f"{i}:{s}:{r}", 150) for i, s, r in EVALS]
              + [("contrast", "contrast", "-", 120), ("summary", "summary", "-", 30)])
# Amendment 1 (declared after `prepare` exited 1 at 3.4 s): subsample does not create its output directory; prepare now creates
# it. prepare.a1 replaces prepare (its directory is used by every later unit); the frozen order is otherwise unchanged.
MANIFEST_A1 = [("prepare.a1", "prepare", "-", 90)] + [m for m in MANIFEST_P if m[0] != "prepare"]
# Amendment 2 (declared after `contrast` and `summary` exited 1 at 827.1 s): the analysis read the timed-out runs' missing outputs;
# missing outputs are now absent results. contrast.a2 and summary.a2 replace them.
MANIFEST_A2 = [("contrast.a2", "contrast", "-", 120), ("summary.a2", "summary", "-", 30)]
MANIFESTS = {"P": MANIFEST_P, "A1": MANIFEST_A1, "A2": MANIFEST_A2}
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


def donor(i):
    return [l.split("\t")[1] for l in open(f"{OUT}/runs/prepare.a1/prepare.tsv") if l.startswith("donor_")][int(i[1]) - 1]


def unit_name(kind, target):   # concrete donor ids in run directories
    i, s, r = target.split(":")
    return f"{kind}.{donor(i)}.{s}.{r}"


def command(kind, target, d):
    ok = lambda u: next((r["exit"] for r in ledger_rows() if r["unit"] == u), None) == "0"
    P = f"{OUT}/runs/prepare.a1"
    if kind == "prepare":
        return [sys.executable, os.path.join(HERE, "prepare.py"), d], "prepare (donors, panels, reads, profiles)"
    if kind == "call":
        assert ok("prepare.a1")
        i, s, r = target.split(":")
        dn = donor(i)
        reads = f"{P}/reads/{dn}/cov5"
        return (["nice", f"{BUILD}/panvar", "genotype", "--profile", f"{P}/profile.{r}.tsv", "--threads", "4", "--gfa", ACOT_GFA, "--reference-path", ACOT_REF,
                 "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq", "--library", f"{reads}/library.txt", "--out-prefix", f"{d}/call", "--quiet",
                 "--panel-curation", f"{P}/curation.{s}.{dn}.txt", "--dosage-mode", "repair"], "caller")
    if kind == "checks":
        return [sys.executable, os.path.join(HERE, "analyse.py"), "checks", f"{OUT}/runs", d], "arm checks"
    if kind == "eval":
        i, s, r = target.split(":")
        dn = donor(i)
        return ([f"{BUILD}/genotype_probe", "evaluate", "--scheme", s.lower(), "--donor", dn, "--reference", ACOT_REF, "--run-prefix", f"{OUT}/runs/call.{dn}.{s}.{r}/call",
                 "--model-gfa", ACOT_GFA, "--truth-gfa", ACOT_GFA, "--out", f"{d}/eval", "--panel-curation", f"{P}/curation.{s}.{dn}.txt"], "truth evaluation")
    if kind in ("contrast", "summary"):
        return [sys.executable, os.path.join(HERE, "analyse.py"), kind, f"{OUT}/runs", d], kind
    raise ValueError(kind)


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in reversed(list(MANIFESTS.items())) for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and unit not in done()
    rev, dirty = state()
    assert dirty == 0, "the tree must be committed and clean"
    used = spent()
    verify()
    if used + timeout + (RESERVE if kind not in ("contrast", "summary") else 0) > BUDGET:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {used:.1f} s + {timeout} s > {BUDGET} s")
        print(f"{unit}: NOT RUN")
        return 3
    d = f"{OUT}/runs/{unit_name(kind, target) if kind in ('call', 'eval') else unit}"
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
