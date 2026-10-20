#!/usr/bin/env python3
"""M10-D runner (DESIGN.md; frozen manifest; one process at a time). Budget 900 s = disclosed pre-freeze development + spent
+ the declared timeout of the next unit. Every unit runs under a frozen watchdog: its declared wall-clock timeout and a 6 GiB
memory envelope on the summed resident size of its whole process tree (polled every 0.2 s; exceeding it kills the tree and
is recorded as ENVELOPE_EXCEEDED). Failed, killed or unsupported units stay in the ledger.

  m10d_run.py manifest P | next P | <unit>
"""
import csv
import os
import signal
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m10d  # noqa: E402

WT, OUT = m10d.WT, m10d.OUT
BUDGET = 900
ENVELOPE = 6 << 30
PRE_FREEZE = 10.0   # disclosed: the --distance mode build (2.0 s) and configure checks, rounded up; no fixture or case run
TRUTHS = list(m10d.m10c.TRUTH)
MANIFESTS = {"P": [("build", "build", "-", 120), ("gate.contract", "py", "gate", 120)] +
             [(f"truth.{c}", "truth", c, 120) for c in TRUTHS] + [("cases", "py", "cases", 400)],
             "A1": [("cases.a1", "py", "cases", 400)],   # Amendment 1: provenance-column harness defect; same cases, same truth units
             "A2": [("cases.a2", "py", "cases", 400)]}   # Amendment 2: framed-prefix rows for interior-unsupported cases


def ledger_rows():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def spent():
    return PRE_FREEZE + sum(float(r["wall_seconds"]) for r in ledger_rows())


def done():
    return {r["unit"] for r in ledger_rows()}


def exit_of(u):
    return next((r["exit"] for r in ledger_rows() if r["unit"] == u), None)


def tree_rss(root):
    out = subprocess.run(["ps", "-A", "-o", "pid=,ppid=,rss="], capture_output=True, text=True).stdout.split("\n")
    kids, rss = {}, {}
    for line in out:
        f = line.split()
        if len(f) == 3:
            kids.setdefault(int(f[1]), []).append(int(f[0]))
            rss[int(f[0])] = int(f[2]) * 1024
    total, stack = 0, [root]
    while stack:
        p = stack.pop()
        total += rss.get(p, 0)
        stack += kids.get(p, [])
    return total


def run_watched(cmd, log, timeout):
    p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
    t0, peak, verdict = time.time(), 0, "ok"
    while p.poll() is None:
        peak = max(peak, tree_rss(p.pid))
        if peak > ENVELOPE:
            os.killpg(p.pid, signal.SIGKILL)
            verdict = "ENVELOPE_EXCEEDED"
            break
        if time.time() - t0 > timeout:
            os.killpg(p.pid, signal.SIGTERM)
            time.sleep(3)
            if p.poll() is None:
                os.killpg(p.pid, signal.SIGKILL)
            verdict = "TIMEOUT"
            break
        time.sleep(0.2)
    p.wait()
    return p.returncode, verdict, peak, time.time() - t0


def record(unit, rev, dirty, timeout, wall, peak, code, verdict, note):
    p = f"{OUT}/ledger.tsv"
    if not os.path.exists(p):
        open(p, "w").write("unit\trevision\tdirty_tree\tdeclared_timeout_seconds\twall_seconds\tpeak_tree_rss_gib\texit\twatchdog\tnote\n")
    with open(p, "a") as f:
        f.write(f"{unit}\t{rev}\t{dirty}\t{timeout}\t{wall:.2f}\t{peak / 2**30:.3f}\t{code}\t{verdict}\t{note}\n")


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in MANIFESTS.items() for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and unit not in done()
    rev = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = len([l for l in subprocess.run(["git", "-C", WT, "status", "--porcelain"], capture_output=True, text=True).stdout.split("\n") if l.strip()])
    assert dirty == 0, "the tree must be committed and clean"
    why = None
    if kind != "build" and exit_of("build") != "0":
        why = "the build failed"
    elif kind in ("truth", "py") and unit != "gate.contract" and exit_of("gate.contract") != "0":
        why = "the contract gate failed: stop before any case"
    elif unit == "cases" and any(exit_of(f"truth.{c}") != "0" for c in TRUTHS):
        why = None   # cases run anyway: a missing truth unit makes its comparisons UNSUPPORTED (kept in the table)
    if why is None and spent() + timeout > BUDGET:
        why = f"{spent():.1f} + {timeout} > {BUDGET} s"
    if why:
        record(unit, rev, dirty, timeout, 0.0, 0, "not_run", "-", f"NOT RUN: {why}")
        print(f"{unit}: NOT RUN ({why})")
        return 0
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d)
    os.makedirs(d)
    py = [sys.executable, os.path.join(HERE, "m10d.py")]
    cmd = {"build": ["cmake", "--build", f"{WT}/build", "--target", "m10c_align_eval", "genotype_probe", "-j", "8"],
           "truth": py + ["truth", target, d], "py": py + [target, d]}[kind]
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s; envelope 6 GiB) " + " ".join(cmd) + "\n")
        log.flush()
        code, verdict, peak, wall = run_watched(cmd, log, timeout)
    record(unit, rev, dirty, timeout, wall, peak, code, verdict, f"{kind} {target}")
    text = open(f"{d}/run.log").read()
    print("\n".join(l for l in text.splitlines() if not l.startswith("$"))[-4000:])
    print(f"{unit}: exit {code} ({verdict}), {wall:.1f} s, peak tree RSS {peak / 2**30:.2f} GiB (spent {spent():.1f}/{BUDGET})")
    return 0


if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        os.makedirs(OUT, exist_ok=True)
        text = "unit\tkind\ttarget\tdeclared_timeout_seconds\n" + "".join("\t".join(map(str, m)) + "\n" for m in MANIFESTS[sys.argv[2]])
        p = f"{OUT}/MANIFEST_{sys.argv[2]}.tsv"
        assert not os.path.exists(p) or open(p).read() == text
        open(p, "w").write(text)
    elif sys.argv[1] == "next":
        todo = [m[0] for m in MANIFESTS[sys.argv[2]] if m[0] not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(sys.argv[1]))
