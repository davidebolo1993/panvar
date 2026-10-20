#!/usr/bin/env python3
"""M10-B runner (DESIGN.md; frozen manifest; one process at a time; 1,200 s = disclosed pre-freeze development + spent +
the declared timeout of the next unit).

  m10b_run.py manifest P | next P | <unit>
"""
import csv
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
sys.path.insert(0, HERE)
import p1_run  # noqa: E402
import m10b  # noqa: E402

WT = m10b.WT
OUT = m10b.OUT
BUDGET = 1200
PRE_FREEZE = 45.0   # disclosed: m10b_walk_dump development build (22.9 s), four exploratory dumps (7.2 s), walk/skeleton
#                     invariance probes and fixture runs (about 10 s), rounded up
MANIFESTS = {"P": [("build", "build", "-", 300), ("gate.fixtures", "gate", "-", 30)] +
             [(f"dump.{n}", "dump", n, 120 if n.startswith("acot") else 60) for n in m10b.MODELS] +
             [("coverage", "py", "coverage", 300), ("cases", "py", "cases", 300)]}


def ledger_rows():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def spent():
    return PRE_FREEZE + sum(float(r["wall_seconds"]) for r in ledger_rows())


def done():
    return {r["unit"] for r in ledger_rows()}


def exit_of(u):
    return next((r["exit"] for r in ledger_rows() if r["unit"] == u), None)


def record(unit, rev, dirty, timeout, wall, code, timed_out, note):
    p = f"{OUT}/ledger.tsv"
    if not os.path.exists(p):
        open(p, "w").write("unit\trevision\tdirty_tree\tdeclared_timeout_seconds\twall_seconds\texit\ttimed_out\tnote\n")
    with open(p, "a") as f:
        f.write(f"{unit}\t{rev}\t{dirty}\t{timeout}\t{wall:.2f}\t{code}\t{int(bool(timed_out))}\t{note}\n")


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in MANIFESTS.items() for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and unit not in done()
    rev = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = len([l for l in subprocess.run(["git", "-C", WT, "status", "--porcelain"], capture_output=True, text=True).stdout.split("\n") if l.strip()])
    assert dirty == 0, "the tree must be committed and clean"
    why = None
    if kind != "build" and exit_of("build") != "0":
        why = "the build failed"
    elif kind in ("dump", "py") and exit_of("gate.fixtures") != "0":
        why = "the fixtures failed: stop before any real case"
    elif kind == "py" and any(exit_of(f"dump.{n}") != "0" for n in m10b.MODELS):
        why = "a model dump failed"
    elif spent() + timeout > BUDGET:
        why = f"{spent():.1f} + {timeout} > {BUDGET} s"
    if why:
        record(unit, rev, dirty, timeout, 0.0, "not_run", 0, f"NOT RUN: {why}")
        print(f"{unit}: NOT RUN ({why})")
        return 0
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d)
    os.makedirs(d)
    py = [sys.executable]
    cmd = {"build": ["cmake", "--build", f"{WT}/build", "--target", "m10b_walk_dump", "-j", "8"],
           "gate": py + [os.path.join(HERE, "test_project.py")],
           "dump": py + [os.path.join(HERE, "m10b.py"), "dump", target, d],
           "py": py + [os.path.join(HERE, "m10b.py"), target, d]}[kind]
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s) " + " ".join(cmd) + "\n")
        log.flush()
        code, timed_out = p1_run.run_timed(["/usr/bin/time", "-l", *cmd], log, timeout)
    wall = time.time() - t
    record(unit, rev, dirty, timeout, wall, code, timed_out, f"{kind} {target}")
    text = open(f"{d}/run.log").read()
    print("\n".join(l for l in text.splitlines() if not l.startswith(" ") and not l.startswith("$"))[-3000:])
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s (spent {spent():.1f}/{BUDGET})")
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
