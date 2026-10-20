#!/usr/bin/env python3
"""M10-C runner (DESIGN.md; frozen manifest; one process at a time; 900 s = disclosed pre-freeze development + spent +
the declared timeout of the next unit).

  m10c_run.py manifest P | next P | <unit>
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
import m10c  # noqa: E402


WT = m10c.WT
OUT = m10c.OUT
BUDGET = 900
PRE_FREEZE = 30.0   # disclosed: m10c_align_eval development builds (22.9 + 1.3 s) and configure, rounded up; no fixture or case was run
TRUTHS = list(m10c.TRUTH)
MANIFESTS = {"P": [("build", "build", "-", 300), ("gate.aligner", "gate", "-", 60), ("fixtures", "py", "fixtures", 120)] +
             [(f"truth.{c}", "truth", c, 120) for c in TRUTHS] + [("cases", "py", "cases", 300)],
             # Amendment 1 (after `fixtures`: method A cannot locate events): method B on the fixtures (+F10c); cases only if
             # method B meets every structural and concordance fixture (m10c.B_REQUIRED)
             "A1": [("build.a1", "build", "-", 300), ("fixtures_b.a1", "py", "fixtures_b", 120)] +
             [(f"truth.{c}.a1", "truth", c, 120) for c in TRUTHS] + [("cases.a1", "py", "cases", 300)]}
MANIFESTS["A2"] = [("build.a2", "build", "-", 300), ("fixtures_c.a2", "py", "fixtures_c", 120)] + \
    [(f"truth.{c}.a2", "truth", c, 120) for c in TRUTHS] + [("cases.a2", "py", "cases", 300)]   # Amendment 2: method C
AMENDMENT_DEV = 35.0   # disclosed: Amendments 1-2 development (minimap2-mode builds incl. two failed links 26 s; abPOA-mode build 1.2 s), rounded up


def ledger_rows():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def spent():
    return PRE_FREEZE + AMENDMENT_DEV + sum(float(r["wall_seconds"]) for r in ledger_rows())


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
    elif kind in ("truth", "py") and exit_of("gate.aligner") != "0":
        why = "the aligner gate failed"
    elif unit == "cases":
        why = "superseded by Amendment 1: method A cannot locate events (fixtures), so cases run only under method B (cases.a1)"
    elif unit.startswith("truth.") and not unit.endswith(".a1"):
        why = "superseded by Amendment 1 (truth extraction runs as .a1, only if method B is supported)"
    elif unit.endswith(".a1") and (unit.startswith("truth.") or unit == "cases.a1") and \
            open(f"{OUT}/runs/fixtures_b.a1/method_b_supported.txt").read().strip() != "True":
        why = "method B does not meet every structural and concordance fixture: no saved case is inspected"
    elif unit == "cases.a1" and any(exit_of(f"truth.{c}.a1") != "0" for c in TRUTHS):
        why = "a truth extraction failed"
    elif unit.endswith(".a2") and (unit.startswith("truth.") or unit == "cases.a2") and \
            open(f"{OUT}/runs/fixtures_c.a2/method_c_supported.txt").read().strip() == "NONE":
        why = "method C does not meet the structural and concordance fixtures (not even inversion-blind): no saved case is inspected"
    elif unit == "cases.a2" and any(exit_of(f"truth.{c}.a2") != "0" for c in TRUTHS):
        why = "a truth extraction failed"
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
    cmd = {"build": ["cmake", "--build", f"{WT}/build", "--target", "m10c_align_eval", "genotype_probe", "-j", "8"],
           "gate": py + [os.path.join(HERE, "m10c.py"), "gate", d],
           "truth": py + [os.path.join(HERE, "m10c.py"), "truth", target, d],
           "py": py + [os.path.join(HERE, "m10c.py"), target, d]}[kind]
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
