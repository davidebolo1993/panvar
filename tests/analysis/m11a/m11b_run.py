#!/usr/bin/env python3
"""M11-B runner: two units (manifest, audit), one process at a time, under M10-D's watchdog (declared timeout, 6 GiB
process-tree envelope). Budget 300 s = disclosed development + spent + the next unit's timeout.
  m11b_run.py <unit>
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m10d"))
import m10d_run  # noqa: E402  (watchdog only)

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = "/Users/davide.bolognini/panvar/tests/results/m11a_confidence_design_20261020/m11b"
BUDGET, PRE = 300, 15.0   # disclosed: M11-A inventory run and path inspection, rounded up
UNITS = {"manifest": 60, "audit": 120}


def main(unit):
    os.makedirs(f"{OUT}/runs", exist_ok=True)
    led = f"{OUT}/ledger.tsv"
    rows = [l.split("\t") for l in open(led).read().splitlines()[1:]] if os.path.exists(led) else []
    assert unit not in {r[0] for r in rows}
    spent = PRE + sum(float(r[3]) for r in rows)
    rev = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", WT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
    assert not dirty, "commit first"
    if unit == "audit" and not any(r[0] == "manifest" and r[4] == "0" for r in rows):
        sys.exit("manifest missing")
    if spent + UNITS[unit] > BUDGET:
        sys.exit(f"budget: {spent} + {UNITS[unit]} > {BUDGET}")
    d = f"{OUT}/runs/{unit}"
    os.makedirs(d)
    with open(f"{d}/run.log", "w") as log:
        code, verdict, peak, wall = m10d_run.run_watched([sys.executable, os.path.join(HERE, "m11b.py"), unit, d], log, UNITS[unit])
    if not os.path.exists(led):
        open(led, "w").write("unit\trevision\tdeclared_timeout\twall_seconds\texit\twatchdog\tpeak_gib\n")
    open(led, "a").write(f"{unit}\t{rev}\t{UNITS[unit]}\t{wall:.2f}\t{code}\t{verdict}\t{peak / 2**30:.3f}\n")
    print(open(f"{d}/run.log").read()[-6000:])
    print(f"{unit}: exit {code} ({verdict}) {wall:.1f} s; spent {spent + wall:.1f}/{BUDGET}")


if __name__ == "__main__":
    main(sys.argv[1])
