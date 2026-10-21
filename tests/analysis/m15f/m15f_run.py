#!/usr/bin/env python3
"""M15-F runner: budget 600 s; read-only units over existing code facts, synthetic facts and saved development outputs.
Units: fidelity (the Python port of headline_filters_v2 reproduces every recorded headline), synthetic (the decision
table), saved (current -> proposed on representative saved runs; originals hashed before and after).
Usage: m15f_run.py next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import runner  # noqa: E402
import m11d  # noqa: E402

WT = m11d.WT
OUT = f"{m11d.R}/m15f_filter_contract_20261021"
# disclosed development: one scan of saved outputs for demonstration cases (15 s), a few file inspections, three dev
# runs of the synthetic table (<1 s each); the two code inventories were read-only code reading (no execution)
BUDGET, PRE = 600, 30.0
UNITS = [("fidelity", 240), ("synthetic", 30), ("saved", 120)]


def command(unit, d):
    return [sys.executable, os.path.join(HERE, "contract.py"), unit, d]


R = runner.Runner(OUT, BUDGET, PRE, UNITS, command, lambda u: [], revision=lambda: runner.git_revision(WT),
                  require_clean=runner.git_clean(WT))

if __name__ == "__main__":
    sys.exit(R.next() if sys.argv[1] == "next" else R.run_unit(sys.argv[1]))
