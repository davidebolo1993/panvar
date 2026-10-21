#!/usr/bin/env python3
"""Regression for the M11-D runner defect (M11-E gate.runner).

M11-D: `outcomes` exited 1 under watchdog verdict "ok"; the dependent `fit` was admitted because its precondition never
checked the outcomes unit, and the operator loop (which tested the verdict, not the exit code) did not stop.
This test shows, on throwaway units in <dir>, that the M11-E runner (a) records the failure, (b) returns a non-zero
status from `next`, (c) refuses every dependent unit while the failed unit's watchdog verdict is "ok", (d) treats a
watchdog timeout the same way, and (e) still runs independent units. It also reads the sealed M11-D ledger to show the
old precondition admitted `fit` after the failed `outcomes` (read-only)."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import runner  # noqa: E402

PY = sys.executable


def main(d):
    out = os.path.join(d, "sandbox")
    os.makedirs(out)
    cmds = {"analysis": [PY, "-c", "import sys; print('analysis fails'); sys.exit(1)"],
            "dependent": [PY, "-c", "print('must not run')"],
            "independent": [PY, "-c", "print('independent ok')"],
            "slow": [PY, "-c", "import time; time.sleep(30)"],
            "after_slow": [PY, "-c", "print('must not run')"]}
    units = [("analysis", 20), ("dependent", 20), ("independent", 20), ("slow", 2), ("after_slow", 20)]
    needs = {"dependent": ["analysis"], "after_slow": ["slow"]}
    r = runner.Runner(out, 600, 0.0, units, lambda u, _d: cmds[u], lambda u: needs.get(u, []))
    checks = []
    s1 = r.next()
    row = r.ledger()[0]
    checks.append(("failed analysis: exit recorded 1, watchdog verdict ok", row["exit"] == "1" and row["watchdog"] == "ok"))
    checks.append(("failed analysis: next returns non-zero (stops a `while next` loop)", s1 == 1))
    s2 = r.next()
    row = r.ledger()[1]
    checks.append(("dependent of a failed unit is NOT RUN, naming it", s2 == runner.NOT_RUN and row["exit"] == "not_run"
                   and "analysis" in row["note"] and not os.path.exists(f"{out}/runs/dependent")))
    s3 = r.next()
    checks.append(("independent unit still runs and passes", s3 == 0 and r.passed("independent")))
    s4 = r.next()
    row = r.ledger()[3]
    checks.append(("watchdog timeout: verdict TIMEOUT, next non-zero", row["watchdog"] == "TIMEOUT" and s4 != 0))
    s5 = r.next()
    checks.append(("dependent of a timed-out unit is NOT RUN", s5 == runner.NOT_RUN and "slow" in r.ledger()[4]["note"]))
    checks.append(("all recorded: next returns DONE (2), also stopping the loop", r.next() == runner.DONE))
    # the sealed M11-D defect, read-only: the old precondition admitted fit although outcomes exited 1
    sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
    import m11d_run  # noqa: E402
    old = {row["unit"]: row for row in m11d_run.ledger()}
    checks.append(("M11-D ledger: outcomes exit 1 with verdict ok, yet the old precondition admitted fit",
                   old["outcomes"]["exit"] == "1" and old["outcomes"]["watchdog"] == "ok" and m11d_run.precondition("fit") is None
                   and old["fit"]["exit"] == "0"))
    for name, ok in checks:
        print(("PASS " if ok else "FAIL ") + name)
    return 0 if all(ok for _, ok in checks) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
