"""Budgeted unit runner with dependency gating and exit-code propagation (M11-E).

Corrects the M11-D runner defect: a unit counts as passed only when its process exited 0 AND the watchdog verdict is
"ok". A dependent unit runs only when every unit it needs has passed; otherwise it is recorded NOT RUN, naming the
failed dependency. `next` returns the unit's status so a shell loop (`while python3 run.py next; do :; done`) stops on
the first failure:
  0 = passed, 1.. = the unit's own non-zero exit, 3 = NOT RUN, 4 = watchdog verdict not ok with exit 0, 2 = all recorded.
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "m10d"))
import m10d_run  # noqa: E402

DONE, NOT_RUN, VERDICT_FAILED = 2, 3, 4
HEADER = "unit\trevision\tdeclared_timeout\twall_seconds\texit\twatchdog\tpeak_gib\tnote\n"


class Runner:
    def __init__(self, out, budget, pre, units, command, needs, revision=None, require_clean=None):
        """units: [(name, timeout)] in order; command(unit, dir) -> argv; needs(unit) -> [unit names];
        revision() -> str; require_clean() -> None or raises."""
        self.out, self.budget, self.pre, self.units = out, budget, pre, units
        self.command, self.needs = command, needs
        self.revision = revision or (lambda: "-")
        self.require_clean = require_clean or (lambda: None)

    def ledger(self):
        p = f"{self.out}/ledger.tsv"
        if not os.path.exists(p):
            return []
        rows = [l.rstrip("\n").split("\t") for l in open(p)][1:]
        return [dict(zip(HEADER.strip().split("\t"), r)) for r in rows]

    def passed(self, unit):
        r = next((r for r in self.ledger() if r["unit"] == unit), None)
        return r is not None and r["exit"] == "0" and r["watchdog"] == "ok"

    def spent(self):
        return self.pre + sum(float(r["wall_seconds"]) for r in self.ledger())

    def precondition(self, unit):
        failed = [u for u in self.needs(unit) if not self.passed(u)]
        return f"needs {','.join(failed)} (not passed)" if failed else None

    def _append(self, line):
        p = f"{self.out}/ledger.tsv"
        if not os.path.exists(p):
            open(p, "w").write(HEADER)
        open(p, "a").write(line)

    def run_unit(self, unit):
        timeout = dict(self.units)[unit]
        assert unit not in {r["unit"] for r in self.ledger()}, f"{unit} already recorded"
        self.require_clean()
        rev, spent = self.revision(), self.spent()
        why = self.precondition(unit) or (f"{spent:.1f} + {timeout} > {self.budget} s" if spent + timeout > self.budget else None)
        if why:
            self._append(f"{unit}\t{rev}\t{timeout}\t0.00\tnot_run\t-\t0\tNOT RUN: {why}\n")
            print(f"{unit}: NOT RUN ({why})")
            return NOT_RUN
        d = f"{self.out}/runs/{unit}"
        os.makedirs(d)
        cmd = self.command(unit, d)
        with open(f"{d}/run.log", "w") as log:
            log.write(f"$ (declared timeout {timeout} s; envelope 6 GiB) " + " ".join(cmd) + "\n")
            log.flush()
            code, verdict, peak, wall = m10d_run.run_watched(cmd, log, timeout)
        self._append(f"{unit}\t{rev}\t{timeout}\t{wall:.2f}\t{code}\t{verdict}\t{peak / 2**30:.3f}\t-\n")
        tail = [l for l in open(f"{d}/run.log").read().splitlines() if not l.startswith((" ", "$"))]
        print("\n".join(tail[-12:]))
        print(f"{unit}: exit {code} ({verdict}) {wall:.1f} s, peak {peak / 2**30:.2f} GiB; spent {spent + wall:.1f}/{self.budget}")
        if code != 0:
            return code if 0 < code < 256 and code not in (DONE, NOT_RUN, VERDICT_FAILED) else 1
        return 0 if verdict == "ok" else VERDICT_FAILED

    def next(self):
        todo = [u for u, _ in self.units if u not in {r["unit"] for r in self.ledger()}]
        if not todo:
            print("all units recorded")
            return DONE
        return self.run_unit(todo[0])


def git_revision(wt):
    return subprocess.run(["git", "-C", wt, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()


def git_clean(wt):
    def check():
        dirty = subprocess.run(["git", "-C", wt, "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
        assert not dirty, "commit first"
    return check
