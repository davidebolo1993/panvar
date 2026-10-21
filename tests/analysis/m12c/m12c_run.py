#!/usr/bin/env python3
"""M12-C runner: budget 900 s; 6 GiB; one scoring process at a time. Score units reuse M12-B's frozen candidate lists and its
exact score-set command (m12b_run.score_set); each is admitted only if spent + its M12-B timeout + the evaluate reserve fits.
Usage: m12c_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "m12b"))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
import m12c  # noqa: E402
import m12b_run  # noqa: E402
import runner  # noqa: E402

m12b = m12c.m12b
WT, OUT = m12b.WT, m12c.OUT
BUDGET = 900
PRE = 60.0    # disclosed: binding checks by hand, script writing (no scoring)
EVALUATE_RESERVE = 300
PLAN = {r["bid"]: r for r in m12b.T(f"{m12c.M12B}/runs/select/score_plan.tsv")}
UNITS = [("bind", 60)] + [(f"score.{b}", int(PLAN[b]["timeout_s"])) for b in m12c.CONTINUE] + [("evaluate", EVALUATE_RESERVE)]


def needs(unit):
    if unit.startswith("score."):
        return ["bind"]
    if unit == "evaluate":
        return ["bind"]
    return []


def command(unit, d):
    py = [sys.executable, os.path.join(HERE, "m12c.py")]
    if unit in ("bind", "evaluate"):
        return py + [unit, d]
    r = PLAN[unit.split(".", 1)[1]]
    return m12b_run.score_set(r["cell"], r["block"], f"{m12c.M12B}/runs/select/{r['bid']}.candidates.tsv", d)


class ReservingRunner(runner.Runner):
    def precondition(self, unit):
        why = super().precondition(unit)
        if why or not unit.startswith("score."):
            return why
        spent, timeout = self.spent(), dict(self.units)[unit]
        if spent + timeout + EVALUATE_RESERVE > self.budget:
            return f"budget: {spent:.1f} + {timeout} + evaluate reserve {EVALUATE_RESERVE} > {self.budget} s"
        return None


R = ReservingRunner(OUT, BUDGET, PRE, UNITS, command, needs, revision=lambda: runner.git_revision(WT), require_clean=runner.git_clean(WT))


def manifest():
    os.makedirs(f"{OUT}/runs", exist_ok=True)
    lines = [f"# executable\t{b}\t{m12b.m11d.sha(f'{WT}/build/{b}')}" for b in ("genotype_probe", "m10c_align_eval")]
    for b in m12c.CONTINUE:
        p = f"{m12c.M12B}/runs/select/{b}.candidates.tsv"
        lines.append(f"# continue\t{b}\t{PLAN[b]['cell']}\t{PLAN[b]['block']}\tcandidates_sha256 {m12c.sha(p)}\testimate {PLAN[b]['estimate_s']}\ttimeout {PLAN[b]['timeout_s']}")
    lines += ["unit\tdeclared_timeout_seconds\tneeds"] + [f"{u}\t{t}\t{','.join(needs(u)) or '-'}" for u, t in UNITS]
    text = "\n".join(lines) + "\n"
    p = f"{OUT}/MANIFEST.tsv"
    assert not os.path.exists(p) or open(p).read() == text
    open(p, "w").write(text)
    print(text)


if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        manifest()
    elif sys.argv[1] == "next":
        assert os.path.exists(f"{OUT}/MANIFEST.tsv"), "manifest first"
        sys.exit(R.next())
    else:
        sys.exit(R.run_unit(sys.argv[1]))
