#!/usr/bin/env python3
"""M11-E runner: budget 900 s including gates (disclosed development + time spent + the next unit's declared timeout).
One scoring process at a time, 4 threads, 6 GiB envelope. Uses runner.Runner (exit-code propagation, dependency gating).
Usage: m11e_run.py manifest | next | <unit>.  `next` returns non-zero on any failure or NOT RUN, and 2 when all units
are recorded, so `while python3 m11e_run.py next; do :; done` stops at the first event needing review."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m11e  # noqa: E402
import runner  # noqa: E402

m11d = m11e.m11d
WT, OUT = m11d.WT, m11e.OUT
BUDGET = 900
PRE = 120.0   # disclosed: trace/table inspection, control hash selection, dry runs of fixes/prepare (no scoring)
TIMEOUT = {"C4": 120, "CYP2D6": 120, "ACOT": 240}
UNITS = [("build", 300), ("gate.runner", 60), ("fixes", 120), ("prepare", 120)]
UNITS += [(f"score.{c['case']}", TIMEOUT[m11d.cell(c["cell"])["locus"]]) for c in m11e.CASES]
UNITS += [("trace", 120)]
NEEDS = {"gate.runner": ["build"], "fixes": ["gate.runner"], "prepare": ["gate.runner"], "trace": ["prepare"]}


def needs(unit):
    if unit.startswith("score."):
        return ["build", "prepare"]
    return NEEDS.get(unit, [])


def score_cmd(c, d):
    cell = m11d.cell(c["cell"])
    i = m11d.inputs(cell)
    args = [f"{WT}/build/genotype_probe", "score-set", "--out", f"{d}/ss", "--reference", i["ref"], "--donor", cell["donor"],
            "--run", f"{m11e.M11D}/runs/call.{c['cell']}/call", "--model-gfa", i["model"]]
    if i["cur"]:
        args += ["--panel-curation", i["cur"]]
    args += ["--library", f"{i['reads']}/library.txt", "--reads1", f"{i['reads']}/reads_1.fq", "--reads2", f"{i['reads']}/reads_2.fq",
             "--threads", "4", "--candidates", f"{OUT}/runs/prepare/{c['case']}.candidates.tsv", "--rank-blocks", c["block"],
             "--decompose", "assign_A,incumbent;assign_B,incumbent", "--decompose-block", c["block"]]
    return ["nice", *args]


def command(unit, d):
    py = [sys.executable, os.path.join(HERE, "m11e.py")]
    if unit == "build":
        return ["cmake", "--build", f"{WT}/build", "--target", "genotype_probe", "-j", "8"]
    if unit == "gate.runner":
        return [sys.executable, os.path.join(HERE, "test_runner.py"), d]
    if unit.startswith("score."):
        return score_cmd(m11e.case(unit.split(".", 1)[1]), d)
    return py + [unit, d]


R = runner.Runner(OUT, BUDGET, PRE, UNITS, command, needs, revision=lambda: runner.git_revision(WT), require_clean=runner.git_clean(WT))


def manifest():
    os.makedirs(OUT, exist_ok=True)
    lines = [f"# executable\tgenotype_probe\t{m11d.sha(f'{WT}/build/genotype_probe')}"]
    for c in m11e.CASES:
        cell = m11d.cell(c["cell"])
        i = m11d.inputs(cell)
        lines.append(f"# case\t{c['case']}\t{c['role']}\t{c['cell']}\t{c['block']}\treads1_sha256 {m11d.sha(i['reads'] + '/reads_1.fq')[:16]}"
                     f"\tmodel {i['model']}\tcall {m11e.M11D}/runs/call.{c['cell']}/call")
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
