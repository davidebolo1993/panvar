#!/usr/bin/env python3
"""M12-B runner: budget 1,800 s including gates and disclosed development; 6 GiB; one scoring process at a time; 4 threads.
Uses the M11-E runner (exit-code propagation, dependency gating). The score.<block> units come from select's frozen plan
(score_plan.tsv: estimate and declared timeout per admitted block). Usage: m12b_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
import m12b  # noqa: E402
import runner  # noqa: E402

m11d = m12b.m11d
WT, OUT = m11d.WT, m12b.OUT
BUDGET = 1800
PRE = 200.0   # disclosed: code reading, two probe builds (about 50 s compute), selftest iterations, manifest checks (no scoring)
BLOCKS = m12b.blocks_manifest()
STATIC = [("build", 300), ("gate.selftest", 60), ("preflight_candidates", 60),
          ("preflight.CYP2D6", 180), ("preflight.C4", 180), ("preflight.ACOT", 300), ("gate.order_reversed", 180), ("gate.order", 60)]
STATIC += [(f"count.{b['bid']}", 300 if b["locus"] == "ACOT" else 120) for b in BLOCKS]
STATIC += [("gate.features", 120), ("select", 120)]
# amendment A1 (AMENDMENT_A1.md): gate.features was dispatched as mode "gate.features" (the mode is "gate_features") and failed
# before running; re-run as gate.features.a1, which the score units now need. The failed unit stays recorded.
STATIC.insert(STATIC.index(("select", 120)), ("gate.features.a1", 120))


def plan():
    p = f"{OUT}/runs/select/score_plan.tsv"
    return m12b.T(p) if os.path.exists(p) else []


def units():
    return STATIC + [(r["unit"], int(r["timeout_s"])) for r in plan()] + [("evaluate", 300)]


def needs(unit):
    if unit in ("gate.selftest", "preflight_candidates"):
        return []
    if unit.startswith("preflight.") or unit == "gate.order_reversed":
        return ["build", "preflight_candidates"]
    if unit == "gate.order":
        return ["preflight.CYP2D6", "gate.order_reversed"]
    if unit.startswith("count."):
        return ["build"]
    if unit in ("gate.features", "gate.features.a1", "select"):
        return ["gate.selftest"]
    if unit.startswith("score."):
        return ["build", "gate.selftest", "gate.order", "gate.features.a1", "select"]
    if unit == "evaluate":
        return ["select"]
    return []


def score_set(cell_id, block, candidates, d):
    cell = m11d.cell(cell_id)
    i = m11d.inputs(cell)
    args = [f"{WT}/build/genotype_probe", "score-set", "--out", f"{d}/ss", "--reference", i["ref"], "--donor", cell["donor"],
            "--run", f"{m12b.M11D}/runs/call.{cell_id}/call", "--model-gfa", i["model"]]
    if i["cur"]:
        args += ["--panel-curation", i["cur"]]
    args += ["--library", f"{i['reads']}/library.txt", "--reads1", f"{i['reads']}/reads_1.fq", "--reads2", f"{i['reads']}/reads_2.fq",
             "--threads", "4", "--candidates", candidates]
    return ["nice", *args]


def command(unit, d):
    py = [sys.executable, os.path.join(HERE, "m12b.py")]
    if unit == "build":
        return ["cmake", "--build", f"{WT}/build", "--target", "genotype_probe", "m10c_align_eval", "-j", "8"]
    if unit == "gate.selftest":
        return py + ["selftest"]
    if unit in ("preflight_candidates", "gate.features", "select", "evaluate"):
        return py + [unit, d]
    if unit == "gate.order":
        return py + ["gate_order", d]
    if unit == "gate.features.a1":
        return py + ["gate_features", d]
    pc = f"{OUT}/runs/preflight_candidates"
    if unit.startswith("preflight."):
        locus = unit.split(".", 1)[1]
        b = m12b.json.loads(open(f"{pc}/{locus}.block").read())
        return score_set(b["cell"], b["block"], f"{pc}/{locus}.candidates.tsv", d)
    if unit == "gate.order_reversed":
        b = m12b.json.loads(open(f"{pc}/CYP2D6.block").read())
        return score_set(b["cell"], b["block"], f"{pc}/CYP2D6.reversed.tsv", d)
    if unit.startswith("count."):
        b = next(x for x in BLOCKS if x["bid"] == unit.split(".", 1)[1])
        cell = m11d.cell(b["cell"])
        i = m11d.inputs(cell)
        win = m12b.winner(b["cell"])
        k = [x["block_id"] for x in m12b.m11e.eval_rows(b["cell"])].index(b["block"])
        h1, h2 = m12b.m11e.vec(win["hap1_alleles"]), m12b.m11e.vec(win["hap2_alleles"])
        open(f"{d}/incumbent.tsv", "w").write(f"incumbent\talleles\t{h1[k]}\t{h2[k]}\t0\n")
        args = [f"{WT}/build/genotype_probe", "count-attrib", "--out", f"{d}/ca", "--run", f"{m12b.M11D}/runs/call.{b['cell']}/call",
                "--block", b["block"], "--model-gfa", i["model"], "--reference", i["ref"]]
        if i["cur"]:
            args += ["--panel-curation", i["cur"]]
        args += ["--reads1", f"{i['reads']}/reads_1.fq", "--reads2", f"{i['reads']}/reads_2.fq", "--candidates", f"{d}/incumbent.tsv",
                 "--features", "1"]
        return ["nice", *args]
    if unit.startswith("score."):
        r = next(x for x in plan() if x["unit"] == unit)
        return score_set(r["cell"], r["block"], f"{OUT}/runs/select/{r['bid']}.candidates.tsv", d)
    raise ValueError(unit)


def make_runner():
    return runner.Runner(OUT, BUDGET, PRE, units(), command, needs, revision=lambda: runner.git_revision(WT), require_clean=runner.git_clean(WT))


def manifest():
    os.makedirs(OUT, exist_ok=True)
    lines = [f"# executable\t{b}\t{m11d.sha(f'{WT}/build/{b}')}" for b in ("genotype_probe", "m10c_align_eval")]
    for b in BLOCKS:
        i = m11d.inputs(m11d.cell(b["cell"]))
        lines.append(f"# block\t{b['bid']}\t{b['cell']}\t{b['block']}\tM11-D excess {b['orig_excess']:g}\treads1_sha256 {m11d.sha(i['reads'] + '/reads_1.fq')[:16]}")
    lines += ["unit\tdeclared_timeout_seconds\tneeds"] + [f"{u}\t{t}\t{','.join(needs(u)) or '-'}" for u, t in STATIC + [("evaluate", 300)]]
    lines.append("# score.<block> units: from runs/select/score_plan.tsv (frozen rule: admitted iff estimate <= 240 s; timeout = ceil(1.5 x estimate + 60))")
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
        sys.exit(make_runner().next())
    else:
        sys.exit(make_runner().run_unit(sys.argv[1]))
