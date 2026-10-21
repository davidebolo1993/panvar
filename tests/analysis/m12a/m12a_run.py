#!/usr/bin/env python3
"""M12-A runner: budget 900 s (disclosed development + time spent + the next unit's declared timeout); 6 GiB; one scoring
process at a time; 4 threads. Uses the M11-E runner (exit-code propagation, dependency gating).
Usage: m12a_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
import m12a  # noqa: E402
import runner  # noqa: E402

m11d = m12a.m11d
WT, OUT = m11d.WT, m12a.OUT
BUDGET = 900
PRE = 150.0   # disclosed: code tracing, provenance reads, probe and fixture builds (about 45 s compute), read-length scans
LOCUS_TIMEOUT = {"C4": 120, "CYP2D6": 120, "ACOT": 240}
UNITS = [("build", 300), ("gate.fixtures", 120), ("contract", 120), ("prepare", 60)]
UNITS += [(f"lik.{c['case']}", LOCUS_TIMEOUT[m11d.cell(c["cell"])["locus"]]) for c in m12a.CASES]
UNITS += [(f"count.{c['case']}", 240 if m11d.cell(c["cell"])["locus"] == "ACOT" else 180) for c in m12a.CASES if c["count"] == "run"]
UNITS += [("analyse", 120)]


def needs(unit):
    if unit == "gate.fixtures":
        return ["build"]
    if unit in ("contract", "prepare"):
        return []
    if unit.startswith(("lik.", "count.")):
        return ["build", "prepare"]
    if unit == "analyse":
        return ["contract", "prepare"]
    return []


def lik_cmd(c, d):
    cell = m11d.cell(c["cell"])
    i = m11d.inputs(cell)
    p = f"{OUT}/runs/prepare/{c['case']}"
    args = [f"{WT}/build/genotype_probe", "score-set", "--out", f"{d}/ss", "--reference", i["ref"], "--donor", cell["donor"],
            "--run", f"{m12a.M11D}/runs/call.{c['cell']}/call", "--model-gfa", i["model"]]
    if i["cur"]:
        args += ["--panel-curation", i["cur"]]
    args += ["--library", f"{i['reads']}/library.txt", "--reads1", f"{i['reads']}/reads_1.fq", "--reads2", f"{i['reads']}/reads_2.fq",
             "--threads", "4", "--candidates", f"{p}.candidates.tsv", "--explicit", f"{p}.explicit.tsv",
             "--decompose", "assign_A,incumbent;truth_A,incumbent;truth_B,incumbent", "--decompose-block", c["block"]]
    return ["nice", *args]


def count_cmd(c, d):
    cell = m11d.cell(c["cell"])
    i = m11d.inputs(cell)
    p = f"{OUT}/runs/prepare/{c['case']}"
    args = [f"{WT}/build/genotype_probe", "count-attrib", "--out", f"{d}/ca", "--run", f"{m12a.M11D}/runs/call.{c['cell']}/call",
            "--block", c["block"], "--model-gfa", i["model"], "--reference", i["ref"]]
    if i["cur"]:
        args += ["--panel-curation", i["cur"]]
    args += ["--reads1", f"{i['reads']}/reads_1.fq", "--reads2", f"{i['reads']}/reads_2.fq", "--candidates", f"{p}.count_candidates.tsv",
             "--block-fasta", f"{p}.block.fa", "--contrasts", "best_A,incumbent;truth_A,incumbent", "--truth-full", f"{p}.truth_full.fa",
             "--skip-catalogue", "1"]
    return ["nice", *args]


def command(unit, d):
    py = [sys.executable, os.path.join(HERE, "m12a.py")]
    if unit == "build":
        return ["cmake", "--build", f"{WT}/build", "--target", "genotype_probe", "m12a_copy_balance", "-j", "8"]
    if unit == "gate.fixtures":
        return [f"{WT}/build/m12a_copy_balance"]
    if unit.startswith("lik."):
        return lik_cmd(m12a.case(unit.split(".", 1)[1]), d)
    if unit.startswith("count."):
        return count_cmd(m12a.case(unit.split(".", 1)[1]), d)
    return py + [unit, d]


R = runner.Runner(OUT, BUDGET, PRE, UNITS, command, needs, revision=lambda: runner.git_revision(WT), require_clean=runner.git_clean(WT))


def manifest():
    os.makedirs(OUT, exist_ok=True)
    lines = [f"# executable\t{b}\t{m11d.sha(f'{WT}/build/{b}')}" for b in ("genotype_probe", "m12a_copy_balance")]
    for c in m12a.CASES:
        cell = m11d.cell(c["cell"])
        p = m12a.provenance(c["cell"])
        lines.append(f"# case\t{c['case']}\t{c['role']}\t{c['cell']}\t{c['block']}\treads1_sha256 {m11d.sha(p['reads'] + '/reads_1.fq')[:16]}"
                     f"\tsimulation {p['sim']}\tcount {c['count']}")
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
