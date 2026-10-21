#!/usr/bin/env python3
"""M14-A runner: budget 2,400 s including development, preprocessing, failures and evaluation; 4 threads; one process at a
time; 6 GiB. Caller and evaluator commands are M11-D's (m11d_run.call_cmd / its evaluate arguments), with only the active-panel
curation changed for reduced settings. Call and eval units are admitted only if spent + timeout + the analyse reserve fits.
Usage: m14a_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
import m14a  # noqa: E402
import m11d_run  # noqa: E402
import runner  # noqa: E402

m11d = m14a.m11d
WT, OUT = m14a.WT, m14a.OUT
BUDGET = 2400
PRE = 60.0    # disclosed: probe build (16 s), three development sketches (9 s), one distance timing (1.5 s), reading
ANALYSE_RESERVE = 600
UNITS = [("build", 300)] + [(f"sketch.{c['id']}", 60) for c in m14a.CELLS] + [("cluster", 120)]
for c in m14a.CELLS:
    for s in m14a.SETTINGS:
        UNITS += [(f"call.{c['id']}.{s}", m11d.WATCHDOG[c["locus"]]), (f"eval.{c['id']}.{s}", 120)]
UNITS += [("analyse", ANALYSE_RESERVE)]


def needs(unit):
    if unit.startswith("sketch."):
        return ["build"]
    if unit == "cluster":
        return [f"sketch.{c['id']}" for c in m14a.CELLS]
    if unit.startswith("call."):
        return ["build", "cluster"]
    if unit.startswith("eval."):
        return ["call." + unit.split(".", 1)[1]]
    if unit == "analyse":
        return ["cluster"]
    return []


def cell_of(unit):
    _, loc, donor, *rest = unit.split(".")
    return m11d.cell(f"{loc}.{donor}"), (rest[0] if rest else None)


def command(unit, d):
    py = [sys.executable, os.path.join(HERE, "m14a.py")]
    if unit == "build":
        return ["cmake", "--build", f"{WT}/build", "--target", "panvar", "genotype_probe", "m10c_align_eval", "-j", "8"]
    if unit in ("cluster", "analyse"):
        return py + [unit, d]
    c, setting = cell_of(unit)
    i = m11d.inputs(c)
    if unit.startswith("sketch."):
        cmd = [f"{WT}/build/genotype_probe", "panel-sketch", "--out", f"{d}/panel", "--model-gfa", i["model"], "--reference", i["ref"],
               "--k", str(m14a.K_MER), "--s", str(m14a.SKETCH)]
        return cmd + (["--panel-curation", i["cur"]] if i["cur"] else [])
    cur = m14a.curation_of(c, setting)
    if unit.startswith("call."):
        cmd = m11d_run.call_cmd(c, d)            # M11-D's exact command (unified configuration, diagnostics profile)
        if "--panel-curation" in cmd:
            j = cmd.index("--panel-curation")
            del cmd[j:j + 2]
        if cur:
            j = cmd.index("--reads1")
            cmd[j:j] = ["--panel-curation", cur]
        return cmd
    cmd = [f"{WT}/build/genotype_probe", "evaluate", "--scheme", "loo", "--donor", c["donor"], "--reference", i["ref"],
           "--run-prefix", f"{OUT}/runs/call.{c['id']}.{setting}/call", "--model-gfa", i["model"], "--truth-gfa", i["truth"], "--out", f"{d}/eval"]
    return cmd + (["--panel-curation", cur] if cur else [])


class ReservingRunner(runner.Runner):
    def precondition(self, unit):
        why = super().precondition(unit)
        if why or not unit.startswith(("call.", "eval.")):
            return why
        spent, timeout = self.spent(), dict(self.units)[unit]
        if spent + timeout + ANALYSE_RESERVE > self.budget:
            return f"budget: {spent:.1f} + {timeout} + analyse reserve {ANALYSE_RESERVE} > {self.budget} s"
        return None


R = ReservingRunner(OUT, BUDGET, PRE, UNITS, command, needs, revision=lambda: runner.git_revision(WT), require_clean=runner.git_clean(WT))


def manifest():
    os.makedirs(f"{OUT}/runs", exist_ok=True)
    lines = [f"# executable\t{b}\t{m11d.sha(f'{WT}/build/{b}')}" for b in ("panvar", "genotype_probe", "m10c_align_eval")]
    for c in m14a.CELLS:
        i = m11d.inputs(c)
        lines.append(f"# cell\t{c['id']}\t{c['locus']}\tmodel {i['model']}\tbase curation {i['cur'] or 'none'}\treads1_sha256 {m11d.sha(i['reads'] + '/reads_1.fq')[:16]}")
    lines += ["unit\tdeclared_timeout_seconds\tneeds"] + [f"{u}\t{t}\t{','.join(needs(u)) or '-'}" for u, t in UNITS]
    text = "\n".join(lines) + "\n"
    p = f"{OUT}/MANIFEST.tsv"
    assert not os.path.exists(p) or open(p).read() == text
    open(p, "w").write(text)
    print(text[:2000])


if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        manifest()
    elif sys.argv[1] == "next":
        assert os.path.exists(f"{OUT}/MANIFEST.tsv"), "manifest first"
        sys.exit(R.next())
    else:
        sys.exit(R.run_unit(sys.argv[1]))
