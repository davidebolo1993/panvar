#!/usr/bin/env python3
"""M11-D runner (DESIGN.md; frozen manifest; one caller process at a time, 4 threads, M10-D watchdog = declared timeout +
6 GiB process-tree envelope; the caller's time guard = watchdog - 20 s). Budget 3,600 s = disclosed development + spent +
the next unit's timeout; a unit that does not fit is recorded NOT RUN (watchdogs are never shortened).
  m11d_run.py manifest | next | <unit>
"""
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("m10d", "m10a", "p1"):
    sys.path.insert(0, os.path.join(HERE, "..", p))
sys.path.insert(0, HERE)
import m10d_run  # noqa: E402
import compose_m10a  # noqa: E402
import m11d  # noqa: E402

WT, OUT = m11d.WT, m11d.OUT
BUDGET = 3600
PRE = 60.0   # disclosed: inventory reconciliation, hash selection and path probes, rounded up (no caller run)
PERF = "/Users/davide.bolognini/panvar/tests/results/genotype_perf_bounded_20261015/inputs/profiles"
UNITS = [("build", 300), ("gate.selftest", 60)]
for c in m11d.CELLS:
    UNITS += [(f"call.{c['id']}", m11d.WATCHDOG[c["locus"]]), (f"eval.{c['id']}", 120), (f"truth.{c['id']}", 60)]
UNITS += [("inventory", 120), ("outcomes", 300), ("fit", 120)]


def ledger():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def exit_of(u):
    return next((r["exit"] for r in ledger() if r["unit"] == u), None)


def call_cmd(c, d):
    i = m11d.inputs(c)
    w = m11d.WATCHDOG[c["locus"]]
    prof = compose_m10a.compose(f"{d}/profile.tsv", f"m11d-D-{c['id']}", "D", f"{PERF}/{i['base']}.tsv", "none",
                                [f"time_limit_seconds={w - 20}:declared watchdog minus 20 s (the caller records TIME_LIMIT first)"])
    diag = f"{d}/profile.block_diagnostics.tsv"
    open(diag, "w").write("profile\tm11d-block-diagnostics\nblock_diagnostics\tdiagnostic_v1\n")
    args = ["--preset", "short-read", "--profile", prof, "--threads", "4", "--graph", i["model"], "--reference-path", i["ref"]]
    if i["cur"]:
        args += ["--panel-curation", i["cur"]]
    args += ["--reads1", f"{i['reads']}/reads_1.fq", "--reads2", f"{i['reads']}/reads_2.fq", "--library", f"{i['reads']}/library.txt",
             "--dosage-mode", "audit", "--out", f"{d}/call", "--quiet", "--profile", diag]
    return ["nice", f"{WT}/build/panvar", "genotype", *args]


def eval_cmd(c, d):
    i = m11d.inputs(c)
    cmd = [f"{WT}/build/genotype_probe", "evaluate", "--scheme", "loo", "--donor", c["donor"], "--reference", i["ref"],
           "--run-prefix", f"{OUT}/runs/call.{c['id']}/call", "--model-gfa", i["model"], "--truth-gfa", i["truth"], "--out", f"{d}/eval"]
    return cmd + (["--panel-curation", i["cur"]] if i["cur"] else [])


def command(unit, d):
    py = [sys.executable, os.path.join(HERE, "m11d.py")]
    if unit == "build":
        return ["cmake", "--build", f"{WT}/build", "--target", "panvar", "genotype_probe", "m10c_align_eval", "cyp2d6_block_alleles", "-j", "8"]
    if unit == "gate.selftest":
        return py + ["selftest", d]
    kind, _, cid = unit.partition(".")
    if kind == "call":
        return call_cmd(m11d.cell(cid), d)
    if kind == "eval":
        return eval_cmd(m11d.cell(cid), d)
    if kind == "truth":
        return py + ["truth", cid, d]
    return py + [unit, d]


def precondition(unit):
    if unit != "build" and exit_of("build") != "0":
        return "the build failed"
    if unit != "gate.selftest" and unit != "build" and exit_of("gate.selftest") != "0":
        return "the fitting gate failed"
    kind, _, cid = unit.partition(".")
    if kind == "eval" and not os.path.exists(f"{OUT}/runs/call.{cid}/call.sequence_locus.tsv"):
        return "no call outputs to evaluate"
    if kind == "truth" and exit_of(f"eval.{cid}") != "0":
        return "no evaluation"
    return None


def run_unit(unit):
    timeout = dict(UNITS)[unit]
    assert os.path.exists(f"{OUT}/MANIFEST.tsv") and unit not in {r["unit"] for r in ledger()}
    rev = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    assert not subprocess.run(["git", "-C", WT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip(), "commit first"
    spent = PRE + sum(float(r["wall_seconds"]) for r in ledger())
    why = precondition(unit) or (f"{spent:.1f} + {timeout} > {BUDGET} s" if spent + timeout > BUDGET else None)
    p = f"{OUT}/ledger.tsv"
    if not os.path.exists(p):
        open(p, "w").write("unit\trevision\tdeclared_timeout\twall_seconds\texit\twatchdog\tpeak_gib\tnote\n")
    if why:
        open(p, "a").write(f"{unit}\t{rev}\t{timeout}\t0.00\tnot_run\t-\t0\tNOT RUN: {why}\n")
        print(f"{unit}: NOT RUN ({why})")
        return 0
    d = f"{OUT}/runs/{unit}"
    os.makedirs(d)
    cmd = command(unit, d)
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s; envelope 6 GiB) " + " ".join(cmd) + "\n")
        log.flush()
        code, verdict, peak, wall = m10d_run.run_watched(cmd, log, timeout)
    open(p, "a").write(f"{unit}\t{rev}\t{timeout}\t{wall:.2f}\t{code}\t{verdict}\t{peak / 2**30:.3f}\t-\n")
    tail = [l for l in open(f"{d}/run.log").read().splitlines() if not l.startswith((" ", "$"))]
    print("\n".join(tail[-12:]))
    print(f"{unit}: exit {code} ({verdict}) {wall:.1f} s, peak {peak / 2**30:.2f} GiB; spent {spent + wall:.1f}/{BUDGET}")
    return 0


def manifest():
    os.makedirs(OUT, exist_ok=True)
    lines = ["unit\tdeclared_timeout_seconds"] + [f"{u}\t{t}" for u, t in UNITS]
    bins = [f"# executable\t{b}\t{m11d.sha(f'{WT}/build/{b}')}" for b in ("panvar", "genotype_probe", "m10c_align_eval", "cyp2d6_block_alleles")]
    cells = []
    for c in m11d.CELLS:
        i = m11d.inputs(c)
        files = [i["model"], i["truth"], f"{i['reads']}/reads_1.fq", f"{i['reads']}/reads_2.fq", f"{i['reads']}/library.txt"] + ([i["cur"]] if i["cur"] else [])
        missing = [f for f in files if not os.path.exists(f)]
        assert not missing, f"{c['id']}: missing inputs {missing}"
        cells.append(f"# cell\t{c['id']}\t{c['locus']}\t{c['donor']}\tfamily {m11d.FAMILY[c['donor']]}\treads1_sha256 {m11d.sha(files[2])[:16]}\tmodel {i['model']}")
    text = "\n".join(bins + cells + lines) + "\n"
    p = f"{OUT}/MANIFEST.tsv"
    assert not os.path.exists(p) or open(p).read() == text
    open(p, "w").write(text)
    print(text)


if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        manifest()
    elif sys.argv[1] == "next":
        todo = [u for u, _ in UNITS if u not in {r["unit"] for r in ledger()}]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(sys.argv[1]))
