#!/usr/bin/env python3
"""M11-C runner (DESIGN.md; frozen manifest; one process at a time; M10-D watchdog: declared timeout + 6 GiB process-tree
envelope). Budget 900 s = disclosed development + spent + the next unit's timeout.
  m11c_run.py manifest | next | <unit>
"""
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
for p in ("m10d", "m10a", "p1"):
    sys.path.insert(0, os.path.join(HERE, "..", p))
import m10d_run  # noqa: E402  (watchdog)
import m10a_run  # noqa: E402  (M10-A command builder: genotype_args, cells, profile composition)
import compose_m10a  # noqa: E402
import m10a  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = "/Users/davide.bolognini/panvar/tests/results/m11c_confidence_inputs_20261020"
BUDGET = 900
PRE = 150.0   # disclosed: development builds of the diagnostic (about 12 s failed + the full rebuild), path probes, rounded up
M1B_CYP5 = "/Users/davide.bolognini/panvar/tests/results/m1b_structural_20261019/runs/call.cyp_5"
UNITS = [("build", 300), ("gate.junction_unit", 30), ("gate.compound_unit", 30), ("gate.sequence_report", 30), ("gate.e2e_sequence_first", 240),
         ("smoke.c4.off", 90), ("smoke.c4.on", 90), ("smoke.c4.on.r", 90), ("smoke.c4.dosageoff.on", 90),
         ("smoke.acot.off", 200), ("smoke.acot.on", 200), ("smoke.cyp.off", 150), ("smoke.cyp.on", 150),
         ("cmp.c4", 60), ("cmp.acot", 60), ("cmp.cyp", 60), ("diag_check", 30), ("truth", 300), ("qvfrac", 300)]


def ledger():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def exit_of(u):
    return next((r["exit"] for r in ledger() if r["unit"] == u), None)


def diag_profile(d):
    p = f"{d}/profile.block_diagnostics.tsv"
    open(p, "w").write("profile\tm11c-block-diagnostics\nblock_diagnostics\tdiagnostic_v1\n")
    return p


def d_config(cell_id, d, on, dosage_off=False):
    cell = m10a.cells()[cell_id]
    t = m10a_run.TIMEOUT["D"][cell_id]
    prof = compose_m10a.compose(f"{d}/profile.{cell_id}.tsv", f"m10a-D-{cell_id}", "D", cell["base_profile"],
                                f"{m10a.OUT}/runs/declare/{cell_id.split('_')[0]}.declarations.tsv",
                                [f"time_limit_seconds={t - 20}:declared unit timeout minus 20 s (records TIME_LIMIT before the kill)"])
    args = m10a_run.genotype_args(cell, prof, f"{d}/call")
    if dosage_off:
        args[args.index("--dosage-mode") + 1] = "off"
    if on:
        args += ["--profile", diag_profile(d)]
    return ["nice", f"{WT}/build/panvar", "genotype", *args]


def m1b_cyp(d, on):
    first = open(f"{M1B_CYP5}/run.log").readline().split()
    cmd = first[first.index("-l") + 1:]
    cmd[cmd.index("--out") + 1] = f"{d}/call"
    if on:
        cmd += ["--profile", diag_profile(d)]
    return cmd


def command(unit, d):
    py = [sys.executable, os.path.join(HERE, "m11c.py")]
    runs = f"{OUT}/runs"
    return {
        "build": ["cmake", "--build", f"{WT}/build", "--target", "panvar", "genotype_probe", "junction_proposals_test", "compound_annotation_test",
                  "genotype_sequence_report_test", "genotype_command_e2e_test", "m10c_align_eval", "-j", "8"],
        "gate.junction_unit": [f"{WT}/build/junction_proposals_test"], "gate.compound_unit": [f"{WT}/build/compound_annotation_test"],
        "gate.sequence_report": [f"{WT}/build/genotype_sequence_report_test", d],
        "gate.e2e_sequence_first": [f"{WT}/build/genotype_command_e2e_test", d, WT],
        "smoke.c4.off": lambda: d_config("c4_5", d, False), "smoke.c4.on": lambda: d_config("c4_5", d, True),
        "smoke.c4.on.r": lambda: d_config("c4_5", d, True), "smoke.c4.dosageoff.on": lambda: d_config("c4_5", d, True, dosage_off=True),
        "smoke.acot.off": lambda: d_config("acot_5", d, False), "smoke.acot.on": lambda: d_config("acot_5", d, True),
        "smoke.cyp.off": lambda: m1b_cyp(d, False), "smoke.cyp.on": lambda: m1b_cyp(d, True),
        "cmp.c4": py + ["compare", f"{runs}/smoke.c4.off", f"{runs}/smoke.c4.on", d],
        "cmp.acot": py + ["compare", f"{runs}/smoke.acot.off", f"{runs}/smoke.acot.on", d],
        "cmp.cyp": py + ["compare", f"{runs}/smoke.cyp.off", f"{runs}/smoke.cyp.on", d],
        "diag_check": py + ["diag_check", d], "truth": py + ["truth", d], "qvfrac": py + ["qvfrac", d]}[unit]


def run_unit(unit):
    timeout = dict(UNITS)[unit]
    assert os.path.exists(f"{OUT}/MANIFEST.tsv") and unit not in {r["unit"] for r in ledger()}
    rev = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    assert not subprocess.run(["git", "-C", WT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip(), "commit first"
    spent = PRE + sum(float(r["wall_seconds"]) for r in ledger())
    why = None
    if unit != "build" and exit_of("build") != "0":
        why = "the build failed"
    elif unit.startswith("smoke.") and any(exit_of(g) != "0" for g in ("gate.junction_unit", "gate.compound_unit", "gate.sequence_report", "gate.e2e_sequence_first")):
        why = "a gate failed: stop before smoke runs"
    elif spent + timeout > BUDGET:
        why = f"{spent:.1f} + {timeout} > {BUDGET} s"
    p = f"{OUT}/ledger.tsv"
    if not os.path.exists(p):
        open(p, "w").write("unit\trevision\tdeclared_timeout\twall_seconds\texit\twatchdog\tpeak_gib\tnote\n")
    if why:
        open(p, "a").write(f"{unit}\t{rev}\t{timeout}\t0.00\tnot_run\t-\t0\tNOT RUN: {why}\n")
        print(f"{unit}: NOT RUN ({why})")
        return 0
    d = f"{OUT}/runs/{unit}"
    os.makedirs(d)
    c = command(unit, d)
    cmd = c() if callable(c) else c
    env = {"PANVAR_E2E_ONLY": "sequence_first"} if unit == "gate.e2e_sequence_first" else None
    if env:
        os.environ.update(env)
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s; envelope 6 GiB) " + " ".join(cmd) + "\n")
        log.flush()
        code, verdict, peak, wall = m10d_run.run_watched(cmd, log, timeout)
    for k in (env or {}):
        os.environ.pop(k, None)
    open(p, "a").write(f"{unit}\t{rev}\t{timeout}\t{wall:.2f}\t{code}\t{verdict}\t{peak / 2**30:.3f}\t-\n")
    print("\n".join(l for l in open(f"{d}/run.log").read().splitlines() if not l.startswith((" ", "$")))[-2500:])
    print(f"{unit}: exit {code} ({verdict}) {wall:.1f} s, peak {peak / 2**30:.2f} GiB; spent {spent + wall:.1f}/{BUDGET}")
    return 0


if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        os.makedirs(OUT, exist_ok=True)
        text = "unit\tdeclared_timeout_seconds\n" + "".join(f"{u}\t{t}\n" for u, t in UNITS)
        assert not os.path.exists(f"{OUT}/MANIFEST.tsv") or open(f"{OUT}/MANIFEST.tsv").read() == text
        open(f"{OUT}/MANIFEST.tsv", "w").write(text)
    elif sys.argv[1] == "next":
        todo = [u for u, _ in UNITS if u not in {r["unit"] for r in ledger()}]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(sys.argv[1]))
