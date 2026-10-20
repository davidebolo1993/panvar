#!/usr/bin/env python3
"""M9-G runner -- default-off caller integration of junction-aware proposals: gates, then paired off/on development runs
(DESIGN.md; frozen manifest; one process at a time; 2,400 s counted as the disclosed pre-freeze development plus time spent
plus the declared timeout of the next unit plus, for builds/gates/calls, the declared timeouts of every pending
freeze/eval/analyse unit; refused units are recorded NOT RUN).

  m9g_run.py manifest P | next P | <unit>
"""
import csv
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
import p1_run  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/m9g_junction_integration_20261020"
M9C = "/Users/davide.bolognini/panvar/tests/results/m9c_acot_reliability_20261020"
R = "/Users/davide.bolognini/panvar/tests/results"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
BUDGET = 2400
PRE_FREEZE = 45.0   # disclosed: development compile (34.0 s + 1.5 s) and one unit-test run, rounded up
F8 = ["LOO.d1", "LOO.d3", "LOO.d6", "LZO.d1", "LZO.d2", "LOO.d2", "LOO.d4", "LOO.d5"]
CELLS = ["LOO.d1", "LOO.d3", "LOO.d6", "LZO.d1", "LZO.d2"]
WATCHDOG = 930   # caller guard 900 s + 30 s
MANIFEST_P = ([("build", "build", "-", 300), ("prepare", "prepare", "-", 10), ("gate.unit", "gate_unit", "-", 30),
               ("gate.sequence_report", "gate_seq", "-", 30), ("gate.e2e_sequence_first", "gate_e2e", "-", 240), ("gate.refusal", "refusal", "-", 60)]
              + [(f"parity.{c}", "parity_run", c, 60) for c in F8] + [("gate.parity", "parity", "-", 30)]
              + [(f"call.off.{c}", "call", f"off:{c}", WATCHDOG) for c in CELLS] + [("gate.offcheck", "offcheck", "-", 30)]
              + [(f"call.on.{c}", "call", f"on:{c}", WATCHDOG) for c in CELLS] + [("freeze", "freeze", "-", 10)]
              + [(f"eval.{a}.{c}", "eval", f"{a}:{c}", 50) for a in ("off", "on") for c in CELLS] + [("analyse", "analyse", "-", 30)])
MANIFEST_A2 = [("gate.refusal.a2", "refusal2", "-", 60)]   # Amendment 2: refusal tested without output_schema v2
MANIFEST_A3 = [("gate.offcheck.a3", "offcheck", "-", 30)]   # Amendment 3: dosage_status timing rows excluded
MANIFESTS = {"P": MANIFEST_P, "A2": MANIFEST_A2, "A3": MANIFEST_A3}


def ledger_rows():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def spent():
    return PRE_FREEZE + sum(float(r["wall_seconds"]) for r in ledger_rows())


def done():
    return {r["unit"] for r in ledger_rows()}


def exit_of(u):
    return next((r["exit"] for r in ledger_rows() if r["unit"] == u), None)


def state():
    rev = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", WT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
    return rev, len([line for line in dirty.split("\n") if line.strip()])


def write_manifest(part):
    os.makedirs(OUT, exist_ok=True)
    text = "unit\tkind\ttarget\tdeclared_timeout_seconds\n" + "".join("\t".join(str(x) for x in m) + "\n" for m in MANIFESTS[part])
    p = f"{OUT}/MANIFEST_{part}.tsv"
    if os.path.exists(p):
        assert open(p).read() == text, "a different manifest is already frozen"
    else:
        open(p, "w").write(text)
    print(text, end="")


def record(unit, rev, dirty, timeout, wall, rss, code, timed_out, note):
    p = f"{OUT}/ledger.tsv"
    if not os.path.exists(p):
        open(p, "w").write("unit\trevision\tdirty_tree\tdeclared_timeout_seconds\twall_seconds\tpeak_rss_gib\texit\ttimed_out\tnote\n")
    with open(p, "a") as f:
        f.write(f"{unit}\t{rev}\t{dirty}\t{timeout}\t{wall:.2f}\t{rss:.3f}\t{code}\t{int(bool(timed_out))}\t{note}\n")


def donor(cell):
    prep = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in open(f"{M9C}/runs/prepare/prepare.tsv")}
    return prep[f"donor_{cell.split('.d')[1]}"][0]


GATES = ["gate.unit", "gate.sequence_report", "gate.e2e_sequence_first", "gate.refusal.a2", "gate.parity"]   # Amendment 2


def precondition(kind, target):
    if kind in ("gate_unit", "gate_seq", "gate_e2e", "parity_run", "prepare") and exit_of("build") != "0":
        return "the build did not succeed"
    if kind in ("refusal", "refusal2") and exit_of("prepare") != "0":
        return "prepare did not succeed"
    if kind == "call":
        if any(exit_of(g) != "0" for g in GATES):
            return "a gate failed: stop before real caller runs"
        if target.startswith("on:") and exit_of("gate.offcheck.a3") != "0":   # Amendment 3
            return "the off arm does not reproduce M9-C: stop before on runs"
    if kind == "offcheck" and any(f"call.off.{c}" not in done() for c in CELLS):
        return "not every off call has an outcome"
    if kind == "freeze" and any(f"call.{a}.{c}" not in done() for a in ("off", "on") for c in CELLS):
        return "not every call has an outcome"
    if kind == "eval":
        arm, cell = target.split(":")
        if exit_of("freeze") != "0":
            return "predictions are not frozen"
        if not os.path.exists(f"{OUT}/runs/call.{arm}.{cell}/call.mosaic_status.tsv"):
            return "no caller status to evaluate"
    return None


def command(kind, target, d):
    P = f"{M9C}/runs/prepare"
    py = lambda *a: [sys.executable, os.path.join(HERE, "m9g.py"), *a]
    if kind == "build":
        return ["cmake", "--build", BUILD, "--target", "panvar", "genotype_probe", "junction_proposals_test", "genotype_sequence_report_test",
                "genotype_command_e2e_test", "-j", "8"], None, "build every touched target"
    if kind == "refusal2":
        return py("refusal", d, "legacy"), None, "refusal without output_schema v2"
    if kind in ("prepare", "refusal", "parity", "offcheck", "freeze", "analyse"):
        return py(kind, d), None, kind
    if kind == "gate_unit":
        return [f"{BUILD}/junction_proposals_test"], None, "JP1.1 unit gates (oracle, invariances, ties, cache)"
    if kind == "gate_seq":
        return [f"{BUILD}/genotype_sequence_report_test", d], None, "sequence report fixtures"
    if kind == "gate_e2e":
        return [f"{BUILD}/genotype_command_e2e_test", d, WT], {"PANVAR_E2E_ONLY": "sequence_first"}, "caller-level sequence_first"
    if kind == "parity_run":
        dn = donor(target)
        reads = f"{P}/reads/{dn}/cov5"
        cur = f"{P}/curation.LOO.txt" if target.startswith("LOO") else f"{P}/curation.LZO.{dn}.txt"
        return [f"{BUILD}/genotype_probe", "junction-proposals", "--out", f"{d}/jp", "--model-gfa", ACOT_GFA, "--reference", ACOT_REF,
                "--panel-curation", cur, "--run", f"{M9C}/runs/call.{target}/call", "--reads1", f"{reads}/reads_1.fq", "--reads2",
                f"{reads}/reads_2.fq", "--block", "sv:5"], None, f"JP1.1 C++ ranking {target}"
    arm, cell = target.split(":")
    dn = donor(cell)
    reads = f"{P}/reads/{dn}/cov5"
    cur = f"{P}/curation.LOO.txt" if cell.startswith("LOO") else f"{P}/curation.LZO.{dn}.txt"
    if kind == "call":
        prof = f"{P}/profile.tsv" if arm == "off" else f"{OUT}/runs/prepare/profile.on.tsv"
        return (["nice", f"{BUILD}/panvar", "genotype", "--profile", prof, "--threads", "4", "--gfa", ACOT_GFA, "--reference-path", ACOT_REF,
                 "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq", "--library", f"{reads}/library.txt", "--out-prefix", f"{d}/call",
                 "--quiet", "--panel-curation", cur, "--dosage-mode", "audit"], None, f"caller {arm} {cell} ({dn})")
    if kind == "eval":
        return ([f"{BUILD}/genotype_probe", "evaluate", "--scheme", cell.split(".")[0].lower(), "--donor", dn, "--reference", ACOT_REF,
                 "--run-prefix", f"{OUT}/runs/call.{arm}.{cell}/call", "--model-gfa", ACOT_GFA, "--truth-gfa", ACOT_GFA, "--out", f"{d}/eval",
                 "--panel-curation", cur], None, f"truth evaluation {arm} {cell}")
    raise ValueError(kind)


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in MANIFESTS.items() for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and unit not in done()
    rev, dirty = state()
    assert dirty == 0, "the tree must be committed and clean"
    why = precondition(kind, target)
    if why:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {why}")
        print(f"{unit}: NOT RUN ({why})")
        return 0
    used = spent()
    pending = sum(m[3] for m in MANIFESTS[part] if m[1] in ("freeze", "eval", "analyse") and m[0] not in done() and m[0] != unit)
    reserve = pending if kind not in ("freeze", "eval", "analyse") else (30 if kind != "analyse" else 0)
    if used + timeout + reserve > BUDGET:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {used:.1f} + {timeout} + reserve {reserve} > {BUDGET} s")
        print(f"{unit}: NOT RUN (budget)")
        return 0
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d), f"{d} exists; a run is never repeated into the same place"
    os.makedirs(d)
    cmd, env, note = command(kind, target, d)
    cmd = ["/usr/bin/time", "-l", *cmd]
    if env:
        os.environ.update(env)
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s){' ' + str(env) if env else ''} " + " ".join(cmd) + "\n")
        log.flush()
        code, timed_out = p1_run.run_timed(cmd, log, timeout)
    for k in (env or {}):
        os.environ.pop(k, None)
    wall = time.time() - t
    text = open(f"{d}/run.log").read()
    rss = next((int(line.split()[0]) for line in text.splitlines() if "maximum resident set size" in line), 0)
    record(unit, rev, dirty, timeout, wall, rss / 2**30, code, timed_out, note)
    print("\n".join(line for line in text.splitlines() if not line.startswith(" ") and not line.startswith("$"))[-2500:])
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s (spent {spent():.1f}/{BUDGET})")
    return 0


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest(sys.argv[2])
    elif arg == "next":
        todo = [m[0] for m in MANIFESTS[sys.argv[2]] if m[0] not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(arg))
