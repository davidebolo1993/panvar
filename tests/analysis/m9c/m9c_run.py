#!/usr/bin/env python3
"""M9-C runner -- bounded ACOT reliability pilot under sequence_first_v2 (DESIGN.md; frozen manifest; one process at a time;
2,400 s counted as time spent plus the declared timeout of the next unit plus the declared timeouts of every pending
non-call unit; refused units are recorded NOT RUN and never substituted).

  m9c_run.py manifest P | next P | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m9c_acot_reliability_20261020"
R = "/Users/davide.bolognini/panvar/tests/results"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
BUDGET = 2400
WATCHDOG = 930   # caller guard 900 s (profile time_limit_seconds) + 30 s for shutdown and output finalisation
CELLS = [f"LOO.d{i}" for i in range(1, 7)] + ["LZO.d1", "LZO.d2"]
MANIFEST_P = ([("build", "build", "-", 400), ("prepare", "prepare", "-", 120), ("gate.sequence_report", "gate_unit", "-", 60),
               ("gate.e2e_sequence_first", "gate_e2e", "-", 240)]
              + [(f"call.{c}", "call", c, WATCHDOG) for c in CELLS] + [("freeze", "freeze", "-", 10)]
              + [(f"eval.{c}", "eval", c, 60) for c in CELLS] + [("analyse", "analyse", "-", 60)])
MANIFESTS = {"P": MANIFEST_P}


def ledger_rows():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def spent():
    return sum(float(r["wall_seconds"]) for r in ledger_rows())


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
    prep = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in open(f"{OUT}/runs/prepare/prepare.tsv")}
    return prep[f"donor_{cell.split('.d')[1]}"][0]


def precondition(kind, target):
    """None if the unit may run; otherwise the reason it is refused (recorded NOT RUN)."""
    if kind in ("gate_unit", "gate_e2e") and exit_of("build") != "0":
        return "the test build did not succeed"
    if kind == "call":
        if exit_of("prepare") != "0":
            return "prepare did not succeed"
        if exit_of("gate.sequence_report") != "0" or exit_of("gate.e2e_sequence_first") != "0":
            return "a contract fixture gate failed: stop before real-data calls"
    if kind == "freeze" and any(f"call.{c}" not in done() for c in CELLS):
        return "not every call has a recorded outcome"
    if kind == "eval":
        if exit_of("freeze") != "0":
            return "predictions are not frozen"
        if not os.path.exists(f"{OUT}/runs/call.{target}/call.mosaic_status.tsv"):
            return "no caller status to evaluate (resource failure or refused call)"
    return None


def command(kind, target, d):
    P = f"{OUT}/runs/prepare"
    if kind == "build":
        return ["cmake", "--build", BUILD, "--target", "genotype_sequence_report_test", "genotype_command_e2e_test", "-j", "4"], None, "relink the focused test binaries"
    if kind == "prepare":
        return [sys.executable, os.path.join(HERE, "prepare.py"), d], None, "selection, panels, reads, profile, hashes"
    if kind == "gate_unit":
        return [f"{BUILD}/genotype_sequence_report_test", d], None, "sequence report fixtures (coupled label ambiguity, bounds, panel membership/FASTA)"
    if kind == "gate_e2e":
        return [f"{BUILD}/genotype_command_e2e_test", d, WT], {"PANVAR_E2E_ONLY": "sequence_first"}, "caller-level sequence_first (v1/v2, FASTA policy)"
    if kind == "call":
        dn = donor(target)
        reads = f"{P}/reads/{dn}/cov5"
        cur = f"{P}/curation.LOO.txt" if target.startswith("LOO") else f"{P}/curation.LZO.{dn}.txt"
        return (["nice", f"{BUILD}/panvar", "genotype", "--profile", f"{P}/profile.tsv", "--threads", "4", "--gfa", ACOT_GFA, "--reference-path", ACOT_REF,
                 "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq", "--library", f"{reads}/library.txt", "--out-prefix", f"{d}/call",
                 "--quiet", "--panel-curation", cur, "--dosage-mode", "audit"], None, f"caller {target} ({dn})")
    if kind == "eval":
        dn = donor(target)
        cur = f"{P}/curation.LOO.txt" if target.startswith("LOO") else f"{P}/curation.LZO.{dn}.txt"
        return ([f"{BUILD}/genotype_probe", "evaluate", "--scheme", target.split(".")[0].lower(), "--donor", dn, "--reference", ACOT_REF,
                 "--run-prefix", f"{OUT}/runs/call.{target}/call", "--model-gfa", ACOT_GFA, "--truth-gfa", ACOT_GFA, "--out", f"{d}/eval",
                 "--panel-curation", cur], None, f"truth evaluation {target} ({dn})")
    if kind in ("freeze", "analyse"):
        return [sys.executable, os.path.join(HERE, "m9c.py"), kind, d], None, kind
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
    pending = sum(m[3] for m in MANIFESTS[part] if m[1] not in ("call", "build", "prepare", "gate_unit", "gate_e2e") and m[0] not in done() and m[0] != unit)
    reserve = pending if kind in ("build", "prepare", "gate_unit", "gate_e2e", "call") else 0
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
