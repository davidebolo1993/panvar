#!/usr/bin/env python3
"""M9-H runner -- bounded paired validation of JP1.1 on four new ACOT families (DESIGN.md; frozen manifest; one process at a
time; 3,600 s counted as time spent plus the declared timeout of the next unit plus, for preparation/gates/calls, the declared
timeouts of every pending freeze/eval/truth/intervals/analyse unit; refused units are recorded NOT RUN, never replaced).

  m9h_run.py manifest P | next P | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m9h_jp11_validation_20261020"
M9C = "/Users/davide.bolognini/panvar/tests/results/m9c_acot_reliability_20261020"
R = "/Users/davide.bolognini/panvar/tests/results"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
BUDGET = 3600
PRE_FREEZE = 2.0   # disclosed: read-only inspections while designing (rounded up)
CELLS = ["L1", "L2", "L3", "L4", "Z1", "Z2"]
WATCHDOG = 930   # caller guard 900 s + 30 s
ORDER = [(c, arms) for c, arms in zip(CELLS, [("off", "on"), ("on", "off")] * 3)]   # deterministic alternation
MANIFEST_P = ([("correction", "correction", "-", 10), ("prepare", "prepare", "-", 120), ("models", "models", "-", 60), ("gate.unit", "gate_unit", "-", 30)]
              + [(f"call.{a}.{c}", "call", f"{a}:{c}", WATCHDOG) for c, arms in ORDER for a in arms] + [("freeze", "freeze", "-", 10)]
              + [(f"eval.{a}.{c}", "eval", f"{a}:{c}", 50) for c in CELLS for a in ("off", "on")]
              + [(f"truth.{c}", "truth", c, 30) for c in CELLS] + [("intervals", "intervals", "-", 120), ("analyse", "analyse", "-", 30)])
MANIFESTS = {"P": MANIFEST_P}


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


def cell_info(cell):
    p = {}
    for line in open(f"{OUT}/runs/prepare/prepare.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0].startswith("family_"):
            p[int(f[0].split("_")[1])] = f[2].split()[1]
    dn = p[int(cell[1])]
    P = f"{OUT}/runs/prepare"
    cur = f"{P}/curation.LOO.txt" if cell[0] == "L" else f"{P}/curation.LZO.{dn}.txt"
    return dn, cur, ("loo" if cell[0] == "L" else "lzo")


def precondition(kind, target):
    if kind in ("models", "gate_unit") and exit_of("prepare") != "0":
        return "prepare did not succeed"
    if kind == "call" and any(exit_of(u) != "0" for u in ("prepare", "models", "gate.unit")):
        return "preparation or a gate failed: stop before caller runs"
    if kind == "freeze" and any(f"call.{a}.{c}" not in done() for c in CELLS for a in ("off", "on")):
        return "not every call has an outcome"
    if kind == "eval":
        arm, cell = target.split(":")
        if exit_of("freeze") != "0":
            return "predictions are not frozen"
        if not os.path.exists(f"{OUT}/runs/call.{arm}.{cell}/call.mosaic_status.tsv"):
            return "no caller status to evaluate"
    if kind in ("truth", "intervals", "analyse") and exit_of("freeze") != "0":
        return "predictions are not frozen"
    return None


def command(kind, target, d):
    P = f"{OUT}/runs/prepare"
    py = lambda *a: [sys.executable, os.path.join(HERE, "m9h.py"), *a]
    if kind == "prepare":
        return [sys.executable, os.path.join(HERE, "prepare.py"), d], None, "selection, panels, reads, lineage, profiles"
    if kind in ("correction", "models", "freeze", "intervals", "analyse"):
        return py(kind, d), None, kind
    if kind == "gate_unit":
        return [f"{BUILD}/junction_proposals_test"], None, "JP1.1 unit gates (unchanged binaries since M9-G)"
    if kind == "truth":
        dn, cur, _ = cell_info(target)
        return [f"{BUILD}/genotype_probe", "block-sequences", "--out", f"{d}/bs", "--reference", ACOT_REF, "--model-gfa", ACOT_GFA,
                "--panel-curation", cur, "--blocks", "sv:5,ctx:6", "--truth-gfa", ACOT_GFA, "--donor", dn], None, f"truth bytes {target}"
    arm, cell = target.split(":")
    dn, cur, scheme = cell_info(cell)
    reads = f"{P}/reads/{dn}/cov5"
    if kind == "call":
        return (["nice", f"{BUILD}/panvar", "genotype", "--profile", f"{P}/profile.{arm}.tsv", "--threads", "4", "--gfa", ACOT_GFA, "--reference-path", ACOT_REF,
                 "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq", "--library", f"{reads}/library.txt", "--out-prefix", f"{d}/call",
                 "--quiet", "--panel-curation", cur, "--dosage-mode", "audit"], None, f"caller {arm} {cell} ({dn})")
    if kind == "eval":
        return ([f"{BUILD}/genotype_probe", "evaluate", "--scheme", scheme, "--donor", dn, "--reference", ACOT_REF,
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
    pending = sum(m[3] for m in MANIFESTS[part] if m[1] in ("freeze", "eval", "truth", "intervals", "analyse") and m[0] not in done() and m[0] != unit)
    reserve = pending if kind not in ("freeze", "eval", "truth", "intervals", "analyse") else (30 if kind != "analyse" else 0)
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
