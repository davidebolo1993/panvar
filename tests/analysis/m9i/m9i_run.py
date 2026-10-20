#!/usr/bin/env python3
"""M9-I runner -- representation-aware compound targets, read-only on saved outputs (DESIGN.md; frozen manifest; one process
at a time; 600 s counted as the disclosed pre-freeze inspection plus time spent plus the declared timeout of the next unit plus
the analysis).

  m9i_run.py manifest P | next P | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m9i_compound_targets_20261020"
M9C = "/Users/davide.bolognini/panvar/tests/results/m9c_acot_reliability_20261020"
R = "/Users/davide.bolognini/panvar/tests/results"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
BUDGET = 600
PRE_FREEZE = 3.0   # disclosed: read-only table inspections while designing (rounded up)
R0 = "/Users/davide.bolognini/panvar/tests/results"
CP, HP = f"{R0}/m9c_acot_reliability_20261020/runs/prepare", f"{R0}/m9h_jp11_validation_20261020/runs/prepare"
CYPG = f"{R0}/cyp2d6_shortread_pilot_v1/data/graphs"
CYPREF = "grch38#1#chr22:42031864-42245566"
CYPCUR = "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"
TRUTH = ([(f"truth.ACOT.LOO.{dn}", ACOT_GFA, ACOT_REF, f"{CP}/curation.LOO.txt", ACOT_GFA, dn, "sv:5,ctx:6")
          for dn in ("NA19909", "NA18747", "HG02178", "HG00639", "NA21106", "HG03831", "HG02723", "HG02074", "NA18971", "HG02083")]
         + [("truth.ACOT.CLZO.d1.NA19909", ACOT_GFA, ACOT_REF, f"{CP}/curation.LZO.NA19909.txt", ACOT_GFA, "NA19909", "sv:5,ctx:6"),
            ("truth.ACOT.CLZO.d2.NA18747", ACOT_GFA, ACOT_REF, f"{CP}/curation.LZO.NA18747.txt", ACOT_GFA, "NA18747", "sv:5,ctx:6"),
            ("truth.ACOT.HZ1.HG02723", ACOT_GFA, ACOT_REF, f"{HP}/curation.LZO.HG02723.txt", ACOT_GFA, "HG02723", "sv:5,ctx:6"),
            ("truth.ACOT.HZ2.HG02074", ACOT_GFA, ACOT_REF, f"{HP}/curation.LZO.HG02074.txt", ACOT_GFA, "HG02074", "sv:5,ctx:6"),
            ("truth.CYP.original.HG02769", f"{CYPG}/cyp2d6.gfa", CYPREF, CYPCUR, f"{CYPG}/cyp2d6.gfa", "HG02769", "sv:3,ctx:4,sv:5"),
            ("truth.CYP.loo.HG02953", f"{CYPG}/loo.HG02953.gfa", CYPREF, CYPCUR, f"{CYPG}/cyp2d6.gfa", "HG02953", "sv:3,ctx:4,sv:5"),
            ("truth.CYP.loo.NA19650", f"{CYPG}/loo.NA19650.gfa", CYPREF, CYPCUR, f"{CYPG}/cyp2d6.gfa", "NA19650", "sv:3,ctx:4,sv:5")])
MANIFEST_P = ([("gates", "gates", "-", 10), ("inventory", "inventory", "-", 60)] + [(t[0], "truth", t[0], 30) for t in TRUTH]
              + [("analyse", "analyse", "-", 120)])
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


def precondition(kind, target):
    if kind in ("truth", "inventory", "analyse") and exit_of("gates") != "0":
        return "the gates failed"
    if kind == "analyse" and any(t[0] not in done() for t in TRUTH):
        return "not every truth unit has an outcome"
    return None


def command(kind, target, d):
    py = lambda *a: [sys.executable, os.path.join(HERE, "m9i.py"), *a]
    if kind in ("gates", "inventory", "analyse"):
        return py(kind, d), None, kind
    if kind == "truth":
        name, gfa, ref, cur, truth, dn, blocks = next(t for t in TRUTH if t[0] == target)
        return [f"{BUILD}/genotype_probe", "block-sequences", "--out", f"{d}/bs", "--reference", ref, "--model-gfa", gfa, "--panel-curation", cur,
                "--blocks", blocks, "--truth-gfa", truth, "--donor", dn], None, f"truth bytes {target}"
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
    reserve = 120 if kind != "analyse" else 0   # the analysis
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
