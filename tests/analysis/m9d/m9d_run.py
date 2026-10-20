#!/usr/bin/env python3
"""M9-D runner -- isolated attribution of the M9-C ACOT high-confidence errors (DESIGN.md; frozen manifest; one process at a
time; 900 s counted as time spent plus the declared timeout of the next unit; refused units are recorded NOT RUN).

  m9d_run.py manifest P | next P | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m9d_attribution_20261020"
M9C = "/Users/davide.bolognini/panvar/tests/results/m9c_acot_reliability_20261020"
R = "/Users/davide.bolognini/panvar/tests/results"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
BUDGET = 900
CASES = {"E1": ("LOO.d1", "sv:5"), "E2": ("LOO.d2", "sv:5"), "E3": ("LOO.d3", "sv:5"), "E4": ("LOO.d3", "sv:9"), "E5": ("LOO.d4", "sv:5"),
         "E6": ("LOO.d6", "sv:5"), "C1": ("LZO.d1", "sv:5"), "C2": ("LZO.d2", "sv:5")}
MANIFEST_P = ([("build", "build", "-", 120), ("models", "models", "-", 120), ("construct", "construct", "-", 30)]
              + [(f"score.{c}", "score", c, 90) for c in CASES] + [("degeneracy", "degeneracy", "-", 120), ("analyse", "analyse", "-", 30)])
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
    prep = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in open(f"{M9C}/runs/prepare/prepare.tsv")}
    return prep[f"donor_{cell.split('.d')[1]}"][0]


def precondition(kind, target):
    if kind == "models" and exit_of("build") != "0":
        return "the dump tool did not build"
    if kind == "degeneracy" and exit_of("models") != "0":
        return "no model dumps"
    if kind == "construct" and exit_of("models") != "0":
        return "model dumps or fingerprints failed"
    if kind == "score" and exit_of("construct") != "0":
        return "construction or spelling gates failed"
    if kind == "analyse" and not any(exit_of(f"score.{c}") == "0" for c in CASES):
        return "no score-set completed"
    return None


def command(kind, target, d):
    P = f"{M9C}/runs/prepare"
    py = lambda mode: ([sys.executable, os.path.join(HERE, "m9d.py"), mode, d], None, mode)
    if kind == "build":
        return ["cmake", "--build", BUILD, "--target", "cyp2d6_block_alleles", "-j", "4"], None, "relink the model dump tool"
    if kind in ("models", "construct", "degeneracy", "analyse"):
        return py(kind)
    if kind == "score":
        cell, block = CASES[target]
        dn = donor(cell)
        reads = f"{P}/reads/{dn}/cov5"
        cur = f"{P}/curation.LOO.txt" if cell.startswith("LOO") else f"{P}/curation.LZO.{dn}.txt"
        con = f"{OUT}/runs/construct"
        dec = open(f"{con}/{target}.decompose").read().strip()
        return ([f"{BUILD}/genotype_probe", "score-set", "--out", f"{d}/score", "--model-gfa", ACOT_GFA, "--reference", ACOT_REF, "--donor", dn,
                 "--run", f"{M9C}/runs/call.{cell}/call", "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq", "--library", f"{reads}/library.txt",
                 "--threads", "4", "--switch-probability", "0.01", "--candidates", f"{con}/{target}.candidates.tsv", "--panel-curation", cur]
                + (["--decompose", dec, "--decompose-block", block] if dec else []), None, f"score-set {target} ({cell} {block})")
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
    reserve = 30 if kind != "analyse" else 0   # the analysis
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
