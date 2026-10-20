#!/usr/bin/env python3
"""M9-E runner -- bounded reachability audit of the three missed improving neighbours (DESIGN.md; frozen manifest; one process at
a time; 300 s counted as the disclosed pre-freeze execution plus time spent plus the declared timeout of the next unit).

  m9e_run.py manifest P | next P | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m9e_reachability_20261020"
M9C = "/Users/davide.bolognini/panvar/tests/results/m9c_acot_reliability_20261020"
R = "/Users/davide.bolognini/panvar/tests/results"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
BUDGET = 300
PRE_FREEZE = 2.0   # disclosed: small read-only table inspections while designing (rounded up)
CELLS = ["LOO.d1", "LOO.d3", "LOO.d6"]
RANK_BLOCKS = "sv:5,ctx:2,ctx:4,ctx:6,ctx:8,ctx:12,ctx:16,ctx:18,sv:13,sv:17"   # sv:5 plus every block of a recorded trunk move
MANIFEST_P = [("fixture", "fixture", "-", 10)] + [(f"ranks.{c}", "ranks", c, 60) for c in CELLS] + [("analyse", "analyse", "-", 30)]
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


def donor(cell):
    prep = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in open(f"{M9C}/runs/prepare/prepare.tsv")}
    return prep[f"donor_{cell.split('.d')[1]}"][0]


def precondition(kind, target):
    if kind == "analyse" and (exit_of("fixture") != "0" or any(exit_of(f"ranks.{c}") != "0" for c in CELLS)):
        return "the fixture or a rank reproduction did not succeed"
    return None


def command(kind, target, d):
    P = f"{M9C}/runs/prepare"
    if kind in ("fixture", "analyse"):
        return [sys.executable, os.path.join(HERE, "m9e.py"), kind, d], None, kind
    if kind == "ranks":
        dn = donor(target)
        reads = f"{P}/reads/{dn}/cov5"
        lg = [l.rstrip("\n").split("\t") for l in open(f"{M9C}/runs/call.{target}/call.locus_genotype.tsv") if not l.startswith("#")]
        row = dict(zip(lg[0], lg[1]))
        open(f"{d}/winner.candidates.tsv", "w").write(f"winner\t{row['hap1_block_path']}\t{row['hap2_block_path']}\t{row['total_score']}\n")
        return ([f"{BUILD}/genotype_probe", "score-set", "--out", f"{d}/score", "--model-gfa", ACOT_GFA, "--reference", ACOT_REF, "--donor", dn,
                 "--run", f"{M9C}/runs/call.{target}/call", "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq", "--library", f"{reads}/library.txt",
                 "--threads", "4", "--switch-probability", "0.01", "--candidates", f"{d}/winner.candidates.tsv", "--panel-curation", f"{P}/curation.LOO.txt",
                 "--rank-blocks", RANK_BLOCKS, "--emissions-only", "1"], None, f"local-pair ranks {target} ({dn})")
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
