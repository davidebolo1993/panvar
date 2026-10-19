#!/usr/bin/env python3
"""M4-A runner -- distance-aware copying-prior contract: invariance gates and model geometry (DESIGN.md; frozen manifest; one
process at a time; 600 s budget counted as time spent plus the declared timeout of the next unit). No caller change, no
reads, no genotyping.

  m4a_run.py manifest A | next A | <unit> | dev <name> <note> -- <command...>
"""
import csv
import hashlib
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
import p1_run  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/m4a_distance_prior_20261019"
BUDGET = 600
C4 = "/Users/davide.bolognini/panvar/tests/results/genotype_c4_production_11fe681"
REF = "HG00096#1#haplotype1-0000024:31848049-32074186"
GR = f"{C4}/graphs"
CYP = "/Users/davide.bolognini/panvar/tests/results/cyp2d6_shortread_pilot_v1/data/graphs"
CYP_REF, CYP_CUR = "grch38#1#chr22:42031864-42245566", "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"
TS1 = f"{C4}/mosaic_inputs/three_switches_1/staged.gfa"
MCUR = f"{C4}/mosaic_inputs/mosaic_exclude.tsv"
ACOT, ACOT_REF = "/Users/davide.bolognini/panvar/tests/results/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", "GRCh38#0#chr14:73416151-73668669"
GEOM = {"c4": (f"{GR}/c4.gfa", REF, None), "c4_loo_NA19650": (f"{GR}/NA19650.loo.gfa", REF, None),
        "ts1_curated": (TS1, REF, MCUR), "ts1_uncurated": (TS1, REF, None),
        "cyp_curated": (f"{CYP}/cyp2d6.gfa", CYP_REF, CYP_CUR), "cyp_uncurated": (f"{CYP}/cyp2d6.gfa", CYP_REF, None),
        "cyp_loo_NA19650": (f"{CYP}/loo.NA19650.gfa", CYP_REF, CYP_CUR), "acot": (ACOT, ACOT_REF, None)}
PAIRS = ["ts1_curated,ts1_uncurated", "cyp_curated,cyp_uncurated", "c4,c4_loo_NA19650", "cyp_curated,cyp_loo_NA19650"]
MANIFEST_A = ([("gates.m4a", "gates", "m4a_distance_prior", 60), ("gates.m3a_regress", "gates", "m3a_prior_representation", 60)] +
              [(f"geometry.{k}", "geometry", k, 60) for k in GEOM] + [("geometry.compare", "compare", "-", 30)])
# Amendment 1: gates.m4a failed one CONTROL row (G9 fixture: the removed donor carried no core-defining allele); fixed fixture
MANIFEST_A1 = [("gates.m4a.a1", "gates", "m4a_distance_prior", 60)]
MANIFESTS = {"A": MANIFEST_A, "A1": MANIFEST_A1}
INPUTS = list(csv.DictReader(open(os.path.join(HERE, "INPUTS.tsv")), delimiter="\t"))


def ledger_rows():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def spent():
    return sum(float(r["wall_seconds"]) for r in ledger_rows())


def done():
    return {r["unit"] for r in ledger_rows()}


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
    os.makedirs(OUT, exist_ok=True)
    p = f"{OUT}/ledger.tsv"
    if not os.path.exists(p):
        open(p, "w").write("unit\trevision\tdirty_tree\tdeclared_timeout_seconds\twall_seconds\tpeak_rss_gib\texit\ttimed_out\tnote\n")
    with open(p, "a") as f:
        f.write(f"{unit}\t{rev}\t{dirty}\t{timeout}\t{wall:.2f}\t{rss:.3f}\t{code}\t{int(bool(timed_out))}\t{note}\n")


def verify():
    for r in INPUTS:
        assert hashlib.sha256(open(r["file"], "rb").read()).hexdigest() == r["sha256"], f"{r['file']} differs from the frozen input"


def command(kind, target, d):
    verify()
    if kind == "gates":
        return [f"{BUILD}/{target}"], f"invariance gates ({target})"
    if kind == "geometry":
        gfa, ref, cur = GEOM[target]
        cmd = ["nice", f"{BUILD}/genotype_probe", "prior-geometry", "--out", f"{d}/geom", "--model-gfa", gfa, "--reference", ref]
        return (cmd + (["--panel-curation", cur] if cur else [])), "model geometry (truth-, read- and candidate-free)"
    if kind == "compare":
        return [sys.executable, os.path.join(HERE, "geometry_compare.py"), f"{d}/geometry",
                *[f"{k}={OUT}/runs/geometry.{k}/geom" for k in GEOM], "--", *PAIRS], "geometry comparison"
    raise ValueError(kind)


def dev(name, note, cmd):
    """A development run: counted against the budget, recorded, never a result."""
    used = spent()
    assert used + 60 <= BUDGET, f"development run refused: {used:.1f} s + 60 s > {BUDGET} s"
    d = f"{OUT}/dev/{name}"
    assert not os.path.exists(d), f"{d} exists"
    os.makedirs(d)
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        code, timed_out = p1_run.run_timed(cmd, log, 60)
    wall = time.time() - t
    record(name, "development", "-", "-", wall, 0.0, code, timed_out, f"development run (counted): {note}")
    print(open(f"{d}/run.log").read()[-3000:])
    print(f"{name}: exit {code}, {wall:.1f} s (spent {spent():.1f}/{BUDGET})")
    return code


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in MANIFESTS.items() for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and unit not in done()
    rev, dirty = state()
    assert dirty == 0, "the tree must be committed and clean"
    used = spent()
    if used + timeout > BUDGET:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {used:.1f} s + {timeout} s > {BUDGET} s")
        print(f"{unit}: NOT RUN")
        return 3
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d), f"{d} exists; a run is never repeated into the same place"
    os.makedirs(d)
    cmd, note = command(kind, target, d)
    cmd = ["/usr/bin/time", "-l", *cmd]
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s) " + " ".join(cmd) + "\n")
        log.flush()
        code, timed_out = p1_run.run_timed(cmd, log, timeout)
    wall = time.time() - t
    text = open(f"{d}/run.log").read()
    rss = next((int(l.split()[0]) for l in text.splitlines() if "maximum resident set size" in l), 0)
    record(unit, rev, dirty, timeout, wall, rss / 2**30, code, timed_out, note)
    print("\n".join(l for l in text.splitlines() if not l.startswith(" ") and not l.startswith("$"))[-4000:])
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s (spent {spent():.1f}/{BUDGET})")
    return 0


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest(sys.argv[2])
    elif arg == "dev":
        sep = sys.argv.index("--")
        sys.exit(dev(sys.argv[2], sys.argv[3], sys.argv[sep + 1:]))
    elif arg == "next":
        todo = [m[0] for m in MANIFESTS[sys.argv[2]] if m[0] not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(arg))
