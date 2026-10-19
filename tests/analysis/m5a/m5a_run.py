#!/usr/bin/env python3
"""M5-A runner -- bounded pooled-winner refinement pilot (diagnostic; DESIGN.md; frozen manifests; one process at a time;
4 threads; 600 s counted as time spent plus the declared timeout of the next unit). No caller change: the legs are fresh
diagnostic replays through genotype_probe repair-search (the ordinary single-block engine) from saved M4-B outputs.

  m5a_run.py manifest X|P | next X|P | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m5a_pooled_winner_20261019"
BUDGET = 600
M4B = "/Users/davide.bolognini/panvar/tests/results/m4b_distance_prior_20261019/runs"
C4 = "/Users/davide.bolognini/panvar/tests/results/genotype_c4_production_11fe681"
REF_C4 = "HG00096#1#haplotype1-0000024:31848049-32074186"
LAMBDA_C4 = "1872698.6117028969"   # M4-B LAMBDA.tsv, C4 cells (cross-score only)
CELL = {
    "c4_lzo_hg00171_5": dict(gfa=f"{C4}/graphs/c4.gfa", ref=REF_C4, curation=None, reads=f"{C4}/reads/HG00171/cov5", donor="HG00171",
                             footprint=4294967296),
    "cyp_lzo_hg02769_5": dict(gfa="/Users/davide.bolognini/panvar/tests/results/cyp2d6_shortread_pilot_v1/data/graphs/cyp2d6.gfa",
                              ref="grch38#1#chr22:42031864-42245566", curation="/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv",
                              reads="/Users/davide.bolognini/panvar/tests/results/cyp2d6_shortread_pilot_v1/data/reads/HG02769/cov5", donor="HG02769",
                              footprint=6442450944),
    "c4_three_switches_1_2": dict(gfa=f"{C4}/mosaic_inputs/three_switches_1/staged.gfa", ref=REF_C4, curation=f"{C4}/mosaic_inputs/mosaic_exclude.tsv",
                                  reads=f"{C4}/mosaic_inputs/three_switches_1/cov2", donor="MOSAIC", footprint=4294967296),
}
ORDER = ["c4_lzo_hg00171_5", "cyp_lzo_hg02769_5", "c4_three_switches_1_2"]   # controls first
ROOTS = list(csv.DictReader(open(os.path.join(HERE, "ROOTS.tsv")), delimiter="\t"))
LEG_TIMEOUT = {"c4_lzo_hg00171_5": 120, "cyp_lzo_hg02769_5": 120, "c4_three_switches_1_2": 90}
MANIFEST_X = [("xscore.const", "xscore", "const", 60), ("xscore.dist", "xscore", "dist", 60), ("xscore.table", "xtable", "-", 30)]
MANIFEST_P = ([(f"leg.{r['cell']}.{r['root']}", "leg", f"{r['cell']}:{r['root']}", LEG_TIMEOUT[r["cell"]]) for c in ORDER for r in ROOTS if r["cell"] == c] +
              [("combine", "combine", "-", 30)] + [(f"eval.{c}", "eval", c, 90) for c in ORDER] + [("summary", "summary", "-", 30)])
MANIFESTS = {"X": MANIFEST_X, "P": MANIFEST_P}
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


def kv(p):
    return {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l}


def model_args(c):
    return ["--model-gfa", c["gfa"], "--reference", c["ref"], *(["--panel-curation", c["curation"]] if c["curation"] else [])]


def read_args(c):
    return ["--reads1", f"{c['reads']}/reads_1.fq", "--reads2", f"{c['reads']}/reads_2.fq", "--library", f"{c['reads']}/library.txt", "--threads", "4"]


def command(kind, target, d, timeout):
    verify()
    if kind == "xscore":
        c = CELL["c4_three_switches_1_2"]
        run = f"{M4B}/run.c4_three_switches_1_2.{target}/call"
        extra = [] if target == "const" else ["--copying-transition", "distance_ref_core_v1_experimental", "--copying-length-bp", LAMBDA_C4]
        return ["nice", f"{BUILD}/genotype_probe", "score-set", "--out", f"{d}/xs", *model_args(c), "--donor", "NA", "--run", run, *read_args(c),
                "--candidates", os.path.join(HERE, f"xscore_{target}.tsv"), *extra], f"cross-score under the {target} prior (truth-free)"
    if kind == "xtable":
        return [sys.executable, os.path.join(HERE, "xscore_table.py"), f"{OUT}/runs/xscore.const/xs.scores_truth_free.tsv",
                f"{OUT}/runs/xscore.dist/xs.scores_truth_free.tsv", f"{d}/xscore"], "cross-score table"
    if kind == "leg":
        cell, root = target.split(":")
        c = CELL[cell]
        r = next(x for x in ROOTS if x["cell"] == cell and x["root"] == root)
        st = kv(f"{M4B}/run.{cell}.const/call.mosaic_status.tsv")
        return ["nice", f"{BUILD}/genotype_probe", "repair-search", "--out", f"{d}/leg", *model_args(c), "--donor", "NA",
                "--run", f"{M4B}/run.{cell}.const/call", *read_args(c), "--local-pairs", st["repair_local_pairs"], "--rounds", "32",
                "--context", "1", "--delta", "1", "--origin-cache", "1", "--origin-cache-policy", "planned", "--origin-cache-bytes", st["repair_origin_cache_bytes"],
                "--max-neighbors", st["repair_max_neighbours"], "--max-physical-origins", st["repair_max_physical_origins"],
                "--max-state-classes", st["repair_max_state_classes"], "--max-piece-hits", st["repair_max_piece_hits"],
                "--max-latent-pairs", st["repair_max_latent_pairs"], "--switch-probability", st["switch_probability"],
                "--start", f"{r['hap1']}|{r['hap2']}", "--start-score", r["score"], "--time-limit-seconds", str(timeout - 10),
                "--footprint-limit-bytes", str(c["footprint"])], f"one ordinary continuation leg from root {root} (constant mode)"
    if kind == "combine":
        return [sys.executable, os.path.join(HERE, "combine.py"), f"{OUT}/runs", os.path.join(HERE, "ROOTS.tsv"), f"{d}/m5a"], \
            "truth-free combination + prediction freeze"
    if kind == "eval":
        assert os.path.exists(f"{OUT}/runs/combine/m5a.selection.sha256"), "evaluation needs the frozen diagnostic predictions"
        c = CELL[target]
        cands = f"{OUT}/runs/combine/m5a.candidates.{target}.tsv"
        return ["nice", f"{BUILD}/genotype_probe", "score-set", "--out", f"{d}/ev", *model_args(c), "--donor", c["donor"],
                "--run", f"{M4B}/run.{target}.const/call", *read_args(c), "--candidates", cands, "--truth-gfa", c["gfa"]], \
            "truth evaluation of the original roots and the leg terminals (after the freeze)"
    if kind == "summary":
        return [sys.executable, os.path.join(HERE, "summary.py"), f"{OUT}/runs", f"{OUT}/ledger.tsv", f"{d}/m5a", *ORDER], "summary"
    raise ValueError(kind)


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
    cmd, note = command(kind, target, d, timeout)
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
    print("\n".join(l for l in text.splitlines() if not l.startswith(" ") and not l.startswith("$"))[-3000:])
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
