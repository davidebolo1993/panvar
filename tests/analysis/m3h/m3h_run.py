#!/usr/bin/env python3
"""M3-H runner -- bounded evidence decomposition of the ctx:0 matched pair (DESIGN.md; frozen manifest; one process at a
time; 4 threads; 300 s budget counted as time spent plus the declared timeout of the next unit; existing resource settings).
No model or caller change; nothing here feeds a reported prediction.

  m3h_run.py manifest H | next H | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m3h_evidence_decomposition_20261019"
BUDGET = 300
C4 = "/Users/davide.bolognini/panvar/tests/results/genotype_c4_production_11fe681"
REF = "HG00096#1#haplotype1-0000024:31848049-32074186"
PLAN = f"{C4}/mosaic_inputs/three_switches_1"
GFA, CUR = f"{PLAN}/staged.gfa", f"{C4}/mosaic_inputs/mosaic_exclude.tsv"
M3F = "/Users/davide.bolognini/panvar/tests/results/m3f_validation_20261019/runs"
RUN = {"2x": f"{M3F}/run.ts1_2.off/call", "5x": f"{M3F}/run.ts1_5.off/call"}   # saved M3-F off arms: settings + model fingerprint
READS = {"2x": f"{PLAN}/cov2", "5x": f"{PLAN}/cov5"}
MODEL = ["--model-gfa", GFA, "--reference", REF, "--panel-curation", CUR]
MANIFEST_H = [("score.2x", "score", "2x", 60), ("score.5x", "score", "5x", 90), ("spell", "spell", "-", 40), ("panel", "panel", "-", 40),
              ("prior.attrib", "attrib", "-", 20), ("prior.pool", "pool", "-", 40), ("reads.detail", "reads", "-", 60)]
MANIFESTS = {"H": MANIFEST_H}
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
    R = f"{OUT}/runs"
    if kind == "score":
        rd = READS[target]
        return ["nice", f"{BUILD}/genotype_probe", "score-set", "--out", f"{d}/score", *MODEL, "--donor", "MOSAIC", "--run", RUN[target],
                "--reads1", f"{rd}/reads_1.fq", "--reads2", f"{rd}/reads_2.fq", "--library", f"{rd}/library.txt", "--threads", "4",
                "--candidates", os.path.join(HERE, f"pair_{target}.tsv"), "--decompose", "R_t0,R_c0", "--decompose-block", "ctx:0",
                "--contrast", "R_t0,R_c0"], f"matched pair scored against the {target} reads (its own frozen library), per-pair decomposition"
    if kind == "spell":
        return ["nice", f"{BUILD}/genotype_probe", "spell", "--out", f"{d}/spell", *MODEL, "--block", "ctx:0", "--candidates",
                os.path.join(HERE, "spell.tsv"), "--truth-gfa", GFA, "--donor", "MOSAIC"], "the pair's bytes + truth bytes (identity only)"
    if kind == "panel":
        return ["nice", f"{BUILD}/genotype_probe", "panel-alleles", "--gfa", GFA, "--reference", REF, "--out", f"{d}/panel.tsv"], \
            "panel allele matrix (MOSAIC rows dropped downstream = the curation)"
    if kind == "attrib":
        return [sys.executable, os.path.join(HERE, "prior_attrib.py"), f"{R}/panel/panel.tsv", f"{R}/score.2x/score.scores_truth_free.tsv",
                f"{d}/prior"], "sequence-level prior attribution (gated)"
    if kind == "pool":
        return ["nice", f"{BUILD}/genotype_probe", "prior-pool", "--out", f"{d}/pool", "--cell", "ts1_2", "--run", RUN["2x"], *MODEL,
                "--counterfactual", os.path.join(HERE, "counterfactual.tsv")], "M3-A prior-pool: saved-prior gate + sequence-level counterfactual priors"
    if kind == "reads":
        return [sys.executable, os.path.join(HERE, "read_detail.py"), f"{R}/score.2x", f"{R}/score.5x", READS["2x"], READS["5x"],
                f"{R}/spell/spell", f"{d}/reads"], "read-level decomposition"
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
    elif arg == "next":
        todo = [m[0] for m in MANIFESTS[sys.argv[2]] if m[0] not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    else:
        sys.exit(run_unit(arg))
