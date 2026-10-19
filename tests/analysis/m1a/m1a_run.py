#!/usr/bin/env python3
"""M1-A runner (frozen manifest; one process at a time; <= 4 threads; <= 300 s per unit; 1200 s budget counted as time
spent in this batch -- development runs included -- plus the declared timeout of the next unit; 6 GiB envelope).
Isolated evaluator units only -- no production caller run, no truth in any caller:
  * score.c4_mosaic_2x, score.cyp_loo_sv9: genotype_probe score-set (the unchanged reference scorer, built from the
    baseline revision) on declared allele-path diplotypes (tests/analysis/m1a/candidates/*.tsv). Every saved candidate's
    expected total must reproduce or the tool refuses before reporting. Truth is not loaded (no --truth-gfa).
  * check.acot_truth: whole-sequence comparison of the ACOT simulation truth with the saved S2-L calls.

  m1a_run.py manifest B | next B | <unit> | dev <name> <note> -- <command...>
"""
import csv
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
import p1_run  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/m1a_attribution_20261019"
BUDGET = 1200
C4 = "/Users/davide.bolognini/panvar/tests/results/genotype_c4_production_11fe681"
CYP = "/Users/davide.bolognini/panvar/tests/results/cyp2d6_shortread_pilot_v1"
SCORE = {
    "c4_mosaic_2x": dict(run=f"{C4}/mosaic/one_homologue_once_1/call.2", gfa=f"{C4}/mosaic_inputs/one_homologue_once_1/staged.gfa",
                         curation=f"{C4}/mosaic_inputs/mosaic_exclude.tsv", reference="HG00096#1#haplotype1-0000024:31848049-32074186",
                         reads=f"{C4}/mosaic_inputs/one_homologue_once_1/cov2", donor="MOSAIC"),
    "cyp_loo_sv9": dict(run=f"{CYP}/data/runs_a5/loo.NA19650.5", gfa=f"{CYP}/data/graphs/loo.NA19650.gfa",
                        curation="/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv",
                        reference="grch38#1#chr22:42031864-42245566", reads=f"{CYP}/data/reads/NA19650/cov5", donor="NA19650"),
}
MANIFEST_B = [("score.c4_mosaic_2x", "score", "c4_mosaic_2x", 300),
              ("score.cyp_loo_sv9", "score", "cyp_loo_sv9", 300),
              ("check.acot_truth", "acot", "-", 30)]
MANIFESTS = {"B": MANIFEST_B}


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
    rows = MANIFESTS[part]
    assert rows, f"manifest {part} is empty"
    text = "unit\tkind\ttarget\tdeclared_timeout_seconds\n" + "".join("\t".join(str(x) for x in m) + "\n" for m in rows)
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


def timed(cmd, d, timeout):
    cmd = ["/usr/bin/time", "-l", *cmd]
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s) " + " ".join(cmd) + "\n")
        log.flush()
        code, timed_out = p1_run.run_timed(cmd, log, timeout)
    wall = time.time() - t
    text = open(f"{d}/run.log").read()
    rss = 0
    for line in text.splitlines():
        if "maximum resident set size" in line:
            rss = int(line.split()[0])
    return code, timed_out, wall, rss, text


def tsv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path, errors="replace"):
            k, _, v = line.rstrip("\n").partition("\t")
            out[k] = v
    return out


def command(unit, kind, target, timeout, d):
    if kind == "acot":
        return [sys.executable, os.path.join(HERE, "acot_truth_check.py"), f"{d}/acot_truth.tsv"], "-"
    c = SCORE[target]
    return ["nice", f"{BUILD}/genotype_probe", "score-set", "--out", f"{d}/score", "--reference", c["reference"], "--donor", c["donor"],
            "--run", c["run"], "--model-gfa", c["gfa"], "--panel-curation", c["curation"],
            "--candidates", os.path.join(HERE, "candidates", f"{target}.tsv"), "--library", f"{c['reads']}/library.txt",
            "--reads1", f"{c['reads']}/reads_1.fq", "--reads2", f"{c['reads']}/reads_2.fq", "--threads", "4"], "isolated evaluator"


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in MANIFESTS.items() for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and timeout <= 300 and unit not in done()
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
    cmd, note = command(unit, kind, target, timeout, d)
    if cmd is None:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, note)
        print(f"{unit}: {note}")
        return 3
    code, timed_out, wall, rss, text = timed(cmd, d, timeout)
    if kind == "ctest":
        m = re.search(r"Test +#\d+: \S+ \.+\s*\**\s*(\w[\w ]*?)\s+([\d.]+) sec", text)
        note = f"ctest: {m.group(1)} {m.group(2)} s" if m else "ctest: no result line"
    record(unit, rev, dirty, timeout, wall, rss / 2**30, code, timed_out, note)
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s, {note} (spent {spent():.1f}/{BUDGET})")
    return 0


def dev(name, note, cmd):
    """A development run: counted against the budget, recorded, never a result."""
    used = spent()
    assert used + 300 <= BUDGET, f"development run refused: {used:.1f} s + 300 s > {BUDGET} s"
    d = f"{OUT}/dev/{name}"
    assert not os.path.exists(d), f"{d} exists"
    os.makedirs(d)
    code, timed_out, wall, rss, text = timed(cmd, d, 300)
    record(name, "development", "-", "-", wall, rss / 2**30, code, timed_out, f"development run (counted): {note}")
    print(text[-4000:])
    print(f"{name}: exit {code}, {wall:.1f} s (spent {spent():.1f}/{BUDGET})")
    return code


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest(sys.argv[2])
    elif arg == "next":
        rows = MANIFESTS[sys.argv[2]]
        todo = [m[0] for m in rows if m[0] not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    elif arg == "dev":
        sep = sys.argv.index("--")
        sys.exit(dev(sys.argv[2], sys.argv[3], sys.argv[sep + 1:]))
    else:
        sys.exit(run_unit(arg))
