#!/usr/bin/env python3
"""M6-A runner -- bounded 31-mer count-model comparator (diagnostic prototype; DESIGN.md; frozen manifests; one process at a
time; 4 threads where used; 600 s counted as time spent plus the declared timeout of the next unit; 6 GiB envelope).

  m6a_run.py manifest G|C|E | next G|C|E | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m6a_count_model_20261019"
BUDGET = 600
C4 = "/Users/davide.bolognini/panvar/tests/results/genotype_c4_production_11fe681"
CYPD = "/Users/davide.bolognini/panvar/tests/results/cyp2d6_shortread_pilot_v1/data"
M4B = "/Users/davide.bolognini/panvar/tests/results/m4b_distance_prior_20261019/runs"
M1B = "/Users/davide.bolognini/panvar/tests/results/m1b_structural_20261019/runs"
REF_C4, REF_CYP = "HG00096#1#haplotype1-0000024:31848049-32074186", "grch38#1#chr22:42031864-42245566"
CYP_CUR = "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"
CELL = {   # controls first
    "c4_lzo_hg00171_5": dict(gfa=f"{C4}/graphs/c4.gfa", ref=REF_C4, cur=None, reads=f"{C4}/reads/HG00171/cov5", run=f"{M4B}/run.c4_lzo_hg00171_5.const/call",
                             truth=f"{C4}/graphs/c4.gfa", donor="HG00171", block="sv:7"),
    "cyp_lzo_hg02769_5": dict(gfa=f"{CYPD}/graphs/cyp2d6.gfa", ref=REF_CYP, cur=CYP_CUR, reads=f"{CYPD}/reads/HG02769/cov5", run=f"{M4B}/run.cyp_lzo_hg02769_5.const/call",
                              truth=f"{CYPD}/graphs/cyp2d6.gfa", donor="HG02769", block="sv:9"),
    "c4_loo_na19650_10": dict(gfa=f"{C4}/graphs/NA19650.loo.gfa", ref=REF_C4, cur=None, reads=f"{C4}/reads/NA19650/cov10", run=f"{M1B}/call.c4_10/call",
                              truth=f"{C4}/graphs/c4.gfa", donor="NA19650", block="sv:7"),
    "cyp_loo_na19650_5": dict(gfa=f"{CYPD}/graphs/loo.NA19650.gfa", ref=REF_CYP, cur=CYP_CUR, reads=f"{CYPD}/reads/NA19650/cov5", run=f"{M1B}/call.cyp_5/call",
                              truth=f"{CYPD}/graphs/cyp2d6.gfa", donor="NA19650", block="sv:9"),
}
ORDER = list(CELL)
TO = {"c4_lzo_hg00171_5": 120, "cyp_lzo_hg02769_5": 90, "c4_loo_na19650_10": 150, "cyp_loo_na19650_5": 90}
MANIFEST_G = [("gate.count_fit", "exe", "count_fit_gates", 60)]
MANIFEST_C = [(f"countfit.{c}", "countfit", c, TO[c]) for c in ORDER] + [("freeze", "freeze", "-", 30)]
MANIFEST_E = ([(f"blockseq.{c}", "blockseq", c, 60) for c in ORDER] + [(f"units.{c}", "units", c, 60) for c in ORDER] +
              [(f"fragscore.{c}", "fragscore", c, 120) for c in ORDER] + [("summary", "summary", "-", 30)])
MANIFESTS = {"G": MANIFEST_G, "C": MANIFEST_C, "E": MANIFEST_E}
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


def model(c):
    return ["--model-gfa", c["gfa"], "--reference", c["ref"], *(["--panel-curation", c["cur"]] if c["cur"] else [])]


def reads(c):
    return ["--reads1", f"{c['reads']}/reads_1.fq", "--reads2", f"{c['reads']}/reads_2.fq"]


def command(kind, target, d):
    verify()
    R = f"{OUT}/runs"
    if kind == "exe":
        return [f"{BUILD}/{target}"], "synthetic count-fit gates"
    if kind == "countfit":
        c = CELL[target]
        return ["nice", f"{BUILD}/genotype_probe", "count-fit", "--out", f"{d}/cf", *model(c), "--run", c["run"], "--block", c["block"], *reads(c),
                "--max-candidates", "40000"], "count-fit comparator (truth-free)"
    if kind == "freeze":
        return [sys.executable, os.path.join(HERE, "select_eval.py"), R, d, *ORDER], "count-model freeze + evaluation selection (truth-free)"
    if kind == "blockseq":
        assert os.path.exists(f"{R}/freeze/count_selection.sha256"), "truth needs the frozen count-model selections"
        c = CELL[target]
        return ["nice", f"{BUILD}/genotype_probe", "block-sequences", "--out", f"{d}/seq", *model(c), "--blocks", c["block"], "--truth-gfa", c["truth"],
                "--donor", c["donor"]], "block alleles + truth block bytes (evaluation only)"
    if kind == "units":
        return [sys.executable, os.path.join(HERE, "units.py"), f"{R}/blockseq.{target}/seq", CELL[target]["block"], f"{d}/units.tsv"], \
            "repeat-unit copies by the M1-B rule (evaluation only)"
    if kind == "fragscore":
        c = CELL[target]
        return ["nice", f"{BUILD}/genotype_probe", "score-set", "--out", f"{d}/fs", *model(c), "--donor", c["donor"], "--run", c["run"], *reads(c),
                "--library", f"{c['reads']}/library.txt", "--threads", "4", "--candidates", f"{R}/freeze/eval_candidates.{target}.tsv",
                "--truth-gfa", c["truth"]], "fragment-model scores + sequence excess of the frozen selection (evaluation)"
    if kind == "summary":
        return [sys.executable, os.path.join(HERE, "summary.py"), R, f"{OUT}/ledger.tsv", f"{d}/m6a", *[f"{c}:{CELL[c]['block']}" for c in ORDER]], "summary"
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
    print("\n".join(l for l in text.splitlines() if not l.startswith(" ") and not l.startswith("$"))[-3500:])
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
