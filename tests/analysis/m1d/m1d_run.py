#!/usr/bin/env python3
"""M1-D runner -- can active-panel units recover better structural approximations? (frozen manifest; one process at a
time; <= 4 threads; <= 300 s per unit; 1200 s budget counted as time spent in this batch -- development runs included --
plus the declared timeout of the next unit; 6 GiB envelope). Diagnostic evaluator only: no caller run, no production change.

Order: panel-only block sequences (no truth) -> panel-only unit inventory -> [inventory frozen and committed:
tests/analysis/m1d/INVENTORY_FROZEN.tsv] -> truth partition + exhaustive distance pairs (refuses an inventory that differs
from its freeze record) -> edit distances -> nearest components and oracle compositions -> candidate construction per
coverage (reported flanks fixed) -> candidate block distances -> explicit-sequence scoring (no origin recording) -> summary.

  m1d_run.py manifest B | next B | <unit> | dev <name> <note> -- <command...>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m1d_unit_reachability_20261019"
M1B = "/Users/davide.bolognini/panvar/tests/results/m1b_structural_20261019"
M1C = "/Users/davide.bolognini/panvar/tests/results/m1c_allocation_20261019"
TRUTH_BLOCKS = f"{M1B}/runs/seq.c4/seq.sv:7.truth.fa"
FROZEN = os.path.join(HERE, "INVENTORY_FROZEN.tsv")
BUDGET = 1200
CASES = {r["case"]: r for r in csv.DictReader(open(os.path.join(HERE, "..", "m1b", "AUDIT_SET.tsv")), delimiter="\t")}
C4 = CASES["c4_30"]
MANIFEST_B = [("seq.panel", "seqpanel", "-", 30), ("inventory", "inventory", "-", 120),
              ("reach.pairs", "reachpairs", "-", 60), ("reach.dist", "dist", "reach", 240), ("reach.select", "reachselect", "-", 60),
              ("construct.c4_10", "construct", "c4_10", 60), ("construct.c4_30", "construct", "c4_30", 60),
              ("cand.dist.c4_10", "dist", "c4_10", 120), ("cand.dist.c4_30", "dist", "c4_30", 120),
              ("score.c4_10", "score", "c4_10", 120), ("score.c4_30", "score", "c4_30", 240),
              ("summary", "summary", "-", 60)]
# Amendment 1 (declared before further execution): `inventory` completed with boundaries_reliable False -- one carried
# allele (allele_75) of 117 gets different unit counts from the two anchor phases (recorded, unchanged). inventory.a1 reruns
# the same rule with --exclude-disagreeing: such alleles are excluded from the inventory and reported (positions, length,
# carriers), and boundaries are reliable for the rest if at most 5 % are excluded and every spacing is in [20, 40] kb.
# Every later unit uses inventory.a1.
MANIFEST_A1 = [("inventory.a1", "inventory", "a1", 120)]
INVENTORY_UNIT = "inventory.a1"
# Amendment 2 (declared before further execution): construct.c4_10/c4_30 failed (KeyError 'whole.truth.1') on a naming
# defect -- the reachability records are named after the M1-B truth file (truth1/truth2), construct.py expected the spell
# naming (truth.1/truth.2) -- and cand.dist.c4_10/c4_30 then failed on the missing construction (all four recorded,
# unchanged). construct.py now uses the reachability names and asserts that both truth sources carry identical sv:7 bytes;
# the selection rule, candidates and every later unit are unchanged. New units rerun construction and distances; the
# score and summary units of manifest B read the .a2 outputs.
MANIFEST_A2 = [("construct.c4_10.a2", "construct", "c4_10", 60), ("construct.c4_30.a2", "construct", "c4_30", 60),
               ("cand.dist.c4_10.a2", "dist", "c4_10", 120), ("cand.dist.c4_30.a2", "dist", "c4_30", 120)]
MANIFESTS = {"B": MANIFEST_B, "A1": MANIFEST_A1, "A2": MANIFEST_A2}


def ledger_rows():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def spent():
    return sum(float(r["wall_seconds"]) for r in ledger_rows())


def done():
    return {r["unit"] for r in ledger_rows()}


def exit_of(unit):
    return next((r["exit"] for r in ledger_rows() if r["unit"] == unit), None)


def state():
    rev = subprocess.run(["git", "-C", WT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", WT, "status", "--porcelain"], capture_output=True, text=True).stdout.strip()
    return rev, len([line for line in dirty.split("\n") if line.strip()])


def write_manifest(part):
    os.makedirs(OUT, exist_ok=True)
    rows = MANIFESTS[part]
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


def reported_total(case):
    rows = list(csv.DictReader(open(f"{M1B}/runs/call.{case}/call.equivalence_set.tsv"), delimiter="\t"))
    return next(r for r in rows if r["rank"] == "1")["score"]


def comparator_total(case):
    for line in open(f"{M1B}/runs/alts.{case}/{case}.candidates.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] == "saved_2+2":
            return f[3]
    raise KeyError("saved_2+2")


def command(unit, kind, target, timeout, d):
    py = sys.executable
    if kind == "seqpanel":
        return ["nice", f"{BUILD}/genotype_probe", "block-sequences", "--out", f"{d}/seq", "--reference", C4["reference"],
                "--model-gfa", C4["model_graph"], "--blocks", C4["target_block"]], "panel only (no truth)"
    if kind == "inventory":
        return [py, os.path.join(HERE, "inventory.py"), "--sequences", f"{OUT}/runs/seq.panel/seq", "--block",
                C4["target_block"], "--out", f"{d}/inventory"] + (["--exclude-disagreeing"] if target == "a1" else []), \
            "panel only (no truth)"
    inv = f"{OUT}/runs/{INVENTORY_UNIT}/inventory"
    alleles = f"{OUT}/runs/seq.panel/seq.{C4['target_block']}.alleles.fa"
    if kind == "reachpairs":
        if not os.path.exists(FROZEN):
            return None, "NOT RUN (condition): the inventory has no committed freeze record"
        return [py, os.path.join(HERE, "reach.py"), "pairs", "--inventory", inv, "--frozen", FROZEN, "--truth", TRUTH_BLOCKS,
                "--alleles", alleles, "--out", f"{d}/reach"], "diagnostic evaluator (truth)"
    if kind == "dist":
        pairs = f"{OUT}/runs/reach.pairs/reach.pairs.tsv" if target == "reach" else f"{OUT}/runs/construct.{target}.a2/candidates.pairs.tsv"
        return ["nice", f"{BUILD}/genotype_probe", "edit-distances", "--pairs", pairs, "--threads", "4", "--out",
                f"{d}/distances.tsv"], "diagnostic evaluator (truth)"
    if kind == "reachselect":
        return [py, os.path.join(HERE, "reach.py"), "select", "--inventory", inv, "--truth", TRUTH_BLOCKS, "--alleles", alleles,
                "--distances", f"{OUT}/runs/reach.dist/distances.tsv", "--components",
                f"{OUT}/runs/reach.pairs/reach.truth_components.fa", "--out", f"{d}/reach"], "diagnostic evaluator (truth)"
    if kind == "construct":
        return [py, os.path.join(HERE, "construct.py"), "--spell", f"{M1C}/runs/spell.{target}/spell", "--inventory", inv,
                "--reach", f"{OUT}/runs/reach.select/reach", "--truth-blocks", TRUTH_BLOCKS, "--out", f"{d}/candidates",
                "--expected-reported",
                reported_total(target), "--expected-comparator", comparator_total(target)], "oracle construction (truth)"
    if kind == "score":
        c = CASES[target]
        return ["nice", f"{BUILD}/genotype_probe", "score-set", "--out", f"{d}/score", "--reference", c["reference"],
                "--donor", c["sample"], "--run", f"{M1B}/runs/call.{target}/call", "--model-gfa", c["model_graph"],
                "--library", c["library"], "--reads1", c["reads1"], "--reads2", c["reads2"], "--threads", "4",
                "--explicit", f"{OUT}/runs/construct.{target}.a2/candidates.explicit.tsv"], "explicit-sequence scoring"
    if kind == "summary":
        return [py, os.path.join(HERE, "summary.py"), f"{OUT}/runs", f"{d}/m1d"], "-"
    raise ValueError(kind)


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
