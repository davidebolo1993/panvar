#!/usr/bin/env python3
"""M1-E runner -- residual-sequence sensitivity and allocation robustness (frozen manifest; one process at a time; <= 4
threads; <= 180 s per unit; 600 s budget counted as time spent in this batch -- development runs included -- plus the
declared timeout of the next unit; 6 GiB envelope). Diagnostic evaluator only: no caller run, no production change.

Order: construction per coverage (oracle restoration patch sets x {2+2, 1+3} + controls) -> explicit-sequence scoring, one
batch per coverage, no origin recording, compact --contrast summaries -> content -> summary.

  m1e_run.py manifest B | next B | <unit> | dev <name> <note> -- <command...>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m1e_residual_sensitivity_20261019"
M1B = "/Users/davide.bolognini/panvar/tests/results/m1b_structural_20261019"
M1C = "/Users/davide.bolognini/panvar/tests/results/m1c_allocation_20261019"
M1D = "/Users/davide.bolognini/panvar/tests/results/m1d_unit_reachability_20261019"
BUDGET = 600
CASES = {r["case"]: r for r in csv.DictReader(open(os.path.join(HERE, "..", "m1b", "AUDIT_SET.tsv")), delimiter="\t")}
PATCHES = ("U1", "T1", "U2", "T2", "all")
CONTRASTS = ([f"a22_{p},a22_none" for p in PATCHES] + [f"e13_{p},e13_none" for p in PATCHES]
             + [f"e13_{p},a22_{p}" for p in ("none",) + PATCHES] + ["a22_none,reported", "truthblock,reported"])
MANIFEST_B = [("construct.c4_10", "construct", "c4_10", 30), ("construct.c4_30", "construct", "c4_30", 30),
              ("score.c4_10", "score", "c4_10", 90), ("score.c4_30", "score", "c4_30", 180),
              ("content", "content", "-", 60), ("summary", "summary", "-", 30)]
MANIFESTS = {"B": MANIFEST_B}


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


def command(unit, kind, target, timeout, d):
    py = sys.executable
    if kind == "construct":
        return [py, os.path.join(HERE, "construct.py"), "--inventory", f"{M1D}/runs/inventory.a1/inventory", "--structure",
                f"{M1D}/runs/reach.select/reach.composed_structure.tsv", "--truth-components",
                f"{M1D}/runs/reach.pairs/reach.truth_components.fa", "--truth-blocks", f"{M1B}/runs/seq.c4/seq.sv:7.truth.fa",
                "--composed", f"{M1D}/runs/reach.select/reach.composed.fa", "--spell", f"{M1C}/runs/spell.{target}/spell",
                "--expected-reported", reported_total(target), "--out", f"{d}/candidates"], "oracle restoration construction"
    if kind == "score":
        c = CASES[target]
        return ["nice", f"{BUILD}/genotype_probe", "score-set", "--out", f"{d}/score", "--reference", c["reference"],
                "--donor", c["sample"], "--run", f"{M1B}/runs/call.{target}/call", "--model-gfa", c["model_graph"],
                "--library", c["library"], "--reads1", c["reads1"], "--reads2", c["reads2"], "--threads", "4",
                "--explicit", f"{OUT}/runs/construct.{target}/candidates.explicit.tsv", "--contrast", ";".join(CONTRASTS)], \
            "explicit-sequence scoring"
    if kind == "content":
        return [py, os.path.join(HERE, "content.py"), f"{OUT}/runs/construct.c4_30/candidates.blocks.fa", f"{d}/content.tsv"], "-"
    if kind == "summary":
        return [py, os.path.join(HERE, "summary.py"), f"{OUT}/runs", f"{d}/m1e"], "-"
    raise ValueError(kind)


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in MANIFESTS.items() for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and timeout <= 180 and unit not in done()
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
    assert used + 180 <= BUDGET, f"development run refused: {used:.1f} s + 180 s > {BUDGET} s"
    d = f"{OUT}/dev/{name}"
    assert not os.path.exists(d), f"{d} exists"
    os.makedirs(d)
    code, timed_out, wall, rss, text = timed(cmd, d, 180)
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
