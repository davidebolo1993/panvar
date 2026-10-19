#!/usr/bin/env python3
"""M6-E runner -- panel-only C4 sv:7 sequence-variation measurements (no reads, no fitting) (DESIGN.md; frozen manifests; one process
at a time; 4 threads; 300 s counted as time spent plus the declared timeout of the next unit).

  m6e_run.py manifest P | next P | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m6e_panel_variation_20261019"
M1D = "/Users/davide.bolognini/panvar/tests/results/m1d_unit_reachability_20261019/runs"
BUDGET = 300
MANIFEST_P = [("prepare", "prepare", "-", 30), ("align", "align", "-", 120), ("analyse", "analyse", "-", 90)]
# Amendment 1 (declared after `analyse` exited 1 at 3.5 s spent): the input guard refused the script's read-back of its own
# summary for the final print (all outputs were already written); reads under the unit's own directory are allowed. analyse.a1.
MANIFEST_A1 = [("analyse.a1", "analyse", "-", 90)]
# Amendment 2 (declared after analyse.a1, 5.2 s spent): panel-only measurements for the state-count estimate and structural
# layer (length classes, arrangements, local distinct patterns per window, non-simple interval signatures). extra.a2.
MANIFEST_A2 = [("extra.a2", "extra", "-", 60)]
MANIFESTS = {"P": MANIFEST_P, "A1": MANIFEST_A1, "A2": MANIFEST_A2}
INPUTS = [   # panel-only inputs (frozen digests); nothing else is read
    {"file": f"{M1D}/inventory.a1/inventory.tsv", "sha256": "e890464ef2a1d6b91cb15b2d53a1577ddcdc032ec85f698033a1cc57a99b1039"},
    {"file": f"{M1D}/inventory.a1/inventory.fa", "sha256": "8f425a37e811428616d1d40b707fe6813e0f710b79ef36cabcb44215deb470c4"},
    {"file": f"{M1D}/inventory.a1/inventory.anchors.tsv", "sha256": "8b9272d23584a6012ba1e8e7b27ee6f03ccb7915b4f38628fcc9a6fa5974b8ed"},
    {"file": f"{M1D}/inventory.a1/inventory.log.tsv", "sha256": "1fdfc9515bfc6caf7372bb88e0535c22c2bb02325497b16bc686734e6c3e80f2"},
    {"file": f"{M1D}/seq.panel/seq.sv:7.alleles.fa", "sha256": "45a12a8c99e7c861d63e5df5c1f6ffe120253e5344dfd5f41ddeb9a467fcf662"},
    {"file": f"{M1D}/seq.panel/seq.block_sequences.tsv", "sha256": "7eb86610677a445059283f7056d8e3e19afae6bd893dcff0eb34f294e5961d1c"},
    {"file": f"{M1D}/seq.panel/run.log", "sha256": "a278fc12d616ecfa4c9c5db459414ef9b86f19cb2f6f59799f9dc6be2a979893"},
]


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
    ok = lambda u: next((r["exit"] for r in ledger_rows() if r["unit"] == u), None) == "0"
    if kind == "prepare":
        return [sys.executable, os.path.join(HERE, "prepare.py"), d], "verify bindings, component table, references, pairs"
    if kind == "align":
        assert ok("prepare"), "align runs only after prepare passed"
        return [f"{BUILD}/unit_correspondence", "--pairs", f"{OUT}/runs/prepare/pairs.tsv", "--out", f"{d}/events.tsv", "--threads", "4"], "edlib alignment events"
    if kind == "analyse":
        assert ok("align"), "analyse runs only after align passed"
        return [sys.executable, os.path.join(HERE, "analyse.py"), f"{OUT}/runs/prepare", f"{OUT}/runs/align/events.tsv", f"{d}/m6e"], "measurements"
    if kind == "extra":
        assert ok("analyse.a1"), "extra runs only after analyse.a1 passed"
        return [sys.executable, os.path.join(HERE, "extra.py"), f"{OUT}/runs/prepare", f"{OUT}/runs/align/events.tsv", f"{OUT}/runs/analyse.a1/m6e", f"{d}/m6e_extra"], "state-count measurements"
    raise ValueError(kind)


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in reversed(list(MANIFESTS.items())) for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and unit not in done()
    rev, dirty = state()
    assert dirty == 0, "the tree must be committed and clean"
    used = spent()
    verify()
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
