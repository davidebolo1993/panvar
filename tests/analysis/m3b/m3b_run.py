#!/usr/bin/env python3
"""M3-B runner -- bounded two-block search audit: synthetic gate only (frozen manifest; one process at a time; <= 4
threads; <= 180 s per unit; 600 s budget counted as time spent in this batch plus the declared timeout of the next unit;
6 GiB envelope). No real-data execution.

  m3b_run.py manifest B | next B | <unit> | dev <name> <note> -- <command...>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m3b_two_block_20261019"
BUDGET = 600
C4 = "/Users/davide.bolognini/panvar/tests/results/genotype_c4_production_11fe681"
REF = "HG00096#1#haplotype1-0000024:31848049-32074186"
MOS = f"{C4}/mosaic_inputs/one_homologue_once_1"
CELL = {
    "mosaic_2x": dict(run=f"{C4}/mosaic/one_homologue_once_1/call.2", gfa=f"{MOS}/staged.gfa", curation=f"{C4}/mosaic_inputs/mosaic_exclude.tsv",
                      reads=f"{MOS}/cov2", donor="MOSAIC", truth=["--truth-tsv", f"{MOS}/truth.tsv"]),
    "mosaic_5x": dict(run=f"{C4}/mosaic/one_homologue_once_1/call.5", gfa=f"{MOS}/staged.gfa", curation=f"{C4}/mosaic_inputs/mosaic_exclude.tsv",
                      reads=f"{MOS}/cov5", donor="MOSAIC", truth=["--truth-tsv", f"{MOS}/truth.tsv"]),
    "lzo_5x": dict(run=f"{C4}/lzo/HG00096/call.5", gfa=f"{C4}/graphs/c4.gfa", curation="none", reads=f"{C4}/reads/HG00096/cov5", donor="HG00096",
                   truth=["--truth-templates", "HG00096#1#haplotype1-0000024:31848049-32074186,HG00096#2#haplotype2-0000117:31856177-32049825"]),
}
CELLS = list(CELL)
MANIFEST_B = ([(f"cat.{c}", "cat", c, 60) for c in CELLS] + [(f"score.{c}", "score", c, 180) for c in CELLS] +
              [("rank", "rank", "-", 30)] + [(f"eval.{c}", "eval", c, 60) for c in CELLS])
# Amendment 1: 2x scoring took 126 s for 1,849 candidates, so a 5x cell projects to ~300 s (> the 180-s unit cap). Each 5x
# catalogue is split deterministically into two parts (labels in order, halves), each carrying the incumbent and its
# expected total; the budget then admits one complete 5x control: the mosaic 5x twin. lzo_5x scoring is NOT run (gap).
MANIFEST_A1 = [("split.mosaic_5x", "split", "mosaic_5x", 30), ("score.mosaic_5x.part1", "score_part", "mosaic_5x:1", 180),
               ("score.mosaic_5x.part2", "score_part", "mosaic_5x:2", 180), ("rank.a1", "rank", "-", 30),
               ("eval.a1.mosaic_2x", "eval", "mosaic_2x", 60), ("eval.a1.mosaic_5x", "eval", "mosaic_5x", 60)]
MANIFESTS = {"B": MANIFEST_B, "A1": MANIFEST_A1}
CARRIED = 0.0
INPUTS = {r["cell"] + ":" + r["file"]: r["sha256"] for r in csv.DictReader(open(os.path.join(HERE, "INPUTS.tsv")), delimiter="\t")}


def ledger_rows():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def spent():
    return CARRIED + sum(float(r["wall_seconds"]) for r in ledger_rows())


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


def verify(cell):
    for key, want in INPUTS.items():
        c, path = key.split(":", 1)
        if c == cell:
            assert hashlib.sha256(open(path, "rb").read()).hexdigest() == want, f"{cell}: {path} differs from the frozen input"


def common(c):
    return ["--model-gfa", c["gfa"], "--reference", REF, *([] if c["curation"] == "none" else ["--panel-curation", c["curation"]])]


def score_files(cell):
    whole = f"{OUT}/runs/score.{cell}/score.scores_truth_free.tsv"
    parts = [f"{OUT}/runs/score.{cell}.part{i}/score.scores_truth_free.tsv" for i in (1, 2)]
    if os.path.exists(whole):
        return [whole]
    return parts if all(os.path.exists(p) for p in parts) else []


def command(unit, kind, target, timeout, d):
    if kind == "split":
        verify(target)
        return [sys.executable, os.path.join(HERE, "split_catalogue.py"), f"{OUT}/runs/cat.{target}/cat.catalogue.tsv", f"{d}/part"], \
            "deterministic two-part split (Amendment 1)"
    if kind == "score_part":
        cell, part = target.split(":")
        verify(cell)
        c = CELL[cell]
        return ["nice", f"{BUILD}/genotype_probe", "score-set", "--out", f"{d}/score", "--donor", c["donor"], "--run", c["run"], *common(c),
                "--candidates", f"{OUT}/runs/split.{cell}/part{part}.tsv", "--library", f"{c['reads']}/library.txt",
                "--reads1", f"{c['reads']}/reads_1.fq", "--reads2", f"{c['reads']}/reads_2.fq", "--threads", "4"], \
            f"unchanged reference scorer, catalogue part {part} of 2 (Amendment 1)"
    if kind in ("cat", "score", "eval"):
        verify(target)
        c = CELL[target]
    if kind == "cat":
        return ["nice", f"{BUILD}/genotype_probe", "two-block-catalogue", "--out", f"{d}/cat", "--run", c["run"], *common(c),
                "--reads1", f"{c['reads']}/reads_1.fq", "--reads2", f"{c['reads']}/reads_2.fq", "--library", f"{c['reads']}/library.txt",
                "--top", "4", "--limit", "2048"], "two-block catalogue (truth-free; top 4; limit 2,048)"
    if kind == "score":
        cat = f"{OUT}/runs/cat.{target}/cat.catalogue.tsv"
        if not os.path.exists(cat):
            return None, "NOT RUN: no frozen catalogue (preflight stopped or failed)"
        return ["nice", f"{BUILD}/genotype_probe", "score-set", "--out", f"{d}/score", "--donor", c["donor"], "--run", c["run"], *common(c),
                "--candidates", cat, "--library", f"{c['reads']}/library.txt", "--reads1", f"{c['reads']}/reads_1.fq",
                "--reads2", f"{c['reads']}/reads_2.fq", "--threads", "4"], "unchanged reference scorer (score-set) on the frozen catalogue"
    if kind == "rank":
        done_cells = [c for c in CELLS if score_files(c)]
        return [sys.executable, os.path.join(HERE, "rank.py"), f"{OUT}/runs", f"{d}/m3b", *done_cells], "truth-free selection (frozen before evaluation)"
    if kind == "eval":
        rank_dir = f"{OUT}/runs/rank.a1" if unit.startswith("eval.a1.") else f"{OUT}/runs/rank"
        assert os.path.exists(f"{rank_dir}/m3b.selection.sha256"), "evaluation needs the frozen selection"
        sc = f"{rank_dir}/merged.{target}.tsv"
        if not os.path.exists(sc):
            return None, "NOT RUN: not scored"
        labels = set(["incumbent"])
        for r in csv.DictReader(open(f"{rank_dir}/m3b.selection.tsv"), delimiter="\t"):
            if r["cell"] == target:
                labels.add(r["winner"])
        for r in csv.DictReader(open(f"{rank_dir}/m3b.decomposition.tsv"), delimiter="\t"):
            if r["cell"] == target and r["label"] != "absent":
                labels.add(r["label"])
        for r in csv.DictReader(open(f"{rank_dir}/m3b.top.tsv"), delimiter="\t"):
            if r["cell"] == target:
                labels.add(r["label"])
        return ["nice", f"{BUILD}/genotype_probe", "two-block-eval", "--out", f"{d}/eval", *common(c), "--scores", sc, *c["truth"],
                "--ranks", f"{OUT}/runs/cat.{target}/cat.emission_ranks.tsv", "--labels", ",".join(sorted(labels))], "evaluation (truth; after selection)"
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
