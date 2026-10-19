#!/usr/bin/env python3
"""M1-B runner -- structural-error attribution on the frozen current caller (frozen manifest; one process at a time;
4 threads; <= 300 s per unit; 1200 s budget counted as time spent in this batch -- development runs included -- plus the
declared timeout of the next unit; 6 GiB envelope).

Caller units: the current frozen caller (233cb65 code; panvar binary = the M1-A baseline binary) on the validated combined
route, with each case's OWN inputs (graph, reads, library.txt, curation) and the frozen baseline composition for its
locus: base profile = the P1 perf base of the matching baseline cell (C4: c4_30, whose science equals the historical v6
campaign profile apart from max_physical_origins, see DESIGN.md; CYP2D6: cyp_5 / cyp_30, = the pilot's a4 profile apart
from the guards), composed in mode `full` with profiles/genotype_dev_performance_v1.tsv; overrides time_limit_seconds =
declared timeout - 20 s and footprint_limit_bytes = 6 GiB; plus early_incumbent experimental_v1, union_fragment_blocks
fine_core_4096_experimental, repair_round1_origin_admission reuse_only_experimental; --preset short-read; --threads 4;
--output-schema sequence_first_v2; no --dosage-mode. No other sample's profile; no library parameter changed.

Evaluator units (truth only here): genotype_probe evaluate (unchanged), block-sequences (diagnostic; new in this batch),
structure.py (unit copies), outcome.py (historical vs current).

  m1b_run.py manifest B | next B | <unit> | dev <name> <note> -- <command...>
"""
import csv
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
import compose_profile  # noqa: E402
import p1_run  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/m1b_structural_20261019"
BASELINE = os.path.join(HERE, "..", "m1a", "baseline")
BUDGET = 1200
GIB6 = 6442450944
CASES = {r["case"]: r for r in csv.DictReader(open(os.path.join(HERE, "AUDIT_SET.tsv")), delimiter="\t")}
MANIFEST_B = [("check.inputs", "inputs", "-", 60),
              ("call.c4_10", "caller", "c4_10", 120), ("eval.c4_10", "eval", "c4_10", 60),
              ("call.c4_30", "caller", "c4_30", 300), ("eval.c4_30", "eval", "c4_30", 60),
              ("call.cyp_5", "caller", "cyp_5", 90), ("eval.cyp_5", "eval", "cyp_5", 60),
              ("call.cyp_30", "caller", "cyp_30", 300), ("eval.cyp_30", "eval", "cyp_30", 60),
              ("seq.c4", "blockseq", "c4_30", 60), ("seq.cyp", "blockseq", "cyp_5", 60),
              ("structure.c4_10", "structure", "c4_10", 30), ("structure.c4_30", "structure", "c4_30", 30),
              ("structure.cyp_5", "structure", "cyp_5", 30), ("structure.cyp_30", "structure", "cyp_30", 30),
              ("outcome", "outcome", "-", 30)]
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


def curation_args(c):
    return [] if c["curation"] == "none" else ["--panel-curation", c["curation"]]


def caller_command(c, timeout, d):
    prof = compose_profile.compose(f"{d}/profile.{c['case']}.tsv", f"m1b-{c['case']}", c["baseline_base_profile"], "full", [
        f"time_limit_seconds={timeout - 20}:declared unit timeout minus 20 s (records TIME_LIMIT before the kill)",
        f"footprint_limit_bytes={GIB6}:6 GiB local envelope (M1-B)"])
    extra = f"{d}/profile.m1b.tsv"
    open(extra, "w").write("profile\tm1b-combined\nearly_incumbent\texperimental_v1\n"
                           "union_fragment_blocks\tfine_core_4096_experimental\n"
                           "repair_round1_origin_admission\treuse_only_experimental\n")
    return ["nice", f"{BUILD}/panvar", "genotype", "--preset", "short-read", "--profile", prof, "--profile", extra,
            "--threads", "4", "--graph", c["model_graph"], "--reference-path", c["reference"], *curation_args(c),
            "--reads1", c["reads1"], "--reads2", c["reads2"], "--library", c["library"], "--out", f"{d}/call", "--quiet",
            "--output-schema", "sequence_first_v2"]


def command(unit, kind, target, timeout, d):
    if kind == "inputs":
        return [sys.executable, os.path.join(HERE, "freeze_audit_set.py"), d, HERE, f"{BUILD}/panvar",
                os.path.join(BASELINE, "baseline_manifest.tsv")], "-"
    if kind == "outcome":
        return [sys.executable, os.path.join(HERE, "outcome.py"), os.path.join(HERE, "AUDIT_SET.tsv"), f"{OUT}/runs",
                os.path.join(BASELINE, "baseline_settings.tsv"), f"{d}/m1b"], "-"
    c = CASES[target]
    if kind == "caller":
        return caller_command(c, timeout, d), "current frozen caller, combined route"
    if kind == "eval":
        run = f"{OUT}/runs/call.{target}/call"
        if not os.path.exists(run + ".mosaic_status.tsv"):
            return None, "NOT RUN (condition): the caller wrote no status"
        return ["nice", f"{BUILD}/genotype_probe", "evaluate", "--scheme", "loo", "--donor", c["sample"], "--reference",
                c["reference"], "--run-prefix", run, "--model-gfa", c["model_graph"], "--truth-gfa", c["truth_graph"],
                *curation_args(c), "--out", f"{d}/eval"], "diagnostic evaluator (truth)"
    if kind == "blockseq":
        return ["nice", f"{BUILD}/genotype_probe", "block-sequences", "--out", f"{d}/seq", "--reference", c["reference"],
                "--model-gfa", c["model_graph"], *curation_args(c), "--blocks", c["target_block"], "--truth-gfa",
                c["truth_graph"], "--donor", c["sample"]], "diagnostic evaluator (truth)"
    if kind == "structure":
        locus = "c4" if target.startswith("c4") else "cyp"
        fresh = f"{OUT}/runs/eval.{target}/eval.blocks.tsv"
        return [sys.executable, os.path.join(HERE, "structure.py"), "--case", target, "--block", c["target_block"],
                "--sequences", f"{OUT}/runs/seq.{locus}/seq", "--hist-eval", c["historical_eval"] + ".blocks.tsv",
                "--fresh-eval", fresh, "--out", f"{d}/{target}"], "diagnostic evaluator (truth)"
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
