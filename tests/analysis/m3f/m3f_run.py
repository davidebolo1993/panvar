#!/usr/bin/env python3
"""M3-F runner -- bounded validation of joint search (off vs repair_joint_search experimental_v1, support v2 in both arms)
beyond the original recovery (SELECTION.md; frozen manifests; one process at a time; 4 threads; 1,200 s budget counted as
time spent in this batch plus the declared timeout of the next unit; 6 GiB envelope). No algorithm change.

  m3f_run.py manifest C|E|D | next C|E|D | <unit> | dev <name> <note> -- <command...>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m3f_validation_20261019"
BUDGET = 1200
C4 = "/Users/davide.bolognini/panvar/tests/results/genotype_c4_production_11fe681"
REF = "HG00096#1#haplotype1-0000024:31848049-32074186"
CUR = f"{C4}/mosaic_inputs/mosaic_exclude.tsv"
M1B_SEQ = "/Users/davide.bolognini/panvar/tests/results/m1b_structural_20261019/runs/seq.c4/seq"


def mosaic(plan, cov, role):
    m = f"{C4}/mosaic_inputs/{plan}"
    return dict(campaign=f"{C4}/mosaic/{plan}/call.{cov}", gfa=f"{m}/staged.gfa", curation=CUR, reads=f"{m}/cov{cov}", scheme="mosaic",
                donor="MOSAIC", truth_gfa=f"{m}/staged.gfa", truth_tsv=f"{m}/truth.tsv", role=role)


CELL = {
    "lzo_hg00171_5": dict(campaign=f"{C4}/lzo/HG00171/call.5", gfa=f"{C4}/graphs/c4.gfa", curation="none", reads=f"{C4}/reads/HG00171/cov5",
                          scheme="lzo", donor="HG00171", truth_gfa=f"{C4}/graphs/c4.gfa", truth_tsv=None, role="control: ordinary LZO, exact"),
    "ts1_5": mosaic("three_switches_1", 5, "control: GE twin of ts1_2 (no exact-sequence twin exists)"),
    "ds1_10": mosaic("dense_switching_1", 10, "control: GE twin of ds1_2 (no exact-sequence twin exists)"),
    "ds1_2": mosaic("dense_switching_1", 2, "error: low-coverage mosaic, positive excess"),
    "ts1_2": mosaic("three_switches_1", 2, "error: low-coverage mosaic, positive excess"),
    "loo_na19650_10": dict(campaign=f"{C4}/loo/NA19650/call.10", gfa=f"{C4}/graphs/NA19650.loo.gfa", curation="none",
                           reads=f"{C4}/reads/NA19650/cov10", scheme="loo", donor="NA19650", truth_gfa=f"{C4}/graphs/c4.gfa", truth_tsv=None,
                           role="error: LOO structural (M1-B attribution set)"),
}
CONTROLS = ["lzo_hg00171_5", "ts1_5", "ds1_10"]
ERRORS = ["ds1_2", "ts1_2", "loo_na19650_10"]
PROFILE = {"off": os.path.join(HERE, "profiles", "support_v2_off.tsv"), "on": os.path.join(HERE, "profiles", "support_v2_on.tsv")}
TIMEOUT = {"lzo_hg00171_5": (60, 120), "ts1_5": (100, 160), "ds1_10": (180, 240),
           "ds1_2": (70, 120), "ts1_2": (70, 120), "loo_na19650_10": (120, 180)}
MANIFEST_C = ([(f"run.{c}.{arm}", "run", f"{c}:{arm}", TIMEOUT[c][i]) for c in CONTROLS for i, arm in enumerate(("off", "on"))] +
              [("control.check", "cmp", "control", 30)] +
              [(f"eval.{c}.{arm}", "eval", f"{c}:{arm}", 40) for c in CONTROLS for arm in ("off", "on")] +
              [("control.gate", "acc", "control", 30)])
MANIFEST_E = ([(f"run.{c}.{arm}", "run", f"{c}:{arm}", TIMEOUT[c][i]) for c in ERRORS for i, arm in enumerate(("off", "on"))] +
              [("compare", "cmp", "all", 30)] +
              [(f"eval.{c}.{arm}", "eval", f"{c}:{arm}", 40) for c in ERRORS for arm in ("off", "on")] +
              [("accuracy", "acc", "all", 30)])
MANIFEST_D = ([(f"lik.{c}", "lik", c, 120) for c in ("ds1_2", "ts1_2")] +
              [(f"tbeval.{c}", "tbeval", c, 40) for c in ("ds1_2", "ts1_2")] +
              [("structure.loo_na19650_10", "structure", "loo_na19650_10", 30), ("attribution", "attr", "-", 30)])
MANIFESTS = {"C": MANIFEST_C, "E": MANIFEST_E, "D": MANIFEST_D}
CARRIED = 0.0
INPUTS = [r for r in csv.DictReader(open(os.path.join(HERE, "INPUTS.tsv")), delimiter="\t")]


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
    for r in INPUTS:
        if r["cell"] == cell:
            assert hashlib.sha256(open(r["file"], "rb").read()).hexdigest() == r["sha256"], f"{cell}: {r['file']} differs from the frozen input"


def run_dir(cell, arm):
    return f"{OUT}/runs/run.{cell}.{arm}"


def eval_value(cell, arm, key):
    p = f"{OUT}/runs/eval.{cell}.{arm}/eval.summary.tsv"
    if not os.path.exists(p):
        return None
    return next((line.rstrip("\n").split("\t")[1] for line in open(p) if line.split("\t")[0] == key), None)


def status_value(cell, arm, key):
    p = f"{run_dir(cell, arm)}/call.mosaic_status.tsv"
    if not os.path.exists(p):
        return None
    return next((line.rstrip("\n").split("\t")[1] for line in open(p) if line.split("\t")[0] == key), None)


def command(unit, kind, target, timeout, d):
    if kind == "run":
        cell, arm = target.split(":")
        if cell in ERRORS:
            assert exit_of("control.gate") == "0", "error cells run only after the control gate passes (stop rule)"
        verify(cell)
        verify("profiles")
        c = CELL[cell]
        cmd = ["nice", f"{BUILD}/panvar", "genotype", "--profile", PROFILE[arm], "--threads", "4", "--gfa", c["gfa"], "--reference-path", REF,
               "--reads1", f"{c['reads']}/reads_1.fq", "--reads2", f"{c['reads']}/reads_2.fq", "--library", f"{c['reads']}/library.txt",
               "--out-prefix", f"{d}/call", "--quiet"]
        if c["curation"] != "none":
            cmd += ["--panel-curation", c["curation"]]
        return cmd, f"caller ({arm})"
    if kind == "cmp":
        cells = CONTROLS if target == "control" else CONTROLS + ERRORS
        return [sys.executable, os.path.join(HERE, "compare_m3f.py"), f"{OUT}/runs", f"{OUT}/ledger.tsv", f"{d}/m3f", *cells], \
            "truth-free matched comparison + prediction freeze"
    if kind == "eval":
        cell, arm = target.split(":")
        freeze = "control.check" if cell in CONTROLS else "compare"
        assert os.path.exists(f"{OUT}/runs/{freeze}/m3f.selection.sha256"), "evaluation needs the frozen predictions"
        if not os.path.exists(f"{run_dir(cell, arm)}/call.mosaic_status.tsv"):
            return None, "NOT RUN (condition): the caller wrote no status"
        c = CELL[cell]
        cmd = ["nice", f"{BUILD}/genotype_probe", "evaluate", "--scheme", c["scheme"], "--donor", c["donor"], "--reference", REF,
               "--run-prefix", f"{run_dir(cell, arm)}/call", "--model-gfa", c["gfa"], "--truth-gfa", c["truth_gfa"], "--out", f"{d}/eval"]
        if c["curation"] != "none":
            cmd += ["--panel-curation", c["curation"]]
        return cmd, "truth evaluation (after the frozen predictions)"
    if kind == "acc":
        cells = CONTROLS if target == "control" else CONTROLS + ERRORS
        return [sys.executable, os.path.join(HERE, "accuracy_m3f.py"), f"{OUT}/runs", target, f"{d}/m3f", *cells], \
            "accuracy per arm (truth); control stop rule" if target == "control" else "accuracy per arm (truth)"
    if kind == "lik":
        # diagnostic only, after the freeze: the frozen shortlist rule (k = 4, 2,048 bounds, repair configuration) from the OFF
        # run's incumbent, for an error cell whose on arm still has positive excess and whose joint pass completed
        assert os.path.exists(f"{OUT}/runs/compare/m3f.selection.sha256")
        excess, js = eval_value(target, "on", "excess_over_floor_total"), status_value(target, "on", "joint_search_status") or "-"
        if excess is None or excess in ("NA", "0"):
            return None, f"NOT RUN (condition): on-arm excess {excess}"
        if not js.startswith("complete"):
            return None, f"NOT RUN (condition): joint pass {js} (classified as pass skipped or incomplete)"
        c = CELL[target]
        return ["nice", f"{BUILD}/genotype_probe", "two-block-lik", "--out", f"{d}/lik", "--run", f"{run_dir(target, 'off')}/call",
                "--model-gfa", c["gfa"], "--reference", REF, "--panel-curation", c["curation"],
                "--reads1", f"{c['reads']}/reads_1.fq", "--reads2", f"{c['reads']}/reads_2.fq", "--library", f"{c['reads']}/library.txt",
                "--top", "4", "--single-limit", "2048", "--limit", "2048", "--threads", "4", "--config", "repair"], \
            "diagnostic: shortlist + joint catalogue from the off incumbent (truth-free)"
    if kind == "tbeval":
        sc = f"{OUT}/runs/lik.{target}/lik.scores.tsv"
        if not os.path.exists(sc):
            return None, "NOT RUN (condition): no diagnostic catalogue"
        c = CELL[target]
        return ["nice", f"{BUILD}/genotype_probe", "two-block-eval", "--out", f"{d}/eval", "--model-gfa", c["gfa"], "--reference", REF,
                "--panel-curation", c["curation"], "--scores", sc, "--truth-tsv", c["truth_tsv"], "--ranks", f"{OUT}/runs/lik.{target}/lik.singles.tsv",
                "--labels", "incumbent"], "diagnostic (truth): is the truth in the catalogue, and where"
    if kind == "structure":
        offe, one = (f"{OUT}/runs/eval.{target}.{arm}/eval.blocks.tsv" for arm in ("off", "on"))
        if not (os.path.exists(offe) and os.path.exists(one)):
            return None, "NOT RUN (condition): an arm was not evaluated"
        verify("m1b_seq")
        # structure.py labels its two tables 'historical' and 'fresh': here historical = OFF arm, fresh = ON arm
        return [sys.executable, os.path.join(WT, "tests", "analysis", "m1b", "structure.py"), "--case", target, "--block", "sv:7",
                "--sequences", M1B_SEQ, "--hist-eval", offe, "--fresh-eval", one, "--out", f"{d}/{target}"], \
            "diagnostic (truth): sv:7 copy allocation, off (historical column) vs on (fresh column)"
    if kind == "attr":
        return [sys.executable, os.path.join(HERE, "attribute_m3f.py"), f"{OUT}/runs", f"{d}/m3f", *ERRORS], "attribution of unchanged errors (diagnostic)"
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
    cmd, note = command(unit, kind, target, timeout, d)
    if cmd is None:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, note)
        print(f"{unit}: {note}")
        return 3
    code, timed_out, wall, rss, text = timed(cmd, d, timeout)
    record(unit, rev, dirty, timeout, wall, rss / 2**30, code, timed_out, note)
    print(text[-2500:] if kind in ("cmp", "acc", "attr") else "", end="")
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s, {note} (spent {spent():.1f}/{BUDGET})")
    return 0


def dev(name, note, cmd):
    """A development run: counted against the budget, recorded, never a result."""
    used = spent()
    assert used + 120 <= BUDGET, f"development run refused: {used:.1f} s + 120 s > {BUDGET} s"
    d = f"{OUT}/dev/{name}"
    assert not os.path.exists(d), f"{d} exists"
    os.makedirs(d)
    code, timed_out, wall, rss, text = timed(cmd, d, 120)
    record(name, "development", "-", "-", wall, rss / 2**30, code, timed_out, f"development run (counted): {note}")
    print(text[-4000:])
    print(f"{name}: exit {code}, {wall:.1f} s (spent {spent():.1f}/{BUDGET})")
    return code


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest(sys.argv[2])
    elif arg == "next":
        todo = [m[0] for m in MANIFESTS[sys.argv[2]] if m[0] not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    elif arg == "dev":
        sep = sys.argv.index("--")
        sys.exit(dev(sys.argv[2], sys.argv[3], sys.argv[sep + 1:]))
    else:
        sys.exit(run_unit(arg))
