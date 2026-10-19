#!/usr/bin/env python3
"""M3-D runner -- bounded default-off integration of likelihood-guided joint-block search: synthetic gate only (frozen manifest; one process at a time; <= 4
threads; <= 180 s per unit; 600 s budget counted as time spent in this batch plus the declared timeout of the next unit;
6 GiB envelope). No real-data execution.

  m3d_run.py manifest B | next B | <unit> | dev <name> <note> -- <command...>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m3d_joint_search_20261019"
BUDGET = 900
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
ORDER = ["lzo_5x", "mosaic_5x", "mosaic_2x"]   # controls first
V6 = "/Users/davide.bolognini/panvar/profiles/genotype_c4_v6_campaign_local.tsv"
JOINT_PROFILE = os.path.join(HERE, "profiles", "genotype_c4_v6_joint_experimental.tsv")
BASELINE_BIN = f"{OUT}/bin/panvar_baseline"
CYP = "/Users/davide.bolognini/panvar/tests/results/cyp2d6_shortread_pilot_v1/data/graphs"
ELIGIBLE = {"cyp2d6_lzo": (f"{CYP}/cyp2d6.gfa", "grch38#1#chr22:42031864-42245566", "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"),
            "cyp2d6_loo_NA19650": (f"{CYP}/loo.NA19650.gfa", "grch38#1#chr22:42031864-42245566", "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"),
            "acot": ("/Users/davide.bolognini/panvar/tests/results/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", "GRCh38#0#chr14:73416151-73668669", "none")}
TIMEOUT = {"lzo_5x": 120, "mosaic_5x": 150, "mosaic_2x": 120}
MANIFEST_G = [("gate.library", "gate_lib", "-", 60), ("gate.e2e", "gate_e2e", "-", 60), ("gate.regress", "gate_regress", "-", 120)]
MANIFEST_R = ([(f"run.{c}.{arm}", "run", f"{c}:{arm}", TIMEOUT[c] + (30 if arm == "on" else 0)) for c in ORDER for arm in ("baseline", "off", "on")] +
              [(f"eligible.{k}", "eligible", k, 60) for k in ELIGIBLE] + [("compare", "cmp", "-", 30)] +
              [(f"eval.{c}.{arm}", "eval", f"{c}:{arm}", 60) for c in ORDER for arm in ("off", "on")])
# Amendment 1: a repeat of the baseline binary on mosaic_2x shows whether the cache work counters vary run to run; compare.a1
# then treats exactly those seven keys as volatile (compare stays recorded with exit 1)
MANIFEST_A1 = [("run.mosaic_2x.baseline_repeat", "run", "mosaic_2x:baseline_repeat", 120), ("repeat.check", "repeat", "-", 30),
               ("compare.a1", "cmp", "-", 30)] + [(f"eval.a1.{c}.{arm}", "eval", f"{c}:{arm}", 60) for c in ORDER for arm in ("off", "on")]
MANIFESTS = {"G": MANIFEST_G, "R": MANIFEST_R, "A1": MANIFEST_A1}
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


def verify(cell):
    for key, want in INPUTS.items():
        c, path = key.split(":", 1)
        if c == cell:
            assert hashlib.sha256(open(path, "rb").read()).hexdigest() == want, f"{cell}: {path} differs from the frozen input"


def command(unit, kind, target, timeout, d):
    if kind == "gate_lib":
        return [f"{BUILD}/genotype_joint_search_test"], "library gates (joint search)"
    if kind == "gate_e2e":
        return [f"{BUILD}/genotype_joint_search_e2e_test"], "caller gates (joint search, synthetic valley)"
    if kind == "gate_regress":
        return ["/bin/sh", "-c", f"{BUILD}/genotype_mosaic_core_test && {BUILD}/genotype_block_delta_test && {BUILD}/genotype_scoring_context_test"], \
            "existing exactness regressions (mosaic core, block delta, scoring context)"
    if kind == "run":
        cell, arm = target.split(":")
        verify(cell)
        c = CELL[cell]
        binary = BASELINE_BIN if arm.startswith("baseline") else f"{BUILD}/panvar"
        profile = JOINT_PROFILE if arm == "on" else V6
        cmd = ["nice", binary, "genotype", "--profile", profile, "--threads", "4", "--gfa", c["gfa"], "--reference-path", REF,
               "--reads1", f"{c['reads']}/reads_1.fq", "--reads2", f"{c['reads']}/reads_2.fq", "--library", f"{c['reads']}/library.txt",
               "--out-prefix", f"{d}/call"]
        if c["curation"] != "none":
            cmd += ["--panel-curation", c["curation"]]
        return cmd, f"caller ({arm}; {'baseline binary' if arm.startswith('baseline') else 'M3-D binary'})"
    if kind == "eligible":
        gfa, ref, cur = ELIGIBLE[target]
        cmd = [f"{BUILD}/genotype_probe", "joint-bounds", "--out", f"{d}/bounds", "--model-gfa", gfa, "--reference", ref]
        if cur != "none":
            cmd += ["--panel-curation", cur]
        return cmd, "model-only joint-search eligibility"
    if kind == "repeat":
        a, b = f"{OUT}/runs/run.mosaic_2x.baseline/call", f"{OUT}/runs/run.mosaic_2x.baseline_repeat/call"
        return ["/bin/bash", "-c", f"diff <(sed 's#{a}#P#g' {a}.mosaic_status.tsv) <(sed 's#{b}#P#g' {b}.mosaic_status.tsv) | grep '^[<>]' | "
                f"grep -vE 'seconds_|threads_|resources_|footprint|rss|_bytes_end|memory|allocator' > {d}/baseline_vs_repeat.txt; cat {d}/baseline_vs_repeat.txt"], \
            "baseline vs baseline repeat: non-volatile status differences"
    if kind == "cmp":
        if unit == "compare.a1":
            os.environ["M3D_CACHE_COUNTERS_VOLATILE"] = "1"
        return [sys.executable, os.path.join(HERE, "compare_m3d.py"), f"{OUT}/runs", f"{OUT}/ledger.tsv", f"{d}/m3d", *ORDER], \
            "truth-free matched comparison + selection freeze"
    if kind == "eval":
        cell, arm = target.split(":")
        assert os.path.exists(f"{OUT}/runs/compare.a1/m3d.selection.sha256"), "evaluation needs the frozen selection (compare.a1)"
        c = CELL[cell]
        scheme = "lzo" if cell == "lzo_5x" else "mosaic"
        truth_gfa = c["gfa"]
        cmd = ["nice", f"{BUILD}/genotype_probe", "evaluate", "--scheme", scheme, "--donor", c["donor"], "--reference", REF,
               "--run-prefix", f"{OUT}/runs/run.{cell}.{arm}/call", "--model-gfa", c["gfa"], "--truth-gfa", truth_gfa, "--out", f"{d}/eval"]
        if c["curation"] != "none":
            cmd += ["--panel-curation", c["curation"]]
        return cmd, "truth evaluation (after the frozen selection)"
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
