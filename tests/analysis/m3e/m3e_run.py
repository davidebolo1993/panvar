#!/usr/bin/env python3
"""M3-E runner -- separate genotype-data support from phase support (experimental annotation): synthetic gate only (frozen manifest; one process at a time; <= 4
threads; <= 180 s per unit; 600 s budget counted as time spent in this batch plus the declared timeout of the next unit;
6 GiB envelope). No real-data execution.

  m3e_run.py manifest B | next B | <unit> | dev <name> <note> -- <command...>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m3e_support_semantics_20261019"
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
ORDER = ["lzo_5x", "mosaic_5x", "mosaic_2x"]   # controls first
V6 = "/Users/davide.bolognini/panvar/profiles/genotype_c4_v6_campaign_local.tsv"
JOINT_PROFILE = os.path.join(HERE, "profiles", "genotype_c4_v6_joint_experimental.tsv")
BASELINE_BIN = f"{OUT}/bin/panvar_baseline"
CYP = "/Users/davide.bolognini/panvar/tests/results/cyp2d6_shortread_pilot_v1/data/graphs"
ELIGIBLE = {"cyp2d6_lzo": (f"{CYP}/cyp2d6.gfa", "grch38#1#chr22:42031864-42245566", "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"),
            "cyp2d6_loo_NA19650": (f"{CYP}/loo.NA19650.gfa", "grch38#1#chr22:42031864-42245566", "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"),
            "acot": ("/Users/davide.bolognini/panvar/tests/results/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", "GRCh38#0#chr14:73416151-73668669", "none")}
TIMEOUT = {"lzo_5x": 120, "mosaic_5x": 150, "mosaic_2x": 120}
M3D = "/Users/davide.bolognini/panvar/tests/results/m3d_joint_search_20261019/runs"
SUPPORT_LINE = "block_support_rule\tfull_data_likelihood_v1_experimental\n"
MANIFEST_A = [("audit.saved", "audit", "-", 60)]
MANIFEST_G = [("gate.support", "gate_support", "-", 60)]
MANIFEST_R = ([(f"run.{c}.{arm}", "run", f"{c}:{arm}", TIMEOUT[c] + (30 if arm == "on" else 0)) for c in ORDER for arm in ("off", "on")] +
              [("compare", "cmp", "-", 30)])
# gate.support (manifest G) failed on gate fixtures, recorded; the fixtures were corrected (commits 0a9044b..) and the final
# gate units below run on committed code
MANIFEST_G2 = [("gate.support.final", "gate_support", "-", 60), ("gate.e2e.final", "gate_e2e", "-", 60)]
MANIFEST_A1 = [("compare.a1", "cmp", "-", 30)]   # Amendment 1: block table by column name; profile identity normalised
MANIFESTS = {"A": MANIFEST_A, "G": MANIFEST_G, "G2": MANIFEST_G2, "R": MANIFEST_R, "A1": MANIFEST_A1}
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
    if kind == "audit":
        runs = [f"{c}.{arm}={M3D}/run.{c}.{arm}/call" for c in ORDER for arm in ("off", "on")]
        return [sys.executable, os.path.join(HERE, "audit_saved.py"), f"{d}/audit", *runs], "saved-evidence audit (independent Python, truth-free)"
    if kind == "gate_e2e":
        return [f"{BUILD}/genotype_joint_search_e2e_test"], "caller gates incl. the annotation-only and thread checks"
    if kind == "gate_support":
        return [f"{BUILD}/genotype_block_support_test"], "support-rule gates (shared function)"
    if kind == "run":
        cell, arm = target.split(":")
        verify(cell)
        c = CELL[cell]
        profile = os.path.join(HERE, "profiles", "support_off.tsv" if arm == "off" else "support_on.tsv")
        cmd = ["nice", f"{BUILD}/panvar", "genotype", "--profile", profile, "--threads", "4", "--gfa", c["gfa"], "--reference-path", REF,
               "--reads1", f"{c['reads']}/reads_1.fq", "--reads2", f"{c['reads']}/reads_2.fq", "--library", f"{c['reads']}/library.txt",
               "--out-prefix", f"{d}/call"]
        if c["curation"] != "none":
            cmd += ["--panel-curation", c["curation"]]
        return cmd, f"caller with the experimental support annotation (joint search {arm})"
    if kind == "cmp":
        return [sys.executable, os.path.join(HERE, "compare_m3e.py"), f"{OUT}/runs", M3D, f"{OUT}/runs/audit.saved/audit.blocks.tsv", f"{d}/m3e", *ORDER], \
            "annotation-only check vs M3-D + cross-check vs the saved audit"
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
