#!/usr/bin/env python3
"""M3-E.1 runner -- numerical certification and mixed-projection corrections (experimental support rule v2): synthetic gate only (frozen manifest; one process at a time; <= 4
threads; <= 180 s per unit; 600 s budget counted as time spent in this batch plus the declared timeout of the next unit;
6 GiB envelope). No real-data execution.

  m3e1_run.py manifest B | next B | <unit> | dev <name> <note> -- <command...>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m3e1_certification_20261019"
BUDGET = 300
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
M3E = "/Users/davide.bolognini/panvar/tests/results/m3e_support_semantics_20261019/runs"
REPRO_BIN = f"{OUT}/bin/repro_v1"
MANIFEST_P = [("repro.v1", "repro", "-", 30)]
RUNS6 = [f"{c}.{arm}" for c in ("lzo_5x", "mosaic_5x", "mosaic_2x") for arm in ("off", "on")]
MANIFEST_V = [("gate.support.v2", "gate_support", "-", 60), ("gate.e2e.v2", "gate_e2e", "-", 60), ("repro.v1.after", "repro_after", "-", 30),
              ("audit.v2", "audit_v2", "-", 60), ("audit.cpp_v2", "audit_cpp", "-", 60), ("crosscheck.v2", "cross", "-", 30)]
MANIFESTS = {"P": MANIFEST_P, "V": MANIFEST_V}
CARRIED = 0.0
INPUTS = {}


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
    if kind == "gate_support":
        return [f"{BUILD}/genotype_block_support_test"], "support gates v1 (unchanged) + v2"
    if kind == "gate_e2e":
        return [f"{BUILD}/genotype_joint_search_e2e_test"], "focused caller checks (v1 and v2 annotation-only, threads, absent never PASS)"
    if kind == "repro_after":
        return [f"{OUT}/bin/repro_v1_after"], "Sol's reproducer rebuilt against the M3-E.1 library: v1 must behave as before"
    if kind == "audit_v2":
        runs = [f"{r}={M3D}/run.{r}/call" for r in RUNS6]
        return [sys.executable, os.path.join(HERE, "audit_saved_v2.py"), f"{d}/audit_v2", f"{M3E}/audit.saved/audit.blocks.tsv", *runs], \
            "independent Python v2 audit of the saved M3-D tables"
    if kind == "audit_cpp":
        cmds = " && ".join(f"{BUILD}/genotype_probe support-audit --run {M3D}/run.{r}/call --out {d}/{r} > /dev/null" for r in RUNS6)
        return ["/bin/bash", "-c", cmds], "C++ block_data_support_v2 on the same saved tables"
    if kind == "cross":
        return [sys.executable, os.path.join(HERE, "cross_v2.py"), f"{OUT}/runs/audit.v2/audit_v2.blocks.tsv", f"{OUT}/runs/audit.cpp_v2",
                f"{d}/cross.tsv"], "C++ v2 vs Python v2"
    if kind == "repro":
        return [REPRO_BIN], "Sol's reproducer (unchanged source) against the current shared function (v1)"
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
