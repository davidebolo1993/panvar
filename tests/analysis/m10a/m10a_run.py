#!/usr/bin/env python3
"""M10-A validation runner (DESIGN.md; frozen manifest; one process at a time). PREPARED ONLY in M10-A: the campaign is
executed in a later batch, after review. Budget: spent + declared timeout of the next unit <= BUDGET.

  m10a_run.py manifest P | next P | <unit>
"""
import csv
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
sys.path.insert(0, HERE)
import p1_run  # noqa: E402
import compose_m10a  # noqa: E402
import m10a  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = m10a.OUT
BUDGET = 2400
PRE_FREEZE = 0.0
GIB6 = str(6 << 30)
LOCI = {"c4": ("c4_5", "c4_30"), "cyp": ("cyp_5", "cyp_30"), "acot": ("acot_5", "acot_30")}
TIMEOUT = {   # declared per arm and cell; the run's time limit is the timeout minus 20 s (an operational guard override)
    "REF": {"c4_5": 150, "cyp_5": 150, "acot_5": 300},
    "O": {"c4_5": 90, "cyp_5": 90, "acot_5": 150, "c4_30": 150, "cyp_30": 150},
    "D0": {"c4_5": 90, "cyp_5": 90, "acot_5": 150, "c4_30": 150, "cyp_30": 150, "acot_30": 300},
    "D": {"c4_5": 90, "cyp_5": 90, "acot_5": 180, "c4_30": 180, "cyp_30": 180, "acot_30": 360},
}
CMP = {"backend": ("REF", "O"), "route": ("O", "D0"), "jp": ("D0", "D"), "prepared": ("D", "P"), "repeat": ("D", "D.r")}


def manifest():
    m = [("build", "build", "-", 300), ("gate.compose", "gate", "compose", 30), ("gate.compound_unit", "gate", "compound", 30),
         ("gate.junction_unit", "gate", "junction", 30), ("gate.sequence_report", "gate", "seq", 30),
         ("gate.e2e_sequence_first", "gate", "e2e:sequence_first", 240), ("gate.e2e_prepared_locus", "gate", "e2e:prepared_locus", 300),
         ("controls", "py", "controls", 120), ("declare", "py", "declare", 300)]
    for cell in ("c4_5", "cyp_5", "acot_5", "c4_30", "cyp_30", "acot_30"):
        for arm in ("REF", "O", "D0", "D"):
            if cell in TIMEOUT[arm]:
                m.append((f"{arm}.{cell}", "cli", f"{arm}:{cell}", TIMEOUT[arm][cell]))
        for mode in ("backend", "route", "jp"):
            x, y = CMP[mode]
            if cell in TIMEOUT[x] and cell in TIMEOUT[y]:
                m.append((f"cmp.{mode}.{cell}", "cmp", f"{mode}:{x}.{cell}:{y}.{cell}", 60))
    for locus, (lo, hi) in LOCI.items():
        m.append((f"write.{locus}", "write", locus, 120 if locus == "acot" else 60))
        m.append((f"loaded.{locus}", "loaded", locus, TIMEOUT["D"][lo] + TIMEOUT["D"][hi]))
        for cell in (lo, hi):
            m.append((f"cmp.prepared.{cell}", "cmp", f"prepared:D.{cell}:loaded.{locus}/{cell}", 60))
    m += [("loaded.acot.single30", "loaded1", "acot", TIMEOUT["D"]["acot_30"]),   # conditional (DESIGN.md 4.4)
          ("cmp.prepared.acot_30.single", "cmp", "prepared:D.acot_30:loaded.acot.single30/acot_30", 60),
          ("D.r.acot_30", "cli", "D:acot_30", TIMEOUT["D"]["acot_30"]),
          ("cmp.repeat.acot_30", "cmp", "repeat:D.acot_30:D.r.acot_30", 60), ("summary", "py", "summary", 60)]
    return m


MANIFESTS = {"P": manifest()}
STOPPING = ("cmp.backend.", "cmp.route.", "cmp.prepared.", "cmp.repeat.")


def ledger_rows():
    p = f"{OUT}/ledger.tsv"
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []


def spent():
    return PRE_FREEZE + sum(float(r["wall_seconds"]) for r in ledger_rows())


def done():
    return {r["unit"] for r in ledger_rows()}


def exit_of(u):
    return next((r["exit"] for r in ledger_rows() if r["unit"] == u), None)


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
    p = f"{OUT}/ledger.tsv"
    if not os.path.exists(p):
        open(p, "w").write("unit\trevision\tdirty_tree\tdeclared_timeout_seconds\twall_seconds\tpeak_rss_gib\texit\ttimed_out\tnote\n")
    with open(p, "a") as f:
        f.write(f"{unit}\t{rev}\t{dirty}\t{timeout}\t{wall:.2f}\t{rss:.3f}\t{code}\t{int(bool(timed_out))}\t{note}\n")


def run_status(d):
    return m10a.KV(f"{d}/call.mosaic_status.tsv").get("status", "missing")


def precondition(unit, kind, target):
    if kind != "build" and exit_of("build") != "0":
        return "the build failed"
    gates = [m[0] for m in MANIFESTS["P"] if m[1] == "gate"] + ["controls", "declare"]
    if kind in ("cli", "write", "loaded", "loaded1") and any(exit_of(g) != "0" for g in gates):
        return "a gate, the comparator controls or the declarations failed: stop before any caller run"
    if kind in ("cli", "write", "loaded", "loaded1") and any(r["unit"].startswith(STOPPING) and r["exit"] == "1" for r in ledger_rows()):
        return "an equivalence comparison failed: stop (an amendment is required before any further caller run)"
    if kind == "loaded" and exit_of(f"write.{target}") != "0":
        return "the index was not written"
    if kind == "loaded1":   # only if the two-sample process did not deliver a completed acot_30 sample (resource outcome)
        if exit_of("loaded.acot") is None:
            return "loaded.acot has not run"
        if run_status(f"{OUT}/runs/loaded.acot/acot_30") in p1_compare_completed():
            return "condition not met: loaded.acot delivered a completed acot_30 sample"
    if kind == "cmp":
        mode, a, b = target.split(":")
        for side in (a, b):
            if not os.path.exists(f"{OUT}/runs/{side}"):
                return f"{side} has no run directory"
    return None


def p1_compare_completed():
    return {"COMPLETE", "CONDITIONAL_GENOTYPE_EQUIVALENT", "BLOCK_CONSENSUS"}


def genotype_args(cell, profile, prefix):
    args = ["--preset", "short-read", "--profile", profile, "--threads", "4", "--graph", cell["graph"], "--reference-path", cell["reference"],
            "--reads1", cell["reads1"], "--reads2", cell["reads2"], "--library", cell["library"], "--dosage-mode", "audit", "--out", prefix, "--quiet"]
    if cell["curation"] != "none":
        args[args.index("--reads1"):args.index("--reads1")] = ["--panel-curation", cell["curation"]]
    return args


def profile_for(arm, cell_id, d, tag):
    cell = m10a.cells()[cell_id]
    locus = cell_id.split("_")[0]
    timeout = TIMEOUT["D" if arm == "P" else arm][cell_id]
    return compose_m10a.compose(f"{d}/profile.{tag}.tsv", f"m10a-{arm}-{cell_id}", "D" if arm == "P" else arm, cell["base_profile"],
                                f"{OUT}/runs/declare/{locus}.declarations.tsv",
                                [f"time_limit_seconds={timeout - 20}:declared timeout of the cell's D unit minus 20 s" if arm == "P" else
                                 f"time_limit_seconds={timeout - 20}:declared unit timeout minus 20 s (records TIME_LIMIT before the kill)"])


def command(unit, kind, target, d):
    py = lambda *a: [sys.executable, *a]  # noqa: E731
    if kind == "build":
        return ["cmake", "--build", BUILD, "--target", "panvar", "genotype_probe", "compound_annotation_test", "junction_proposals_test",
                "genotype_sequence_report_test", "genotype_command_e2e_test", "genotype_prepared_locus_driver", "cyp2d6_block_alleles",
                "-j", "8"], None
    if kind == "gate":
        if target == "compose":
            return py(os.path.join(HERE, "test_compose_m10a.py")), None
        if target.startswith("e2e:"):
            return [f"{BUILD}/genotype_command_e2e_test", d, WT], {"PANVAR_E2E_ONLY": target[4:]}
        exe = {"compound": "compound_annotation_test", "junction": "junction_proposals_test", "seq": "genotype_sequence_report_test"}[target]
        return [f"{BUILD}/{exe}"] + ([d] if target == "seq" else []), None
    if kind == "py":
        return py(os.path.join(HERE, "m10a.py"), target, d), None
    if kind == "cli":
        arm, cell_id = target.split(":")
        prof = profile_for(arm, cell_id, d, cell_id)
        return ["nice", f"{BUILD}/panvar", "genotype", *genotype_args(m10a.cells()[cell_id], prof, f"{d}/call")], None
    if kind == "cmp":
        mode, a, b = target.split(":")
        return py(os.path.join(HERE, "m10a_compare.py"), mode, f"{OUT}/runs/{a}", f"{OUT}/runs/{b}", f"{d}/compare.tsv"), None
    # persisted-index route (genotype_prepared_locus_driver): write-panel from the 5x args, then load and analyse in ONE process
    lo, hi = LOCI[target]
    samples = (lo,) if kind == "write" else (hi,) if kind == "loaded1" else (lo, hi)
    files = []
    for cell_id in samples:
        os.makedirs(f"{d}/{cell_id}")
        prof = profile_for("P", cell_id, f"{d}/{cell_id}", cell_id)
        files.append(f"{d}/{cell_id}/args.txt")
        open(files[-1], "w").write("".join(a + "\n" for a in genotype_args(m10a.cells()[cell_id], prof, f"{d}/{cell_id}/call")))
    index = f"{OUT}/index/{target}.panel.pvpl"
    driver = f"{BUILD}/genotype_prepared_locus_driver"
    if kind == "write":
        assert not os.path.exists(index), "an index is never overwritten"
        os.makedirs(os.path.dirname(index), exist_ok=True)
        return ["nice", driver, "write-panel", GIB6, "0", f"{d}/driver", index, *files], None
    return ["nice", driver, "loaded", GIB6, "1", f"{d}/driver", index, *files], None


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in MANIFESTS.items() for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and unit not in done()
    rev, dirty = state()
    assert dirty == 0, "the tree must be committed and clean"
    why = precondition(unit, kind, target)
    if why:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {why}")
        print(f"{unit}: NOT RUN ({why})")
        return 0
    used = spent()
    if used + timeout > BUDGET:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {used:.1f} + {timeout} > {BUDGET} s")
        print(f"{unit}: NOT RUN (budget)")
        return 0
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d), f"{d} exists; a run is never repeated into the same place"
    os.makedirs(d)
    cmd, env = command(unit, kind, target, d)
    cmd = ["/usr/bin/time", "-l", *cmd]
    if env:
        os.environ.update(env)
    t = time.time()
    with open(f"{d}/run.log", "w") as log:
        log.write(f"$ (declared timeout {timeout} s){' ' + str(env) if env else ''} " + " ".join(cmd) + "\n")
        log.flush()
        code, timed_out = p1_run.run_timed(cmd, log, timeout)
    for k in (env or {}):
        os.environ.pop(k, None)
    wall = time.time() - t
    text = open(f"{d}/run.log").read()
    rss = next((int(line.split()[0]) for line in text.splitlines() if "maximum resident set size" in line), 0)
    record(unit, rev, dirty, timeout, wall, rss / 2**30, code, timed_out, f"{kind} {target}")
    print("\n".join(line for line in text.splitlines() if not line.startswith(" ") and not line.startswith("$"))[-2500:])
    print(f"{unit}: exit {code}{' TIMEOUT' if timed_out else ''}, {wall:.1f} s (spent {spent():.1f}/{BUDGET})")
    return 0


if __name__ == "__main__":
    arg = sys.argv[1]
    if arg == "manifest":
        write_manifest(sys.argv[2])
    elif arg == "next":
        todo = [m[0] for m in MANIFESTS[sys.argv[2]] if m[0] not in done()]
        sys.exit(run_unit(todo[0]) if todo else print("all units recorded") or 0)
    elif arg == "print":
        print("".join("\t".join(str(x) for x in m) + "\n" for m in MANIFESTS["P"]), end="")
    else:
        sys.exit(run_unit(arg))
