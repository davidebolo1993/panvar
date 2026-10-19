#!/usr/bin/env python3
"""M4-B runner -- default-off distance-prior integration: gates, hazard-matched lambda freeze, six matched fresh-search
comparisons (controls first, stop rule), evaluation after the frozen predictions (DESIGN.md; one process at a time; 4
threads; 1,200 s counted as time spent plus the declared timeout of the next unit; 6 GiB execution envelope).

  m4b_run.py manifest G|L|C|E|D | next G|L|C|E|D | <unit> | dev <name> <note> -- <command...>
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
import compose_profile  # noqa: E402

WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
BUILD = f"{WT}/build"
OUT = "/Users/davide.bolognini/panvar/tests/results/m4b_distance_prior_20261019"
BASE_BIN = f"{OUT}/bin/panvar_baseline"   # a728eb71 (the M3/M4-A caller), copied before any M4-B library change
BUDGET = 1200
C4 = "/Users/davide.bolognini/panvar/tests/results/genotype_c4_production_11fe681"
V6 = "/Users/davide.bolognini/panvar/profiles/genotype_c4_v6_campaign_local.tsv"
REF_C4 = "HG00096#1#haplotype1-0000024:31848049-32074186"
MCUR = f"{C4}/mosaic_inputs/mosaic_exclude.tsv"
PERF = "/Users/davide.bolognini/panvar/tests/results/genotype_perf_bounded_20261015"
CYP_CUR = "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"


def mosaic(plan, cov):
    m = f"{C4}/mosaic_inputs/{plan}"
    return dict(family="c4", gfa=f"{m}/staged.gfa", ref=REF_C4, curation=MCUR, reads=f"{m}/cov{cov}", scheme="mosaic", donor="MOSAIC",
                truth=f"{m}/staged.gfa", profile=V6, dosage=None)


PERF_CELLS = {r["id"]: r for r in csv.DictReader(open(f"{PERF}/cells.tsv"), delimiter="\t")}
CELL = {
    "c4_lzo_hg00171_5": dict(family="c4", gfa=f"{C4}/graphs/c4.gfa", ref=REF_C4, curation=None, reads=f"{C4}/reads/HG00171/cov5", scheme="lzo",
                             donor="HG00171", truth=f"{C4}/graphs/c4.gfa", profile=V6, dosage=None),
    "c4_panel_control_1_2": mosaic("panel_control_1", 2),
    "cyp_lzo_hg02769_5": dict(family="p1", gfa=PERF_CELLS["cyp_5"]["graph"], ref=PERF_CELLS["cyp_5"]["reference"], curation=CYP_CUR,
                              reads=os.path.dirname(PERF_CELLS["cyp_5"]["reads1"]), scheme="lzo", donor="HG02769", truth=PERF_CELLS["cyp_5"]["graph"],
                              profile=PERF_CELLS["cyp_5"]["base_profile"], dosage="repair"),
    "acot_hg00097_5": dict(family="p1", gfa=PERF_CELLS["acot_5"]["graph"], ref=PERF_CELLS["acot_5"]["reference"], curation=None,
                           reads=os.path.dirname(PERF_CELLS["acot_5"]["reads1"]), scheme="lzo", donor="HG00097", truth=PERF_CELLS["acot_5"]["graph"],
                           profile=PERF_CELLS["acot_5"]["base_profile"], dosage="repair"),
    "c4_both_homologues_1_5": mosaic("both_homologues_1", 5),
    "c4_three_switches_1_2": mosaic("three_switches_1", 2),
}
CONTROLS = ["c4_lzo_hg00171_5", "c4_panel_control_1_2", "cyp_lzo_hg02769_5", "acot_hg00097_5", "c4_both_homologues_1_5"]
ERRORS = ["c4_three_switches_1_2"]
TIMEOUT = {"c4_lzo_hg00171_5": 90, "c4_panel_control_1_2": 90, "cyp_lzo_hg02769_5": 120, "acot_hg00097_5": 180,
           "c4_both_homologues_1_5": 120, "c4_three_switches_1_2": 90}
REGRESS = ("^panvar_genotype_(mosaic_core|joint_search|block_support|joint_search_e2e|locus_adapter|cli_contract|curation_contract|dosage|"
           "headline|sequence_report|prepared_locus_file|v6_profile|mosaic_proposal|grouped_beam|review_regressions|source_tree_clean)$")

MANIFEST_G = [("gate.lib", "exe", "genotype_copying_transitions_test", 60), ("gate.e2e", "exe", "genotype_copying_transitions_e2e_test", 240),
              ("gate.m4a_prototype", "exe", "m4a_distance_prior", 30), ("gate.regress", "ctest", REGRESS, 420),
              ("gate.cli.switch_conflict", "refuse", "switch_flag", 30), ("gate.cli.profile_conflict", "refuse", "switch_profile", 30),
              ("gate.cli.length_without_mode", "refuse", "length_only", 30)]
MANIFEST_L = [(f"geometry.{c}", "geometry", c, 60) for c in CELL] + [("lambda.freeze", "freeze", "-", 30)]
MANIFEST_C = ([("identity.c4_lzo_hg00171_5.baseline", "run", "c4_lzo_hg00171_5:baseline", 90)] +
              [(f"run.{c}.{arm}", "run", f"{c}:{arm}", TIMEOUT[c]) for c in CONTROLS for arm in ("const", "dist")] +
              [("identity.check.c4_lzo_hg00171_5", "identity", "c4_lzo_hg00171_5", 30),
               ("agreement.c4_lzo_hg00171_5", "agree", "c4_lzo_hg00171_5", 60),
               ("incompatible.c4_lzo_hg00171_5", "incompat", "c4_lzo_hg00171_5", 60),
               ("control.check", "cmp", "control", 30)] +
              [(f"eval.{c}.{arm}", "eval", f"{c}:{arm}", 60) for c in CONTROLS for arm in ("const", "dist")] +
              [("control.gate", "acc", "control", 30)])
MANIFEST_E = ([("identity.c4_three_switches_1_2.baseline", "run", "c4_three_switches_1_2:baseline", 90)] +
              [(f"run.{c}.{arm}", "run", f"{c}:{arm}", TIMEOUT[c]) for c in ERRORS for arm in ("const", "dist")] +
              [("identity.check.c4_three_switches_1_2", "identity", "c4_three_switches_1_2", 30), ("compare", "cmp", "all", 30)] +
              [(f"eval.{c}.{arm}", "eval", f"{c}:{arm}", 60) for c in ERRORS for arm in ("const", "dist")] + [("accuracy", "acc", "all", 30)])
MANIFEST_D = [(f"decomp.{c}", "decomp", c, 90) for c in CONTROLS + ERRORS]
# Amendment 1: gate.lib failed one row on a test-assertion defect (G4 carved case); corrected assertion, re-run
MANIFEST_A1 = [("gate.lib.a1", "exe", "genotype_copying_transitions_test", 60)]
MANIFESTS = {"A1": MANIFEST_A1, "G": MANIFEST_G, "L": MANIFEST_L, "C": MANIFEST_C, "E": MANIFEST_E, "D": MANIFEST_D}
INPUTS = list(csv.DictReader(open(os.path.join(HERE, "INPUTS.tsv")), delimiter="\t"))
LAMBDA = os.path.join(HERE, "LAMBDA.tsv")


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


def lam(cell):
    rows = {r["cell"]: r for r in csv.DictReader(open(LAMBDA), delimiter="\t")}
    return rows[cell]


def profile_for(cell, arm, d):
    """const: the cell's validated settings; dist: the same minus any scalar switch setting plus the frozen distance keys."""
    c = CELL[cell]
    if c["family"] == "c4":
        base = c["profile"]
    else:   # the P1 'opt' composition of the cell's base profile, its own guards unchanged (no operational overrides)
        base = f"{d}/profile.const.tsv"
        compose_profile.compose(base, f"m4b-const-{cell}", c["profile"], "full", [])
    if arm in ("const", "baseline"):
        return base
    lines = [l for l in open(base).read().split("\n") if not l.startswith("switch_probability\t")]
    lines = [l if not l.startswith("profile\t") else f"profile\tm4b-dist-{cell}" for l in lines]
    row = lam(cell)
    out = f"{d}/profile.dist.tsv"
    open(out, "w").write("\n".join(lines).rstrip("\n") + f"\n# M4-B: hazard-matched control length frozen in LAMBDA.tsv (geometry {row['geometry_digest']},"
                         f" model {row['model_fingerprint']})\ncopying_transition\tdistance_ref_core_v1_experimental\ncopying_length_bp\t{row['lambda_bp']}\n")
    return out


def caller(cell, arm, d, binary):
    c = CELL[cell]
    prof = profile_for(cell, arm, d)
    cmd = ["nice", binary, "genotype", "--profile", prof, "--threads", "4", "--gfa", c["gfa"], "--reference-path", c["ref"],
           "--reads1", f"{c['reads']}/reads_1.fq", "--reads2", f"{c['reads']}/reads_2.fq", "--library", f"{c['reads']}/library.txt",
           "--out-prefix", f"{d}/call", "--quiet"]
    if c["curation"]:
        cmd += ["--panel-curation", c["curation"]]
    if c["dosage"]:
        cmd += ["--dosage-mode", c["dosage"]]
    return cmd


def run_dir(cell, arm):
    return f"{OUT}/runs/run.{cell}.{arm}" if arm != "baseline" else f"{OUT}/runs/identity.{cell}.baseline"


def kv(p):
    return {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l} if os.path.exists(p) else {}


def command(kind, target, d):
    verify()
    if kind == "exe":
        return [f"{BUILD}/{target}"], f"gate executable {target}"
    if kind == "ctest":
        return ["nice", "ctest", "--test-dir", BUILD, "-R", target, "--output-on-failure", "-j1"], "existing regressions (constant mode)"
    if kind == "refuse":
        c = CELL["c4_lzo_hg00171_5"]
        base = [f"{BUILD}/panvar", "genotype", "--threads", "4", "--gfa", c["gfa"], "--reference-path", c["ref"], "--reads1", f"{c['reads']}/reads_1.fq",
                "--reads2", f"{c['reads']}/reads_2.fq", "--library", f"{c['reads']}/library.txt", "--out-prefix", f"{d}/call", "--quiet"]
        if target == "switch_flag":
            open(f"{d}/p.tsv", "w").write("profile\tm4b-refuse\ncopying_transition\tdistance_ref_core_v1_experimental\ncopying_length_bp\t1000\n")
            return ["/bin/sh", "-c", " ".join([*base, "--profile", f"{d}/p.tsv", "--switch-probability", "0.01"]) +
                    " 2>&1 | tee " + f"{d}/out.txt" + "; grep -q 'conflicts with switch_probability' " + f"{d}/out.txt"], "expect refusal: distance + --switch-probability"
        if target == "switch_profile":
            open(f"{d}/p.tsv", "w").write("profile\tm4b-refuse\nswitch_probability\t0.01\ncopying_transition\tdistance_ref_core_v1_experimental\ncopying_length_bp\t1000\n")
            return ["/bin/sh", "-c", " ".join([*base, "--profile", f"{d}/p.tsv"]) + " 2>&1 | tee " + f"{d}/out.txt" +
                    "; grep -q 'conflicts with switch_probability' " + f"{d}/out.txt"], "expect refusal: distance + profile switch_probability"
        open(f"{d}/p.tsv", "w").write("profile\tm4b-refuse\ncopying_length_bp\t1000\n")
        return ["/bin/sh", "-c", " ".join([*base, "--profile", f"{d}/p.tsv"]) + " 2>&1 | tee " + f"{d}/out.txt" +
                "; grep -q 'refused with constant_v1' " + f"{d}/out.txt"], "expect refusal: copying_length_bp without the distance mode"
    if kind == "geometry":
        c = CELL[target]
        cmd = ["nice", f"{BUILD}/genotype_probe", "prior-geometry", "--out", f"{d}/geom", "--model-gfa", c["gfa"], "--reference", c["ref"], "--constant-r", "0.01"]
        return cmd + (["--panel-curation", c["curation"]] if c["curation"] else []), "adapter geometry + hazard-matched control length"
    if kind == "freeze":
        return [sys.executable, os.path.join(HERE, "freeze_lambda.py"), f"{OUT}/runs", f"{d}/LAMBDA.tsv", *CELL], "freeze lambda_HM per cell"
    if kind == "run":
        cell, arm = target.split(":")
        if cell in ERRORS:
            assert exit_of("control.gate") == "0", "the error cell runs only after the control gate passes (stop rule)"
        if arm == "dist":
            assert os.path.exists(LAMBDA), "the distance arm needs the frozen LAMBDA.tsv"
        return caller(cell, arm, d, BASE_BIN if arm == "baseline" else f"{BUILD}/panvar"), f"caller ({arm})"
    if kind == "identity":
        return [sys.executable, os.path.join(HERE, "identity.py"), run_dir(target, "baseline") + "/call", run_dir(target, "const") + "/call",
                f"{d}/identity.tsv"], "constant-mode identity: new binary vs baseline binary"
    if kind in ("agree", "incompat"):
        c = CELL[target]
        return [sys.executable, os.path.join(HERE, "agree.py"), kind, f"{BUILD}/genotype_probe", run_dir(target, "const") + "/call",
                run_dir(target, "dist") + "/call", c["gfa"], c["ref"], c["curation"] or "none", c["reads"], lam(target)["lambda_bp"], d], \
            "caller vs score-set priors" if kind == "agree" else "incompatible specification refused"
    if kind == "cmp":
        cells = CONTROLS if target == "control" else CONTROLS + ERRORS
        return [sys.executable, os.path.join(HERE, "compare_m4b.py"), f"{OUT}/runs", f"{OUT}/ledger.tsv", LAMBDA, f"{d}/m4b", *cells], \
            "truth-free matched comparison + prediction freeze"
    if kind == "eval":
        cell, arm = target.split(":")
        freeze = "control.check" if cell in CONTROLS else "compare"
        assert os.path.exists(f"{OUT}/runs/{freeze}/m4b.selection.sha256"), "evaluation needs the frozen predictions"
        if not os.path.exists(f"{run_dir(cell, arm)}/call.mosaic_status.tsv"):
            return None, "NOT RUN (condition): the caller wrote no status"
        c = CELL[cell]
        cmd = ["nice", f"{BUILD}/genotype_probe", "evaluate", "--scheme", c["scheme"], "--donor", c["donor"], "--reference", c["ref"],
               "--run-prefix", f"{run_dir(cell, arm)}/call", "--model-gfa", c["gfa"], "--truth-gfa", c["truth"], "--out", f"{d}/eval"]
        return cmd + (["--panel-curation", c["curation"]] if c["curation"] else []), "truth evaluation (after the frozen predictions)"
    if kind == "acc":
        cells = CONTROLS if target == "control" else CONTROLS + ERRORS
        return [sys.executable, os.path.join(HERE, "accuracy_m4b.py"), f"{OUT}/runs", target, f"{d}/m4b", *cells], \
            "accuracy per arm (truth); control stop rule" if target == "control" else "accuracy per arm (truth)"
    if kind == "decomp":
        cm = kv(f"{OUT}/runs/compare/m4b.matched.tsv")
        c = CELL[target]
        return [sys.executable, os.path.join(HERE, "decomp.py"), f"{BUILD}/genotype_probe", target, run_dir(target, "const") + "/call",
                run_dir(target, "dist") + "/call", c["gfa"], c["ref"], c["curation"] or "none", c["reads"], lam(target)["lambda_bp"], d], \
            "changed-call decomposition: likelihood, prior, pool (runs only for a changed call)"
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
    cmd, note = command(kind, target, d)
    if cmd is None:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, note)
        print(f"{unit}: {note}")
        return 3
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
    print("\n".join(l for l in text.splitlines() if not l.startswith(" ") and not l.startswith("$"))[-3000:])
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
