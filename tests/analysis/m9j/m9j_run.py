#!/usr/bin/env python3
"""M9-J runner -- opt-in compound annotations: build, gates, replay, smoke off/on pairs (DESIGN.md; frozen manifest; one process
at a time; 1,800 s counted as the disclosed pre-freeze development plus time spent plus the declared timeout of the next unit plus,
for builds/gates/smoke calls, the smoke check).

  m9j_run.py manifest P | next P | <unit>
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
OUT = "/Users/davide.bolognini/panvar/tests/results/m9j_compound_integration_20261020"
M9C = "/Users/davide.bolognini/panvar/tests/results/m9c_acot_reliability_20261020"
R = "/Users/davide.bolognini/panvar/tests/results"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
BUDGET = 1800
PRE_FREEZE = 75.0   # disclosed: development compiles (48 + 5 + 2 s, including two failed attempts) and unit-test runs, rounded up
SMOKE_SRC = {"acot": ("/Users/davide.bolognini/panvar/tests/results/m9h_jp11_validation_20261020/runs/call.on.L4", 930),
             "cyp": ("/Users/davide.bolognini/panvar/tests/results/m1b_structural_20261019/runs/call.cyp_30", 310),
             "c4": ("/Users/davide.bolognini/panvar/tests/results/m1b_structural_20261019/runs/call.c4_10", 130)}
MANIFEST_P = ([("build", "build", "-", 300), ("gate.compound_unit", "gate_compound", "-", 30), ("gate.junction_unit", "gate_junction", "-", 30),
               ("gate.sequence_report", "gate_seq", "-", 30), ("gate.e2e_sequence_first", "gate_e2e", "-", 240), ("c4model", "c4model", "-", 150),
               ("declare", "declare", "-", 30), ("gate.replay", "replay", "-", 150),
               ("smoke.acot.off", "smoke", "acot:off", 930), ("smoke.acot.on", "smoke", "acot:on", 930), ("smoke.acot.on.t1", "smoke", "acot:on:t1", 930),
               ("smoke.cyp.off", "smoke", "cyp:off", 310), ("smoke.cyp.on", "smoke", "cyp:on", 310),
               ("smoke.c4.off", "smoke", "c4:off", 130), ("smoke.c4.on", "smoke", "c4:on", 130), ("smoke_check", "smoke_check", "-", 120)])
MANIFEST_A1 = [("build.a1", "build", "-", 300), ("gate.compound_unit.a1", "gate_compound", "-", 30), ("gate.replay.a1", "replay", "-", 150)]   # Amendment 1
MANIFESTS = {"P": MANIFEST_P, "A1": MANIFEST_A1}


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


GATES = ["build.a1", "gate.compound_unit.a1", "gate.junction_unit", "gate.sequence_report", "gate.e2e_sequence_first", "declare", "gate.replay.a1"]   # Amendment 1


def precondition(kind, target):
    if kind not in ("build",) and exit_of("build.a1" if "build.a1" in done() else "build") != "0":
        return "the build failed"
    if kind == "declare" and exit_of("c4model") != "0":
        return "the C4 model inventory failed"
    if kind == "smoke" and any(exit_of(g) != "0" for g in GATES):
        return "a gate failed: stop before caller smoke runs"
    if kind == "smoke_check" and not any(m[0] in done() for m in MANIFEST_P if m[1] == "smoke"):
        return "no smoke run"
    return None


def command(kind, target, d):
    py = lambda *a: [sys.executable, os.path.join(HERE, "m9j.py"), *a]
    if kind == "build":
        return ["cmake", "--build", BUILD, "--target", "panvar", "genotype_probe", "compound_annotation_test", "junction_proposals_test",
                "genotype_sequence_report_test", "genotype_command_e2e_test", "cyp2d6_block_alleles", "-j", "8"], None, "build every touched target"
    if kind == "gate_compound":
        return [f"{BUILD}/compound_annotation_test"], None, "compound annotation gates"
    if kind == "gate_junction":
        return [f"{BUILD}/junction_proposals_test"], None, "JP1.1 gates (unchanged)"
    if kind == "gate_seq":
        return [f"{BUILD}/genotype_sequence_report_test", d], None, "sequence report fixtures (publication/failure handling)"
    if kind == "gate_e2e":
        return [f"{BUILD}/genotype_command_e2e_test", d, WT], {"PANVAR_E2E_ONLY": "sequence_first"}, "caller-level sequence_first lifecycle"
    if kind in ("c4model", "declare", "replay", "smoke_check"):
        return py(kind, d), None, kind
    if kind == "smoke":
        parts = target.split(":")
        name, arm = parts[0], parts[1]
        src, _ = SMOKE_SRC[name]
        first = open(f"{src}/run.log").readline().split()
        cmd = first[first.index("-l") + 1:]
        cmd[cmd.index("--out-prefix") + 1] = f"{d}/call"
        if arm == "on":
            cmd += ["--profile", f"{OUT}/runs/declare/{name}.compound_profile.tsv"]
        if len(parts) > 2 and parts[2] == "t1":
            cmd[cmd.index("--threads") + 1] = "1"
        return cmd, None, f"smoke {target} (from {src})"
    raise ValueError(kind)


def run_unit(unit):
    part, (u, kind, target, timeout) = next((p, m) for p, rows in MANIFESTS.items() for m in rows if m[0] == unit)
    assert os.path.exists(f"{OUT}/MANIFEST_{part}.tsv") and unit not in done()
    rev, dirty = state()
    assert dirty == 0, "the tree must be committed and clean"
    why = precondition(kind, target)
    if why:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {why}")
        print(f"{unit}: NOT RUN ({why})")
        return 0
    used = spent()
    reserve = 120 if kind != "smoke_check" else 0   # the smoke check
    if used + timeout + reserve > BUDGET:
        record(unit, rev, dirty, timeout, 0.0, 0.0, "not_run", 0, f"NOT RUN: {used:.1f} + {timeout} + reserve {reserve} > {BUDGET} s")
        print(f"{unit}: NOT RUN (budget)")
        return 0
    d = f"{OUT}/runs/{unit}"
    assert not os.path.exists(d), f"{d} exists; a run is never repeated into the same place"
    os.makedirs(d)
    cmd, env, note = command(kind, target, d)
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
    record(unit, rev, dirty, timeout, wall, rss / 2**30, code, timed_out, note)
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
    else:
        sys.exit(run_unit(arg))
