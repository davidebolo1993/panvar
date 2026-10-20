#!/usr/bin/env python3
"""M9-J.1 closure checks (DESIGN.md): Sol's reproducer against the fixed library, the 34-case replay compared with the saved M9-J
replay, and the saved M9-J smoke annotations re-assessed from their own saved tables. No genotyping.

  m9j1.py sol_repro <out> | replay_compare <out> | saved_smoke <out>
"""
import csv
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
M9J = "/Users/davide.bolognini/panvar/tests/results/m9j_compound_integration_20261020"
OUT = "/Users/davide.bolognini/panvar/tests/results/m9j1_closure_20261020"
REPRO = "/private/tmp/panvar-m9j-review.lm2ynx/repro.cpp"
mode, out = sys.argv[1], sys.argv[2]
NEW = ("decision_set_groups", "decision_set_missing_groups")


def T(p):
    return list(csv.DictReader((line for line in open(p) if not line.startswith("#")), delimiter="\t"))


def write(name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(rows)


if mode == "sol_repro":
    shutil.copy(REPRO, f"{out}/repro.cpp")   # the reviewer's file, kept verbatim with the unit
    subprocess.run(["c++", "-std=c++17", "-O1", f"-I{WT}/include", f"{out}/repro.cpp", f"{WT}/src/cli_utils.cpp", f"{WT}/build/libpanvarlib.a",
                    f"{WT}/build/libedlib.a", f"{WT}/external/abpoa/lib/libabpoa.a",   # Amendment 1: the test target's link line
                    f"{WT}/external/minigraph/libminigraph.a", "-lz", "-lpthread", "-o", f"{out}/repro"], check=True)
    lines = subprocess.run([f"{out}/repro"], capture_output=True, text=True, check=True).stdout.splitlines()
    print("\n".join(lines))
    want = ["INCOMPLETE_DECISION_SET", "INCOMPLETE_DECISION_SET"]
    got = [line.split(" / ")[1] for line in lines]
    ok = got == want and all("the_resolved_target" not in line for line in lines)
    print(f"sol_repro: {'PASS' if ok else 'FAIL'} (expected {want}, no the_resolved_target)")
    sys.exit(0 if ok else 1)

if mode == "replay_compare":
    # (1) the unchanged M9-J replay against the M9-I.1 Python reference, into this unit's directory
    code = subprocess.run([sys.executable, os.path.join(HERE, "..", "m9j", "m9j.py"), "replay", out]).returncode
    # (2) every column of the saved M9-J replay must be identical; the two new columns must show complete decision sets
    old = {r["case"]: r for r in T(f"{M9J}/runs/gate.replay.a2/cpp.tsv")}
    rows = []
    for r in T(f"{out}/cpp.tsv"):
        o = old[r["case"]]
        diff = [k for k in o if o[k] != r.get(k)]
        rows.append(dict(case=r["case"], identical_saved_columns=not diff, differing=";".join(diff) or "-", target_resolution=r["target_resolution"],
                         decision_set_groups=r["decision_set_groups"], decision_set_missing_groups=r["decision_set_missing_groups"],
                         ok=not diff and r["decision_set_missing_groups"] == "0" and int(r["decision_set_groups"]) > 0))
    write("saved_replay_comparison", rows)
    n = sum(r["ok"] for r in rows)
    print(f"replay_compare: m9j replay exit {code}; {n}/{len(rows)} cases identical to the saved M9-J replay with complete decision sets")
    sys.exit(0 if code == 0 and len(rows) == 34 == len(old) and n == 34 else 1)

if mode == "saved_smoke":
    # the saved M9-J caller annotations, re-assessed by the fixed helper from each run's own saved tables (manifests from smoke_check.a4)
    rows = []
    for name, run in (("acot", "smoke.acot.on"), ("acot_t1", "smoke.acot.on.t1"), ("cyp", "smoke.cyp.on.a3"), ("c4", "smoke.c4.on.a3")):
        saved = T(f"{M9J}/runs/{run}/call.sequence_compound.tsv")
        src = f"{M9J}/runs/smoke_check.a4/replay_{name.split('_')[0]}.manifest.tsv"
        man = f"{out}/replay_{name}.manifest.tsv"
        open(man, "w").write(open(src).read().replace(f"{M9J}/runs/smoke.acot.on/", f"{M9J}/runs/{run}/") if name == "acot_t1" else open(src).read())
        subprocess.run([f"{WT}/build/compound_annotation_test", "replay", man, f"{out}/replay_{name}.tsv"], check=True)
        rep = {r["case"]: r for r in T(f"{out}/replay_{name}.tsv")}
        if not saved:
            rows.append(dict(run=run, interval="(no rows)", saved_resolution="-", fixed_resolution="-", saved_mass_reference="-",
                             decision_set_groups="-", decision_set_missing_groups="-", ok=not rep))
        for s in saved:
            r = rep[s["interval"]]
            rows.append(dict(run=run, interval=s["interval"], saved_resolution=s["target_resolution"], fixed_resolution=r["target_resolution"],
                             saved_mass_reference=s["mass_reference"], decision_set_groups=r["decision_set_groups"],
                             decision_set_missing_groups=r["decision_set_missing_groups"],
                             ok=s["target_resolution"] == r["target_resolution"] and r["decision_set_missing_groups"] == "0"
                             and (s["mass_reference"] == "the_resolved_target") == (r["target_resolution"] == "RESOLVED")))
    write("saved_smoke", rows)
    print(open(f"{out}/saved_smoke.tsv").read())
    sys.exit(0 if rows and all(r["ok"] for r in rows) else 1)
raise SystemExit(f"unknown mode {mode}")
