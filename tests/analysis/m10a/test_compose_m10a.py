#!/usr/bin/env python3
"""M10-A gate: the composer's layer separation (pure Python; no genotyping)."""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compose_m10a  # noqa: E402

BASES = "/Users/davide.bolognini/panvar/tests/results/genotype_perf_bounded_20261015/inputs/profiles"
fails = 0


def check(cond, what):
    global fails
    print(("  PASS " if cond else "  FAIL ") + what)
    fails += 0 if cond else 1


def keys(p):
    return {l.split("\t", 1)[0]: l.split("\t", 1)[1] for l in open(p).read().split("\n") if "\t" in l and not l.startswith("#")}


with tempfile.TemporaryDirectory() as t:
    decl = os.path.join(t, "d.tsv")
    open(decl, "w").write("# compound_interval_declarations_v1\nmodel_fingerprint\tx\nblock_ids\ta\n")
    ov = ["time_limit_seconds=130:declared timeout minus 20 s"]
    got = {arm: keys(compose_m10a.compose(os.path.join(t, f"{arm}.tsv"), f"t-{arm}", arm, f"{BASES}/acot_5.tsv", decl, ov))
           for arm in compose_m10a.ARMS}
    d = got["D"]
    check(d["repair_junction_proposals"] == "experimental_v1" and d["union_fragment_blocks"] == "fine_core_4096_experimental"
          and d["compound_interval_declarations"] == decl and d["union_scoring_backend"] == "indexed_fine_experimental"
          and d["repair_max_neighbours"] == "1004" and d["time_limit_seconds"] == "130", "D = all layers, override applied")
    check("repair_junction_proposals" not in got["D0"] and set(d) - set(got["D0"]) == {"repair_junction_proposals"},
          "D0 = D without JP1.1 only")
    check(set(got["D0"]) - set(got["O"]) == {"early_incumbent", "union_fragment_blocks", "repair_round1_origin_admission"},
          "O = D0 without the retained route only")
    check("union_scoring_backend" not in got["REF"] and got["REF"]["output_schema"] == "sequence_first_v2"
          and got["REF"]["compound_annotation"] == "contract_v2", "REF: reference backends, same contract and annotation")
    check(all(got[a][k] == got["D"][k] for a in got for k in ("output_schema", "compound_annotation", "repair_joint_search",
                                                                 "max_state_classes", "repair_max_neighbours")),
          "contract, annotation and resource layers equal in every arm")
    bad = os.path.join(t, "bad.tsv")
    open(bad, "w").write("profile\tx\ntime_limit_seconds\t9\nswitch_probability\t0.1\n")
    try:
        compose_m10a.compose(os.path.join(t, "x.tsv"), "x", "D", bad, "none", [])
        check(False, "a scientific key in layer 3 is refused")
    except ValueError:
        check(True, "a scientific key in layer 3 is refused")
    dup = os.path.join(t, "dup.tsv")
    open(dup, "w").write("profile\tx\ntime_limit_seconds\t9\nrepair_max_neighbours\t5\n")
    try:
        compose_m10a.compose(os.path.join(t, "y.tsv"), "y", "D", dup, "none", ["output_schema=v1:no"])
        check(False, "a non-guard override is refused")
    except ValueError:
        check(True, "a non-guard override is refused")
    try:
        compose_m10a.compose(os.path.join(t, "z.tsv"), "z", "X", f"{BASES}/c4_5.tsv", "none", [])
        check(False, "an unknown arm is refused")
    except ValueError:
        check(True, "an unknown arm is refused")
print("compose_m10a: all passed" if not fails else "compose_m10a: FAILURES")
sys.exit(1 if fails else 0)
