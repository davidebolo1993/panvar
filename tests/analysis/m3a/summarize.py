#!/usr/bin/env python3
"""M3-A summary (frozen before execution). Reads the unit outputs and writes:
  <prefix>.sensitivity.tsv   per cell x prior model: winner, call rank, margins, intact-vs-mosaic, ambiguity (truth-free)
  <prefix>.answers.tsv       per cell: sensitive to prior strength / boundary placement / zero redraw (frozen rule below)
  <prefix>.flank.tsv         prior-side flank association (pool contrasts) + fragment-side classes from EXISTING M1-B
                             decompositions (descriptive; M1-C later reconciled the mass reduction they flagged)
  <prefix>.evaluator_only.tsv  counterfactual rows (truth or truth-informed selections), never mixed into the above
  <prefix>.representation.tsv  the representation-test rows
Frozen rule: a cell is SENSITIVE to a factor when the pool winner (argmax of saved fragment + exposure + recomputed prior)
differs from the const_0.01 winner under that factor's models: strength = const_{0.001,0.003,0.03,0.1}; placement =
dist_matched (equal total hazard); zero redraw = zero (diagnostic only). Margins are reported either way.

  summarize.py <runs dir> <out prefix>
"""
import csv
import os
import sys

runs, prefix = sys.argv[1:3]
CELLS = ["c4_control", "c4_mosaic_2x", "c4_sv7_10", "cyp_label", "cyp_loo_sv9", "acot_5", "acot_30"]
M1B = "/Users/davide.bolognini/panvar/tests/results/m1b_structural_20261019/runs"
DECOMP = {"c4_sv7_10": f"{M1B}/score.c4_10/score.decomposition_summary.tsv", "cyp_loo_sv9": f"{M1B}/score.cyp_5/score.decomposition_summary.tsv"}


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.exists(path) else []


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path):
            k, _, v = line.rstrip("\n").partition("\t")
            out[k] = v
    return out


sens, answers, flank, evaluator = [], [], [], []
for cell in CELLS:
    d = f"{runs}/pool.{cell}"
    s = {r["model"]: r for r in rows(f"{d}/pool.summary.tsv")}
    gate = kv(f"{d}/pool.gate.tsv")
    if not s:
        answers.append({"cell": cell, "status": "NO OUTPUT"})
        continue
    for m, r in s.items():
        sens.append({"cell": cell, **r})
    base = s["const_0.01"]["winner_seq_id"]
    strength = [m for m in ("const_0.001", "const_0.003", "const_0.03", "const_0.1") if s[m]["winner_seq_id"] != base]
    placement = s["dist_matched"]["winner_seq_id"] != base
    zero = s["zero"]["winner_seq_id"] != base
    answers.append({
        "cell": cell, "status": "OK", "pool": gate.get(f"{cell}_pool_candidates"), "unique_haplotypes": gate.get(f"{cell}_unique_haplotypes"),
        "saved_prior_reproduced_max_abs": gate.get(f"{cell}_current_prior_max_abs_diff"),
        "current_winner_is_saved_call": s["const_0.01"]["winner_is_call"],
        "current_margin_to_second": s["const_0.01"]["gap_to_second"], "current_winner_probability": s["const_0.01"]["winner_probability"],
        "sensitive_strength": "yes: " + ",".join(strength) if strength else "no",
        "sensitive_placement": "yes" if placement else "no",
        "placement_winner_minus_call": s["dist_matched"]["winner_minus_call_total"],
        "sensitive_zero_redraw": "yes" if zero else "no",
        "zero_finite": s["zero"]["finite"], "zero_winner_intact": s["zero"]["winner_intact"],
        "best_intact_minus_best_mosaic(const_0.01)": s["const_0.01"]["best_intact_minus_best_mosaic"],
        "best_intact_minus_best_mosaic(dist_matched)": s["dist_matched"]["best_intact_minus_best_mosaic"],
        "within_1nat(const_0.001..0.1)": "/".join(s[m]["within_1nat"] for m in ("const_0.001", "const_0.003", "const_0.01", "const_0.03", "const_0.1")),
    })
    for r in rows(f"{d}/pool.flank_prior.tsv"):
        if r["model"] not in ("const_0.01", "dist_matched", "const_0.001", "const_0.1"):
            continue
        flank.append({"cell": cell, "source": "prior (pool contrast)", **r})
    for r in rows(f"{d}/pool.counterfactual_evaluator_only.tsv"):
        evaluator.append(r)
    if cell in DECOMP and os.path.exists(DECOMP[cell]):
        k = kv(DECOMP[cell])
        flank.append({"cell": cell, "source": "fragment (existing M1-B decomposition, descriptive)", "contrast": f"{k['label_A']} vs {k['label_B']}",
                      "model": "-", "target_blocks": k["block"],
                      "delta_prior_full_path": "fragment_delta=" + k["fragment_delta"],
                      "delta_prior_target_marginal": f"interior={k['class_block_interior_on_both_delta']} ({k['class_block_interior_on_both_fragments']} frags)",
                      "delta_flank_association": f"junction-spanning={k['class_block_junction_on_both_delta']} ({k['class_block_junction_on_both_fragments']} frags); "
                                                 f"outside={k['class_outside_block_on_all_delta']} ({k['class_outside_block_on_all_fragments']} frags); "
                                                 f"exposure A-B={float(k['exposure_term_A']) - float(k['exposure_term_B']):.3f}"})


def write(name, items):
    if not items:
        open(f"{prefix}.{name}.tsv", "w").write("")
        return
    keys = []
    for it in items:
        for k in it:
            if k not in keys:
                keys.append(k)
    with open(f"{prefix}.{name}.tsv", "w") as f:
        f.write("\t".join(keys) + "\n")
        for it in items:
            f.write("\t".join(str(it.get(k, "")) for k in keys) + "\n")


write("sensitivity", sens)
write("answers", answers)
write("flank", flank)
write("evaluator_only", evaluator)
rep = [line for line in open(f"{runs}/repr/run.log") if line.count("\t") >= 2 or line.startswith("m3a_prior_representation")]
open(f"{prefix}.representation.tsv", "w").write("".join(l for l in rep if not l.startswith("$")))
for a in answers:
    print("\t".join(f"{k}={v}" for k, v in a.items()))
