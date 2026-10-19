#!/usr/bin/env python3
"""M5-A cross-score of the M4-B ts1 2x terminal incumbents (truth-free). Both terminals are scored under BOTH priors by
score-set (each run's own terminal declared with its saved total, which must reproduce; fragment and exposure are
spec-independent and must agree across the two scorings). Classification (frozen):
  under the DISTANCE prior, C_term (constant-arm terminal) vs D_term (distance-arm terminal):
    D_term >= C_term - 1e-6  -> CHANGED_OBJECTIVE: the distance trunk is not below the constant trunk under its own objective
    C_term >  D_term + 1e-6  -> DIFFERENT_BASIN: the distance trunk stopped below a point its own objective ranks higher
  (and the same comparison under the constant prior is reported). Never 'worse search' from cross-prior score comparisons.

  xscore_table.py <const scores_truth_free> <dist scores_truth_free> <out prefix>
"""
import csv
import sys

c = {r["label"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
d = {r["label"]: r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
out = sys.argv[3]
F = lambda r: float(r["fragment_term"])
E = lambda r: float(r["exposure_term"])
lines = ["label\tfragment\texposure\tlikelihood\tprior_constant\tprior_distance\ttotal_constant\ttotal_distance\tlikelihood_identical_across_specs"]
for k in ("C_term", "D_term"):
    lines.append("\t".join(map(str, [k, f"{F(c[k]):.6f}", f"{E(c[k]):.6f}", f"{F(c[k]) + E(c[k]):.6f}", f"{float(c[k]['prior']):.6f}",
                                     f"{float(d[k]['prior']):.6f}", f"{float(c[k]['total']):.6f}", f"{float(d[k]['total']):.6f}",
                                     abs(F(c[k]) - F(d[k])) < 1e-9 and abs(E(c[k]) - E(d[k])) < 1e-9])))
dc = float(c["C_term"]["total"]) - float(c["D_term"]["total"])
dd = float(d["C_term"]["total"]) - float(d["D_term"]["total"])
cls = "DIFFERENT_BASIN" if dd > 1e-6 else "CHANGED_OBJECTIVE"
reproduced = c["C_term"]["expected_total"] != "NA" and abs(float(c["C_term"]["total"]) - float(c["C_term"]["expected_total"])) <= 1e-6 and \
             abs(float(d["D_term"]["total"]) - float(d["D_term"]["expected_total"])) <= 1e-6
lines += ["", f"C_term - D_term under the constant prior\t{dc:.6f}", f"C_term - D_term under the distance prior\t{dd:.6f}",
          f"C_term - D_term likelihood (spec-independent)\t{(F(c['C_term']) + E(c['C_term'])) - (F(c['D_term']) + E(c['D_term'])):.6f}",
          f"saved terminal totals reproduced\t{reproduced}", f"CLASS\t{cls}"]
open(out + ".tsv", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
sys.exit(0 if reproduced else 1)
