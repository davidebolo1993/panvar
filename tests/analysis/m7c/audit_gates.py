#!/usr/bin/env python3
"""M7-C maintenance (added at M8-A): explicit numerical gates and perturbation controls for the attempt reconstruction.

M7-C's audit.py computed reconstruction errors but its exit condition did not enforce tolerances. This script recomputes
them from the archived M7-B.1 outputs (no genotyping; no historical package altered) and FAILS unless every run passes:
  G_mean      |mean of the 20 stage-2 attempt `lik` - res.json lik_mean|                 <= 1e-6
  G_var       |sample variance(attempt lik) * ln 10 - res.json lik_sd|                     <= 1e-4
  G_weighted  max over attempts |0.15 aln_lik + 1.85 depth_lik - lik|                      <= 1e-6
  G_attempts  every declared genotype has exactly attempts 1..20 at stage 2
Perturbation controls (each must make its gate FAIL, in every run):
  P1 +1e-3 on one attempt's lik            -> G_mean
  P2 res lik_sd + 1e-2                     -> G_var
  P3 population instead of sample variance -> G_var
  P4 +1e-3 on one attempt's aln_lik        -> G_weighted
  P5 depth weight 1.84 instead of 1.85     -> G_weighted
  P6 one attempt removed                   -> G_attempts

  audit_gates.py <out dir>
"""
import csv
import gzip
import json
import math
import os
import statistics
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from brotli_members import decode_members  # noqa: E402

out = sys.argv[1]
B1 = "/Users/davide.bolognini/panvar/tests/results/m7b1_fixed_pool_20261020/runs"
B0 = "/Users/davide.bolognini/panvar/tests/results/m7b_fixed_pool_20261020/runs"
LN10 = math.log(10)
TOL = {"mean": 1e-6, "var": 1e-4, "weighted": 1e-6}
canon = lambda g: ",".join(sorted(g.split(",")))
LOCUS = {"c4": "C4", "cyp": "CYP2D6"}


def gates(rows, opts, declared, wa=0.15, wd=1.85, var=statistics.variance):
    st2 = defaultdict(list)
    for r in rows:
        if r["stage"] == "2":
            st2[canon(r["genotype"])].append(r)
    att_ok = set(st2) == declared and all(sorted(int(r["attempt"]) for r in rs) == list(range(1, 21)) for rs in st2.values())
    lik = {g: [float(r["lik"]) for r in rs] for g, rs in st2.items()}
    e_mean = max(abs(statistics.mean(v) - opts[g]["lik_mean"]) for g, v in lik.items() if g in opts)
    e_var = max(abs(var(v) * LN10 - opts[g]["lik_sd"]) for g, v in lik.items() if g in opts and len(v) > 1)
    e_w = max(abs(wa * float(r["aln_lik"]) + wd * float(r["depth_lik"]) - float(r["lik"])) for rs in st2.values() for r in rs)
    return dict(G_attempts=att_ok, G_mean=e_mean <= TOL["mean"], G_var=e_var <= TOL["var"], G_weighted=e_w <= TOL["weighted"],
                e_mean=e_mean, e_var=e_var, e_weighted=e_w)


rows_out, ok = [], True
for u in sorted(x for x in os.listdir(B1) if x.startswith("loc.") and os.path.isdir(f"{B1}/{x}/gt")):
    cell, pool = u.split(".")[1], u.split(".")[2]
    g = f"{B1}/{u}/gt/loci/{LOCUS[cell.split('_')[0]]}"
    txt, members = decode_members(open(f"{g}/sol_ext..csv.br", "rb").read())
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    opts = {canon(o["genotype"]): o for o in json.load(gzip.open(f"{g}/res.json.gz", "rt"))["options"]}
    declared = {l.split("\t")[1] for l in open(f"{B0}/pools/{cell}.priors.{pool}.tsv")}
    base = gates(rows, opts, declared)
    i = next(k for k, r in enumerate(rows) if r["stage"] == "2")
    p1 = [dict(r) for r in rows]; p1[i]["lik"] = str(float(p1[i]["lik"]) + 1e-3)
    o2 = {k: dict(v) for k, v in opts.items()}; first = next(iter(o2)); o2[first]["lik_sd"] += 1e-2
    p4 = [dict(r) for r in rows]; p4[i]["aln_lik"] = str(float(p4[i]["aln_lik"]) + 1e-3)
    p6 = [r for k, r in enumerate(rows) if k != i]
    ctrl = {"P1_lik+1e-3": not gates(p1, opts, declared)["G_mean"],
            "P2_lik_sd+1e-2": not gates(rows, o2, declared)["G_var"],
            "P3_population_variance": not gates(rows, opts, declared, var=statistics.pvariance)["G_var"],
            "P4_aln+1e-3": not gates(p4, opts, declared)["G_weighted"],
            "P5_depth_weight_1.84": not gates(rows, opts, declared, wd=1.84)["G_weighted"],
            "P6_attempt_removed": not gates(p6, opts, declared)["G_attempts"]}
    run_ok = all(base[k] for k in ("G_attempts", "G_mean", "G_var", "G_weighted")) and all(ctrl.values())
    ok &= run_ok
    rows_out.append(dict(run=u, streams=len(members), **{k: base[k] for k in ("G_attempts", "G_mean", "G_var", "G_weighted")},
                         e_mean=f"{base['e_mean']:.2e}", e_var=f"{base['e_var']:.2e}", e_weighted=f"{base['e_weighted']:.2e}",
                         **{k: ("detected" if v else "NOT DETECTED") for k, v in ctrl.items()}, verdict="PASS" if run_ok else "FAIL"))
os.makedirs(out, exist_ok=True)
keys = list(rows_out[0])
with open(f"{out}/audit_gates.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows_out)
print(open(f"{out}/audit_gates.tsv").read())
sys.exit(0 if ok else 1)
