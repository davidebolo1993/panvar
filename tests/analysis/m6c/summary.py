#!/usr/bin/env python3
"""M6-C summary (synthetic; truth is the declared fixture truth).

Per run: the structural posterior over the declared catalogue under a UNIFORM structural weighting (no structural or
copying prior), from the nuisance-marginalised likelihood; the same from fixed-sequence scoring (base content only); the
set of structures within 1e-6 of the best (ties); and the best structure's decomposition (fragment / exposure / log prior
of its best nuisance state, flips, effective states).
Frozen classification (marginal and fixed separately):
  CORRECT        truth's posterior >= 0.95
  UNCERTAIN      truth is in the tied top set together with other structures (no structure >= 0.95)
  WRONG          another structure has posterior >= 0.95
  SPLIT          otherwise (report the posterior of truth and of the top structure)
Attribution per fixture (frozen reading):
  CONTENT_ACCOMMODATION  fixed WRONG/SPLIT and marginal CORRECT in a fixture with missing variants (F2*)
  ALLOCATION_EVIDENCE    marginal CORRECT in an adjacency fixture (F4*) with truth separated from same-total structures
  PRIOR_DRIVEN           the marginal winner and the runner-up have equal best-state likelihood terms (|dF|,|dE| < 1e-6)
                         and differ only through the nuisance prior
  UNRESOLVED / MISLEADING UNCERTAIN where structure is unobservable (expected for F3) / WRONG with posterior >= 0.95

  summary.py <runs dir> <ledger> <out prefix>
"""
import csv
import glob
import math
import os
import sys

runs, ledger, prefix = sys.argv[1:4]
T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l}
led = {r["unit"]: r for r in T(ledger)}


def classify(post, truth, logs):
    # Amendment 3: ties are decided on the log marginals (a posterior can underflow to 0; log(0) failed)
    top = max(logs, key=logs.get)
    tied = [k for k in logs if abs(logs[k] - logs[top]) <= 1e-6]
    if post.get(truth, 0) >= 0.95:
        return "CORRECT", tied
    if truth in tied and len(tied) > 1:
        return "UNCERTAIN", tied
    if post[top] >= 0.95:
        return "WRONG", tied
    return "SPLIT", tied


def posterior(vals):
    m = max(vals.values())
    z = sum(math.exp(v - m) for v in vals.values())
    return {k: math.exp(v - m) / z for k, v in vals.items()}


lines = ["fixture\tseed\ttruth\tfragments\tmarginal_class\tmarginal_top\tmarginal_tied\tP_truth_marginal\tfixed_class\tfixed_top\tP_truth_fixed\t"
         "winner_best_state_flips\twinner_log_prior\twinner_vs_runnerup_dF\twinner_vs_runnerup_dE\twinner_vs_runnerup_dlogprior\tattribution\tseconds\tpeak_rss_gib"]
for d in sorted(glob.glob(f"{runs}/run.*")):
    if not os.path.exists(f"{d}/n.run.tsv"):
        lines.append(f"{os.path.basename(d)}\tNOT_COMPLETED")
        continue
    info = KV(f"{d}/n.run.tsv")
    rows = T(f"{d}/n.scaffolds.tsv")
    truth = info["truth"]
    marg = {r["scaffold"]: float(r["log_marginal"]) for r in rows}
    fixd = {r["scaffold"]: float(r["log_fixed"]) for r in rows}
    pm, pf = posterior(marg), posterior(fixd)
    cm, tm = classify(pm, truth, marg)
    cf, tf = classify(pf, truth, fixd)
    order = sorted(rows, key=lambda r: -float(r["log_marginal"]))
    w, ru = order[0], order[1]
    dF = float(w["best_state_fragment"]) - float(ru["best_state_fragment"])
    dE = float(w["best_state_exposure"]) - float(ru["best_state_exposure"])
    dP = float(w["best_state_log_prior"]) - float(ru["best_state_log_prior"])
    f = info["fixture"]
    if cm == "WRONG":
        attr = "MISLEADING"
    elif f.startswith("F2") and cf != "CORRECT" and cm == "CORRECT":
        attr = "CONTENT_ACCOMMODATION"
    elif abs(dF) < 1e-6 and abs(dE) < 1e-6 and abs(dP) > 1e-9:
        attr = "PRIOR_DRIVEN"
    elif f.startswith("F4") and cm == "CORRECT":
        attr = "ALLOCATION_EVIDENCE"
    elif cm == "UNCERTAIN":
        attr = "UNRESOLVED(uncertainty retained)"
    elif cm == "CORRECT":
        attr = "CORRECT"
    else:
        attr = "SPLIT"
    u = led.get(os.path.basename(d), {})
    lines.append("\t".join(map(str, [f, info["seed"], truth, info["fragments"], cm, max(pm, key=pm.get), ",".join(sorted(tm)), f"{pm.get(truth, 0):.4f}",
                                     cf, max(pf, key=pf.get), f"{pf.get(truth, 0):.4f}", w["best_state_flips"], w["best_state_log_prior"],
                                     f"{dF:.4f}", f"{dE:.4f}", f"{dP:.4f}", attr, info["seconds"], info["peak_rss_gib"]])))
open(prefix + ".tsv", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
