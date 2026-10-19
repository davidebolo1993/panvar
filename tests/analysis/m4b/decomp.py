#!/usr/bin/env python3
"""M4-B changed-call decomposition (truth-free; runs only when the two arms selected different unordered sequence pairs).

W_c = the constant arm's winning representative; W_d = the distance arm's. Both are scored under BOTH specifications by
score-set (the constant spec against the constant run, the distance spec against the distance run; each run's own winner
declared with its saved total, which must reproduce). Reported: fragment, exposure and likelihood (spec-independent) and the
prior under each spec, so a changed call separates into
  likelihood   L(W_d) - L(W_c)
  prior        P_const(W_d) - P_const(W_c) and P_dist(W_d) - P_dist(W_c)
  search pool  whether W_d's unordered allele paths are in the constant arm's scored pool, and W_c's in the distance arm's.
Classification (frozen):
  PRIOR_RESCORING  W_d was scored by the constant arm (and W_c by the distance arm): the same pool, re-ranked by the prior
  POOL_CHANGE      a winner was never scored by the other arm: the search reached a different pool
Neither is read as an improvement or a deterioration; accuracy comes from the evaluator.

  decomp.py <probe> <cell> <const prefix> <dist prefix> <gfa> <reference> <curation|none> <reads dir> <lambda> <out dir>
"""
import csv
import subprocess
import sys

probe, cell, cpre, dpre, gfa, ref, cur, reads, lam, out = sys.argv[1:11]


def rep(prefix):
    r = next(r for r in csv.DictReader(open(prefix + ".equivalence_set.tsv"), delimiter="\t") if r["rank"] == "1" and r["representative"] == "true")
    return r["hap1_alleles"], r["hap2_alleles"], r["score"]


def unordered(a, b):
    return tuple(sorted((a, b)))


def pool(prefix):
    return {unordered(r["hap1_alleles"], r["hap2_alleles"]) for r in csv.DictReader(open(prefix + ".mosaic_candidates.tsv"), delimiter="\t")}


def scored(prefix, rows, tag, extra):
    cands = f"{out}/{tag}.tsv"
    with open(cands, "w") as f:
        for label, h1, h2, exp in rows:
            f.write(f"{label}\t{h1}\t{h2}\t{exp}\n")
    cmd = [probe, "score-set", "--out", f"{out}/{tag}", "--model-gfa", gfa, "--reference", ref, "--donor", "NA", "--run", prefix,
           "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq", "--library", f"{reads}/library.txt", "--threads", "4",
           "--candidates", cands] + ([] if cur == "none" else ["--panel-curation", cur]) + extra
    p = subprocess.run(cmd, capture_output=True, text=True)
    open(f"{out}/{tag}.log", "w").write(" ".join(cmd) + "\n" + p.stdout + p.stderr)
    if p.returncode != 0:
        raise SystemExit(f"score-set {tag} failed (exit {p.returncode}); see {out}/{tag}.log")
    return {r["label"]: r for r in csv.DictReader(open(f"{out}/{tag}.scores_truth_free.tsv"), delimiter="\t")}


try:
    c1, c2, cs = rep(cpre)
    d1, d2, ds = rep(dpre)
except (FileNotFoundError, StopIteration):
    open(f"{out}/decomp.txt", "w").write(f"{cell}\tNOT_ASSESSED\tan arm has no winning representative (refusal or no-call)\n")
    print(f"{cell}: NOT_ASSESSED (an arm has no winning representative)")
    sys.exit(0)
if unordered(c1, c2) == unordered(d1, d2):
    open(f"{out}/decomp.txt", "w").write(f"{cell}\tUNCHANGED\tthe two arms selected the same unordered allele paths\n")
    print(f"{cell}: UNCHANGED (same unordered allele paths); nothing to decompose")
    sys.exit(0)
sc = scored(cpre, [("W_c", c1, c2, cs), ("W_d", d1, d2, "NA")], "const_spec", [])
sd = scored(dpre, [("W_c", c1, c2, "NA"), ("W_d", d1, d2, ds)], "dist_spec",
            ["--copying-transition", "distance_ref_core_v1_experimental", "--copying-length-bp", lam])
L = lambda r: float(r["fragment_term"]) + float(r["exposure_term"])
in_const, in_dist = unordered(d1, d2) in pool(cpre), unordered(c1, c2) in pool(dpre)
cls = "PRIOR_RESCORING" if in_const and in_dist else "POOL_CHANGE"
lines = [f"cell\t{cell}", f"class\t{cls}", f"W_d_in_constant_pool\t{in_const}", f"W_c_in_distance_pool\t{in_dist}",
         f"delta_fragment(W_d-W_c)\t{float(sc['W_d']['fragment_term']) - float(sc['W_c']['fragment_term']):.6f}",
         f"delta_exposure(W_d-W_c)\t{float(sc['W_d']['exposure_term']) - float(sc['W_c']['exposure_term']):.6f}",
         f"delta_likelihood(W_d-W_c)\t{L(sc['W_d']) - L(sc['W_c']):.6f}",
         f"delta_prior_constant(W_d-W_c)\t{float(sc['W_d']['prior']) - float(sc['W_c']['prior']):.6f}",
         f"delta_prior_distance(W_d-W_c)\t{float(sd['W_d']['prior']) - float(sd['W_c']['prior']):.6f}",
         f"likelihood_spec_independent\t{abs(L(sc['W_d']) - L(sd['W_d'])) < 1e-9 and abs(L(sc['W_c']) - L(sd['W_c'])) < 1e-9}"]
open(f"{out}/decomp.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
