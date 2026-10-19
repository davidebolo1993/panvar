#!/usr/bin/env python3
"""M1-E summary (saved outputs only): per coverage, the candidate table (fragment, exposure, likelihood, prior support), the
gain of each single-component correction and of all four over no correction (for each arrangement), the interaction
(combined gain minus the sum of single gains), the allocation contrast (1+3 minus 2+2) per patch set, the compact
fragment-concentration summaries, and the content differences. Nothing is assumed additive.

  summary.py RUNS_DIR OUT_PREFIX
"""
import csv
import glob
import os
import sys


def kv(path):
    return dict(line.rstrip("\n").split("\t", 1) for line in open(path) if "\t" in line)


def main(runs, out):
    content = {}
    for r in csv.DictReader(open(f"{runs}/content/content.tsv"), delimiter="\t"):
        content[(r["label"], r["k"])] = r
    t, g, al, cn = [], [], [], []
    for cov in ("c4_10", "c4_30"):
        rows = {r["label"]: r for r in csv.DictReader(open(f"{runs}/score.{cov}/score.scores_truth_free.tsv"), delimiter="\t")}
        L = {k: float(v["likelihood"]) for k, v in rows.items()}
        for k, r in rows.items():
            c31, c151 = content.get((k, "31"), {}), content.get((k, "151"), {})
            t.append("\t".join([cov, k, r["fragment_term"], r["exposure_term"], r["likelihood"], r["prior_support"], r["total"],
                                f"{L[k] - L['reported']:.3f}", c31.get("pooled_vs_truth", "NA"), c151.get("pooled_vs_truth", "NA"),
                                c31.get("vs_matched_arrangement", "NA"), c151.get("vs_matched_arrangement", "NA")]) + "\n")
        for arr in ("a22", "e13"):
            singles = {p: L[f"{arr}_{p}"] - L[f"{arr}_none"] for p in ("U1", "T1", "U2", "T2")}
            combined = L[f"{arr}_all"] - L[f"{arr}_none"]
            for p, v in singles.items():
                g.append(f"{cov}\t{arr}\t{p}\t{v:.3f}\t{v / combined:.3f}\n")
            g.append(f"{cov}\t{arr}\tall\t{combined:.3f}\t1.000\n")
            g.append(f"{cov}\t{arr}\tinteraction (all - sum of singles)\t{combined - sum(singles.values()):.3f}\tNA\n")
        for p in ("none", "U1", "T1", "U2", "T2", "all"):
            x, y = rows[f"e13_{p}"], rows[f"a22_{p}"]
            al.append(f"{cov}\t{p}\t{float(x['likelihood']) - float(y['likelihood']):.3f}\t"
                      f"{float(x['fragment_term']) - float(y['fragment_term']):.3f}\t{float(x['exposure_term']) - float(y['exposure_term']):.3f}\t"
                      f"{x['fragment_term_hex'] == y['fragment_term_hex']}\n")
        for path in sorted(glob.glob(f"{runs}/score.{cov}/score.contrast.*.tsv")):
            s = kv(path)
            cn.append("\t".join([cov, f"{s['label_A']}-{s['label_B']}", s["fragment_delta"], s["exposure_delta"],
                                 s["fragments_favouring_A"], s["fragments_favouring_B"], s["top10_fragments_delta"],
                                 s["top50_fragments_delta"], s["top200_fragments_delta"], s["top1000_fragments_delta"]]) + "\n")
    sections = [("candidates", "coverage\tlabel\tfragment_term\texposure_term\tlikelihood\tprior_support\ttotal\td_likelihood_vs_reported\t"
                 "pooled_vs_truth_k31\tpooled_vs_truth_k151\tvs_matched_arrangement_k31\tvs_matched_arrangement_k151\n", t),
                ("gains", "coverage\tarrangement\tpatch\tgain_over_none\tshare_of_combined\n", g),
                ("allocation", "coverage\tpatch\te13_minus_a22_likelihood\tfragment\texposure\tfragment_bit_equal\n", al),
                ("contrasts", "coverage\tpair\tfragment_delta\texposure_delta\tfavouring_A\tfavouring_B\ttop10\ttop50\ttop200\ttop1000\n", cn)]
    for name, header, rows in sections:
        open(f"{out}.{name}.tsv", "w").write(header + "".join(rows))
        print(f"== {name}\n" + header + "".join(rows), end="")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
