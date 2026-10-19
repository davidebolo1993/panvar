#!/usr/bin/env python3
"""M1-C summary: one table per coverage from saved outputs (score-set explicit runs, content, decompositions, mass
reconciliation, construction). Likelihood (fragment + exposure) and prior are kept separate; a total is shown only where
the panel prior is defined. Deltas are alternative - reported.

  summary.py RUNS_DIR OUT_PREFIX
"""
import csv
import glob
import os
import sys


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path):
            k, _, v = line.rstrip("\n").partition("\t")
            out[k] = v
    return out


def main(runs, out):
    lines = ["coverage\tlabel\tfragment_term\texposure_term\tlikelihood\tprior\tprior_support\ttotal\t"
             "d_likelihood\td_fragment\td_exposure\tcontent_vs_reported_k31\tcontent_vs_reported_k151\t"
             "content_vs_truth_pooled_k31\tcontent_vs_truth_per_hap_k31\tcontent_vs_truth_pooled_k151\n"]
    dec = ["coverage\tpair\tfragment_delta\tfavouring_A\tfavouring_B\ttop10\ttop50\ttop200\ttop1000\tclasses\t"
           "placement_evidence_interpretable\n"]
    for cov in ("c4_10", "c4_30"):
        # Amendment 2: c4_30 is scored in two parts that share `reported`, which must be bit-equal across them.
        parts = [f"{runs}/score.{cov}/score"] if os.path.exists(f"{runs}/score.{cov}/score.scores_truth_free.tsv") else \
            [f"{runs}/score.{cov}.part{i}/score" for i in (1, 2) if os.path.exists(f"{runs}/score.{cov}.part{i}/score.scores_truth_free.tsv")]
        if not parts:
            lines.append(f"{cov}\tNOT RUN\n")
            continue
        rows, seen = [], {}
        for score in parts:
            for r in csv.DictReader(open(score + ".scores_truth_free.tsv"), delimiter="\t"):
                if r["label"] in seen:
                    same = all(seen[r["label"]][k] == r[k] for k in ("fragment_term_hex", "exposure_term_hex", "prior"))
                    print(f"{cov}: {r['label']} across parts bit-equal: {same}")
                    assert same, "a shared candidate differs across parts"
                    continue
                seen[r["label"]] = r
                rows.append(r)
        content = {}
        cp = f"{runs}/content.{cov}/{cov}.content.tsv"
        if os.path.exists(cp):
            for r in csv.DictReader(open(cp), delimiter="\t"):
                content[(r["label"], r["k"])] = r
        rep = next(r for r in rows if r["label"] == "reported")
        f0, e0, l0 = float(rep["fragment_term"]), float(rep["exposure_term"]), float(rep["likelihood"])
        for r in rows:
            c31, c151 = content.get((r["label"], "31"), {}), content.get((r["label"], "151"), {})
            lines.append("\t".join([cov, r["label"], r["fragment_term"], r["exposure_term"], r["likelihood"], r["prior"],
                                    r["prior_support"], r["total"], f"{float(r['likelihood']) - l0:.6f}",
                                    f"{float(r['fragment_term']) - f0:.6f}", f"{float(r['exposure_term']) - e0:.6f}",
                                    c31.get("vs_reported_pooled", "NA"), c151.get("vs_reported_pooled", "NA"),
                                    c31.get("vs_truth_pooled", "NA"), c31.get("vs_truth_per_hap", "NA"),
                                    c151.get("vs_truth_pooled", "NA")]) + "\n")
        for path in sorted(p for score in parts for p in glob.glob(score + ".decomp.*.summary.tsv")):
            s = kv(path)
            classes = ";".join(f"{k[6:-6]}={float(v):.3f}" for k, v in s.items() if k.startswith("class_") and k.endswith("_delta")
                               and not k.endswith("median_delta") and not k.endswith("min_delta") and not k.endswith("max_delta"))
            dec.append("\t".join([cov, f"{s['label_A']}-{s['label_B']}", s["fragment_delta"], s["fragments_favouring_A"],
                                  s["fragments_favouring_B"], s["top10_fragments_delta"], s["top50_fragments_delta"],
                                  s["top200_fragments_delta"], s["top1000_fragments_delta"], classes,
                                  s["placement_evidence_interpretable"]]) + "\n")
    open(out + ".alternatives.tsv", "w").write("".join(lines))
    open(out + ".decompositions.tsv", "w").write("".join(dec))
    print("".join(lines), end="")
    print("".join(dec), end="")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
