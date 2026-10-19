#!/usr/bin/env python3
"""M3-C cross-check against M3-B's saved REFERENCE scores (cold one-shot score-set) -- frozen before execution.
Every M3-C candidate (incumbent included) whose unordered allele-path pair also appears in M3-B's scored catalogue must
agree on fragment, exposure, prior and total within 1e-6 (absolute). Cells without M3-B scores (lzo_5x) report 0 matched.

  compare_m3b.py <m3c runs dir> <m3b runs dir> <out.tsv> cell...
"""
import csv
import os
import sys

m3c, m3b, out, cells = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
REF = {"mosaic_2x": "score.mosaic_2x/score.scores_truth_free.tsv", "mosaic_5x": "rank.a2/merged.mosaic_5x.tsv"}


def key(r):
    return tuple(sorted((r["hap1_alleles"], r["hap2_alleles"])))


lines, failed = [], 0
for cell in cells:
    new = {key(r): r for r in csv.DictReader(open(f"{m3c}/lik.{cell}/lik.scores.tsv"), delimiter="\t")}
    if cell not in REF:
        lines.append((cell, "matched", 0, "no M3-B reference scores (cell not scored in M3-B)"))
        continue
    ref = {key(r): r for r in csv.DictReader(open(f"{m3b}/{REF[cell]}"), delimiter="\t")}
    worst = {"fragment_term": 0.0, "exposure_term": 0.0, "prior": 0.0, "total": 0.0}
    matched = bad = 0
    for k, r in new.items():
        if k not in ref:
            continue
        matched += 1
        ok = True
        for col in worst:
            d = abs(float(r[col]) - float(ref[k][col]))
            worst[col] = max(worst[col], d)
            ok = ok and d <= 1e-6
        bad += 0 if ok else 1
    failed += bad
    lines.append((cell, "matched", matched, f"mismatched {bad}; worst " + ", ".join(f"{c}={v:.3g}" for c, v in worst.items())))
with open(out, "w") as f:
    f.write("cell\tkey\tvalue\tdetail\n")
    for l in lines:
        f.write("\t".join(str(x) for x in l) + "\n")
for l in lines:
    print(*l)
sys.exit(1 if failed else 0)
