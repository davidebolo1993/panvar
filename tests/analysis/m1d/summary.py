#!/usr/bin/env python3
"""M1-D summary (saved outputs only). Per coverage and candidate: fragment, exposure, likelihood (their sum), prior support
(never a neutral prior; totals only where defined), the unordered sv:7 block distance to truth, the pooled block k-mer
discrepancy to truth and to the reported blocks (canonical, k = 31 and 151), the likelihood change relative to the
reported pair, and the share of the truth-block reference's gain recovered. Plus allocation comparisons on the matched
composed inventory and the reachability tables.

  summary.py RUNS_DIR OUT_PREFIX
"""
import csv
import os
import sys
from collections import Counter

COMP = str.maketrans("ACGT", "TGCA")


def fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line)
    return {k: "".join(v) for k, v in out.items()}


def kmers(seq, k):
    c, rc, n = Counter(), seq.translate(COMP)[::-1], len(seq)
    for i in range(n - k + 1):
        f, r = seq[i:i + k], rc[n - k - i:n - i]
        c[hash(f if f < r else r)] += 1
    return c


def diff(a, b):
    return sum(abs(a[x] - b[x]) for x in set(a) | set(b)) / 2


def main(runs, out):
    lines = ["coverage\tlabel\tfragment_term\texposure_term\tlikelihood\tprior_support\ttotal\td_likelihood\td_fragment\t"
             "d_exposure\tshare_of_truthblock_gain\tblock_distance_to_truth\tpooled_vs_truth_k31\tpooled_vs_truth_k151\t"
             "pooled_vs_reported_k31\tpooled_vs_reported_k151\n"]
    alloc = ["coverage\tcomparison\td_likelihood\td_fragment\td_exposure\n"]
    for cov in ("c4_10", "c4_30"):
        sp = f"{runs}/score.{cov}/score.scores_truth_free.tsv"
        if not os.path.exists(sp):
            lines.append(f"{cov}\tNOT RUN\n")
            continue
        rows = {r["label"]: r for r in csv.DictReader(open(sp), delimiter="\t")}
        dist = {}
        for r in csv.DictReader(open(f"{runs}/cand.dist.{cov}.a2/distances.tsv"), delimiter="\t"):
            dist[(r["record_a"], r["record_b"])] = int(r["distance"])
        blocks = fasta(f"{runs}/construct.{cov}.a2/candidates.blocks.fa")
        content = {}
        for k in (31, 151):
            cache = {}

            def cnt(name):
                if name not in cache:
                    cache[name] = kmers(blocks[name], k)
                return cache[name]

            tp = cnt("truth.1") + cnt("truth.2")
            rp = cnt("reported.1") + cnt("reported.2")
            for label in rows:
                pooled = cnt(f"{label}.1") + cnt(f"{label}.2")
                content[(label, k)] = (diff(pooled, tp), diff(pooled, rp))
        rep, tb = rows["reported"], rows["truthblock_a"]
        l0, f0, e0 = float(rep["likelihood"]), float(rep["fragment_term"]), float(rep["exposure_term"])
        gain = float(tb["likelihood"]) - l0
        for label, r in rows.items():
            d = min(dist[(f"{label}.1", "truth.1")] + dist[(f"{label}.2", "truth.2")],
                    dist[(f"{label}.1", "truth.2")] + dist[(f"{label}.2", "truth.1")])
            dl = float(r["likelihood"]) - l0
            lines.append("\t".join([cov, label, r["fragment_term"], r["exposure_term"], r["likelihood"], r["prior_support"],
                                    r["total"], f"{dl:.3f}", f"{float(r['fragment_term']) - f0:.3f}",
                                    f"{float(r['exposure_term']) - e0:.3f}", f"{dl / gain:.3f}", str(d),
                                    f"{content[(label, 31)][0]:g}", f"{content[(label, 151)][0]:g}",
                                    f"{content[(label, 31)][1]:g}", f"{content[(label, 151)][1]:g}"]) + "\n")
        if "oracle_compose_a" in rows:
            base = rows["oracle_compose_a"]
            for other in ("oracle_compose_b", "oracle_compose_c", "oracle_compose_d31", "oracle_compose_e13"):
                o = rows[other]
                alloc.append(f"{cov}\t{other} - oracle_compose_a\t"
                             + "\t".join(f"{float(o[k]) - float(base[k]):.3f}" for k in ("likelihood", "fragment_term", "exposure_term"))
                             + "\n")
    open(out + ".candidates.tsv", "w").write("".join(lines))
    open(out + ".allocation.tsv", "w").write("".join(alloc))
    print("".join(lines), end="")
    print("".join(alloc), end="")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
