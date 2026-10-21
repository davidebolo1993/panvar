#!/usr/bin/env python3
"""M13-A: variation-tolerant structural scoring, bounded prototype (diagnostic). Modes: prepare <out> | evaluate <out>.
prepare: the frozen contrasts and candidate identities (truth-free; both whole-haplotype assignments of each candidate in the
reported winner's fixed context). evaluate: rankings, decomposition, approximation sensitivity, and -- only here -- the truth
comparison (per-block excess with M11-D's evaluator convention and floor)."""
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
import m11e  # noqa: E402

csv.field_size_limit(sys.maxsize)
m11d = m11e.m11d
T, KV, write = m11d.T, m11d.KV, m11d.write
M11D, R_ = m11e.M11D, m11d.R
OUT = f"{R_}/m13a_variation_tolerant_20261021"
M12B, M12C = f"{R_}/m12b_dosage_first_20261021", f"{R_}/m12c_dosage_first_completion_20261021"

# frozen contrasts (DESIGN.md §4). Candidate block pairs are (hap1-context allele, hap2-context allele); "comparator" for the
# controls is the run's own best recorded candidate that differs from the winner only at the target block with a different
# unordered block pair (truth-free).
CONTRASTS = [
    dict(cid="K4", cell="cyp.NA19650", block="sv:9", role="copy-number error", m12=(22, 80)),
    dict(cid="K5", cell="cyp.NA19983", block="sv:9", role="within-state content/allocation", m12=(104, 97)),
    dict(cid="K1", cell="c4.NA24385", block="sv:7", role="allocation after pooled-state correction", m12=(115, 48)),
    dict(cid="R1", cell="cyp.HG02953", block="sv:9", role="correct control", m12=None),
    dict(cid="R2", cell="c4.HG00358", block="sv:7", role="correct control", m12=None),
    dict(cid="R3", cell="acot.HG03209", block="sv:13", role="correct control", m12=None),
]


def winner(cell):
    return next(r for r in T(f"{M11D}/runs/call.{cell}/call.equivalence_set.tsv") if r["rank"] == "1" and r["representative"] == "true")


def comparator(cell, k, h1, h2):
    sc = {r["sequence_candidate_id"]: float(r["score"]) for r in T(f"{M11D}/runs/call.{cell}/call.mosaic_scores.tsv")}
    best = None
    for r in T(f"{M11D}/runs/call.{cell}/call.mosaic_candidates.tsv"):
        a, b = r["hap1_alleles"].split(","), r["hap2_alleles"].split(",")
        for x, y in ((a, b), (b, a)):
            if all(x[i] == h1[i] and y[i] == h2[i] for i in range(len(h1)) if i != k) and sorted((x[k], y[k])) != sorted((h1[k], h2[k])):
                s = sc.get(r["sequence_candidate_id"])
                if s is not None and (best is None or s > best[0]):
                    best = (s, int(x[k]), int(y[k]))
    return best


def prepare(out):
    rows = []
    for c in CONTRASTS:
        w = winner(c["cell"])
        h1, h2 = w["hap1_alleles"].split(","), w["hap2_alleles"].split(",")
        ev = m11e.eval_rows(c["cell"])
        k = [e["block_id"] for e in ev].index(c["block"])
        e = ev[k]
        inc = (int(h1[k]), int(h2[k]))
        named = [("incumbent", inc)]
        if c["m12"] and tuple(sorted(c["m12"])) != tuple(sorted(inc)):
            named.append(("m12_selected", c["m12"]))
        if c["role"] == "correct control":
            cp = comparator(c["cell"], k, h1, h2)
            named.append(("comparator", (cp[1], cp[2])))
        else:
            named.append(("best_available", (int(e["truth1_closest_allele"]), int(e["truth2_closest_allele"]))))
        lines, seen = [], set()
        for name, (a, b) in named:
            for tag, (x, y) in (("A", (a, b)), ("B", (b, a))):
                if (x, y) in seen:
                    continue
                seen.add((x, y))
                p1, p2 = list(h1), list(h2)
                p1[k], p2[k] = str(x), str(y)
                exp = w["score"] if (name == "incumbent" and tag == "A") else "NA"
                lines.append(f"{name}_{tag}\t{','.join(p1)}\t{','.join(p2)}\t{exp}\n")
                rows.append(dict(contrast=c["cid"], cell=c["cell"], block=c["block"], candidate=f"{name}_{tag}", block_pair=f"{x},{y}"))
        open(f"{out}/{c['cid']}.candidates.tsv", "w").writelines(lines)
    write(out, "candidates", rows)
    for r in rows:
        print(r)
    return 0


if __name__ == "__main__":
    mode, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"prepare": prepare}[mode](out))
