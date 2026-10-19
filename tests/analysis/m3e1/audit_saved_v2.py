#!/usr/bin/env python3
"""M3-E.1 independent saved-evidence audit, contract v2 (CONTRACT_V2.md). TRUTH-FREE; no genotyping reruns.
Envelope certification over every relevant group; mixed groups preserved; absent comparison non-affirmative. For each
deciding inequality the slack is reported; a decision whose slack is within ROUNDING (1e-4) of the saved tables'
precision is flagged NEEDS_FULL_PRECISION (not certified from rounded tables).
  audit_saved_v2.py <out prefix> <v1 audit blocks.tsv> <label>=<run prefix> ...
"""
import csv
import sys
from collections import Counter

ROUNDING = 1e-4
NINF = float("-inf")


def rows(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def audit(run):
    groups = {}
    for r in rows(run + ".mosaic_scores.tsv"):
        groups[r["sequence_candidate_id"]] = dict(F=float(r["fragment_term"]), E=float(r["exposure_term"]), P=float(r["sequence_pair_prior"]),
                                                  u=float(r["numeric_error_bound"]), win=r["in_winner_set"] == "true", rows=set())
    for c in rows(run + ".mosaic_candidates.tsv"):
        g = groups.get(c["sequence_candidate_id"])
        if g is not None:
            g["rows"].add(tuple(tuple(sorted(p)) for p in zip(map(int, c["hap1_alleles"].split(",")), map(int, c["hap2_alleles"].split(",")))))
    w = next(g for g in groups.values() if g["win"])
    ids = [r["block_id"] for r in rows(run + ".block_genotypes.tsv")]
    out = []
    for b, bid in enumerate(ids):
        called = {row[b] for row in w["rows"]}
        rec = {"block": bid}
        if len(called) != 1:
            rec.update(label="AMBIGUOUS_PROJECTION", slack=float("nan"))
            out.append(rec)
            continue
        gb = next(iter(called))
        env = {k: [NINF, NINF] for k in "SDM"}
        n = Counter()
        for g in groups.values():
            at = {row[b] for row in g["rows"]}
            cls = "S" if at == {gb} else ("D" if gb not in at else "M")
            n[cls] += 1
            L = g["F"] + g["E"]
            env[cls][0] = max(env[cls][0], L - g["u"])
            env[cls][1] = max(env[cls][1], L + g["u"])
        (s_lo, s_hi), (d_lo, d_hi), (m_lo, m_hi) = env["S"], env["D"], env["M"]
        rec.update(disagreeing=n["D"], mixed=n["M"])
        if n["D"] == 0 and n["M"] == 0:
            rec.update(label="NO_EVALUATED_ALTERNATIVE", slack=float("nan"))
            out.append(rec)
            continue
        top = max(d_hi, m_hi)
        supported = (n["D"] == 0 or s_lo > d_hi) and (n["M"] == 0 or s_lo > m_hi)
        prior_dep = n["D"] > 0 and d_lo > s_hi and (n["M"] == 0 or d_lo > m_hi)
        if supported:
            label, slack = "DATA_SUPPORTED", s_lo - top
        elif prior_dep:
            label, slack = "PRIOR_DEPENDENT", min(d_lo - s_hi, d_lo - m_hi if n["M"] else float("inf"))
        elif n["M"] and m_hi >= s_lo:
            label, slack = "AMBIGUOUS_COMPETITOR", m_hi - s_lo
        else:
            # DATA_TIED: the closest the decision came to either certified side
            label, slack = "DATA_TIED", min(top - s_lo, s_hi - d_lo if n["D"] else float("inf"))
        rec.update(label=label, gap=s_lo - top, slack=slack)
        out.append(rec)
    return out


if __name__ == "__main__":
    prefix, v1_path = sys.argv[1], sys.argv[2]
    v1 = {(r["run"], r["block"]): r["genotype_support_v1x"] for r in rows(v1_path)}
    lines = ["run\tblock\tv1_label\tv2_label\tchanged\tsupport_gap\tdecision_slack\tneeds_full_precision\tdisagreeing\tmixed"]
    summary = ["run\tv2_labels\tchanged_from_v1\tneeds_full_precision"]
    for arg in sys.argv[3:]:
        label, run = arg.split("=", 1)
        recs = audit(run)
        changed = needs = 0
        for r in recs:
            old = v1[(label, r["block"])]
            ch = old != r["label"]
            nfp = r.get("slack") == r.get("slack") and abs(r["slack"]) <= ROUNDING
            changed += ch
            needs += nfp
            lines.append("\t".join(map(str, [label, r["block"], old, r["label"], ch, f"{r.get('gap', float('nan')):.6g}", f"{r['slack']:.6g}",
                                             nfp, r.get("disagreeing", "-"), r.get("mixed", "-")])))
        summary.append(f"{label}\t" + ",".join(f"{k}:{v}" for k, v in sorted(Counter(r['label'] for r in recs).items())) + f"\t{changed}\t{needs}")
    open(prefix + ".blocks.tsv", "w").write("\n".join(lines) + "\n")
    open(prefix + ".summary.tsv", "w").write("\n".join(summary) + "\n")
    print("\n".join(summary))
