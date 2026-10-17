#!/usr/bin/env python3
"""S1-B fixture oracle: what sequence_first_v1 must report for one run, computed independently of the caller's report
from the run's own legacy tables and the model dump. Tiny fixtures only: block projections are COMPLETE (an exhaustive
spelling search over the model's allele tuples), which is exactly what the caller does NOT do (it reports EXPLORED
projections). Where the caller explored every projection of a group the two must agree; where it did not
(projections_unexplored > 0) the oracle shows what EXPLORED coverage misses.

Corrections over S1-A's projection (tests/analysis/s1a/s1a_project.py, kept unchanged as a record):
  * representative: the decision set's canonical member -- every decision-set group is numerically tied with the best
    under the caller's certified bounds, so the smallest sorted pair of full SHA-256 digests of the spelled sequences
    (amendment 1: the first version broke only bitwise-equal 17-digit scores, which sample order perturbs);
  * representation intervals: conservative enclosing spans over every pair of a group's labellings (coupled
    alternatives include the unchanged blocks between them), merged when overlapping; tuples are whole observed
    labellings restricted to the span, as unordered haplotype pairs (never per-block Cartesian combinations);
  * provisional: REPAIR_INCOMPLETE at round_cap or cumulative_cap with an incumbent; NO_PREDICTION otherwise.
"""
import hashlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s1a"))
from s1a_project import Model, kv, rows  # noqa: E402

COMPLETED_CALL = {"COMPLETE", "CONDITIONAL_GENOTYPE_EQUIVALENT"}


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def spans_of(tuples, n):
    spans = []
    for i in range(len(tuples)):
        for j in range(i + 1, len(tuples)):
            d = [b for b in range(n) if tuples[i][0][b] != tuples[j][0][b] or tuples[i][1][b] != tuples[j][1][b]]
            if d:
                spans.append((d[0], d[-1]))
    merged = []
    for s in sorted(spans):
        if merged and s[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], s[1]))
        else:
            merged.append(s)
    return merged


def unordered(a, b):
    return tuple(sorted((tuple(a), tuple(b))))


def oracle(prefix, dump_prefix):
    model = Model(dump_prefix)
    n = len(model.blocks)
    st = kv(prefix + ".mosaic_status.tsv")
    status, reason = st.get("status", "NA"), st.get("no_call_reason", "NA")
    out = {"prediction": "NO_PREDICTION", "sequence_status": "NA", "rep": None, "blocks": [], "intervals": [],
           "decision_set_groups": 0, "projections_unexplored": 0}
    if status == "REPAIR_INCOMPLETE" and reason in ("round_cap", "cumulative_cap") and \
            st.get("repair_final_incumbent_hap1", "NA") not in ("NA", ""):
        s1, s2 = model.spell(st["repair_final_incumbent_hap1"]), model.spell(st["repair_final_incumbent_hap2"])
        out.update(prediction="PROVISIONAL", sequence_status="NOT_DETERMINED", rep=tuple(sorted((sha(s1), sha(s2)))))
        a = [int(x) for x in st["repair_final_incumbent_hap1"].split(",")]
        b = [int(x) for x in st["repair_final_incumbent_hap2"].split(",")]
        out["blocks"] = [(model.ids[k], "PROVISIONAL_LABEL", f"{min(a[k], b[k])}/{max(a[k], b[k])}") for k in range(n)]
        return out
    labellings = {}
    for r in rows(prefix + ".mosaic_candidates.tsv"):
        labellings.setdefault(r["sequence_candidate_id"], set()).add((r["hap1_alleles"], r["hap2_alleles"]))
    exact = {r["sequence_candidate_id"]: float(r["score"]) for r in rows(prefix + ".equivalence_set.tsv")}
    groups = {}
    for r in rows(prefix + ".mosaic_scores.tsv"):
        g = r["sequence_candidate_id"]
        a, b = sorted(labellings[g])[0]
        s1, s2 = model.spell(a), model.spell(b)
        p1, p2 = model.projections(s1), model.projections(s2)
        # tuples oriented by sequence: A spells the SHA-256-smaller sequence (amendment 1: the first version sorted each
        # pair, so two distinct sequences' labellings were compared across haplotypes); unordered only for identical ones
        if sha(s2) < sha(s1):
            p1, p2 = p2, p1
        complete = sorted({(tuple(x), tuple(y)) for x in p1 for y in p2}) if s1 != s2 else \
            sorted({unordered(x, y) for x in p1 for y in p1})
        explored = {unordered([int(v) for v in x.split(",")], [int(v) for v in y.split(",")]) for x, y in labellings[g]}
        out["projections_unexplored"] += len({unordered(*t) for t in complete} - explored)
        geno = [{tuple(sorted((t[0][k], t[1][k]))) for t in complete} for k in range(n)]
        groups[g] = {"score": exact.get(g, float(r["score"])), "prob": float(r["candidate_set_probability"]),
                     "top": r["in_winner_set"] == "true", "sha": tuple(sorted((sha(s1), sha(s2)))), "geno": geno,
                     "tuples": complete}
    if status not in COMPLETED_CALL and not (status == "RESCORE_INCOMPLETE" and reason in ("distinct_tie", "divergent_tie")):
        return out
    top = {g: v for g, v in groups.items() if v["top"]}
    rep = min(top, key=lambda g: top[g]["sha"])
    union = [set().union(*(top[g]["geno"][k] for g in top)) for k in range(n)]
    all_identified = all(len(u) == 1 for u in union)
    out.update(prediction="PREDICTION", rep=top[rep]["sha"], decision_set_groups=len(top),
               sequence_status="RESOLVED" if len(top) == 1 else
               "GENOTYPE_RESOLVED_SEQUENCE_SET" if all_identified else "SEQUENCE_TIED")
    for k in range(n):
        if len(union[k]) == 1:
            g0 = next(iter(union[k]))
            lo = sum(v["prob"] for v in groups.values() if all(x != g0 for x in v["geno"][k]))
            hi = sum(v["prob"] for v in groups.values() if any(x != g0 for x in v["geno"][k]))
            out["blocks"].append((model.ids[k], "IDENTIFIED", f"{g0[0]}/{g0[1]}", lo, hi))
        else:
            out["blocks"].append((model.ids[k], "REPRESENTATION_AMBIGUOUS" if len(top) == 1 else "GENOTYPE_UNRESOLVED",
                                  "NA", None, None))
    for f, l in spans_of(top[rep]["tuples"], n):
        tuples = sorted({unordered(t[0][f:l + 1], t[1][f:l + 1]) for t in top[rep]["tuples"]})
        out["intervals"].append((f"{model.ids[f]}..{model.ids[l]}", tuples))
    return out
