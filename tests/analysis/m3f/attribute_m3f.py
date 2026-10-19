#!/usr/bin/env python3
"""M3-F attribution of errors the joint search leaves in place (DIAGNOSTIC; truth-guided; after the frozen predictions;
never feeds a prediction or a setting).

For each error cell whose ON arm still has positive excess, one class (precedence as listed; frozen):
  PASS_SKIPPED_OR_INCOMPLETE          the on arm's joint_search_status does not start with 'complete'
  EVALUATED_ALTERNATIVE_LOSES         the truth (mosaic: exact spelled pair in the diagnostic catalogue built from the off
                                      incumbent, ranked below first; LOO: a floor-attaining candidate among the caller's scored
                                      candidates) was scored and ranks below the selection under the current model (a mosaic
                                      truth that ranks first in the catalogue while the call stays wrong is UNRESOLVED)
  OUTSIDE_TWO_BLOCK_NEIGHBOURHOOD     mosaic: the truth differs from the off incumbent (orientation aligned as in the probe)
                                      on both homologues or at more than two blocks of one homologue
  ALTERNATIVE_ABSENT_FROM_SHORTLIST   mosaic: within the neighbourhood, but a truth allele was not kept in the top-k
                                      likelihood shortlist
  UNRESOLVED                          otherwise (including the LOO case when the floor is not attained by any scored
                                      candidate: the truth is not spellable and the two-block probe cannot enumerate it)
No class claims exhaustive coverage: each is relative to the frozen k = 4, 2,048-bounded same-homologue catalogue.

  attribute_m3f.py <runs dir> <out prefix> cell...
"""
import csv
import os
import sys

runs, prefix, cells = sys.argv[1], sys.argv[2], sys.argv[3:]


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path):
            k, _, v = line.rstrip("\n").partition("\t")
            out.setdefault(k, v)
    return out


def table(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.exists(path) else []


rows = []
for cell in cells:
    ev_off, ev_on = kv(f"{runs}/eval.{cell}.off/eval.summary.tsv"), kv(f"{runs}/eval.{cell}.on/eval.summary.tsv")
    st_on = kv(f"{runs}/run.{cell}.on/call.mosaic_status.tsv")
    x_on = ev_on.get("excess_over_floor_total")
    rec = {"cell": cell, "excess_off": ev_off.get("excess_over_floor_total", "-"), "excess_on": x_on or "-",
           "joint_search_status": st_on.get("joint_search_status", "-"), "joint_search_decision": st_on.get("joint_search_decision", "-"),
           "pool_floor_attained_on": ev_on.get("floor_attained_in_candidates", "-"),
           "pool_closest_total_gap_on": ev_on.get("closest_candidate_total_score_gap_to_winner", "-")}
    detail = []
    if x_on in (None, "NA", "-"):
        cls = "NOT_ATTRIBUTED (no evaluable on call)"
    elif x_on == "0":
        cls = "NOT_ATTRIBUTED (on arm excess 0)"
    elif not rec["joint_search_status"].startswith("complete"):
        cls = "PASS_SKIPPED_OR_INCOMPLETE"
    elif cell.startswith("loo"):
        if ev_on.get("floor_attained_in_candidates") == "true":
            cls = "EVALUATED_ALTERNATIVE_LOSES"
            detail.append(f"floor-attaining candidate scored; total-score gap to winner {ev_on.get('closest_candidate_total_score_gap_to_winner')}")
        else:
            cls = "UNRESOLVED"
            detail.append(f"availability floor {ev_on.get('availability_floor_total')} (truth not spellable); no scored candidate attains the floor")
        s = kv(f"{runs}/structure.{cell}/{cell}.structure.tsv")
        if s:
            detail.append(f"sv:7 truth units {s.get('truth_total_units')}; off {s.get('historical_called_units')} {s.get('historical_structure_verdict')}; "
                          f"on {s.get('fresh_called_units')} {s.get('fresh_structure_verdict')}")
    else:
        tb = kv(f"{runs}/tbeval.{cell}/eval.eval_summary.tsv")
        lik = {r["label"]: r for r in table(f"{runs}/lik.{cell}/lik.scores.tsv")}
        if not tb or "incumbent" not in lik:
            cls = "UNRESOLVED"
            detail.append("diagnostic catalogue or its evaluation missing")
        else:
            t = [list(map(int, tb["truth_hap1_path"].split(","))), list(map(int, tb["truth_hap2_path"].split(",")))]
            inc = [list(map(int, lik["incumbent"]["hap1_alleles"].split(","))), list(map(int, lik["incumbent"]["hap2_alleles"].split(",")))]
            a = sum((i != j) for i, j in zip(inc[0], t[0])) + sum((i != j) for i, j in zip(inc[1], t[1]))
            c = sum((i != j) for i, j in zip(inc[0], t[1])) + sum((i != j) for i, j in zip(inc[1], t[0]))
            tt = t if a <= c else [t[1], t[0]]
            diff = [[b for b, (i, j) in enumerate(zip(inc[h], tt[h])) if i != j] for h in (0, 1)]
            rec["truth_diff_blocks_h0"], rec["truth_diff_blocks_h1"] = len(diff[0]), len(diff[1])
            label = tb.get("truth_generated_label", "absent")
            rec["truth_in_catalogue"] = label != "absent"
            rec["truth_rank"] = tb.get("truth_rank", "-")
            rec["catalogue_scored"] = tb.get("catalogue_scored", "-")
            ranks = {k: v for k, v in tb.items() if k.startswith("truth_allele_proposal_rank_")}
            if label != "absent":
                # rank counts the incumbent row too: rank 1 means the truth outscores everything in the catalogue, yet the call
                # is wrong (decision, certification or continuation), which is not 'loses under the model'
                cls = "EVALUATED_ALTERNATIVE_LOSES" if tb.get("truth_rank") != "1" else "UNRESOLVED"
                tr, ir = lik[label], lik["incumbent"]
                winner = max(lik.values(), key=lambda r: float(r["total"]))
                detail.append(f"truth {label} rank {tb.get('truth_rank')}/{tb.get('catalogue_scored')}; truth - incumbent: total "
                              f"{float(tr['total']) - float(ir['total']):.4f} (F {float(tr['fragment_term']) - float(ir['fragment_term']):.4f}, "
                              f"E {float(tr['exposure_term']) - float(ir['exposure_term']):.4f}, P {float(tr['prior']) - float(ir['prior']):.4f}); "
                              f"truth - catalogue best ({winner['label']}): {float(tr['total']) - float(winner['total']):.4f}")
            elif (diff[0] and diff[1]) or len(diff[0]) > 2 or len(diff[1]) > 2:
                cls = "OUTSIDE_TWO_BLOCK_NEIGHBOURHOOD"
                detail.append(f"truth differs at {len(diff[0])} block(s) on homologue 0 and {len(diff[1])} on homologue 1")
            elif any("not kept" in v or "not_ranked" in v for v in ranks.values()):
                cls = "ALTERNATIVE_ABSENT_FROM_SHORTLIST"
                detail.append("; ".join(f"{k[len('truth_allele_proposal_rank_'):]} {v}" for k, v in ranks.items()))
            else:
                cls = "UNRESOLVED"
                detail.append("; ".join(f"{k[len('truth_allele_proposal_rank_'):]} {v}" for k, v in ranks.items()) or "no differing allele")
    rec["class"] = cls
    rec["detail"] = " | ".join(detail) or "-"
    rows.append(rec)
cols = ["cell", "class", "excess_off", "excess_on", "joint_search_status", "joint_search_decision", "truth_in_catalogue", "truth_rank",
        "catalogue_scored", "truth_diff_blocks_h0", "truth_diff_blocks_h1", "pool_floor_attained_on", "pool_closest_total_gap_on", "detail"]
with open(f"{prefix}.attribution.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r.get(k, "-")) for k in cols) + "\n")
print(open(f"{prefix}.attribution.tsv").read(), end="")
