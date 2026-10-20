#!/usr/bin/env python3
"""M11-D reporting unit (amendment A3): descriptive tables from outcomes.a2 and fit.a2. Nothing here selects or changes a
predictor; the frozen comparison is fit.a2's performance.tsv. In-sample associations are descriptive only."""
import collections
import os
import statistics
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m11d  # noqa: E402

T, KV, write = m11d.T, m11d.KV, m11d.write
LOCI = ("ALL", "C4", "CYP2D6", "ACOT")


def f(x, nd=4):
    return "NA" if x is None else f"{x:.{nd}g}"


def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * (len(xs) - 1) + 0.5))] if xs else None


def main(out):
    src = os.path.dirname(out)
    blocks, runs = T(f"{src}/outcomes.a2/blocks.tsv"), T(f"{src}/outcomes.a2/runs.tsv")
    preds = {(r["cell"], r["block_id"]): r for r in T(f"{src}/fit.a2/predictions.tsv")}
    sel = lambda rows, L: [r for r in rows if L == "ALL" or r["locus"] == L]  # noqa: E731

    # 1. runs: counts, completeness, prediction coverage, ambiguous-completed vs resource-incomplete
    rr = []
    for L in LOCI:
        R = sel(runs, L)
        rr.append(dict(locus=L, cells=len(R), donors=len({r["donor"] for r in R}), families=len({r["family"] for r in R}),
                       completed_ranking=sum(r["completion"] == "completed_ranking" for r in R),
                       resource_incomplete=sum(r["completion"] != "completed_ranking" for r in R),
                       legacy_status=";".join(f"{k}={v}" for k, v in sorted(collections.Counter(r["legacy_status"] for r in R).items())),
                       sequence_status=";".join(f"{k}={v}" for k, v in sorted(collections.Counter(r["sequence_status"] for r in R).items())),
                       assessable=sum(r["assessable"] == "True" for r in R)))
    write(out, "runs_summary", rr)

    # 2. blocks: continuous excess and categorical outcomes, by locus (all predicted blocks; assessable subset noted)
    br = []
    for L in LOCI:
        B = sel(blocks, L)
        P = [r for r in B if r["predicted"] == "True"]
        ex = [float(r["excess"]) for r in P]
        nx = [float(r["normalized_excess"]) for r in P]
        qf = [float(r["QVfrac"]) for r in P if r["QVfrac"] not in ("", "None")]
        quint = collections.Counter(int(r["quintile"]) for r in P if r["quintile"] not in ("", "None"))
        art = [r for r in P if r["synonymous_artefact"] == "True"]
        br.append(dict(locus=L, blocks=len(B), truth_defined=sum(r["truth_defined"] == "True" for r in B), predicted=len(P),
                       no_prediction=sum(r["truth_defined"] == "True" and r["predicted"] != "True" for r in B),
                       fit_rows=sum(1 for k in preds if (L == "ALL" or preds[k]["locus"] == L)),
                       excess_sum=f(sum(ex), 8), excess_mean=f(statistics.mean(ex)), excess_median=f(statistics.median(ex)),
                       excess_p90=f(q(ex, 0.9)), excess_max=f(max(ex), 8),
                       norm_excess_mean=f(statistics.mean(nx)), norm_excess_median=f(statistics.median(nx)), norm_excess_max=f(max(nx)),
                       floor_attained=f"{sum(r['floor_attained'] == 'True' for r in P)}/{len(P)}",
                       catastrophic_ge_1000=sum(r["catastrophic"] == "True" for r in P),
                       top_QVfrac_ge_0_8=f"{sum(r['top_category'] == 'True' for r in P)}/{len(P)}",
                       quintiles_0_2_4_6_8=" / ".join(str(quint.get(i, 0)) for i in range(5)),
                       QVfrac_min=f(min(qf)), QVfrac_above_1=sum(r["QVfrac_above_1"] == "True" for r in P),
                       artefact_blocks=len(art), artefact_bp=f(sum(float(r["synonymous_artefact_bp"] or 0) for r in art), 8)))
    write(out, "blocks_summary", br)

    # 3. every block outside the top category, and every catastrophic block, with its out-of-family confidence
    cols = ("cell", "block_id", "kind", "excess", "normalized_excess", "QVfrac", "floor_attained", "catastrophic", "gq", "gq_source",
            "read_support", "census_route", "jp_emission_tie", "block_status", "synonymous_artefact", "synonymous_artefact_bp", "assessable")
    low = []
    for r in blocks:
        if r["predicted"] != "True" or (r["top_category"] == "True" and r["catastrophic"] != "True"):
            continue
        p = preds.get((r["cell"], r["block_id"]), {})
        low.append(dict({k: r[k] for k in cols}, p_B0=p.get("p_B0", "not_in_fit"), p_C1=p.get("p_C1", "not_in_fit"),
                        p_C2=p.get("p_C2", "not_in_fit")))
    low.sort(key=lambda r: -float(r["excess"]))
    write(out, "non_top_and_catastrophic", low)

    # 4. QVfrac > 1: the pair edits and alignment lengths behind each
    above = []
    for r in blocks:
        if r["QVfrac_above_1"] != "True":
            continue
        D = {x["id"]: x for x in T(f"{src}/outcomes.a2/{r['cell']}.distance.tsv")}
        pid = f"{r['cell']}.{r['block_id']}".replace(":", "_")
        above.append(dict(cell=r["cell"], block_id=r["block_id"], QVfrac=r["QVfrac"], excess=r["excess"],
                          **{f"{t}_edits_len": f"{D[pid + '.' + t]['edits']}/{D[pid + '.' + t]['aln_len']}"
                             for t in ("p1t1", "p2t2", "p1t2", "p2t1", "b1t1", "b2t2")}))
    write(out, "qvfrac_above_1", above)

    # 5. descriptive, in-sample: direction of association of each input with the label (never used for selection)
    fr = list(preds.values())
    dr = []
    for L in LOCI:
        S = sel(fr, L)
        y = [float(r["y"]) for r in S]
        x = [float(r["x"]) for r in S]
        pr = [float(r["prior"]) for r in S]
        top = [r for r in S if r["y"] == "1.0"]
        non = [r for r in S if r["y"] == "0.0"]
        m1 = m11d.fit_logistic([[a] for a in x], y)
        m2 = m11d.fit_logistic([[a, b] for a, b in zip(x, pr)], y)
        dr.append(dict(locus=L, rows=len(S), non_top=len(non),
                       gq_auc_in_sample=f(m11d.auc(x, y), 3),
                       gq_median_top=f(statistics.median(float(r["x"]) for r in top)) if top else "NA",
                       gq_median_non_top=f(statistics.median(float(r["x"]) for r in non)) if non else "NA",
                       prior_only_rate_top=f"{sum(r['prior'] == '1.0' for r in top)}/{len(top)}",
                       prior_only_rate_non_top=f"{sum(r['prior'] == '1.0' for r in non)}/{len(non)}",
                       c1_std_slope=f(m1["w"][1], 3), c2_std_slope_gq=f(m2["w"][1], 3), c2_std_slope_prior_only=f(m2["w"][2], 3),
                       b0_lofo_prediction_range=f"{f(min(float(r['p_B0']) for r in S))}..{f(max(float(r['p_B0']) for r in S))}"))
    write(out, "descriptive_association", dr)

    # 6. Q2 / Q3 diagnostics (cross-tabulated with the label; diagnostic only)
    qr = []
    for key in ("census_route", "jp_emission_tie", "read_support", "gq_source"):
        src_rows = [r for r in blocks if r["predicted"] == "True"]
        for val, n in sorted(collections.Counter(r[key] for r in src_rows).items()):
            V = [r for r in src_rows if r[key] == val]
            qr.append(dict(input=key, value=val, blocks=n, top=sum(r["top_category"] == "True" for r in V),
                           floor=sum(r["floor_attained"] == "True" for r in V), catastrophic=sum(r["catastrophic"] == "True" for r in V)))
    write(out, "diagnostic_crosstab", qr)

    for name in ("runs_summary", "blocks_summary", "descriptive_association"):
        print(open(f"{out}/{name}.tsv").read())
    print(f"non-top or catastrophic blocks {len(low)}; QVfrac > 1 blocks {len(above)}")
    return 0


if __name__ == "__main__":
    os.makedirs(sys.argv[1], exist_ok=True)
    sys.exit(main(sys.argv[1]))
