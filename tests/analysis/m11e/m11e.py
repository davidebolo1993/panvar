#!/usr/bin/env python3
"""M11-E: reporting/eligibility corrections for M11-D and a bounded attribution of its five catastrophic blocks.

Modes:
  fixes <out>    corrected interval-artefact totals, calibration-eligibility annotations (no refit, no relabelling),
                 the QVfrac per-homologue breakdown of the catastrophic blocks
  prepare <out>  per case: the reported winner, the two single-block replacements by the evaluator's best-available pair,
                 spelled and verified (incumbent bytes = the run's haplotype SHA-256; only the block changes)
  trace <out>    per case: score-set reproduction, probe-byte identity, term differences, proposal / shortlist /
                 evaluation status, dosage audit, and the frozen classification
Inputs are the sealed M11-D package; nothing in it is written.
"""
import functools
import hashlib
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11d"))
import m11d  # noqa: E402

T, KV, write, qv = m11d.T, m11d.KV, m11d.write, m11d.qv
fasta = functools.lru_cache(maxsize=None)(m11d.fasta)
M11D = m11d.OUT
OUT = f"{m11d.R}/m11e_catastrophic_attribution_20261021"
TOL_TOTAL = 1e-6      # score-set's own reproduction tolerance (nats)
TOL_PRINTED = 1e-5    # recorded run scores are printed to ~12 significant digits

# frozen cases (DESIGN.md §2): the five M11-D catastrophic blocks, then one floor-attaining repeated_unit control per locus,
# the lowest sha256("M11-E:control:" + cell + ":" + block) among that locus's floor-attaining repeated_unit blocks
CASES = [dict(case="K1", cell="c4.NA24385", block="sv:7", role="catastrophic"),
         dict(case="K2", cell="acot.NA19909", block="sv:13", role="catastrophic"),
         dict(case="K3", cell="c4.NA19434", block="sv:7", role="catastrophic"),
         dict(case="K4", cell="cyp.NA19650", block="sv:9", role="catastrophic"),
         dict(case="K5", cell="cyp.NA19983", block="sv:9", role="catastrophic"),
         dict(case="R1", cell="cyp.HG02953", block="sv:9", role="control"),
         dict(case="R2", cell="c4.HG00358", block="sv:7", role="control"),
         dict(case="R3", cell="acot.HG03209", block="sv:13", role="control")]


def case(cid):
    return next(c for c in CASES if c["case"] == cid)


def fnv64(s):
    h = 1469598103934665603
    for b in s.encode():
        h = ((h ^ b) * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return f"{h:016x}"


def vec(s):
    return [int(x) for x in s.split(",")]


def eval_rows(cell):
    return T(f"{M11D}/runs/eval.{cell}/eval.blocks.tsv")


# ------------------------------------------------------------------------------------------------ part 1: fixes
def interval_rows():
    """Per evaluated degeneracy interval (M11-D's ARTEFACT_INTERVAL; C4 has none): summed block distance, the
    concatenated-interval distance from the saved alignments, and the artefact bp once per interval."""
    blocks = T(f"{M11D}/runs/outcomes.a2/blocks.tsv")
    out = []
    for c in m11d.CELLS:
        if c["locus"] == "C4":
            continue
        cid = c["id"]
        iv = m11d.ARTEFACT_INTERVAL[c["locus"]]
        ev = eval_rows(cid)
        order = [e["block_id"] for e in ev]
        blk = order[order.index(iv[0]):order.index(iv[-1]) + 1]
        D = {r["id"]: (int(r["edits"]), int(r["aln_len"])) for r in T(f"{M11D}/runs/outcomes.a2/{cid}.distance.tsv")}
        if f"{cid}.iv0.p1t1" not in D:
            out.append(dict(cell=cid, locus=c["locus"], interval=f"{blk[0]}..{blk[-1]}", blocks=",".join(blk), evaluated=False))
            continue
        g = {t: D[f"{cid}.iv0.{t}"][0] for t in ("p1t1", "p2t2", "p1t2", "p2t1")}
        iv_d = min(g["p1t1"] + g["p2t2"], g["p1t2"] + g["p2t1"])
        e = {x["block_id"]: x for x in ev}
        blk_d = sum(float(e[b]["unordered_distance"]) for b in blk)
        flagged = [r for r in blocks if r["cell"] == cid and r["block_id"] in blk and r["synonymous_artefact"] == "True"]
        out.append(dict(cell=cid, locus=c["locus"], interval=f"{blk[0]}..{blk[-1]}", blocks=",".join(blk), evaluated=True,
                        summed_block_distance=blk_d, concatenated_distance=iv_d, artefact_bp=max(0.0, blk_d - iv_d),
                        summed_block_excess=sum(float(e[b]["excess"]) for b in blk), flagged_blocks=len(flagged),
                        report_a4_counted_bp=sum(float(r["synonymous_artefact_bp"]) for r in flagged)))
    return out


def fixes(out):
    ivs = interval_rows()
    write(out, "artefact_intervals", ivs)
    rep = {r["locus"]: r for r in T(f"{M11D}/runs/report.a4/blocks_summary.tsv")}
    corr = []
    for L in ("ALL", "C4", "CYP2D6", "ACOT"):
        S = [r for r in ivs if r.get("evaluated") and (L == "ALL" or r["locus"] == L) and r["artefact_bp"] > 0]
        corr.append(dict(locus=L, report_a4_artefact_bp=rep[L]["artefact_bp"], report_a4_artefact_blocks=rep[L]["artefact_blocks"],
                         flagged_intervals=len(S), corrected_artefact_bp=sum(r["artefact_bp"] for r in S),
                         recount_of_report_a4=sum(r["report_a4_counted_bp"] for r in S),
                         note="report.a4 summed the interval value over each flagged block; corrected = once per interval"))
    write(out, "artefact_correction", corr)

    # calibration-eligibility annotations: labels are copied unchanged; nothing is refitted or relabelled
    by_iv = {}
    for r in ivs:
        for b in r["blocks"].split(","):
            by_iv[(r["cell"], b)] = r
    rows = []
    for r in T(f"{M11D}/runs/outcomes.a2/blocks.tsv"):
        e = next(x for x in eval_rows(r["cell"]) if x["block_id"] == r["block_id"])
        iv = by_iv.get((r["cell"], r["block_id"]))
        all_empty = r["predicted"] == "True" and all(e[k] == "0" for k in ("truth1_bp", "truth2_bp", "called1_bp", "called2_bp"))
        rep_dep = bool(iv and iv.get("evaluated") and iv["artefact_bp"] > 0)
        if r["truth_defined"] != "True":
            status, informative = "UNDEFINED_TRUTH", "no"
        elif r["predicted"] != "True":
            status, informative = "NO_PREDICTION", "no"
        elif all_empty:
            status, informative = "ALL_EMPTY_COUNTED", "no"
        elif rep_dep and iv["concatenated_distance"] == 0:
            # existing interval evidence: the predicted interval sequence equals the truth interval exactly
            status, informative = "REPRESENTATION_DEPENDENT_INTERVAL_EXACT", "interval_level_only"
        elif rep_dep:
            # the interval is not exact and its availability floor was never computed (no new floor search)
            status, informative = "REPRESENTATION_DEPENDENT_UNRESOLVED", "no"
        else:
            status, informative = "BLOCK_OUTCOME", "yes" if r["assessable"] == "True" else "no"
        rows.append(dict(cell=r["cell"], locus=r["locus"], family=r["family"], block_id=r["block_id"], kind=r["kind"],
                         eligibility=status, fit_informative=informative,
                         interval=iv["interval"] if iv else "none",
                         interval_concatenated_distance=iv.get("concatenated_distance", "NA") if iv else "NA",
                         interval_artefact_bp=iv.get("artefact_bp", "NA") if iv else "NA",
                         excess=r.get("excess", ""), floor_attained=r.get("floor_attained", ""),
                         top_category=r.get("top_category", ""), catastrophic=r.get("catastrophic", ""), relabelled="no"))
    write(out, "eligibility", rows)
    summ = []
    for L in ("ALL", "C4", "CYP2D6", "ACOT"):
        S = [r for r in rows if L == "ALL" or r["locus"] == L]
        for st in sorted({r["eligibility"] for r in rows}):
            V = [r for r in S if r["eligibility"] == st]
            summ.append(dict(locus=L, eligibility=st, blocks=len(V), non_top=sum(r["top_category"] == "False" for r in V),
                             catastrophic=sum(r["catastrophic"] == "True" for r in V)))
    write(out, "eligibility_summary", summ)

    # QVfrac per homologue for the catastrophic blocks (descriptive; QVfrac is not a probability and not a cutoff)
    qrows = []
    for r in T(f"{M11D}/runs/outcomes.a2/blocks.tsv"):
        if r.get("catastrophic") != "True":
            continue
        D = {x["id"]: (int(x["edits"]), int(x["aln_len"])) for x in T(f"{M11D}/runs/outcomes.a2/{r['cell']}.distance.tsv")}
        pid = f"{r['cell']}.{r['block_id']}".replace(":", "_")
        g = {t: D[f"{pid}.{t}"] for t in ("p1t1", "p2t2", "p1t2", "p2t1", "b1t1", "b2t2")}
        sel = ("p1t1", "p2t2") if g["p1t1"][0] + g["p2t2"][0] <= g["p1t2"][0] + g["p2t1"][0] else ("p1t2", "p2t1")
        q = {t: qv(*g[t]) for t in g}
        on_t1 = next(t for t in sel if t.endswith("t1"))   # the predicted homologue paired with truth1
        on_t2 = next(t for t in sel if t.endswith("t2"))
        qrows.append(dict(cell=r["cell"], block_id=r["block_id"], excess=r["excess"], QVfrac=r["QVfrac"], top_category=r["top_category"],
                          pairing=",".join(sel),
                          truth1_pred_edits_len=f"{g[on_t1][0]}/{g[on_t1][1]}", QV_pred_t1=round(q[on_t1], 3), QV_best_t1=round(q["b1t1"], 3),
                          truth2_pred_edits_len=f"{g[on_t2][0]}/{g[on_t2][1]}", QV_pred_t2=round(q[on_t2], 3), QV_best_t2=round(q["b2t2"], 3),
                          QV_pred_sum=round(q[on_t1] + q[on_t2], 3), QV_best_sum=round(q["b1t1"] + q["b2t2"], 3),
                          share_of_QV_pred_sum_from_better_homologue=round(max(q[on_t1], q[on_t2]) / (q[on_t1] + q[on_t2]), 3)))
    write(out, "qvfrac_catastrophic", qrows)
    print(open(f"{out}/artefact_correction.tsv").read())
    print(open(f"{out}/qvfrac_catastrophic.tsv").read())
    bad = [c for c in corr if abs(float(c["recount_of_report_a4"]) - float(c["report_a4_artefact_bp"])) > 1e-9]
    print(f"eligibility rows {len(rows)}; report.a4 recount mismatches {len(bad)}")
    return 0 if not bad else 1


# ------------------------------------------------------------------------------------------------ part 2: prepare
def spell(cell, block_ids, path):
    return [fasta(f"{M11D}/runs/truth.{cell}/bs.{bid}.alleles.fa")[f"allele_{a}"] for bid, a in zip(block_ids, path)]


def prepare(out):
    rows, failures = [], 0
    for c in CASES:
        cell, bid = c["cell"], c["block"]
        call = f"{M11D}/runs/call.{cell}/call"
        ev = eval_rows(cell)
        ids = [e["block_id"] for e in ev]
        b = ids.index(bid)
        e = ev[b]
        win = next(r for r in T(call + ".equivalence_set.tsv") if r["rank"] == "1" and r["representative"] == "true")
        h1, h2 = vec(win["hap1_alleles"]), vec(win["hap2_alleles"])
        called = sorted([int(e["called_allele1"]), int(e["called_allele2"])])
        b1, b2 = int(e["truth1_closest_allele"]), int(e["truth2_closest_allele"])
        cands = {"incumbent": (h1, h2)}
        a1, a2 = list(h1), list(h2)
        a1[b], a2[b] = b1, b2
        cands["assign_A"] = (a1, a2)          # hap1 <- best for truth1, hap2 <- best for truth2
        c1, c2 = list(h1), list(h2)
        c1[b], c2[b] = b2, b1
        cands["assign_B"] = (c1, c2)          # the other homologue assignment
        seg = {k: (spell(cell, ids, v[0]), spell(cell, ids, v[1])) for k, v in cands.items()}
        seq = {k: ("".join(s[0]), "".join(s[1])) for k, s in seg.items()}
        sh = {r["haplotype"]: r for r in T(call + ".sequence_haplotypes.tsv")}
        by_path = {sh[h]["block_path"]: sh[h]["sha256"] for h in sh}
        checks = dict(
            winner_matches_evaluator_called=sorted([h1[b], h2[b]]) == called,
            incumbent_hap1_sha256=hashlib.sha256(seq["incumbent"][0].encode()).hexdigest() == by_path.get(win["hap1_alleles"]),
            incumbent_hap2_sha256=hashlib.sha256(seq["incumbent"][1].encode()).hexdigest() == by_path.get(win["hap2_alleles"]))
        for k in ("assign_A", "assign_B"):
            for s in (0, 1):
                inc, alt = seg["incumbent"][s], seg[k][s]
                checks[f"{k}_hap{s + 1}_flanks_unchanged"] = inc[:b] == alt[:b] and inc[b + 1:] == alt[b + 1:]
                checks[f"{k}_hap{s + 1}_block_is_declared_allele"] = alt[b] == fasta(f"{M11D}/runs/truth.{cell}/bs.{bid}.alleles.fa")[
                    f"allele_{cands[k][s][b]}"]
        ok = all(checks.values())
        failures += 0 if ok else 1
        with open(f"{out}/{c['case']}.candidates.tsv", "w") as f:
            f.write("# label\thap1\thap2\texpected_total\n")
            for k, (x, y) in cands.items():
                f.write(f"{k}\t{','.join(map(str, x))}\t{','.join(map(str, y))}\t{win['score'] if k == 'incumbent' else 'NA'}\n")
        for k in cands:
            for s in (0, 1):
                rows.append(dict(case=c["case"], cell=cell, block_id=bid, block_index=b, label=k, haplotype=s + 1,
                                 block_allele=cands[k][s][b], block_bp=len(seg[k][s][b]), haplotype_bp=len(seq[k][s]),
                                 fnv64=fnv64(seq[k][s]), sha256=hashlib.sha256(seq[k][s].encode()).hexdigest()))
        write(out, f"{c['case']}.checks", [dict(case=c["case"], check=k, ok=v) for k, v in checks.items()])
        print(c["case"], cell, bid, "called", called, "best", (b1, b2), "verified" if ok else f"FAILED {[k for k, v in checks.items() if not v]}")
    write(out, "spelled", rows)
    return 0 if failures == 0 else 1


# ------------------------------------------------------------------------------------------------ part 2: trace + classify
def classify(rep_ok, bytes_ok, consistency, d_best, decomp_ok, evaluated_better):
    """Frozen order (DESIGN.md §4)."""
    if not (rep_ok and bytes_ok):
        return "UNRESOLVED_ATTRIBUTION", "harness: incumbent not reproduced or scored bytes not verified"
    if consistency == "INCONSISTENT" or decomp_ok is False:
        return "NUMERICAL_OR_IMPLEMENTATION_INCONSISTENCY", "recorded run score or origin re-summation disagrees with score-set"
    if evaluated_better:
        return "NUMERICAL_OR_IMPLEMENTATION_INCONSISTENCY", "the run scored a higher-total alternative yet reported the incumbent"
    if d_best > TOL_TOTAL:
        return "SEARCH_MISSED_BETTER_ALTERNATIVE", "a single-block replacement scores higher in the winner's own context"
    if d_best < -TOL_TOTAL:
        return "OBJECTIVE_PREFERS_FARTHER_PREDICTION", "locally: both best-available assignments score lower than the winner"
    return "UNRESOLVED_ATTRIBUTION", "tie within tolerance"


def trace(out):
    rows = []
    for c in CASES:
        cell, bid, cid = c["cell"], c["block"], c["case"]
        call = f"{M11D}/runs/call.{cell}/call"
        ss = f"{OUT}/runs/score.{cid}/ss"
        if not os.path.exists(ss + ".scores_truth_free.tsv"):
            rows.append(dict(case=cid, cell=cell, block_id=bid, role=c["role"], classification="UNRESOLVED_ATTRIBUTION",
                             reason="no score-set output"))
            continue
        sc = {r["label"]: r for r in T(ss + ".scores_truth_free.tsv")}
        hp = {r["unique"]: r for r in T(ss + ".haplotypes.tsv")}
        sp = {(r["label"], r["haplotype"]): r for r in T(f"{OUT}/runs/prepare/spelled.tsv") if r["case"] == cid}
        bytes_ok = all(hp[sc[l][f"hap{s}_unique"]]["sequence_fnv64"] == sp[(l, str(s))]["fnv64"]
                       and hp[sc[l][f"hap{s}_unique"]]["bp"] == sp[(l, str(s))]["haplotype_bp"] for l in sc for s in (1, 2))
        inc = sc["incumbent"]
        rep_ok = inc["expected_total"] != "NA" and abs(float(inc["total"]) - float(inc["expected_total"])) <= TOL_TOTAL
        ev = next(x for x in eval_rows(cell) if x["block_id"] == bid)
        b = int(sp[("incumbent", "1")]["block_index"])
        best = sorted([int(ev["truth1_closest_allele"]), int(ev["truth2_closest_allele"])])
        # run record: every union / repair candidate the run scored, keyed by unordered allele-vector pair
        cand = T(call + ".mosaic_candidates.tsv")
        score_of = {r["sequence_candidate_id"]: r["score"] for r in T(call + ".mosaic_scores.tsv")}
        key = lambda x, y: tuple(sorted([x, y]))  # noqa: E731
        rec = {}
        for r in cand:
            rec.setdefault(key(r["hap1_alleles"], r["hap2_alleles"]), r)
        prop = [r for r in cand if sorted([vec(r["hap1_alleles"])[b], vec(r["hap2_alleles"])[b]]) == best]
        lr = {(int(r["allele_x"]), int(r["allele_y"])): r for r in T(ss + ".local_ranks.tsv") if r["block_id"] == bid}
        rank = lr.get(tuple(best), {})
        rounds = T(call + ".repair_rounds.tsv")
        row = dict(case=cid, cell=cell, block_id=bid, role=c["role"], reproduced=rep_ok, scored_bytes_verified=bytes_ok,
                   incumbent_total=inc["total"], called_pair=f"{ev['called_allele1']},{ev['called_allele2']}", best_pair=",".join(map(str, best)),
                   best_pair_local_rank=rank.get("rank", "NA"), best_pair_in_local_list=rank.get("in_local_list", "NA"),
                   best_pair_emission_below_top=rank.get("below_top", "NA"),
                   best_pair_proposed_any_context=len(prop),
                   best_pair_proposed_by=";".join(f"{k}={v}" for k, v in sorted({s: sum(r["scored_in"] == s for r in prop) for s in {r["scored_in"] for r in prop}}.items())) or "none",
                   repair_rounds=len(rounds), repair_last_round_status=rounds[-1]["status"] if rounds else "NA",
                   repair_last_round_moved=rounds[-1]["moved"] if rounds else "NA",
                   repair_blocks_moved=";".join(r["repaired_block"] for r in rounds if r["moved"] == "true") or "none")
        consistency, d_best, decomp_ok, evaluated_better = "NOT_RECORDED", -math.inf, True, False
        for l in ("assign_A", "assign_B"):
            a = sc[l]
            for t in ("fragment_term", "exposure_term", "prior", "total"):
                row[f"{l}_d_{t}"] = float(a[t]) - float(inc[t]) if a[t] not in ("NA", "-inf") else a[t]
            identical = a["sequence_pair"] == inc["sequence_pair"]
            row[f"{l}_identical_to_incumbent"] = identical
            r = rec.get(key(a["hap1_alleles"], a["hap2_alleles"]))
            row[f"{l}_evaluated_by_run"] = r["scored_in"] if r else "not_in_run_record"
            if r and r["sequence_candidate_id"] in score_of:
                diff = abs(float(score_of[r["sequence_candidate_id"]]) - float(a["total"]))
                row[f"{l}_run_recorded_total_abs_diff"] = diff
                consistency = "INCONSISTENT" if diff > TOL_PRINTED else ("CONSISTENT" if consistency != "INCONSISTENT" else consistency)
            if not identical:
                d_best = max(d_best, row[f"{l}_d_total"])
                evaluated_better = evaluated_better or (r is not None and row[f"{l}_d_total"] > TOL_TOTAL)
            dp = f"{ss}.decomp.{l}__incumbent.summary.tsv"
            if os.path.exists(dp) and not identical:
                s = KV(dp)
                ok = s.get("placement_evidence_interpretable") == "true" and all(
                    abs(float(s[k])) <= 1e-6 for k in ("fragment_term_A_from_canonical_minus_scorer", "fragment_term_B_from_canonical_minus_scorer"))
                decomp_ok = decomp_ok and ok
                row[f"{l}_decomp_reconciled"] = ok
                for k in ("fragments_favouring_A", "fragments_favouring_B", "top10_fragments_delta", "top50_fragments_delta"):
                    row[f"{l}_{k.replace('_A', '_alt').replace('_B', '_incumbent')}"] = s.get(k, "NA")
                for cls in ("block_on_A_only", "block_on_B_only", "block_interior_on_both", "block_junction_on_both", "outside_block_on_all"):
                    row[f"{l}_class_{cls.replace('_A_', '_alt_').replace('_B_', '_incumbent_')}"] = \
                        f"{s.get(f'class_{cls}_fragments', '0')}:{s.get(f'class_{cls}_delta', '0')}"
        if d_best == -math.inf:
            d_best = 0.0   # both replacements identical to the incumbent (a control at its best pair)
        row["consistency_with_run_record"] = consistency
        row["best_alternative_d_total"] = d_best
        db = next((r for r in T(call + ".dosage_blocks.tsv") if r["block_id"] == bid), {})
        for k in ("census_route", "dosage_status", "representation_status", "allocation_status", "original_content", "corrected_content",
                  "repair_separation", "repaired_prior_spellable", "repair_class"):
            row[f"dosage_{k}"] = db.get(k, "NA")
        row["classification"], row["reason"] = classify(rep_ok, bytes_ok, consistency, d_best, decomp_ok, evaluated_better)
        if c["role"] == "control":
            row["classification"] = "CONTROL:" + row["classification"]
        rows.append(row)
        print(cid, cell, bid, row["classification"], f"d_best {d_best:+.4f}", "reproduced" if rep_ok else "NOT REPRODUCED",
              "bytes ok" if bytes_ok else "BYTES DIFFER")
    write(out, "attribution", rows)
    return 0


if __name__ == "__main__":
    mode, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"fixes": fixes, "prepare": prepare, "trace": trace}[mode](out))
