#!/usr/bin/env python3
"""M12-B: bounded, offline dosage-first selection experiment (an experimental decision rule; not a posterior, not a
production dosage repair, not a confidence model, no prior change).

Rule (DESIGN.md §2, frozen): at each M11-D repeated_unit block, in the reported winner's fixed context, with existing active-
model alleles only:
  1. copy state of an allele = its copies of every pooled copy-count feature of the block (panvar::dosage::pooled_features,
     the dosage audit's own graph-only census; default settings); state of an ordered pair = the per-feature sum
  2. count loss of every ordered pair = the unchanged M6-A/count-attrib loss (real reads, stable-copy depth, 31-mer dictionary)
  3. per state: its best count loss; select the state(s) with the lowest best loss (exact ties at the printed precision kept)
  4. within the selected state(s): every ordered pair (both homologue assignments) scored by the UNCHANGED score-set objective
     (fragment + exposure + copying prior); the highest total is the selection (exact ties reported, never broken)
Modes: selftest | gate_features <out> | select <out> | evaluate <out>
"""
import collections
import functools
import hashlib
import json
import math
import os
import re
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
import m11e  # noqa: E402

m11d = m11e.m11d
T, KV, write = m11d.T, m11d.KV, m11d.write
fasta = functools.lru_cache(maxsize=None)(m11d.fasta)
M11D = m11e.M11D
OUT = f"{m11d.R}/m12b_dosage_first_20261021"
WT = m11d.WT
SAFETY = 1.5                 # preflight estimate multiplier
BLOCK_CAP_S = 240.0          # a block whose estimated exhaustive scoring exceeds this is a LIMITATION (never shortlisted)
TIE_TOL = 1e-9               # exact-tie tolerance on totals (nats)
LARGE = 1000.0               # catastrophic excess (bp)
NEAR_FLOOR = 10.0            # near-floor control (bp)


def blocks_manifest():
    """The 29 M11-D repeated_unit blocks (predicted, truth-defined), in M11-D cell order."""
    rows = [r for r in T(f"{M11D}/runs/outcomes.a2/blocks.tsv") if r["census_route"] == "repeated_unit" and r["predicted"] == "True"]
    order = {c["id"]: i for i, c in enumerate(m11d.CELLS)}
    rows.sort(key=lambda r: order[r["cell"]])
    return [dict(bid=f"B{i + 1:02d}", cell=r["cell"], locus=r["locus"], block=r["block_id"], orig_excess=float(r["excess"])) for i, r in enumerate(rows)]


def winner(cell):
    return next(r for r in T(f"{M11D}/runs/call.{cell}/call.equivalence_set.tsv") if r["rank"] == "1" and r["representative"] == "true")


# ------------------------------------------------------------------------------------------------ the frozen selection
def select_state(features, losses):
    """features: list of per-allele copy lists (one per pooled feature); losses: {(a, b): loss_string}.
    Returns dict(status, states, selected, ...). Exact ties = equal printed loss."""
    if not features:
        return dict(status="UNSUPPORTED_NO_POOLED_FEATURES")
    if not losses:
        return dict(status="UNSUPPORTED_NO_COUNT_LOSS")
    state = lambda a, b: tuple(f[a] + f[b] for f in features)  # noqa: E731
    best = {}
    for (a, b), l in losses.items():
        s = state(a, b)
        v = float(l)
        if s not in best or v < best[s][0]:
            best[s] = (v, l)
    lowest = min(v for v, _ in best.values())
    lowest_str = {l for v, l in best.values() if v == lowest}
    selected = sorted(s for s, (v, l) in best.items() if l in lowest_str)
    ranked = sorted(best.items(), key=lambda kv: kv[1][0])
    runner = next((v for s, (v, l) in ranked if s not in selected), None)
    in_sel = sorted(float(l) for (a, b), l in losses.items() if state(a, b) in selected)
    return dict(status="SUPPORTED", selected=selected, best={str(s): v for s, (v, _) in best.items()}, lowest=lowest,
                runner_up_state_loss=runner, state_margin=(runner - lowest) if runner is not None else None,
                within_selected_min=in_sel[0], within_selected_median=statistics.median(in_sel),
                content_sensitive=(runner is not None and (runner - lowest) < statistics.median(in_sel) - in_sel[0]),
                candidates=sorted((a, b) for (a, b) in losses if state(a, b) in selected))


def selftest(out):
    checks = []
    f = [[1, 2, 3, 2]]                       # one feature; alleles 0..3 carry 1,2,3,2 copies
    L = {(a, b): f"{10 + abs(f[0][a] + f[0][b] - 4) + 0.1 * a:.6g}" for a in range(4) for b in range(4)}
    r = select_state(f, L)
    checks.append(("total-state grouping: 3+1 and 2+2 share state (4,)", r["selected"] == [(4,)] and (2, 0) in r["candidates"] and (1, 1) in r["candidates"]
                   and (1, 3) in r["candidates"]))
    checks.append(("both homologue assignments are candidates", (0, 2) in r["candidates"] and (2, 0) in r["candidates"]))
    L2 = dict(L)
    best4 = min((L[(a, b)] for a in range(4) for b in range(4) if f[0][a] + f[0][b] == 4), key=float)
    L2[(0, 0)] = best4                       # state (2,) now ties the best loss of state (4,) exactly
    r2 = select_state(f, L2)
    checks.append(("exact ties across states are retained", r2["selected"] == [(2,), (4,)]))
    shuffled = dict(reversed(list(L.items())))
    checks.append(("selection is invariant to candidate order", select_state(f, shuffled) == r))
    checks.append(("no features -> unsupported", select_state([], L)["status"] == "UNSUPPORTED_NO_POOLED_FEATURES"))
    checks.append(("no count loss -> unsupported", select_state(f, {})["status"] == "UNSUPPORTED_NO_COUNT_LOSS"))
    f3 = [[1, 2, 3, 2], [0, 1, 0, 1]]        # a second feature separates alleles 1 and 3
    r3 = select_state(f3, L)
    checks.append(("multi-feature state is the per-feature sum", all(len(s) == 2 for s in r3["selected"])))
    for name, ok in checks:
        print(("PASS " if ok else "FAIL ") + name)
    return 0 if all(ok for _, ok in checks) else 1


# ------------------------------------------------------------------------------------------------ gates on real features
def gate_features(out):
    """The probe's pooled features equal the caller's own dosage_feature_copies.tsv wherever the audit wrote one; and on real
    C4 alleles, (3+1) and (2+2) allocations share one total state (total dosage cannot resolve allocation)."""
    rows, ok_all = [], True
    for b in blocks_manifest():
        fp = f"{OUT}/runs/count.{b['bid']}/ca.features.tsv"
        caller = f"{M11D}/runs/call.{b['cell']}/call.dosage_feature_copies.tsv"
        if not os.path.exists(fp):
            rows.append(dict(bid=b["bid"], cell=b["cell"], probe_features="missing", caller="-", match="-"))
            continue
        probe = {int(r["feature"]): r["copies_by_allele"] for r in T(fp)}
        if os.path.exists(caller):
            cc = {int(r["pooled_index"]): r["copies_by_allele"] for r in T(caller) if r["block_id"] == b["block"]}
            match = all(probe.get(k) == v for k, v in cc.items()) and bool(cc)
            ok_all = ok_all and match
            rows.append(dict(bid=b["bid"], cell=b["cell"], probe_features=len(probe), caller=len(cc), match=match))
        else:
            rows.append(dict(bid=b["bid"], cell=b["cell"], probe_features=len(probe), caller="not_written(audit_failed_or_absent)", match="-"))
    demo = []
    for b in blocks_manifest():
        if b["locus"] != "C4" or not os.path.exists(f"{OUT}/runs/count.{b['bid']}/ca.features.tsv"):
            continue
        feats = [list(map(int, r["copies_by_allele"].split(","))) for r in T(f"{OUT}/runs/count.{b['bid']}/ca.features.tsv")]
        c0 = feats[0]
        a31 = next(((a, x) for a in range(len(c0)) for x in range(len(c0)) if c0[a] == 3 and c0[x] == 1), None)
        a22 = next(((a, x) for a in range(len(c0)) for x in range(len(c0)) if c0[a] == 2 and c0[x] == 2), None)
        if a31 and a22:
            demo.append(dict(bid=b["bid"], pair_3plus1=a31, pair_2plus2=a22, state_feature0=(c0[a31[0]] + c0[a31[1]], c0[a22[0]] + c0[a22[1]])))
            break
    write(out, "feature_gate", rows)
    write(out, "allocation_demo", demo)
    ok_demo = bool(demo) and demo[0]["state_feature0"][0] == demo[0]["state_feature0"][1]
    print(f"feature reproduction: {sum(r['match'] is True for r in rows)} blocks match the caller, "
          f"{sum(r['match'] is False for r in rows)} differ; 3+1 vs 2+2 same feature-0 total: {ok_demo}")
    return 0 if ok_all and ok_demo else 1


# ------------------------------------------------------------------------------------------------ select + preflight plan
def preflight_rates():
    """Per locus: seconds per unique haplotype and fixed overhead, from the preflight units' score-set logs."""
    rates = {}
    for locus in ("CYP2D6", "C4", "ACOT"):
        log = f"{OUT}/runs/preflight.{locus}/run.log"
        led = {r["unit"]: r for r in m12b_ledger()}
        if not os.path.exists(log) or led.get(f"preflight.{locus}", {}).get("exit") != "0":
            continue
        txt = open(log).read()
        m = re.search(r"(\d+) unique haplotypes, .*? ([0-9.]+)s scoring", txt)
        haps, scoring = int(m.group(1)), float(m.group(2))
        wall = float(led[f"preflight.{locus}"]["wall_seconds"])
        rates[locus] = dict(per_hap=scoring / haps, overhead=max(0.0, wall - scoring))
    return rates


def m12b_ledger():
    p = f"{OUT}/ledger.tsv"
    if not os.path.exists(p):
        return []
    lines = [l.rstrip("\n").split("\t") for l in open(p)]
    return [dict(zip(lines[0], r)) for r in lines[1:]]


def select(out):
    rates = preflight_rates()
    plan, rows = [], []
    for b in blocks_manifest():
        cp = f"{OUT}/runs/count.{b['bid']}/ca"
        feats = [list(map(int, r["copies_by_allele"].split(","))) for r in T(cp + ".features.tsv")] if os.path.exists(cp + ".features.tsv") else []
        losses = {}
        if os.path.exists(cp + ".catalogue.tsv"):
            for r in T(cp + ".catalogue.tsv"):
                if r.get("profile") == "real":
                    losses[(int(r["allele_h1"]), int(r["allele_h2"]))] = r["loss"]
        res = select_state(feats, losses)
        win = winner(b["cell"])
        h1, h2 = m11e.vec(win["hap1_alleles"]), m11e.vec(win["hap2_alleles"])
        ev = m11e.eval_rows(b["cell"])
        k = [e["block_id"] for e in ev].index(b["block"])
        row = dict(bid=b["bid"], cell=b["cell"], locus=b["locus"], block=b["block"], status=res["status"], features=len(feats),
                   alleles=len(feats[0]) if feats else "NA", count_pairs=len(losses))
        if res["status"] == "SUPPORTED":
            st = lambda a, c: tuple(f[a] + f[c] for f in feats)  # noqa: E731
            inc_state = st(h1[k], h2[k])
            cands = res["candidates"]
            haps = len({a for a, _ in cands}) + len({c for _, c in cands})
            r = rates.get(b["locus"])
            est = (r["overhead"] + SAFETY * haps * r["per_hap"]) if r else None
            admitted = est is not None and est <= BLOCK_CAP_S
            row.update(original_state=inc_state, original_state_best_loss=res["best"].get(str(inc_state)),
                       incumbent_count_loss=losses.get((h1[k], h2[k])), selected_states=";".join(map(str, res["selected"])),
                       selected_best_loss=res["lowest"], runner_up_state_loss=res["runner_up_state_loss"], state_margin=res["state_margin"],
                       within_selected_min=res["within_selected_min"], within_selected_median=res["within_selected_median"],
                       content_sensitive=res["content_sensitive"], states=len(res["best"]),
                       candidates_in_selected=len(cands), unique_haplotypes=haps,
                       estimate_s=round(est, 1) if est is not None else "NO_PREFLIGHT",
                       admitted=admitted, limitation="" if admitted else ("LIMITATION_NO_PREFLIGHT_RATE" if est is None else "LIMITATION_STATE_TOO_LARGE_TO_SCORE_EXHAUSTIVELY"))
            if admitted:
                with open(f"{out}/{b['bid']}.candidates.tsv", "w") as f:
                    f.write(f"incumbent\t{','.join(map(str, h1))}\t{','.join(map(str, h2))}\t{win['score']}\n")
                    for i, (a, c) in enumerate(cands):
                        x, y = list(h1), list(h2)
                        x[k], y[k] = a, c
                        f.write(f"c{i}_{a}_{c}\t{','.join(map(str, x))}\t{','.join(map(str, y))}\tNA\n")
                plan.append(dict(unit=f"score.{b['bid']}", bid=b["bid"], cell=b["cell"], block=b["block"], estimate_s=round(est, 1),
                                 timeout_s=int(math.ceil(1.5 * est + 60))))
        rows.append(row)
        print(b["bid"], b["cell"], b["block"], row["status"], row.get("selected_states", ""), "incumbent", row.get("original_state", ""),
              "cands", row.get("candidates_in_selected", ""), "haps", row.get("unique_haplotypes", ""), "est", row.get("estimate_s", ""),
              row.get("limitation", ""))
    write(out, "selection", rows)
    write(out, "score_plan", plan)
    return 0


# ------------------------------------------------------------------------------------------------ evaluate (truth opened here)
def canon_codes(seq, k=31):
    enc = {"A": 0, "C": 1, "G": 2, "T": 3}
    mask = (1 << (2 * k)) - 1
    fwd = rev = 0
    filled = 0
    out = collections.Counter()
    for ch in seq:
        b = enc.get(ch, -1)
        if b < 0:
            fwd = rev = filled = 0
            continue
        fwd = ((fwd << 2) | b) & mask
        rev = (rev >> 2) | ((3 - b) << (2 * (k - 1)))
        filled = min(k, filled + 1)
        if filled == k:
            out[min(fwd, rev)] += 1
    return out


def truth_state(bid, seq):
    """Per feature: the median copy count of the feature's k-mers in `seq`; assessable when >= 90% of them sit at that median."""
    vec, ok = [], True
    counts = canon_codes(seq)
    for r in T(f"{OUT}/runs/count.{bid}/ca.feature_kmers.tsv"):
        codes = [int(x) for x in r["kmer_codes"].split(",")]
        c = [counts.get(x, 0) for x in codes]
        med = statistics.median(c)
        vec.append(int(round(med)))
        ok = ok and med == int(med) and sum(v == med for v in c) >= 0.9 * len(c)
    return tuple(vec), ok


def distances(tag, alleles, t1, t2, out):
    os.makedirs(f"{out}/seq", exist_ok=True)
    man = f"{out}/{tag}.manifest.tsv"
    seqs = {"a": alleles[0], "b": alleles[1], "t1": t1, "t2": t2}
    with open(man, "w") as f:
        for x, y in (("a", "t1"), ("b", "t2"), ("a", "t2"), ("b", "t1")):
            for nm, s in ((x, seqs[x]), (y, seqs[y])):
                open(f"{out}/seq/{tag}.{nm}.txt", "w").write(s)
            f.write(f"{tag}.{x}{y}\t{out}/seq/{tag}.{x}.txt\t{out}/seq/{tag}.{y}.txt\t-\n")
    subprocess.run([f"{WT}/build/m10c_align_eval", "--distance", man, f"{out}/{tag}.distance.tsv"], check=True, capture_output=True)
    d = {r["id"].split(".")[-1]: int(r["edits"]) for r in T(f"{out}/{tag}.distance.tsv")}
    return min(d["at1"] + d["bt2"], d["at2"] + d["bt1"])


def evaluate(out):
    sel = {r["bid"]: r for r in T(f"{OUT}/runs/select/selection.tsv")}
    led = {r["unit"]: r for r in m12b_ledger()}
    rows = []
    for b in blocks_manifest():
        s = sel[b["bid"]]
        cell, blk, bid = b["cell"], b["block"], b["bid"]
        e = next(x for x in m11e.eval_rows(cell) if x["block_id"] == blk)
        al = fasta(f"{M11D}/runs/truth.{cell}/bs.{blk}.alleles.fa")
        tr = fasta(f"{M11D}/runs/truth.{cell}/bs.{blk}.truth.fa")
        t1, t2 = tr["truth1"], tr["truth2"]
        floor = float(e["availability_floor"])
        row = dict(bid=bid, cell=cell, locus=b["locus"], block=blk, status=s["status"], limitation=s.get("limitation", ""),
                   original_alleles=f"{e['called_allele1']},{e['called_allele2']}", original_excess=float(e["excess"]), floor=floor,
                   original_state=s.get("original_state", "NA"), selected_states=s.get("selected_states", "NA"),
                   incumbent_count_loss=s.get("incumbent_count_loss", "NA"), selected_best_count_loss=s.get("selected_best_loss", "NA"),
                   state_margin=s.get("state_margin", "NA"), content_sensitive=s.get("content_sensitive", "NA"),
                   candidates_in_selected=s.get("candidates_in_selected", "NA"))
        fp = f"{OUT}/runs/count.{bid}/ca.feature_kmers.tsv"
        if os.path.exists(fp):
            ts1, ok1 = truth_state(bid, t1)
            ts2, ok2 = truth_state(bid, t2)
            row.update(truth_state_hap=f"{ts1}|{ts2}", truth_total_state=str(tuple(x + y for x, y in zip(ts1, ts2))), truth_state_assessable=ok1 and ok2)
        ss = f"{OUT}/runs/score.{bid}/ss.scores_truth_free.tsv"
        if not os.path.exists(ss):
            row["outcome"] = "NOT_SCORED:" + (s.get("limitation") or s["status"] or led.get(f"score.{bid}", {}).get("note", ""))
            rows.append(row)
            continue
        sc = T(ss)
        inc = next(r for r in sc if r["label"] == "incumbent")
        reproduced = abs(float(inc["total"]) - float(inc["expected_total"])) <= 1e-6
        cand = [r for r in sc if r["label"] != "incumbent"]
        top = max(float(r["total"]) for r in cand)
        tied = [r for r in cand if float(r["total"]) >= top - TIE_TOL]
        # the selection(s), each evaluated; distances recomputed with the evaluator-matching aligner
        evals = []
        for r in tied:
            _, a, c = r["label"].split("_")
            a, c = int(a), int(c)
            d = distances(f"{bid}.{a}_{c}", (al[f"allele_{a}"], al[f"allele_{c}"]), t1, t2, out)
            ex = d - floor
            den = max(len(t1) + len(t2), len(al[f"allele_{a}"]) + len(al[f"allele_{c}"]))
            evals.append(dict(pair=(a, c), distance=d, excess=ex, norm=ex / den if den else None, total=float(r["total"]),
                              fragment=float(r["fragment_term"]), exposure=float(r["exposure_term"]), prior=float(r["prior"])))
        inc_d = distances(f"{bid}.incumbent", (al[f"allele_{e['called_allele1']}"], al[f"allele_{e['called_allele2']}"]), t1, t2, out)
        feats = [list(map(int, x["copies_by_allele"].split(","))) for x in T(f"{OUT}/runs/count.{bid}/ca.features.tsv")]
        stt = lambda a, c: tuple(f[a] + f[c] for f in feats)  # noqa: E731
        sel0 = evals[0]
        row.update(reproduced=reproduced, incumbent_distance_matches_evaluator=inc_d == float(e["unordered_distance"]),
                   scored_candidates=len(cand), coverage=f"{len(cand)}/{s['candidates_in_selected']}", tie_size=len(tied),
                   selected_alleles=";".join(f"{x['pair'][0]},{x['pair'][1]}" for x in evals),
                   selected_state=str(stt(*sel0["pair"])), selected_distance=";".join(str(x["distance"]) for x in evals),
                   selected_excess=";".join(f"{x['excess']:g}" for x in evals), selected_excess_max=max(x["excess"] for x in evals),
                   selected_norm_excess=";".join(f"{x['norm']:.4g}" for x in evals),
                   d_total_vs_incumbent=sel0["total"] - float(inc["total"]), d_fragment=sel0["fragment"] - float(inc["fragment_term"]),
                   d_exposure=sel0["exposure"] - float(inc["exposure_term"]), d_prior=sel0["prior"] - float(inc["prior"]),
                   seconds=led.get(f"score.{bid}", {}).get("wall_seconds"), peak_gib=led.get(f"score.{bid}", {}).get("peak_gib"))
        if row.get("truth_state_assessable"):
            ts = [tuple(map(int, re.findall(r"\d+", x))) for x in row["truth_state_hap"].split("|")]
            per = [tuple(f[a] for f in feats) for a in sel0["pair"]]
            row["total_copy_recovered"] = row["selected_state"] == row["truth_total_state"]
            row["original_total_copy_correct"] = row["original_state"] == row["truth_total_state"]
            row["allocation_matches_truth"] = sorted(per) == sorted(ts)
        # gates on this block's scoring: bytes (only the declared block changes), shared candidates, aliases
        ids = [x["block_id"] for x in m11e.eval_rows(cell)]
        k = ids.index(blk)
        win = winner(cell)
        h1, h2 = m11e.vec(win["hap1_alleles"]), m11e.vec(win["hap2_alleles"])
        seg = (m11e.spell(cell, ids, h1), m11e.spell(cell, ids, h2))
        flank_bp = [sum(len(x) for x in seg[i]) - len(seg[i][k]) for i in (0, 1)]
        hp = {r["unique"]: r for r in T(f"{OUT}/runs/score.{bid}/ss.haplotypes.tsv")}
        length_ok, fnv_checked, fnv_ok = True, 0, True
        for j, r in enumerate(sc):
            pa = (h1[k], h2[k]) if r["label"] == "incumbent" else tuple(map(int, r["label"].split("_")[1:]))
            for side in (0, 1):
                want = flank_bp[side] + len(al[f"allele_{pa[side]}"])
                length_ok = length_ok and int(hp[r[f"hap{side + 1}_unique"]]["bp"]) == want
                if j in (0, 1, len(sc) - 1):
                    full = "".join(seg[side][:k]) + al[f"allele_{pa[side]}"] + "".join(seg[side][k + 1:])
                    fnv_ok = fnv_ok and hp[r[f"hap{side + 1}_unique"]]["sequence_fnv64"] == m11e.fnv64(full)
                    fnv_checked += 1
        rec = {}
        for x in T(f"{M11D}/runs/call.{cell}/call.mosaic_candidates.tsv"):
            rec.setdefault(tuple(sorted([x["hap1_alleles"], x["hap2_alleles"]])), x["sequence_candidate_id"])
        rscore = {x["sequence_candidate_id"]: float(x["score"]) for x in T(f"{M11D}/runs/call.{cell}/call.mosaic_scores.tsv")}
        shared, shared_bad = 0, 0
        lines = {x.split("\t")[0]: x.split("\t") for x in open(f"{OUT}/runs/select/{bid}.candidates.tsv")}
        for r in sc:
            key = tuple(sorted([lines[r["label"]][1], lines[r["label"]][2]]))
            if key in rec and rec[key] in rscore:
                shared += 1
                shared_bad += abs(rscore[rec[key]] - float(r["total"])) > 1e-5
        by_seq = collections.defaultdict(list)
        for a_ in range(len(feats[0])):
            by_seq[al[f"allele_{a_}"]].append(a_)
        alias = {a_: min(v) for v in by_seq.values() for a_ in v}
        tot = collections.defaultdict(set)
        for r in cand:
            _, a_, c_ = r["label"].split("_")
            tot[(alias[int(a_)], alias[int(c_)])].add(round(float(r["total"]), 9))
        alias_bad = sum(len(v) > 1 for v in tot.values())
        row.update(gate_bytes_lengths=length_ok, gate_bytes_fnv=f"{fnv_ok}({fnv_checked})", gate_shared=f"{shared - shared_bad}/{shared}",
                   gate_alias_groups_inconsistent=alias_bad, alias_alleles=sum(len(v) - 1 for v in by_seq.values()))
        oe, se = row["original_excess"], row["selected_excess_max"]
        row["outcome"] = ("FIXED_LARGE" if oe >= LARGE and se < LARGE else "RETAINED_LARGE" if oe >= LARGE else
                          "INTRODUCED_LARGE" if se >= LARGE else "UNCHANGED" if se == oe else "IMPROVED" if se < oe else "WORSENED")
        row["near_floor_control"] = oe <= NEAR_FLOOR
        row["near_floor_damaged"] = oe <= NEAR_FLOOR and se > NEAR_FLOOR
        rows.append(row)
        print(bid, cell, blk, row["outcome"], "orig", oe, "sel", row["selected_excess"], "state", row["original_state"], "->", row["selected_state"],
              "truth", row.get("truth_total_state"), "tie", len(tied))
    write(out, "evaluation", rows)
    scored = [r for r in rows if "selected_excess_max" in r]
    cat = [r for r in rows if r["original_excess"] >= LARGE]
    cat_scored = [r for r in cat if "selected_excess_max" in r]
    fixed = sum(r["outcome"] == "FIXED_LARGE" for r in rows)
    introduced = sum(r["outcome"] == "INTRODUCED_LARGE" for r in rows)
    nf_damaged = sum(bool(r.get("near_floor_damaged")) for r in rows)
    if len(cat_scored) >= 3 and fixed >= 3 and introduced == 0 and nf_damaged == 0:
        decision = "MERITS_DEFAULT_OFF_CALLER_EXPERIMENT"
    elif len(cat_scored) >= 3 and (fixed == 0 or introduced >= fixed):
        decision = "FAILS"
    else:
        decision = "INCONCLUSIVE"
    summ = dict(blocks=len(rows), supported=sum(r["status"] == "SUPPORTED" for r in rows), scored=len(scored),
                catastrophic=len(cat), catastrophic_scored=len(cat_scored), fixed_large=fixed, retained_large=sum(r["outcome"] == "RETAINED_LARGE" for r in rows),
                introduced_large=introduced, near_floor_controls_scored=sum(bool(r.get("near_floor_control")) and "selected_excess_max" in r for r in rows),
                near_floor_damaged=nf_damaged, improved=sum(r["outcome"] == "IMPROVED" for r in rows), worsened=sum(r["outcome"] == "WORSENED" for r in rows),
                unchanged=sum(r["outcome"] == "UNCHANGED" for r in rows), decision=decision)
    write(out, "summary", [summ])
    print(json.dumps(summ))
    return 0


def preflight_candidates(out):
    """Timing preflight (no outcome is read): per locus, the manifest's first block, incumbent + the first 24 alleles of the block
    in the hap1 context (hap2 unchanged). Also the same CYP2D6 list in reverse order, for the order-invariance gate."""
    seen = set()
    for b in blocks_manifest():
        if b["locus"] in seen:
            continue
        seen.add(b["locus"])
        win = winner(b["cell"])
        h1, h2 = m11e.vec(win["hap1_alleles"]), m11e.vec(win["hap2_alleles"])
        k = [x["block_id"] for x in m11e.eval_rows(b["cell"])].index(b["block"])
        n = len(fasta(f"{M11D}/runs/truth.{b['cell']}/bs.{b['block']}.alleles.fa"))
        lines = [f"incumbent\t{','.join(map(str, h1))}\t{','.join(map(str, h2))}\t{win['score']}\n"]
        for a in range(min(24, n)):
            x = list(h1)
            x[k] = a
            lines.append(f"p{a}_{a}_{h2[k]}\t{','.join(map(str, x))}\t{','.join(map(str, h2))}\tNA\n")
        open(f"{out}/{b['locus']}.candidates.tsv", "w").writelines(lines)
        if b["locus"] == "CYP2D6":
            open(f"{out}/{b['locus']}.reversed.tsv", "w").writelines([lines[0]] + lines[1:][::-1])
        open(f"{out}/{b['locus']}.block", "w").write(json.dumps(b))
    return 0


def gate_order(out):
    a = {r["label"]: r for r in T(f"{OUT}/runs/preflight.CYP2D6/ss.scores_truth_free.tsv")}
    b = {r["label"]: r for r in T(f"{OUT}/runs/gate.order_reversed/ss.scores_truth_free.tsv")}
    same = a.keys() == b.keys() and all(a[k]["total"] == b[k]["total"] and a[k]["fragment_term_hex"] == b[k]["fragment_term_hex"] for k in a)
    print(f"order invariance: {len(a)} candidates, identical totals and fragment bits in reversed order: {same}")
    return 0 if same else 1


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "selftest":
        sys.exit(selftest(None))
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"gate_features": gate_features, "select": select, "evaluate": evaluate, "preflight_candidates": preflight_candidates,
              "gate_order": gate_order}[mode](out))
