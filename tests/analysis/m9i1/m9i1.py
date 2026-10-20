#!/usr/bin/env python3
"""M9-I.1 hardened compound-annotation bounds (contract v2; frozen with DESIGN.md). M9-I's code is reused unchanged (its
case list, model dumps and projector are loaded by executing m9i.py in a no-op mode) and its tables are not modified.

Contract v2 for a target T (unordered pair of concatenated interval sequences) over a pool of sequence groups:
  inputs      every group mass must be finite and >= 0, else INVALID (no bounds); the original total is recorded, masses are
              normalised by it (math.fsum), and the discrepancy |total - 1| is reported
  empty pool  (no group, or total mass 0): UNASSESSED, never zero error
  projections every recorded labelling of a group is projected; a group WITHOUT recorded labellings is MISSING: its mass enters
              BOTH support_upper and error_upper and neither lower bound (conservative policy), and is counted
  bounds      support_lower = groups whose projections all equal T; support_upper = groups with any projection equal to T, plus
              missing groups; error_lower = groups with labellings and no projection equal to T; error_upper = groups with any
              projection unequal to T, plus missing groups; every bound is a direct math.fsum of normalised masses (never
              1 - something)
  competitor  from projections, not mass: NO_EVALUATED_ALTERNATIVE (every recorded projection equals T and no group is
              missing), EVALUATED_ALTERNATIVE (some recorded projection differs and carries positive mass),
              EVALUATED_ALTERNATIVE_ZERO_MASS (differing projections exist but all their groups have mass exactly 0),
              UNKNOWN (no differing projection recorded but some group is missing). NO_EVALUATED_ALTERNATIVE is not evidence
              of correctness: it is a statement about the explored pool only.
  status      ASSESSED, ASSESSED_WITH_MISSING_PROJECTIONS, UNASSESSED, INVALID
  scope       saved-table masses are printed values (12 significant digits): the bounds are approximate pool-conditional masses
              reconstructed from them, not certified bounds on the original full-precision posterior

  m9i1.py gates <out>     fixtures (independent oracle for the sequence facts)
  m9i1.py reapply <out>   the 34 M9-I cases under contract v2, compared with M9-I's tables
"""
import csv
import math
import os
import sys
import types

mode, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))
R = "/Users/davide.bolognini/panvar/tests/results"
M9I = f"{R}/m9i_compound_targets_20261020/runs"

# M9-I's unchanged module (no-op mode: only definitions run)
m9i = types.ModuleType("m9i")
_argv = sys.argv
sys.argv = ["m9i.py", "__definitions_only__", os.path.join(out, ".m9i_noop")]
exec(compile(open(os.path.join(HERE, "..", "m9i", "m9i.py")).read(), "m9i.py", "exec"), m9i.__dict__)
sys.argv = _argv
Projector, T = m9i.Projector, m9i.T


def write(name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def bounds_v2(groups, proj, target):
    """groups: list of (mass, [(hap1 labels, hap2 labels), ...])."""
    masses = [p for p, _ in groups]
    if any((not isinstance(p, (int, float))) or math.isnan(p) or math.isinf(p) or p < 0 for p in masses):
        return dict(status="INVALID", reason="non-finite or negative mass")
    total = math.fsum(masses)
    if not groups or total == 0.0:
        return dict(status="UNASSESSED", reason="empty pool", original_total=total, groups=len(groups))
    sl, su, el, eu = [], [], [], []
    missing = missing_mass = 0
    disagree_groups = disagree_positive = ambiguous = 0
    for p, members in groups:
        q = p / total
        if not members:
            missing += 1
            missing_mass += q
            su.append(q)
            eu.append(q)
            continue
        projs = {proj.pair(h1, h2) for h1, h2 in members}
        ambiguous += len(projs) > 1
        hit, miss = target in projs, any(x != target for x in projs)
        if hit and not miss:
            sl.append(q)
        if hit:
            su.append(q)
        if not hit:
            el.append(q)
        if miss:
            eu.append(q)
            disagree_groups += 1
            disagree_positive += p > 0
    if disagree_groups and disagree_positive:
        competitor = "EVALUATED_ALTERNATIVE"
    elif disagree_groups:
        competitor = "EVALUATED_ALTERNATIVE_ZERO_MASS"
    elif missing:
        competitor = "UNKNOWN"
    else:
        competitor = "NO_EVALUATED_ALTERNATIVE"
    return dict(status="ASSESSED_WITH_MISSING_PROJECTIONS" if missing else "ASSESSED", original_total=total, normalisation_discrepancy=total - 1.0,
                groups=len(groups), missing_groups=missing, missing_mass=missing_mass, ambiguous_groups=ambiguous,
                support_lower=math.fsum(sl), support_upper=math.fsum(su), error_lower=math.fsum(el), error_upper=math.fsum(eu),
                competitor=competitor, disagreeing_groups=disagree_groups)


def spell(blocks, labels):
    return "".join(blocks[b][a] for b, a in enumerate(labels))


if mode == "gates":
    res = []
    ok = lambda n, c, d="": res.append(dict(gate=n, result="PASS" if c else "FAIL", detail=d))
    blocks = [["TT"], ["A", "AC"], ["CG", "G"], ["CCC"]]
    # M1 missing projections: conservative (both uppers), lowers unchanged, counted
    pj = Projector(blocks, 1, 2)
    tgt = pj.pair([0, 0, 0, 0], [0, 0, 0, 0])
    r = bounds_v2([(0.6, [([0, 0, 0, 0], [0, 0, 0, 0])]), (0.4, [])], pj, tgt)
    ok("M1 missing group: mass in both uppers, neither lower; counted", abs(r["support_lower"] - 0.6) < 1e-15 and abs(r["support_upper"] - 1.0) < 1e-15
       and r["error_lower"] == 0.0 and abs(r["error_upper"] - 0.4) < 1e-15 and r["missing_groups"] == 1 and r["competitor"] == "UNKNOWN"
       and r["status"] == "ASSESSED_WITH_MISSING_PROJECTIONS", str({k: r[k] for k in ("support_lower", "support_upper", "error_lower", "error_upper", "competitor")}))
    # M2 empty pool -> UNASSESSED
    ok("M2 empty pool is UNASSESSED, never zero error", bounds_v2([], pj, tgt)["status"] == "UNASSESSED" and bounds_v2([(0.0, [([0, 0, 0, 0], [0, 0, 0, 0])])], pj, tgt)["status"] == "UNASSESSED")
    # M3 invalid masses
    ok("M3 invalid masses (NaN, inf, negative) are INVALID", all(bounds_v2([(x, [([0, 0, 0, 0], [0, 0, 0, 0])])], pj, tgt)["status"] == "INVALID"
                                                               for x in (float("nan"), float("inf"), -1e-12)))
    # M4 rounded totals slightly above / below one: normalised, discrepancy recorded
    for tot in (1.0000003, 0.9999997):
        r = bounds_v2([(0.7 * tot, [([0, 0, 0, 0], [0, 0, 0, 0])]), (0.3 * tot, [([0, 1, 1, 0], [0, 1, 1, 0])])], pj, tgt)
        ok(f"M4 total {tot}: normalised, discrepancy recorded", abs(r["support_lower"] + r["error_upper"] - 1.0) < 1e-15 and abs(r["normalisation_discrepancy"] - (tot - 1.0)) < 1e-12,
           f"discrepancy {r['normalisation_discrepancy']:.3g}")
    # M5 tiny error masses are kept (direct sums), not lost to 1 - support
    r = bounds_v2([(1.0, [([0, 0, 0, 0], [0, 0, 0, 0])]), (1e-300, [([0, 1, 0, 0], [0, 0, 0, 0])])], pj, tgt)
    ok("M5 a 1e-300 error mass is retained exactly", r["error_lower"] == 1e-300 and r["support_lower"] == 1.0, f"error {r['error_lower']}")
    # M6 zero-mass competitor: evaluated alternative exists even with mass 0
    r = bounds_v2([(1.0, [([0, 0, 0, 0], [0, 0, 0, 0])]), (0.0, [([0, 1, 0, 0], [0, 0, 0, 0])])], pj, tgt)
    ok("M6 zero-mass competitor is EVALUATED_ALTERNATIVE_ZERO_MASS", r["competitor"] == "EVALUATED_ALTERNATIVE_ZERO_MASS" and r["error_upper"] == 0.0)
    r = bounds_v2([(1.0, [([0, 0, 0, 0], [0, 0, 0, 0])])], pj, tgt)
    ok("M6b a wholly supporting pool is NO_EVALUATED_ALTERNATIVE (a pool statement, not evidence of correctness)", r["competitor"] == "NO_EVALUATED_ALTERNATIVE")
    # M7 physically realisable mixed projection: two labellings spell the SAME whole haplotype ("A"+"CG" == "AC"+"G"),
    # but their compensating boundary differs on the chosen single-block interval [1..1]
    x, y = [0, 0, 0, 0], [0, 1, 1, 0]
    same_whole = spell(blocks, x) == spell(blocks, y) == "TTACGCCC"
    p1 = Projector(blocks, 1, 1)
    r = bounds_v2([(0.8, [(x, x), (y, y)]), (0.2, [([0, 1, 0, 0], [0, 1, 0, 0])])], p1, p1.pair(x, x))
    ok("M7 realisable mixed projection: whole sequences identical (independent spell), interval projections differ, mass preserved as an interval",
       same_whole and r["ambiguous_groups"] == 1 and abs(r["support_lower"]) < 1e-15 and abs(r["support_upper"] - 0.8) < 1e-15
       and abs(r["error_lower"] - 0.2) < 1e-15 and abs(r["error_upper"] - 1.0) < 1e-15,
       f"whole identical {same_whole}; support [{r['support_lower']}, {r['support_upper']}], error [{r['error_lower']}, {r['error_upper']}]")
    # the same group on the enclosing compound interval [1..2] has ONE projection (the degeneracy is inside it)
    p12 = Projector(blocks, 1, 2)
    r2 = bounds_v2([(0.8, [(x, x), (y, y)]), (0.2, [([0, 1, 0, 0], [0, 1, 0, 0])])], p12, p12.pair(x, x))
    ok("M7b on the enclosing interval the same group is one projection", r2["ambiguous_groups"] == 0 and abs(r2["support_lower"] - 0.8) < 1e-15)
    # M8 within-interval phase: truth AC/GT, prediction AT/GC; single-position unphased genotypes match, compound pairs do not
    pb = [["A", "G"], ["C", "T"]]
    truth, pred = ([0, 0], [1, 1]), ([0, 1], [1, 0])
    unphased = all(sorted((truth[0][b], truth[1][b])) == sorted((pred[0][b], pred[1][b])) for b in (0, 1))
    cp = Projector(pb, 0, 1)
    ok("M8 phase counterexample: unphased single-position genotypes match, compound sequence pair does not",
       unphased and cp.pair(*truth) != cp.pair(*pred) and sorted([spell(pb, h) for h in truth]) == ["AC", "GT"] and sorted([spell(pb, h) for h in pred]) == ["AT", "GC"])
    write("gates", res)
    print(open(f"{out}/gates.tsv").read())
    sys.exit(0 if all(x["result"] == "PASS" for x in res) else 1)

if mode == "reapply":
    old = {r["case"]: r for r in T(f"{M9I}/analyse/compound.tsv")}
    oldc = {r["case"]: r for r in T(f"{M9I}/analyse/controls.tsv")}
    rows, ctl = [], []
    zero_mass_groups_total = 0
    for cs in m9i.cases():
        blocks, ids = m9i.dump(cs["dump"])
        i, j = ids[m9i.INTERVAL[cs["locus"]][0]], ids[m9i.INTERVAL[cs["locus"]][1]]
        prob = {r["sequence_candidate_id"]: float(r["candidate_set_probability"]) for r in T(f"{cs['prefix']}.mosaic_scores.tsv")}
        members = {}
        for r in T(f"{cs['prefix']}.mosaic_candidates.tsv"):
            members.setdefault(r["sequence_candidate_id"], []).append((m9i.labs(r["hap1_alleles"]), m9i.labs(r["hap2_alleles"])))
        groups = [(p, members.get(g, [])) for g, p in prob.items()]
        zero_mass_groups_total += sum(1 for p, _ in groups if p == 0.0)
        sl = m9i.KV(f"{cs['prefix']}.sequence_locus.tsv")
        if cs["v2"] and sl.get("representative_labelling_hap1"):
            rep = (m9i.labs(sl["representative_labelling_hap1"]), m9i.labs(sl["representative_labelling_hap2"]))
        else:
            eq = T(f"{cs['prefix']}.equivalence_set.tsv")
            rep = (m9i.labs(eq[0]["hap1_alleles"]), m9i.labs(eq[0]["hap2_alleles"]))
        pj = Projector(blocks, i, j)
        b = bounds_v2(groups, pj, pj.pair(*rep))
        o = old[cs["case"]]
        rel = lambda new, k: (abs(new - float(o[k])) / max(abs(new), abs(float(o[k])))) if max(abs(new), abs(float(o[k]))) > 0 else 0.0
        rows.append(dict(case=cs["case"], status=b["status"], groups=b["groups"], missing_groups=b["missing_groups"], ambiguous_groups=b["ambiguous_groups"],
                         zero_mass_groups=sum(1 for p, _ in groups if p == 0.0), original_total=repr(b["original_total"]),
                         normalisation_discrepancy=b["normalisation_discrepancy"], competitor_v2=b["competitor"], competitor_m9i=o["competitor"],
                         support_lower=b["support_lower"], support_upper=b["support_upper"], error_lower=b["error_lower"], error_upper=b["error_upper"],
                         max_rel_change_vs_m9i=max(rel(b[k], f"compound_{k}") for k in ("support_lower", "support_upper", "error_lower", "error_upper")),
                         compound_exact_m9i=o["compound_exact"], compound_distance_m9i=o["compound_distance"]))
        cb = m9i.CONTROL[cs["locus"]]
        pc = Projector(blocks, ids[cb], ids[cb])
        rc = bounds_v2(groups, pc, pc.pair(*rep))
        oc = oldc[cs["case"]]
        ex = {r["block_id"]: r for r in T(f"{cs['prefix']}.sequence_blocks.tsv")}.get(cb, {})
        ctl.append(dict(case=cs["case"], block=cb, error_lower=rc["error_lower"], error_upper=rc["error_upper"],
                        saved_error_mass_lower=ex.get("error_mass_lower", "NA (legacy)"), saved_error_mass_upper=ex.get("error_mass_upper", "NA (legacy)"),
                        agrees=(m9i.agree(rc["error_lower"], ex.get("error_mass_lower")) and m9i.agree(rc["error_upper"], ex.get("error_mass_upper")))
                        if ex.get("error_mass_lower") else "NA", m9i_agrees=oc["agrees"]))
    write("reapply", rows)
    write("controls", ctl)
    summ = dict(cases=len(rows), status_counts=";".join(f"{s}={sum(1 for r in rows if r['status'] == s)}" for s in sorted({r['status'] for r in rows})),
                competitor_changes=sum(1 for r in rows if (r["competitor_v2"] == "NO_EVALUATED_ALTERNATIVE") != (r["competitor_m9i"] == "NO_EVALUATED_ALTERNATIVE")
                                       or r["competitor_v2"] in ("UNKNOWN", "EVALUATED_ALTERNATIVE_ZERO_MASS")),
                competitor_v2_counts=";".join(f"{s}={sum(1 for r in rows if r['competitor_v2'] == s)}" for s in sorted({r['competitor_v2'] for r in rows})),
                max_abs_normalisation_discrepancy=max(abs(r["normalisation_discrepancy"]) for r in rows),
                max_rel_bound_change=max(r["max_rel_change_vs_m9i"] for r in rows), missing_groups=sum(r["missing_groups"] for r in rows),
                ambiguous_groups=sum(r["ambiguous_groups"] for r in rows), zero_mass_groups=zero_mass_groups_total,
                controls_agree=f"{sum(1 for c in ctl if c['agrees'] is True)}/{sum(1 for c in ctl if c['agrees'] != 'NA')}")
    write("summary", [summ])
    print(open(f"{out}/summary.tsv").read())
    sys.exit(0)
