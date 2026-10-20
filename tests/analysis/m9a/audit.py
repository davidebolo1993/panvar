#!/usr/bin/env python3
"""M9-A read-only uncertainty-target audit (frozen with DESIGN.md).

  audit.py gates <out>    synthetic fixtures for the aggregation (swap, phase-only, disagreeing, same sequence under
                          different labels, mixed projection, absent alternatives) and normalisation
  audit.py audit <out>    the eight completed M8-C/M8-D runs: targets A-E from saved scores, candidate projections and
                          equivalence information; reproductions of existing quantities; the M8-D qualification numbers

Hypotheses are the caller's sequence-level groups (one unordered whole-haplotype sequence pair each, sequence_candidate_id).
Each is counted once. Its members are the label pairs (candidate rows) spelling it. For a target with a per-member
predicate, a group is SAME (every member carries the reported state), DISAGREEING (none does) or MIXED (some do):
the existing mixed-projection semantics (genotype_block_support.cpp). Mass is reported as an interval
[SAME, SAME + MIXED]; a mixed group is never assigned to a representative and never discarded. Masses are conditional on
each run's evaluated pool; they are not compared across pools as calibration.
"""
import csv
import math
import os
import sys

mode, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
R = "/Users/davide.bolognini/panvar/tests/results"
C, D = f"{R}/m8c_copying_challenge_20261020/runs", f"{R}/m8d_hg01975_loo_20261020/runs"
RUNS = [(f"{d}.{s}.{r}", f"{C}/call.{d}.{s}.{r}") for d in ("HG03209", "HG01975") for s in ("LZO", "LOO") for r in ("r0.01", "r0.1")
        if not (d == "HG01975" and s == "LOO")] + [(f"HG01975.LOO.{r}", f"{D}/call.HG01975.LOO.{r}") for r in ("r0.01", "r0.1")]
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l and not l.startswith("#")}


def write(name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def lse(xs):
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))


# ---------------------------------------------------------------------------------------------------------------- core
def build(groups):
    """groups: list of dict(id, score, fragment, exposure, prior, members=[(hap1 tuple, hap2 tuple)]).
    Members are deduplicated by unordered whole-haplotype identity (a swapped row is the same label pair); a label pair
    in two groups is inconsistent input. Returns groups with normalised members and log-normalised probabilities."""
    seen = {}
    for g in groups:
        uniq = []
        for h1, h2 in g["members"]:
            key = tuple(sorted((tuple(h1), tuple(h2))))
            if key in seen and seen[key] != g["id"]:
                raise ValueError(f"label pair in two sequence hypotheses: {g['id']} and {seen[key]}")
            if key not in seen:
                seen[key] = g["id"]
                uniq.append(key)
        g["members"] = uniq
        g["G"] = [tuple(tuple(sorted((a, b))) for a, b in zip(m[0], m[1])) for m in uniq]
    ids = [g["id"] for g in groups]
    if len(set(ids)) != len(ids):
        raise ValueError("a sequence hypothesis listed twice")
    z = lse([g["score"] for g in groups])
    for g in groups:
        g["p"] = math.exp(g["score"] - z)
    groups.sort(key=lambda g: -g["score"])
    return groups


geno = lambda m: tuple(tuple(sorted((a, b))) for a, b in zip(m[0], m[1]))


def target(groups, winner, pred, applicable=True):
    """Mass interval of the reported state and its best genuinely disagreeing competitor."""
    if not applicable:
        return dict(applicable=False)
    same = mixed = 0.0
    n = dict(same=0, mixed=0, disagreeing=0)
    best_d = best_s = None
    for g in groups:
        c = [pred(m, gm) for m, gm in zip(g["members"], g["G"])]
        if all(c):
            same += g["p"]
            n["same"] += 1
            if g is not winner and best_s is None:
                best_s = g
        elif not any(c):
            n["disagreeing"] += 1
            if best_d is None:
                best_d = g
        else:
            mixed += g["p"]
            n["mixed"] += 1
    wc = [pred(m, gm) for m, gm in zip(winner["members"], winner["G"])]
    row = dict(applicable=True, reported_state_projection="UNAMBIGUOUS" if all(wc) else "AMBIGUOUS_PROJECTION",
               mass_lower=same, mass_upper=same + mixed, same_groups=n["same"], mixed_groups=n["mixed"], disagreeing_groups=n["disagreeing"])
    if best_d is None:
        row.update(competitor="NO_EVALUATED_ALTERNATIVE")
    else:
        row.update(competitor=best_d["id"], gap_total=winner["score"] - best_d["score"],
                   gap_likelihood=(winner["fragment"] + winner["exposure"]) - (best_d["fragment"] + best_d["exposure"]),
                   gap_fragment=winner["fragment"] - best_d["fragment"], gap_exposure=winner["exposure"] - best_d["exposure"],
                   gap_prior=winner["prior"] - best_d["prior"])
    row.update(best_same_state_alternative=best_s["id"] if best_s else "none",
               gap_to_best_same_state_alternative=winner["score"] - best_s["score"] if best_s else "")
    return row


def phase_vector(m, B, het):
    """Orientation of a member carrying B: the haplotype holding the smaller allele at the first het block is 'A';
    the vector is A's allele at every het block (the label-level phase)."""
    if not het:
        return ()
    b0 = het[0]
    a = m[0] if m[0][b0] == B[b0][0] else m[1]
    return tuple(a[b] for b in het)


def audit_pool(groups, B, reportable):
    """All targets for one pool. B: reported unordered genotype per block; reportable: block flags."""
    winner = groups[0]
    nb = len(B)
    het = [b for b in range(nb) if B[b][0] != B[b][1]]
    rows = {}
    # A: identity is the sequence hypothesis itself, so the predicate is group membership
    a_same = winner["p"]
    best_other = next((g for g in groups if g is not winner), None)
    rows["A exact unordered sequence pair"] = dict(applicable=True, reported_state_projection="UNAMBIGUOUS (by construction)", mass_lower=a_same, mass_upper=a_same,
                                                   same_groups=1, mixed_groups=0, disagreeing_groups=len(groups) - 1,
                                                   **(dict(competitor=best_other["id"], gap_total=winner["score"] - best_other["score"],
                                                           gap_likelihood=(winner["fragment"] + winner["exposure"]) - (best_other["fragment"] + best_other["exposure"]),
                                                           gap_fragment=winner["fragment"] - best_other["fragment"], gap_exposure=winner["exposure"] - best_other["exposure"],
                                                           gap_prior=winner["prior"] - best_other["prior"]) if best_other else dict(competitor="NO_EVALUATED_ALTERNATIVE")),
                                                   best_same_state_alternative="none (the state is one hypothesis)")
    carriesB = lambda m, gm: gm == B
    for g in groups:   # label-level phase of every member carrying B, computed once
        g["PH"] = [phase_vector(m, B, het) if gm == B else None for m, gm in zip(g["members"], g["G"])]
    phase_of = {}
    for g in groups:
        for m, ph in zip(g["members"], g["PH"]):
            phase_of[id(m)] = ph
    rows["B unphased genotype, all blocks"] = target(groups, winner, carriesB)
    sv = [b for b in range(nb) if reportable[b]]
    rows["C reportable-block genotype agreement"] = target(groups, winner, lambda m, gm: all(gm[b] == B[b] for b in sv), bool(sv))
    # E: phase given B (only members carrying B; requires two or more het blocks)
    wph = {ph for ph in winner["PH"] if ph is not None}
    b_unamb = rows["B unphased genotype, all blocks"]["reported_state_projection"] == "UNAMBIGUOUS"
    ph_ok = len(het) >= 2 and len(wph) == 1 and b_unamb
    if ph_ok:
        rph = next(iter(wph))
        joint = target(groups, winner, lambda m, gm: phase_of[id(m)] == rph)
        bmass = rows["B unphased genotype, all blocks"]
        # the genuine phase competitor: best group whose members all carry B with another phase
        comp = next((g for g in groups if g is not winner and all(ph is not None and ph != rph for ph in g["PH"])), None)
        e = dict(applicable=True, reported_state_projection="UNAMBIGUOUS", het_blocks=len(het),
                 mass_lower=joint["mass_lower"] / bmass["mass_upper"], mass_upper=min(1.0, joint["mass_upper"] / bmass["mass_lower"]),
                 joint_mass_lower=joint["mass_lower"], joint_mass_upper=joint["mass_upper"], conditioning="given B (mass of B-carrying hypotheses)")
        if comp is None:
            e.update(competitor="NO_EVALUATED_ALTERNATIVE")
        else:
            e.update(competitor=comp["id"], gap_total=winner["score"] - comp["score"],
                     gap_likelihood=(winner["fragment"] + winner["exposure"]) - (comp["fragment"] + comp["exposure"]),
                     gap_fragment=winner["fragment"] - comp["fragment"], gap_exposure=winner["exposure"] - comp["exposure"], gap_prior=winner["prior"] - comp["prior"])
        rows["E phase given B"] = e
    else:
        rows["E phase given B"] = dict(applicable=False, het_blocks=len(het),
                                       note="fewer than two het blocks" if len(het) < 2 else "reported unphased genotype is an ambiguous projection" if not b_unamb
                                       else "winner members disagree on phase (ambiguous)")
    blocks = []
    for b in range(nb):
        r = target(groups, winner, lambda m, gm, b=b: gm[b] == B[b])
        r["block_index"] = b
        blocks.append(r)
    pairs = []
    if ph_ok:
        bm = rows["B unphased genotype, all blocks"]
        rel = lambda ph, i, j: (ph[i] == B[het[i]][0]) == (ph[j] == B[het[j]][0])
        for i in range(len(het) - 1):
            rr = rel(rph, i, i + 1)
            j = target(groups, winner, lambda m, gm, i=i, rr=rr: phase_of[id(m)] is not None and rel(phase_of[id(m)], i, i + 1) == rr)
            pairs.append(dict(left_block=het[i], right_block=het[i + 1], mass_lower=j["mass_lower"] / bm["mass_upper"],
                              mass_upper=min(1.0, j["mass_upper"] / bm["mass_lower"]), conditioning="given B"))
    return rows, blocks, pairs


def contenders(groups):
    """The caller's fragment contenders (genotype_repair.cpp): rivals the reads alone do not put below the winner."""
    w = groups[0]
    return [g for g in groups[1:] if not (g["fragment"] + g["bound"] < w["fragment"] - w["bound"])]


# ---------------------------------------------------------------------------------------------------------------- gates
def G(i, score, members, fragment=None, exposure=0.0):
    return dict(id=i, score=score, fragment=score if fragment is None else fragment, exposure=exposure, prior=0.0, members=members, bound=0.0)


if mode == "gates":
    res = []
    ok = lambda name, cond, detail="": res.append(dict(gate=name, result="PASS" if cond else "FAIL", detail=detail))
    # F1 whole-haplotype swapping: a swapped row is the same label pair, counted once
    g = build([G(0, 0.0, [((1, 2), (3, 4)), ((3, 4), (1, 2))]), G(1, -1.0, [((1, 4), (3, 2))])])
    B = geno(g[0]["members"][0])
    rows, _, _ = audit_pool(g, B, [True, True])
    ok("F1 swap counted once", len(g[0]["members"]) == 1 and abs(rows["A exact unordered sequence pair"]["mass_lower"] - 1 / (1 + math.e ** -1)) < 1e-12,
       f"members={len(g[0]['members'])} A={rows['A exact unordered sequence pair']['mass_lower']:.6f}")
    # F2 same block genotypes, different phase: B carries both; A and E split; E competitor is the other phase
    ok("F2 phase-only alternative", abs(rows["B unphased genotype, all blocks"]["mass_lower"] - 1) < 1e-12
       and rows["B unphased genotype, all blocks"]["competitor"] == "NO_EVALUATED_ALTERNATIVE"
       and rows["E phase given B"]["competitor"] == 1 and abs(rows["E phase given B"]["mass_lower"] - 1 / (1 + math.e ** -1)) < 1e-12,
       f"B={rows['B unphased genotype, all blocks']['mass_lower']:.3f} E={rows['E phase given B']['mass_lower']:.6f}")
    # F3 genuinely different block genotypes: competitor and gap
    g = build([G(0, 0.0, [((1, 2), (3, 4))], fragment=-1.0, exposure=1.0), G(1, -2.0, [((1, 2), (3, 5))], fragment=-1.5), G(2, -3.0, [((2, 1), (4, 3))])])
    B = geno(g[0]["members"][0])
    rows, blocks, _ = audit_pool(g, B, [False, True])
    ok("F3 disagreeing competitor", rows["B unphased genotype, all blocks"]["competitor"] == 1 and abs(rows["B unphased genotype, all blocks"]["gap_total"] - 2) < 1e-12
       and abs(rows["B unphased genotype, all blocks"]["gap_likelihood"] - 1.5) < 1e-12 and blocks[0]["competitor"] == 2
       and rows["C reportable-block genotype agreement"]["competitor"] == 1, f"gap={rows['B unphased genotype, all blocks']['gap_total']} lik={rows['B unphased genotype, all blocks']['gap_likelihood']}")
    # F4 identical sequence under different labels: one hypothesis, two members, mass once; agreeing members are not mixed
    g = build([G(0, 0.0, [((1, 2), (3, 4)), ((1, 2), (3, 7))]), G(1, 0.0, [((1, 9), (3, 4))])])
    rows, blocks, _ = audit_pool(g, geno(((1, 2), (3, 4))), [True, True])
    tot = sum(x["p"] for x in g)
    ok("F4 same sequence, different labels", abs(tot - 1) < 1e-12 and abs(sum(x["p"] for x in g if x["id"] == 0) - 0.5) < 1e-12
       and blocks[0]["mixed_groups"] == 0 and blocks[1]["mixed_groups"] == 1,
       f"hypotheses={len(g)} mass(id0)={sum(x['p'] for x in g if x['id'] == 0):.3f}")
    # F5 mixed projection: a hypothesis whose label pairs project to different genotypes is MIXED, never assigned or dropped
    res.append(dict(gate="F5a a label pair in two hypotheses is rejected", result="PASS", detail=""))
    try:
        build([G(0, 0.0, [((1, 2), (3, 4))]), G(1, -0.5, [((3, 4), (1, 2))])])
        res[-1]["result"], res[-1]["detail"] = "FAIL", "accepted"
    except ValueError as e:
        res[-1]["detail"] = str(e)
    g = build([G(0, 0.0, [((1, 2), (3, 4))]), G(1, -0.5, [((3, 2), (1, 4)), ((1, 2), (3, 5))]), G(2, -1.0, [((1, 2), (8, 5))])])
    B = geno(((1, 2), (3, 4)))
    rows, blocks, _ = audit_pool(g, B, [True, True])
    p = {x["id"]: x["p"] for x in g}
    r = rows["B unphased genotype, all blocks"]
    ok("F5 mixed projection interval", r["mixed_groups"] == 1 and abs(r["mass_lower"] - p[0]) < 1e-12 and abs(r["mass_upper"] - p[0] - p[1]) < 1e-12
       and r["competitor"] == 2, f"lower={r['mass_lower']:.4f} upper={r['mass_upper']:.4f} competitor={r['competitor']}")
    # F6 absent disagreeing alternatives
    g = build([G(0, 0.0, [((1, 2), (3, 4))])])
    rows, blocks, _ = audit_pool(g, geno(((1, 2), (3, 4))), [True, False])
    ok("F6 no evaluated alternative", all(rows[k].get("competitor") == "NO_EVALUATED_ALTERNATIVE" for k in rows if rows[k].get("applicable"))
       and rows["A exact unordered sequence pair"]["mass_lower"] == 1.0, "single hypothesis")
    # F7 normalisation and the caller's block posterior on a pool without mixed groups
    g = build([G(i, -0.37 * i, [((i, 0), (0, i % 3))]) for i in range(1, 40)])
    B = geno(g[0]["members"][0])
    _, blocks, _ = audit_pool(g, B, [True, True])
    other = [lse([x["score"] for x in g if x["G"][0][b] != B[b]]) - lse([x["score"] for x in g]) for b in range(2)]
    ok("F7 normalisation and block posterior", abs(sum(x["p"] for x in g) - 1) < 1e-12 and all(abs(blocks[b]["mass_lower"] - (-math.expm1(other[b]))) < 1e-12 for b in range(2)),
       f"sum={sum(x['p'] for x in g):.15f}")
    write("gates", res)
    print(open(f"{out}/gates.tsv").read())
    sys.exit(0 if all(r["result"] == "PASS" for r in res) else 1)

# ---------------------------------------------------------------------------------------------------------------- audit
if mode == "audit":
    pools, trows, brows, prows, rep = [], [], [], [], []
    only = os.environ.get("M9A_SMOKE")   # disclosed smoke test on one run before the freeze; never used by the runner
    for name, d in RUNS:
        if only and name != only:
            continue
        sc = T(f"{d}/call.mosaic_scores.tsv")
        members = {}
        for c in T(f"{d}/call.mosaic_candidates.tsv"):
            members.setdefault(c["sequence_candidate_id"], []).append((tuple(map(int, c["hap1_alleles"].split(","))), tuple(map(int, c["hap2_alleles"].split(",")))))
        groups = build([dict(id=int(r["sequence_candidate_id"]), score=float(r["score"]), fragment=float(r["fragment_term"]), exposure=float(r["exposure_term"]),
                             prior=float(r["sequence_pair_prior"]), bound=float(r["numeric_error_bound"]), members=members.pop(r["sequence_candidate_id"]),
                             file_p=float(r["candidate_set_probability"]), file_ambiguous=r["projection_ambiguous"] == "true",
                             file_winner=r["in_winner_set"] == "true", file_carries=r["carries_called_genotype"] == "true") for r in sc])
        assert not members, "candidate rows without a scored hypothesis"
        lg = T(f"{d}/call.locus_genotype.tsv")[0]
        bg = T(f"{d}/call.block_genotypes.tsv")
        eq = T(f"{d}/call.equivalence_set.tsv")
        st = KV(f"{d}/call.mosaic_status.tsv")
        ev = KV(f"{D}/eval.{name.split('.')[-2]}.{name.split('.')[-1]}/eval.summary.tsv" if d.startswith(D) else f"{d.replace('/call.', '/eval.')}/eval.summary.tsv")
        B = tuple(tuple(sorted((int(b["allele1"]), int(b["allele2"])))) for b in bg)
        reportable = [b["reportable"] == "true" for b in bg]
        w = groups[0]
        # the reported state is the caller's: unrepaired, the top group, one unordered genotype
        locus_pair = {tuple(map(int, lg["hap1_block_path"].split(","))), tuple(map(int, lg["hap2_block_path"].split(",")))}
        consistent = (lg["reported_genotype_source"] == "caller" and lg["dosage_content_corrected_blocks"] == "0" and lg["dosage_complete_applied_blocks"] == "0"
                      and locus_pair == set(w["members"][0]) and len(eq) == 1 and int(eq[0]["sequence_candidate_id"]) == w["id"]
                      and sum(g["file_winner"] for g in groups) == 1 and w["file_winner"] and geno(w["members"][0]) == B)
        pools.append(dict(run=name, hypotheses=len(groups), label_pairs=sum(len(g["members"]) for g in groups),
                          multi_label_hypotheses=sum(len(g["members"]) > 1 for g in groups), projection_ambiguous_hypotheses=sum(g["file_ambiguous"] for g in groups),
                          top_set_size=sum(g["file_winner"] for g in groups), candidate_scope=lg["candidate_scope"], search_conditional=lg["search_conditional"],
                          proposal_search_exhaustive=ev.get("proposal_search_exhaustive", ""), reported_state_is_caller_top_group=consistent,
                          het_blocks=sum(1 for b in B if b[0] != b[1]), reportable_blocks=sum(reportable)))
        rows, blocks, pairs = audit_pool(groups, B, reportable)
        for k, r in rows.items():
            trows.append(dict(run=name, target=k, **r))
        sup = {}
        con = contenders(groups)
        for b, r in enumerate(blocks):
            read_sup = all(g["G"][0][b] == B[b] for g in con)   # the caller uses each group's representative genotype
            brows.append(dict(run=name, block=bg[b]["block_id"], kind=bg[b]["kind"], reportable=bg[b]["reportable"], reported=f"{B[b][0]}/{B[b][1]}",
                              **{k: v for k, v in r.items() if k != "block_index"},
                              existing_conditional_genotype_posterior=bg[b]["conditional_genotype_posterior"], existing_GQ=bg[b]["conditional_GQ"],
                              posterior_abs_diff=abs(float(bg[b]["conditional_genotype_posterior"]) - r["mass_lower"]),
                              existing_genotype_support=bg[b]["genotype_support"], reproduced_genotype_support="READS" if read_sup else "PRIOR_ONLY"))
        for p in pairs:
            prows.append(dict(run=name, left=bg[p["left_block"]]["block_id"], right=bg[p["right_block"]]["block_id"], mass_lower=p["mass_lower"],
                              mass_upper=p["mass_upper"], conditioning=p["conditioning"]))
        # reproductions of existing quantities
        sep_file = float(lg["separation"])
        brow = rows["B unphased genotype, all blocks"]
        rep += [
            dict(run=name, quantity="normalisation (sum of recomputed probabilities)", existing=1.0, reproduced=sum(g["p"] for g in groups)),
            dict(run=name, quantity="max |candidate_set_probability - recomputed|", existing=0.0, reproduced=max(abs(g["p"] - g["file_p"]) for g in groups)),
            dict(run=name, quantity="locus genotype separation (caller_genotype_separation_gap)", existing=sep_file,
                 reproduced=brow.get("gap_total", float("inf"))),
            dict(run=name, quantity="mass carrying the reported unphased genotype (carries_called_genotype)", existing=sum(g["p"] for g in groups if g["file_carries"]),
                 reproduced=brow["mass_lower"]),
            dict(run=name, quantity="numeric_separation_gap (rank 1 - rank 2 sequence pair)", existing=float(st["numeric_separation_gap"]),
                 reproduced=rows["A exact unordered sequence pair"].get("gap_total", float("inf"))),
            dict(run=name, quantity="winner_candidate_set_probability", existing=float(st["winner_candidate_set_probability"]), reproduced=w["p"]),
            dict(run=name, quantity="fragment_contenders", existing=int(ev["fragment_contenders"]), reproduced=len(con)),
            dict(run=name, quantity="phase_support", existing=lg["phase_support"], reproduced="READS_WITHIN_CANDIDATES" if not con else "PRIOR_ONLY"),
            dict(run=name, quantity="max |block posterior - existing|", existing=0.0, reproduced=max(x["posterior_abs_diff"] for x in brows if x["run"] == name)),
            dict(run=name, quantity="block genotype_support labels reproduced", existing=len(bg),
                 reproduced=sum(x["existing_genotype_support"] == x["reproduced_genotype_support"] for x in brows if x["run"] == name)),
        ]
    write("pools", pools)
    write("targets", trows)
    write("blocks", brows)
    write("phase_pairs", prows)
    write("reproductions", rep)
    # M8-D qualification (Sol's independent values)
    q = []
    for name, sol_sep, sol_mass in (("HG01975.LOO.r0.01", 0.7021371507, 0.405311), ("HG01975.LOO.r0.1", 0.9018471817, 0.313111)):
        if not any(r["run"] == name for r in trows):
            continue
        d = dict(RUNS)[name]
        lg = T(f"{d}/call.locus_genotype.tsv")[0]
        st = KV(f"{d}/call.mosaic_status.tsv")
        b = next(r for r in trows if r["run"] == name and r["target"].startswith("B"))
        a = next(r for r in trows if r["run"] == name and r["target"].startswith("A"))
        q.append(dict(run=name, sol_locus_separation=sol_sep, archived_locus_separation=lg["separation"], sol_mass_reported_unphased=sol_mass,
                      reproduced_mass_reported_unphased=f"{b['mass_lower']:.6f}", numeric_separation_gap_quoted_in_M8D=st["numeric_separation_gap"],
                      rank2_carries_reported_unphased_genotype=a["competitor"] in {g for g in [b.get("best_same_state_alternative")]},
                      reproduced=abs(float(lg["separation"]) - sol_sep) < 1e-9 and abs(b["mass_lower"] - sol_mass) < 5e-7))
    write("m8d_qualification", q)
    for n in ("pools", "reproductions", "m8d_qualification"):
        print(open(f"{out}/{n}.tsv").read())
    sys.exit(0)
