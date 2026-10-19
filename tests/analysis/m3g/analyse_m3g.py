#!/usr/bin/env python3
"""M3-G decision table (diagnostic; frozen before the probe units run).

Inputs: saved.trace (per-round scored ctx:0 neighbours), rank (emission ranks), score (frozen list scored).
Checks first (a failed check is reported, and the classification below is then UNRESOLVED):
  C1 the saved round moves' local_pair_rank equal the probe's rank of the moved pair at the moved block;
  C2 the set of ctx:0 unordered pairs the saved rounds scored equals the probe's top-repair_local_pairs set at ctx:0;
  C3 every saved total declared in candidates.tsv reproduced (score-set refuses otherwise).
Per comparison (against R, and against the call W1/W2): fragment, exposure, prior and total differences; unordered block
genotype identity; sequence identity (score-set's interned haplotype ids); whether it beats; whether it is one existing
single-block move from R (ordinary refinement) -- true for R_* by construction.
Classification (frozen precedence; one primary class, the context finding reported beside it):
  PROPOSAL_SHORTLIST_OMISSION   target rank > local_pairs at ctx:0, never scored by the repair, and the best target
                                assignment in R beats R on total (the model would have taken the move)
  SEARCH_BOOKKEEPING            target rank <= local_pairs but no round scored it
  PRIOR_REVERSES                the best target assignment in R beats R on likelihood (F+E) but loses on total
  EVALUATED_LOSES               the target assignment loses to R on likelihood and total (whether scored or not)
  UNRESOLVED                    a check failed or the evidence does not fit
  context qualifier             CONTEXT_DEPENDENT when the target beats its context in R but not in W1/W2 (or the
                                reverse), or the called pair's preference over the repair's pair flips between contexts
Also: would the best target candidate have become the final call? (total vs the saved pooled winner, -31870.9193977016)

  analyse_m3g.py <runs dir> <out prefix>
"""
import csv
import sys

runs, prefix = sys.argv[1], sys.argv[2]
T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
LOCAL = 32
WINNER = -31870.9193977016
ranks = T(f"{runs}/rank/rank.local_ranks.tsv")
rank_of = {(r["block_id"], tuple(sorted((int(r["allele_x"]), int(r["allele_y"]))))): int(r["rank"]) for r in ranks}
out, checks = [], []
target = rank_of[("ctx:0", (14, 56))]
called = rank_of[("ctx:0", (4, 52))]
inc = rank_of[("ctx:0", (40, 52))]
emis = {(r["block_id"], int(r["rank"])): r for r in ranks}
out.append(f"ctx:0 emission rank (unordered pairs, repair rule): target (14,56) {target}; called (4,52) {called}; repair's (40,52) {inc}; "
           f"local list = top {LOCAL}; pairs ranked {sum(1 for r in ranks if r['block_id'] == 'ctx:0')}")
tr = next(r for r in ranks if r["block_id"] == "ctx:0" and int(r["rank"]) == target)
top = emis[("ctx:0", 1)]
cut = emis[("ctx:0", LOCAL)]
out.append(f"ctx:0 emission: top {top['emission']}; rank-{LOCAL} {cut['emission']}; target {tr['emission']} (below top {tr['below_top']}; "
           f"below the cutoff by {float(cut['emission']) - float(tr['emission']):.6f})")
# C1: moves
rounds = T(f"{runs}/../../m3f_validation_20261019/runs/run.ts1_2.off/call.repair_rounds.tsv")
ok1 = True
for r in rounds:
    if r["moved"] != "true":
        continue
    p = (r["repaired_block"], tuple(sorted((int(r["allele_x"]), int(r["allele_y"])))))
    probe = rank_of.get(p)
    ok1 &= str(probe) == r["local_pair_rank"]
    out.append(f"C1 round {r['round']} {p[0]} {p[1]}: saved local_pair_rank {r['local_pair_rank']} vs probe {probe}")
checks.append(("C1 move ranks reproduce", ok1))
# C2: scored ctx:0 pairs per round == probe top-LOCAL
probe_top = {p for (b, p), k in rank_of.items() if b == "ctx:0" and k <= LOCAL}
nb = T(f"{runs}/saved.trace/trace.block_neighbours.tsv")
per_round = {}
for r in nb:
    per_round.setdefault(r["round"], set()).add(tuple(sorted(map(int, r["pair"].split(",")))))
ok2 = all(s | {(40, 52)} == probe_top | {(40, 52)} for s in per_round.values()) and len(per_round) == 7
out.append(f"C2 saved scored ctx:0 pairs per round equal the probe top-{LOCAL} (the incumbent pair (40,52) itself only by swap): "
           f"{ok2}; rounds {sorted(per_round)}; target scored in any round: {any((14, 56) in s for s in per_round.values())}")
checks.append(("C2 scored ctx:0 set equals local list", ok2))
# C3 + comparisons
sc = {r["label"]: r for r in T(f"{runs}/score/score.scores_truth_free.tsv")}
checks.append(("C3 saved totals reproduced", all(abs(float(r["total"]) - float(r["expected_total"])) <= 1e-6 for r in sc.values()
                                                  if r["expected_total"] not in ("NA", ""))))
f = lambda r: (float(r["fragment_term"]), float(r["exposure_term"]), float(r["prior"]), float(r["total"]))
geno = lambda r: tuple(tuple(sorted(p)) for p in zip(map(int, r["hap1_alleles"].split(",")), map(int, r["hap2_alleles"].split(","))))
seq = lambda r: tuple(sorted((r["hap1_unique"], r["hap2_unique"])))
rows = ["label\tversus\tdelta_fragment\tdelta_exposure\tdelta_likelihood\tdelta_prior\tdelta_total\tbeats\tsame_unordered_genotype_as_versus\t"
        "same_sequence_pair_as_versus\tsingle_block_move_from_R\tsame_sequence_pair_as_truth\ttotal_minus_pooled_winner"]
truth = sc["R_t0"]   # R_t0 is the truth path by construction (DESIGN.md)
for lab, r in sc.items():
    for vs in ("R", "W1", "W2"):
        if lab == vs or (vs == "R" and not lab.startswith("R")) or (vs != "R" and not lab.startswith(vs)):
            continue
        a, b = f(r), f(sc[vs])
        d = [x - y for x, y in zip(a, b)]
        rows.append("\t".join(map(str, [lab, vs, f"{d[0]:.6f}", f"{d[1]:.6f}", f"{d[0] + d[1]:.6f}", f"{d[2]:.6f}", f"{d[3]:.6f}", d[3] > 0,
                                        geno(r) == geno(sc[vs]), seq(r) == seq(sc[vs]), lab.startswith("R"), seq(r) == seq(truth),
                                        f"{a[3] - WINNER:.6f}"])))
    if lab in ("R", "W1", "W2"):
        rows.append("\t".join(map(str, [lab, "pooled_winner", "", "", "", "", "", "", "", "", lab == "R", seq(r) == seq(truth), f"{f(r)[3] - WINNER:.6f}"])))
bt = max(("R_t0", "R_t1"), key=lambda k: f(sc[k])[3])
dl = (f(sc[bt])[0] + f(sc[bt])[1]) - (f(sc["R"])[0] + f(sc["R"])[1])
dt = f(sc[bt])[3] - f(sc["R"])[3]
ctx = []
for w in ("W1", "W2"):
    btw = max((w + "_t0", w + "_t1"), key=lambda k: f(sc[k])[3])
    ctx.append(f(sc[btw])[3] - f(sc[w])[3])
    rw = max((w + "_r0", w + "_r1"), key=lambda k: f(sc[k])[3])
    ctx.append(-(f(sc[rw])[3] - f(sc[w])[3]))   # called (4,52) over the repair's (40,52) in this context
rc = max(("R_c0", "R_c1"), key=lambda k: f(sc[k])[3])
called_over_repair_in_R = f(sc[rc])[3] - f(sc["R"])[3]
context_dependent = (dt > 0) != (ctx[0] > 0) or (dt > 0) != (ctx[2] > 0) or (called_over_repair_in_R > 0) != (ctx[1] > 0) or \
                    (called_over_repair_in_R > 0) != (ctx[3] > 0)
if not all(ok for _, ok in checks):
    cls = "UNRESOLVED"
elif target > LOCAL and not any((14, 56) in s for s in per_round.values()) and dt > 0:
    cls = "PROPOSAL_SHORTLIST_OMISSION"
elif target <= LOCAL and not any((14, 56) in s for s in per_round.values()):
    cls = "SEARCH_BOOKKEEPING"
elif dl > 0 and dt <= 0:
    cls = "PRIOR_REVERSES"
elif dl <= 0 and dt <= 0:
    cls = "EVALUATED_LOSES"
else:
    cls = "UNRESOLVED"
out.append(f"best target assignment in R: {bt}; total - R {dt:.6f} (likelihood {dl:.6f}); total - pooled winner {f(sc[bt])[3] - WINNER:.6f}")
out.append(f"target over its context: W1 {ctx[0]:.6f}, W2 {ctx[2]:.6f}; called (4,52) over repair's (40,52): in R {called_over_repair_in_R:.6f}, "
           f"in W1 {ctx[1]:.6f}, in W2 {ctx[3]:.6f}")
out.append("checks: " + "; ".join(f"{n} {'PASS' if ok else 'FAIL'}" for n, ok in checks))
out.append(f"CLASS\t{cls}{' + CONTEXT_DEPENDENT' if context_dependent else ''}")
open(prefix + ".comparisons.tsv", "w").write("\n".join(rows) + "\n")
open(prefix + ".decision.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
print("\n".join(rows))
