#!/usr/bin/env python3
"""M6-B attribution summary (oracle diagnostics; nothing here is a prediction).

Tables: candidate losses under the three profiles (real reads/depth; noiseless truth block in the fixed caller context;
noiseless complete true diplotype), binding and out-of-dictionary content; contrast decompositions by location category;
the unchanged active catalogue under each profile (best group's unit structure from M6-A's evaluation units, ranks of the
reported 3+1 and the floor-attaining 2+2, first group holding a 2+2); reused fragment scores (M6-A score-set for the
catalogue members, M1-E score-set for M1-E constructions -- same run, model, reads and library; truth_a/truth_b matched to
M1-E truthblock by bytes).
Frozen causal reading (C4 NA19650 LOO 10x, reported 3+1 vs floor-attaining 2+2):
  SAMPLING_OR_NORMALIZATION  real profile prefers 3+1 but the noiseless in-context truth profile prefers the 2+2
  RESIDUAL_CONTENT_MISMATCH  the noiseless in-context truth profile ALSO prefers the 3+1 over the available 2+2 (the
                             available 2+2 content is further from the truth's counts than the 3+1 is)
  CONTEXT_EFFECT             the complete-true-diplotype profile reverses the in-context conclusion
  COUNT_INDISTINGUISHABLE    for content-matched constructions (M1-E e13 vs a22), the predicted profiles are identical
  (several may hold; each is reported with its evidence)

  summary.py <M6-B runs> <M6-A runs> <M1-E runs> <out prefix>
"""
import csv
import sys

runs, m6a, m1e, prefix = sys.argv[1:5]
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l}
A = f"{runs}/attrib/ca"
cands = {r["label"]: r for r in T(A + ".candidates.tsv")}
contr = T(A + ".contrasts.tsv")
cat = T(A + ".catalogue.tsv")
units = {r["name"]: r["unit_copies"] for r in T(f"{m6a}/units.c4_loo_na19650_10/units.tsv")}
U = lambda a: units.get(f"allele_{a}", "NA")
out = ["== candidates (loss under each profile; oracle = truth-derived diagnostic only)"]
prof = [k for k in next(iter(cands.values())) if k.startswith("loss_")]
out.append("label\toracle\tbound\ttarget_kmers_outside_dictionary\t" + "\t".join(prof))
for k, r in cands.items():
    out.append("\t".join([k, r["oracle"], r["bound"], r["target_kmers_outside_dictionary"]] + [f"{float(r[p]):.3f}" for p in prof]))
out += ["", "== contrasts (A - B; positive = A fits worse)", "contrast\tprofile\tdelta\treconstructed\tdiffering\tinterior\tjunction\tshared\toutside\ttop10\ttop100\ttop1000"]
grp = {}
for r in contr:
    grp.setdefault((r["contrast"], r["profile"]), {})[r["category"]] = r
for (c, p), m in grp.items():
    any_r = next(iter(m.values()))
    out.append("\t".join([c, p, f"{float(any_r['delta_loss']):.4f}", f"{float(any_r['reconstructed']):.4f}", any_r["kmers_differing"]] +
                         [f"{float(m[x]['category_delta']):.3f} ({m[x]['category_kmers']})" for x in ("target_interior_only", "target_junction_only", "shared", "outside_only")] +
                         [f"{float(any_r['top10_share']):.3f}", f"{float(any_r['top100_share']):.3f}", f"{float(any_r['top1000_share']):.3f}"]))
out += ["", "== unchanged active catalogue under each profile", "profile\tbest_group_units\tbest_group_size\treported_3+1_rank\tfloor_2+2_rank\tfirst_2+2_group\tbest_minus_reported"]
verdict_rows = {}
for p in sorted({r["profile"] for r in cat}):
    rows = [r for r in cat if r["profile"] == p]
    best = [r for r in rows if r["group_rank"] == "1"]
    rk = lambda a, b: next(r["group_rank"] for r in rows if r["allele_h1"] == a and r["allele_h2"] == b)
    first22 = next((r["group_rank"] for r in rows if sorted((U(r["allele_h1"]), U(r["allele_h2"]))) == ["2", "2"]), "absent")
    rep = next(r for r in rows if r["is_caller"] == "true")
    verdict_rows[p] = (int(rk("41", "47")), int(rk("37", "46")), first22)
    out.append("\t".join([p, ",".join(sorted({f"{U(r['allele_h1'])}+{U(r['allele_h2'])}" for r in best})), str(len(best)), rk("41", "47"), rk("37", "46"), first22,
                          f"{-float(rep['loss_minus_best']):.3f}"]))
# fragment scores reused
fs = {r["label"]: r for r in T(f"{m6a}/fragscore.c4_loo_na19650_10/fs.scores.tsv")}
es = {r["label"]: r for r in T(f"{m1e}/score.c4_10/score.scores_truth_free.tsv")}
out += ["", "== fragment model (reused; likelihood = fragment + exposure; explicit M1-E sequences have no prior)", "label\tsource\tfragment\texposure\tlikelihood"]
for lab, src in (("g1_m1_caller", "M6-A score-set (reported 3+1)"), ("g4_m1", "M6-A score-set (floor 2+2)")):
    r = fs[lab]
    out.append(f"{lab}\t{src}\t{r['fragment_term']}\t{r['exposure_term']}\t{float(r['fragment_term']) + float(r['exposure_term']):.4f}")
for lab, r in es.items():
    out.append(f"m1e_{lab}\tM1-E score-set\t{r['fragment_term']}\t{r['exposure_term']}\t{float(r['fragment_term']) + float(r['exposure_term']):.4f}")
# truth_a/b vs M1-E truthblock bytes
def fa(p):
    d, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = []
        elif n:
            d[n].append(l)
    return {k: "".join(v) for k, v in d.items()}
tf = fa(f"{m6a}/blockseq.c4_loo_na19650_10/seq.sv:7.truth.fa")
mf = fa(f"{m1e}/construct.c4_10/candidates.fa")
tb = (mf["truthblock.1"][89544:154996], mf["truthblock.2"][89547:154996])
out.append(f"truth_a == M1-E truthblock (truth1 on h1, truth2 on h2): {tb == (tf['truth1'], tf['truth2'])}; truth_b == truthblock: {tb == (tf['truth2'], tf['truth1'])}")
# causal reading
real, o1, o2 = verdict_rows.get("real"), verdict_rows.get("oracle_truth_block_in_context"), verdict_rows.get("oracle_complete_true_diplotype")
d = lambda c, p: float(next(r["delta_loss"] for r in contr if r["contrast"] == c and r["profile"] == p))
reasons = []
r_real = d("reported - floor2x2", "real")
r_o1 = d("reported - floor2x2", "oracle_truth_block_in_context")
r_o2 = d("reported - floor2x2", "oracle_complete_true_diplotype")
if r_real < 0 and r_o1 > 0:
    reasons.append("SAMPLING_OR_NORMALIZATION")
if r_o1 < 0:
    reasons.append("RESIDUAL_CONTENT_MISMATCH")
if (r_o1 < 0) != (r_o2 < 0):
    reasons.append("CONTEXT_EFFECT")
matched = [r for r in contr if r["contrast"] == "m1e_e13_none - m1e_a22_none" and r["profile"] == "real"]
if matched and matched[0]["kmers_differing"] == "0":
    reasons.append("COUNT_INDISTINGUISHABLE(content-matched e13_none vs a22_none)")
out += ["", f"reported - floor2x2 loss: real {r_real:.3f}; noiseless truth-in-context {r_o1:.3f}; noiseless complete truth {r_o2:.3f}",
        "CAUSAL\t" + (";".join(reasons) if reasons else "UNRESOLVED")]
s = KV(A + ".summary.tsv")
out += ["", "== normalization and run"] + [f"{k}\t{v}" for k, v in s.items()]
open(prefix + ".tsv", "w").write("\n".join(out) + "\n")
print("\n".join(out))
