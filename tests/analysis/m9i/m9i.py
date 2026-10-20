#!/usr/bin/env python3
"""M9-I representation-aware compound targets (read-only prototype on saved outputs; frozen with DESIGN.md).

Target T for a run and an interval I = [i..j] of consecutive blocks: the UNORDERED pair of concatenated interval sequences
of the reported labelling (allele bytes of blocks i..j, per homologue). Invariant to homologue exchange and to synonymous
block labellings (different labels spelling the same bytes); sensitive to any genuine byte change and to a different
within-interval arrangement. An EXACT-SEQUENCE target: it tolerates no SNP or indel and is not structural equivalence.

Pool-conditional bounds over the explored sequence groups (mosaic_scores; members = every recorded labelling of the group in
mosaic_candidates), each group's candidate_set_probability counted ONCE:
  support_lower = mass of groups whose recorded projections ALL equal T
  support_upper = mass of groups with AT LEAST ONE recorded projection equal to T
  error_lower = mass of groups with NO projection equal to T      (= 1 - support_upper, summed directly)
  error_upper = mass of groups with SOME projection unequal to T  (= 1 - support_lower, summed directly)
A group with several distinct projections is reported as projection-ambiguous; its mass is never split or assigned.
These bounds cover explored candidates and recorded labellings only; they certify nothing about unsearched alternatives.

  m9i.py gates <out>      synthetic fixtures against an independent oracle
  m9i.py analyse <out>    the saved ACOT (M9-C/G/H) and CYP2D6 cases
"""
import csv
import hashlib
import os
import subprocess
import sys

mode, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
R = "/Users/davide.bolognini/panvar/tests/results"
WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
OUT = f"{R}/m9i_compound_targets_20261020/runs"
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l and not l.startswith("#")} if os.path.exists(p) else {}
labs = lambda s: [int(x) for x in s.split(",")]
H = lambda s: hashlib.sha256(s.encode()).hexdigest()


def write(name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


# ---------------------------------------------------------------------------------------------------------------- core
class Projector:
    """Compound projection of a labelling over interval [i..j], memoised per distinct label tuple."""
    def __init__(self, blocks, i, j):
        self.blocks, self.i, self.j, self.memo = blocks, i, j, {}

    def hap(self, labels):
        key = tuple(labels[self.i:self.j + 1])
        h = self.memo.get(key)
        if h is None:
            h = self.memo[key] = H("".join(self.blocks[b][a] for b, a in zip(range(self.i, self.j + 1), key)))
        return h

    def pair(self, h1, h2):
        a, b = self.hap(h1), self.hap(h2)
        return (a, b) if a <= b else (b, a)


def bounds(groups, proj, target):
    """groups: list of (probability, [ (hap1 labels, hap2 labels), ... ]). Returns the bound record."""
    sl = su = el = eu = total = 0.0
    ambiguous = absent_members = 0
    for p, members in groups:
        total += p
        if not members:
            absent_members += 1
            el += p
            eu += p
            continue
        projs = {proj.pair(h1, h2) for h1, h2 in members}   # duplicates of a projection collapse; mass counted once per group
        ambiguous += len(projs) > 1
        hit = target in projs
        miss = any(x != target for x in projs)
        if hit and not miss:
            sl += p
        if hit:
            su += p
        if not hit:
            el += p
        if miss:
            eu += p
    return dict(total=total, support_lower=sl, support_upper=su, error_lower=el, error_upper=eu, ambiguous_groups=ambiguous,
                groups_without_labellings=absent_members, competitor="NO_EVALUATED_ALTERNATIVE" if eu == 0.0 else "present")


def oracle(groups, blocks, i, j, target_seqs):
    """Independent: spell whole haplotypes, slice the interval by cumulative offsets, compare unordered byte pairs."""
    def interval(labels):
        offs, pos = [], 0
        for b, a in enumerate(labels):
            offs.append(pos)
            pos += len(blocks[b][a])
        full = "".join(blocks[b][a] for b, a in enumerate(labels))
        end = offs[j] + len(blocks[j][labels[j]])
        return full[offs[i]:end]
    tgt = sorted(target_seqs)
    sl = su = 0.0
    for p, members in groups:
        hits = [sorted([interval(h1), interval(h2)]) == tgt for h1, h2 in members]
        if hits and all(hits):
            sl += p
        if any(hits):
            su += p
    return sl, su


# ---------------------------------------------------------------------------------------------------------------- gates
if mode == "gates":
    res = []
    ok = lambda n, c, d="": res.append(dict(gate=n, result="PASS" if c else "FAIL", detail=d))
    # blocks: 0 flank, 1-2 the compound interval (synonymous split: "AC"+"G" == "A"+"CG"), 3 an ordinary block
    blocks = [["TT"], ["AC", "A", "AG"], ["G", "CG", "T"], ["CCC", "CCA"]]
    I, J = 1, 2
    proj = Projector(blocks, I, J)
    tgt = proj.pair([0, 0, 0, 0], [0, 2, 2, 0])          # ACG | AGT
    tseq = ["ACG", "AGT"]
    # G1 two labellings of one group spell the same compound sequence (synonymous split)
    g = [(0.6, [([0, 0, 0, 0], [0, 2, 2, 0]), ([0, 1, 1, 0], [0, 2, 2, 0])]), (0.4, [([0, 2, 1, 0], [0, 2, 2, 0])])]
    r = bounds(g, proj, tgt)
    ok("G1 synonymous labellings are one projection (support 0.6, not ambiguous)", abs(r["support_lower"] - 0.6) < 1e-12 and abs(r["support_upper"] - 0.6) < 1e-12
       and r["ambiguous_groups"] == 0 and (r["support_lower"], r["support_upper"]) == oracle(g, blocks, I, J, tseq), str(r))
    # G2 a real sequence difference is an error
    ok("G2 real difference counted as error", abs(r["error_lower"] - 0.4) < 1e-12 and abs(r["error_upper"] - 0.4) < 1e-12, str(r))
    # G3 homologue exchange and allele-index permutation
    swapped = [(p, [(b, a) for a, b in m]) for p, m in g]
    r3 = bounds(swapped, Projector(blocks, I, J), tgt)
    perm1 = [2, 0, 1]   # new index k holds old allele perm1[k] at block 1
    blocks_p = [blocks[0], [blocks[1][x] for x in perm1], blocks[2], blocks[3]]
    inv = {o: n for n, o in enumerate(perm1)}
    gp = [(p, [([a[0], inv[a[1]], a[2], a[3]], [b[0], inv[b[1]], b[2], b[3]]) for a, b in m]) for p, m in g]
    pp = Projector(blocks_p, I, J)
    r3b = bounds(gp, pp, pp.pair([0, inv[0], 0, 0], [0, inv[2], 2, 0]))
    ok("G3 homologue exchange and allele-index permutation invariance", all(abs(r3[k] - r[k]) < 1e-12 and abs(r3b[k] - r[k]) < 1e-12
                                                                      for k in ("support_lower", "support_upper", "error_lower", "error_upper")))
    # G4 mixed projections: lower < upper, mass not split
    g4 = [(0.5, [([0, 0, 0, 0], [0, 2, 2, 0]), ([0, 0, 2, 0], [0, 2, 2, 0])]), (0.5, [([0, 0, 0, 0], [0, 2, 2, 0])])]
    r4 = bounds(g4, proj, tgt)
    ok("G4 mixed projections give an interval, mass not split", abs(r4["support_lower"] - 0.5) < 1e-12 and abs(r4["support_upper"] - 1.0) < 1e-12
       and r4["ambiguous_groups"] == 1 and (r4["support_lower"], r4["support_upper"]) == oracle(g4, blocks, I, J, tseq), str(r4))
    # G5 duplicate synonymous labellings create no mass
    g5 = [(0.7, [([0, 0, 0, 0], [0, 2, 2, 0])] * 3 + [([0, 1, 1, 0], [0, 2, 2, 0])] * 2), (0.3, [([0, 2, 1, 0], [0, 2, 2, 0])])]
    r5 = bounds(g5, proj, tgt)
    ok("G5 duplicate synonymous labellings do not create mass", abs(r5["support_upper"] - 0.7) < 1e-12 and abs(r5["total"] - 1.0) < 1e-12, str(r5))
    # G6 absent comparison: no group carries the target; and no evaluated alternative when every group supports it
    r6 = bounds([(1.0, [([0, 2, 1, 0], [0, 2, 1, 0])])], proj, tgt)
    r6b = bounds([(1.0, [([0, 0, 0, 0], [0, 2, 2, 0])])], proj, tgt)
    ok("G6 absent target -> support 0; all-supporting pool -> NO_EVALUATED_ALTERNATIVE",
       r6["support_upper"] == 0.0 and r6["error_lower"] == 1.0 and r6b["competitor"] == "NO_EVALUATED_ALTERNATIVE")
    # G7 an ordinary single block reproduces the existing per-block calculation (genotype_sequence_report.cpp: lower = groups
    # all of whose explored genotypes differ, upper = groups with any differing genotype; genotype = unordered allele pair)
    g7 = [(0.5, [([0, 0, 0, 0], [0, 2, 2, 1])]), (0.3, [([0, 0, 0, 1], [0, 2, 2, 0]), ([0, 0, 0, 1], [0, 2, 2, 1])]), (0.2, [([0, 0, 0, 1], [0, 2, 2, 1])])]
    p7 = Projector(blocks, 3, 3)
    r7 = bounds(g7, p7, p7.pair([0, 0, 0, 0], [0, 2, 2, 1]))
    def existing(groups, b, gt):
        lo = up = 0.0
        for p, m in groups:
            gts = {tuple(sorted((h1[b], h2[b]))) for h1, h2 in m}
            if all(x != gt for x in gts):
                lo += p
            if any(x != gt for x in gts):
                up += p
        return lo, up
    ok("G7 single-block target reproduces the existing error-mass bounds", (r7["error_lower"], r7["error_upper"]) == existing(g7, 3, (0, 1)),
       f"{(r7['error_lower'], r7['error_upper'])} vs {existing(g7, 3, (0, 1))}")
    # normalisation and ordering on every fixture
    allr = [r, r3, r3b, r4, r5, r6, r6b, r7]
    ok("normalisation (support_upper + error_lower = total; support_lower + error_upper = total) and lower <= upper",
       all(abs(x["support_upper"] + x["error_lower"] - x["total"]) < 1e-12 and abs(x["support_lower"] + x["error_upper"] - x["total"]) < 1e-12
           and x["support_lower"] <= x["support_upper"] and x["error_lower"] <= x["error_upper"] for x in allr))
    write("gates", res)
    print(open(f"{out}/gates.tsv").read())
    sys.exit(0 if all(x["result"] == "PASS" for x in res) else 1)

# ---------------------------------------------------------------------------------------------------------------- inventory
if mode == "inventory":
    Hm = f"{R}/m9h_jp11_validation_20261020/runs/models"
    code = subprocess.run([sys.executable, f"{WT}/tests/analysis/cyp2d6_pilot/degeneracy_inventory.py", f"{out}/degeneracy_m9h.tsv", f"{out}/degeneracy_m9h_summary.tsv",
                           f"M9H.LOO={Hm}/LOO", f"M9H.Z1={Hm}/Z1", f"M9H.Z2={Hm}/Z2"]).returncode
    rows = []
    for src, path in (("M9-D ACOT models", f"{R}/m9d_attribution_20261020/runs/degeneracy/degeneracy_summary.tsv"),
                      ("M9-I re-run on M9-H ACOT models", f"{out}/degeneracy_m9h_summary.tsv"),
                      ("CYP2D6 pilot (sealed)", f"{R}/cyp2d6_shortread_pilot_v1/degeneracy/degeneracy_summary.tsv")):
        for r in T(path):
            rows.append(dict(source=src, **r))
    write("inventory", rows)
    print(open(f"{out}/inventory.tsv").read())
    sys.exit(code)

# ---------------------------------------------------------------------------------------------------------------- cases
C, G, Hh, D = f"{R}/m9c_acot_reliability_20261020/runs", f"{R}/m9g_junction_integration_20261020/runs", f"{R}/m9h_jp11_validation_20261020/runs", f"{R}/m9d_attribution_20261020/runs"
P, M = f"{R}/cyp2d6_shortread_pilot_v1", f"{R}/m1b_structural_20261019/runs"


def cases():
    out = []
    donors_c = {f"d{i}": l.split("\t")[1] for i in range(1, 7) for l in open(f"{C}/prepare/prepare.tsv") if l.startswith(f"donor_{i}\t")}
    donors_h = {}
    for l in open(f"{Hh}/prepare/prepare.tsv"):
        f = l.rstrip("\n").split("\t")
        if f[0].startswith("family_"):
            donors_h[f[0].split("_")[1]] = f[2].split()[1]
    for c in ["LOO.d1", "LOO.d2", "LOO.d3", "LOO.d4", "LOO.d5", "LOO.d6", "LZO.d1", "LZO.d2"]:
        out.append(dict(case=f"M9-C {c}", locus="ACOT", prefix=f"{C}/call.{c}/call", dump=f"{D}/models/{'LOO' if c.startswith('LOO') else c}.block_alleles.tsv",
                        donor=donors_c[c.split(".")[1]], truth=f"truth.ACOT.{'LOO' if c.startswith('LOO') else 'C' + c}.{donors_c[c.split('.')[1]]}", eval=f"{C}/eval.{c}/eval.blocks.tsv", v2=True))
    for arm in ("off", "on"):
        for c in ["LOO.d1", "LOO.d3", "LOO.d6", "LZO.d1", "LZO.d2"]:
            out.append(dict(case=f"M9-G {arm} {c}", locus="ACOT", prefix=f"{G}/call.{arm}.{c}/call", dump=f"{D}/models/{'LOO' if c.startswith('LOO') else c}.block_alleles.tsv",
                            donor=donors_c[c.split(".")[1]], truth=f"truth.ACOT.{'LOO' if c.startswith('LOO') else 'C' + c}.{donors_c[c.split('.')[1]]}", eval=f"{G}/eval.{arm}.{c}/eval.blocks.tsv", v2=True))
    for arm in ("off", "on"):
        for c in ["L1", "L2", "L3", "L4", "Z1", "Z2"]:
            out.append(dict(case=f"M9-H {arm} {c}", locus="ACOT", prefix=f"{Hh}/call.{arm}.{c}/call", dump=f"{Hh}/models/{'LOO' if c[0] == 'L' else c}.block_alleles.tsv",
                            donor=donors_h[c[1]], truth=f"truth.ACOT.{'LOO' if c[0] == 'L' else 'H' + c}.{donors_h[c[1]]}", eval=f"{Hh}/eval.{arm}.{c}/eval.blocks.tsv", v2=True))
    out += [dict(case="CYP pilot lzo.HG02769.30 (label-ambiguity no-call)", locus="CYP2D6", prefix=f"{P}/data/runs_a5/lzo.HG02769.30", dump=f"{P}/degeneracy/models/original.block_alleles.tsv",
                 donor="HG02769", truth="truth.CYP.original.HG02769", eval="", v2=False),
            dict(case="CYP pilot loo.HG02953.30 (label-ambiguity no-call)", locus="CYP2D6", prefix=f"{P}/data/runs_a5/loo.HG02953.30", dump=f"{P}/degeneracy/models/loo.HG02953.block_alleles.tsv",
                 donor="HG02953", truth="truth.CYP.loo.HG02953", eval="", v2=False),
            dict(case="CYP M1-B cyp_5 (NA19650 LOO 5x)", locus="CYP2D6", prefix=f"{M}/call.cyp_5/call", dump=f"{P}/degeneracy/models/loo.NA19650.block_alleles.tsv",
                 donor="NA19650", truth="truth.CYP.loo.NA19650", eval=f"{M}/eval.cyp_5/eval.blocks.tsv", v2=True),
            dict(case="CYP M1-B cyp_30 (NA19650 LOO 30x)", locus="CYP2D6", prefix=f"{M}/call.cyp_30/call", dump=f"{P}/degeneracy/models/loo.NA19650.block_alleles.tsv",
                 donor="NA19650", truth="truth.CYP.loo.NA19650", eval=f"{M}/eval.cyp_30/eval.blocks.tsv", v2=True)]
    return out


INTERVAL = {"ACOT": ("sv:5", "ctx:6"), "CYP2D6": ("sv:3", "sv:5")}
CONTROL = {"ACOT": "sv:13", "CYP2D6": "sv:9"}


def dump(path):
    blocks = []
    for r in T(path):
        b = int(r["block_index"])
        while len(blocks) <= b:
            blocks.append([])
        blocks[b].append(r["sequence"])
    ids = {}
    for r in T(path):
        ids[r["block_id"]] = int(r["block_index"])
    return blocks, ids


def truth_bytes(name, i_id, j_id, ids):
    seqs = ["", ""]
    for b in range(ids[i_id], ids[j_id] + 1):
        bid = next(k for k, v in ids.items() if v == b)
        recs, cur = {}, None
        for line in open(f"{OUT}/{name}/bs.{bid}.truth.fa"):
            if line.startswith(">"):
                cur = line[1:].split()[0]
                recs[cur] = ""
            else:
                recs[cur] += line.strip()
        for k in (0, 1):
            seqs[k] += recs[f"truth{k + 1}"]
    return seqs


def agree(computed, printed):
    """The saved bounds are printed to 6 significant digits."""
    if printed in (None, "NA"):
        return False
    v = float(printed)
    return abs(computed - v) <= 1e-5 * max(abs(computed), abs(v)) or (v == 0.0 and computed == 0.0)


if mode == "analyse":
    rows, controls, l4, fa, pairs = [], [], [], [], []
    for cs in cases():
        blocks, ids = dump(cs["dump"])
        i, j = ids[INTERVAL[cs["locus"]][0]], ids[INTERVAL[cs["locus"]][1]]
        prob = {r["sequence_candidate_id"]: float(r["candidate_set_probability"]) for r in T(f"{cs['prefix']}.mosaic_scores.tsv")}
        members = {}
        for r in T(f"{cs['prefix']}.mosaic_candidates.tsv"):
            members.setdefault(r["sequence_candidate_id"], []).append((labs(r["hap1_alleles"]), labs(r["hap2_alleles"])))
        groups = [(p, members.get(g, [])) for g, p in prob.items()]
        sl = KV(f"{cs['prefix']}.sequence_locus.tsv")
        if cs["v2"] and sl.get("representative_labelling_hap1"):
            rep = (labs(sl["representative_labelling_hap1"]), labs(sl["representative_labelling_hap2"]))
            top = [r["group"] for r in T(f"{cs['prefix']}.sequence_alternatives.tsv") if r.get("in_decision_set") == "yes"]
            rep_source = "v2 representative labelling"
        else:   # legacy: the rank-1 group's members (the equivalence set); any member labelling
            eq = T(f"{cs['prefix']}.equivalence_set.tsv")
            rep = (labs(eq[0]["hap1_alleles"]), labs(eq[0]["hap2_alleles"]))
            top = sorted({r["sequence_candidate_id"] for r in eq})
            rep_source = "legacy equivalence set (rank-1 group), first recorded labelling"
        proj = Projector(blocks, i, j)
        target = proj.pair(*rep)
        b = bounds(groups, proj, target)
        rep_group = sl.get("representative_group") if cs["v2"] and sl.get("representative_group") in members else next((g for g, m in members.items() if any(key == tuple(map(tuple, rep)) or key == tuple(map(tuple, rep[::-1])) for key in
                                                                 [(tuple(h1), tuple(h2)) for h1, h2 in m])), None)
        rep_projs = {proj.pair(h1, h2) for h1, h2 in members.get(rep_group, [])}
        # block-label view inside the interval: distinct label pairs among the representative group's labellings
        rep_label_variants = {tuple(tuple(sorted((h1[k], h2[k]))) for k in range(i, j + 1)) for h1, h2 in members.get(rep_group, [])}
        # existing block-level bounds inside the interval (v2) or legacy posterior
        sb = {r["block_id"]: r for r in T(f"{cs['prefix']}.sequence_blocks.tsv")}
        bg = {r["block_id"]: r for r in T(f"{cs['prefix']}.block_genotypes.tsv")}
        ev = {r["block_id"]: r for r in T(cs["eval"])} if cs["eval"] else {}
        inside = [k for k, v in ids.items() if i <= v <= j]
        inside.sort(key=lambda k: ids[k])
        blk = []
        for k in inside:
            e = ev.get(k)
            exact = "NA" if not e or e.get("unordered_distance", "NA") == "NA" else e["unordered_distance"] == "0"
            blk.append(f"{k}:{(sb.get(k) or {}).get('block_status', 'legacy')}:{(sb.get(k) or {}).get('error_mass_upper', (bg.get(k) or {}).get('conditional_genotype_posterior', 'NA'))}"
                       f":exact={exact}:d={e.get('unordered_distance', 'NA') if e else 'NA'}")
        # compound truth distance: the representative's interval vs truth (edit distance, homologue exchange)
        rep_seqs = ["".join(blocks[k][h[k]] for k in range(i, j + 1)) for h in rep]
        tb = truth_bytes(cs["truth"], INTERVAL[cs["locus"]][0], INTERVAL[cs["locus"]][1], ids) if os.path.exists(f"{OUT}/{cs['truth']}") else None
        tag = cs["case"].replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_")
        if tb:
            fa += [(f"{tag}_h1", rep_seqs[0]), (f"{tag}_h2", rep_seqs[1]), (f"{tag}_t1", tb[0]), (f"{tag}_t2", tb[1])]
            pairs += [(f"{tag}|{a}|{t}", f"{tag}_{a}", f"{tag}_{t}") for a in ("h1", "h2") for t in ("t1", "t2")]
        rows.append(dict(case=cs["case"], locus=cs["locus"], interval="..".join(INTERVAL[cs["locus"]]), representative_source=rep_source,
                         groups=len(prob), labellings=sum(len(m) for m in members.values()), total_mass=round(b["total"], 9),
                         rep_group_labellings=len(members.get(rep_group, [])), rep_group_block_label_variants=len(rep_label_variants),
                         rep_group_compound_projections=len(rep_projs), projection_ambiguous_groups=b["ambiguous_groups"],
                         compound_support_lower=b["support_lower"], compound_support_upper=b["support_upper"],
                         compound_error_lower=b["error_lower"], compound_error_upper=b["error_upper"], competitor=b["competitor"],
                         compound_exact=(sorted(rep_seqs) == sorted(tb)) if tb else "no truth bytes", block_view=" ; ".join(blk)))
        # control: an ordinary single block -> must reproduce the existing per-block bounds (v2)
        cb = CONTROL[cs["locus"]]
        pc = Projector(blocks, ids[cb], ids[cb])
        rc = bounds(groups, pc, pc.pair(*rep))
        ex = sb.get(cb, {})
        controls.append(dict(case=cs["case"], block=cb, compound_error_lower=rc["error_lower"], compound_error_upper=rc["error_upper"],
                             existing_error_mass_lower=ex.get("error_mass_lower", "NA (legacy)"), existing_error_mass_upper=ex.get("error_mass_upper", "NA (legacy)"),
                             agrees=agree(rc["error_lower"], ex.get("error_mass_lower")) and agree(rc["error_upper"], ex.get("error_mass_upper"))
                             if ex.get("error_mass_lower") not in (None, "NA") else "NA",
                             relative_difference=(abs(rc["error_upper"] - float(ex["error_mass_upper"])) / float(ex["error_mass_upper"]))
                             if ex.get("error_mass_upper") not in (None, "NA") and float(ex["error_mass_upper"]) > 0 else "NA"))
        # M9-H L4: every recorded labelling of every top-set group
        if cs["case"] in ("M9-H on L4", "M9-H off L4") and tb:
            for g in (top or [rep_group]):
                for h1, h2 in members.get(g, []):
                    s = sorted(["".join(blocks[k][h[k]] for k in range(i, j + 1)) for h in (h1, h2)])
                    l4.append(dict(case=cs["case"], group=g, sv5=f"{h1[ids['sv:5']]}/{h2[ids['sv:5']]}", ctx6=f"{h1[ids['ctx:6']]}/{h2[ids['ctx:6']]}",
                                   exact_interval_sequence=s == sorted(tb)))
    with open(f"{out}/intervals.fa", "w") as f:
        for n, s in fa:
            f.write(f">{n}\n{s}\n")
    with open(f"{out}/pairs.tsv", "w") as f:
        for lab, a, b2 in pairs:
            f.write(f"{lab}\t{out}/intervals.fa\t{a}\t{out}/intervals.fa\t{b2}\n")
    subprocess.run([f"{WT}/build/genotype_probe", "edit-distances", "--pairs", f"{out}/pairs.tsv", "--out", f"{out}/distances.tsv", "--threads", "4"], check=True)
    dist = {r["label"]: int(r["distance"]) for r in T(f"{out}/distances.tsv")}
    for r in rows:
        tag = r["case"].replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_")
        k = lambda a, t: dist.get(f"{tag}|{a}|{t}")
        r["compound_distance"] = min(k("h1", "t1") + k("h2", "t2"), k("h1", "t2") + k("h2", "t1")) if k("h1", "t1") is not None else "NA"
    write("compound", rows)
    write("controls", controls)
    write("l4_top_set", l4)
    for n in ("compound", "controls", "l4_top_set"):
        print(f"## {n}\n" + open(f"{out}/{n}.tsv").read())
    sys.exit(0)
