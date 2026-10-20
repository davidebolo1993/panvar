#!/usr/bin/env python3
"""M9-B focused fixtures and saved-evidence extraction (frozen with DESIGN.md). Read-only; no caller, no genotyping.

  m9b.py fixtures <out>   joint vs marginal agreement through M9-A's frozen core aggregation (exec'd from
                          tests/analysis/m9a/audit.py, unchanged): two block marginals of 0.9 with joint agreement 0.8
                          (and 0.9), never the product; plus a block-label vs sequence fixture
  m9b.py evidence <out>   the three illustrative failures from saved outputs only: HG01975 ACOT sv:5 (M8-D r0.01/r0.1),
                          C4 NA19650 LOO sv:7 and CYP2D6 NA19650 LOO sv:9 (M1-B fresh-caller runs); each bound to its graph,
                          reference, curation and export settings; truth enters only from the archived evaluator outputs
"""
import csv
import math
import os
import sys

mode, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))
R = "/Users/davide.bolognini/panvar/tests/results"
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l and not l.startswith("#")}


def write(name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def walks(sources):
    """allele*_sources of .genotypes.tsv: comma-joined source ids, each 'bypass', 'empty' or 'walk:' + comma-joined
    oriented nodes (genotype_command.cpp sources_of; graph_utils.cpp build_walk_signature)."""
    out, cur = [], None
    for tok in sources.split(","):
        if tok.startswith("walk:"):
            cur = [tok[5:]]
            out.append(cur)
        elif tok in ("bypass", "empty"):
            out.append([tok])
            cur = None
        else:
            cur.append(tok)
    return out


def walk_shape(w):
    if w in (["bypass"], ["empty"]):
        return dict(steps=0, max_node_visits=0, reverse_steps=0)
    nodes = [x[:-1] for x in w]
    return dict(steps=len(w), max_node_visits=max(nodes.count(n) for n in set(nodes)), reverse_steps=sum(x.endswith("-") for x in w))


def label_structure(d, block, a1, a2):
    """Per called label: how many distinct source walks spell it, and whether those walks differ in shape."""
    g = next((x for x in T(f"{d}/call.genotypes.tsv") if x["block_id"] == block), None)
    if g is None:
        return dict(label_walks="block not in .genotypes.tsv")
    res = {}
    for k, lab in (("allele1", a1), ("allele2", a2)):
        ws = walks(g[f"{k}_sources"])
        shapes = [walk_shape(w) for w in ws]
        res[f"{k}_label"] = lab
        res[f"{k}_source_walks"] = len(ws)
        res[f"{k}_walk_shapes"] = ";".join(f"{x['steps']}st/{x['max_node_visits']}v/{x['reverse_steps']}r" for x in shapes)
        res[f"{k}_walks_differ_in_nodes"] = len({tuple(w) for w in ws}) > 1
    return res


def m9a_core():
    """M9-A's frozen aggregation (build / target / audit_pool), executed from its committed source without change."""
    src = open(os.path.join(HERE, "..", "m9a", "audit.py")).read()
    core = src[src.index("def lse(xs):"):src.index("\n# ", src.index("def contenders("))]
    ns = {"math": math}
    exec(core, ns)
    return ns


if mode == "fixtures":
    ns = m9a_core()
    build, audit_pool, geno = ns["build"], ns["audit_pool"], ns["geno"]
    G = lambda i, logp, m: dict(id=i, score=logp, fragment=logp, exposure=0.0, prior=0.0, bound=0.0, members=m)
    res = []
    ok = lambda name, cond, detail: res.append(dict(fixture=name, result="PASS" if cond else "FAIL", detail=detail))
    # J1: two blocks, marginals 0.9 each, joint 0.8 (the Frechet lower bound for these marginals); product would be 0.81
    win = ((1, 2), (3, 4))
    B = geno(win)
    pool = build([G(0, math.log(0.8), [win]), G(1, math.log(0.1), [((1, 2), (3, 5))]), G(2, math.log(0.1), [((9, 2), (3, 4))])   # Amendment 1: disagrees at block 0 only])
    rows, blocks, _ = audit_pool(pool, B, [True, True])
    d1, d2, c = blocks[0]["mass_lower"], blocks[1]["mass_lower"], rows["C reportable-block genotype agreement"]["mass_lower"]
    ok("J1 marginals 0.9/0.9, joint 0.8", abs(d1 - 0.9) < 1e-12 and abs(d2 - 0.9) < 1e-12 and abs(c - 0.8) < 1e-12 and abs(c - d1 * d2) > 0.009,
       f"D1={d1:.6f} D2={d2:.6f} C={c:.6f} product={d1 * d2:.6f} Frechet=[{max(0, d1 + d2 - 1):.2f},{min(d1, d2):.2f}]")
    # J2: same marginals, perfectly dependent errors: joint 0.9 (the upper bound)
    pool = build([G(0, math.log(0.9), [win]), G(1, math.log(0.1), [((9, 9), (3, 5))])   # Amendment 1: disagrees at both blocks])
    rows, blocks, _ = audit_pool(pool, B, [True, True])
    d1, d2, c = blocks[0]["mass_lower"], blocks[1]["mass_lower"], rows["C reportable-block genotype agreement"]["mass_lower"]
    ok("J2 marginals 0.9/0.9, joint 0.9", abs(d1 - 0.9) < 1e-12 and abs(d2 - 0.9) < 1e-12 and abs(c - 0.9) < 1e-12,
       f"D1={d1:.6f} D2={d2:.6f} C={c:.6f} product={d1 * d2:.6f}")
    # J3: C over a reportable subset ignores non-reportable blocks; D at a block is unaffected by the others
    pool = build([G(0, math.log(0.6), [((1, 2, 7), (3, 4, 8))]), G(1, math.log(0.3), [((1, 2, 9), (3, 4, 8))]), G(2, math.log(0.1), [((1, 5, 7), (3, 4, 8))])])
    rows, blocks, _ = audit_pool(pool, geno(((1, 2, 7), (3, 4, 8))), [False, True, False])
    c, b = rows["C reportable-block genotype agreement"]["mass_lower"], rows["B unphased genotype, all blocks"]["mass_lower"]
    ok("J3 C covers only reportable blocks", abs(c - 0.9) < 1e-12 and abs(b - 0.6) < 1e-12 and abs(blocks[1]["mass_lower"] - 0.9) < 1e-12,
       f"C={c:.3f} B={b:.3f} D(reportable)={blocks[1]['mass_lower']:.3f}")
    write("fixtures", res)
    print(open(f"{out}/fixtures.tsv").read())
    sys.exit(0 if all(r["result"] == "PASS" for r in res) else 1)

if mode == "evidence":
    rows = []
    # ---- 1. HG01975 ACOT sv:5 (block index 5), M8-D r0.01 and r0.1
    for r in ("r0.01", "r0.1"):
        d = f"{R}/m8d_hg01975_loo_20261020/runs/call.HG01975.LOO.{r}"
        bg = {x["block_id"]: x for x in T(f"{d}/call.block_genotypes.tsv")}
        ev = {x["block_id"]: x for x in T(f"{R}/m8d_hg01975_loo_20261020/runs/eval.{r}/eval.blocks.tsv")}
        st = KV(f"{d}/call.mosaic_status.tsv")
        cmd = open(f"{d}/run.log").readline().split()
        arg = lambda k: cmd[cmd.index(k) + 1] if k in cmd else ""
        mvb = next(x.split("\t")[1] for x in open(f"{d}/call.settings.tsv") if x.startswith("min_variant_bp\t"))
        b = bg["sv:5"]
        e = ev["sv:5"]
        idx = list(bg).index("sv:5")
        truth = tuple(sorted((int(e["truth1_closest_allele"]), int(e["truth2_closest_allele"]))))
        # evaluated hypotheses carrying the truth-closest unordered pair at sv:5 (oracle tier: evaluation only)
        cand = {}
        for c in T(f"{d}/call.mosaic_candidates.tsv"):
            h1, h2 = c["hap1_alleles"].split(","), c["hap2_alleles"].split(",")
            cand[c["sequence_candidate_id"]] = (tuple(sorted((int(h1[idx]), int(h2[idx])))), int(h1[idx]), int(h2[idx]))
        sc = T(f"{d}/call.mosaic_scores.tsv")
        w = max(sc, key=lambda x: float(x["score"]))
        carrying = [x for x in sc if cand[x["sequence_candidate_id"]][0] == truth]
        best = max(carrying, key=lambda x: float(x["score"])) if carrying else None
        lik = lambda x: float(x["fragment_term"]) + float(x["exposure_term"])
        rows.append(dict(case=f"ACOT HG01975 LOO {r} sv:5", graph=arg("--gfa") or arg("--graph"), reference=arg("--reference-path"),
                         curation=arg("--panel-curation"), profile=arg("--profile"), preset=arg("--preset") or "none (profile composition)",
                         min_variant_bp=mvb, model_fingerprint=st["model_fingerprint"],
                         block_kind=b["kind"], reportable=b["reportable"], called=f"{b['allele1']}/{b['allele2']}", called_bp=f"{b['allele1_bp']}/{b['allele2_bp']}",
                         truth_closest=f"{e['truth1_closest_allele']}/{e['truth2_closest_allele']}", truth_bp=f"{e['truth1_bp']}/{e['truth2_bp']}",
                         floor=e["availability_floor"], excess=e["excess"], block_posterior=b["conditional_genotype_posterior"], GQ=b["conditional_GQ"],
                         genotype_support=b["genotype_support"], projection_ambiguous=b["projection_ambiguous"],
                         evaluated_hypotheses_with_truth_pair=len(carrying), truth_pair_mass=sum(float(x["candidate_set_probability"]) for x in carrying),
                         best_truth_pair_gap_total=float(w["score"]) - float(best["score"]) if best else "",
                         best_truth_pair_gap_likelihood=lik(w) - lik(best) if best else "",
                         best_truth_pair_gap_prior=float(w["sequence_pair_prior"]) - float(best["sequence_pair_prior"]) if best else "",
                         **label_structure(d, "sv:5", b["allele1"], b["allele2"])))
    # ---- 2/3. C4 sv:7 and CYP2D6 sv:9, NA19650 LOO (M1-B fresh-caller runs)
    M = f"{R}/m1b_structural_20261019/runs"
    for case, block in (("c4_10", "sv:7"), ("c4_30", "sv:7"), ("cyp_5", "sv:9"), ("cyp_30", "sv:9")):
        d = f"{M}/call.{case}"
        b = next(x for x in T(f"{d}/call.block_genotypes.tsv") if x["block_id"] == block)
        lg = T(f"{d}/call.locus_genotype.tsv")[0]
        st = KV(f"{d}/call.mosaic_status.tsv")
        s = KV(f"{M}/structure.{case}/{case}.structure.tsv")
        mvb = next(x.split("\t")[1] for x in open(f"{d}/call.settings.tsv") if x.startswith("min_variant_bp\t"))
        ds = [x for x in T(f"{d}/call.dosage_blocks.tsv") if x["block_id"] == block] if os.path.exists(f"{d}/call.dosage_blocks.tsv") else []
        cmd = open(f"{d}/run.log").readline().split()
        arg = lambda k: cmd[cmd.index(k) + 1] if k in cmd else ""
        rows.append(dict(case=f"{case} NA19650 LOO {block}", graph=arg("--graph"), reference=arg("--reference-path"), curation=arg("--panel-curation") or "none (LOO graph)",
                         profile=",".join(cmd[i + 1] for i, x in enumerate(cmd) if x == "--profile"), preset=arg("--preset"), min_variant_bp=mvb, model_fingerprint=st["model_fingerprint"],
                         block_kind=b["kind"], reportable=b["reportable"], called=f"{b['allele1']}/{b['allele2']}", called_bp=f"{b['allele1_bp']}/{b['allele2_bp']}",
                         truth_units=f"{s['truth_truth1']} | {s['truth_truth2']}", called_units=s["fresh_called_units"], structure_verdict=s["fresh_structure_verdict"],
                         excess=s["fresh_target_excess"], block_posterior=b["conditional_genotype_posterior"], GQ=b["conditional_GQ"], genotype_support=b["genotype_support"],
                         quality_filter=b["quality_filter"], projection_ambiguous=b["projection_ambiguous"], run_status=lg["run_status"], phase_support=lg["phase_support"],
                         dosage_mode=lg["dosage_mode"], dosage_audit_status=lg["dosage_audit_status"],
                         dosage_status=ds[0]["dosage_status"] if ds else "no dosage block row", allocation_status=ds[0]["allocation_status"] if ds else "",
                         dosage_features=next((v for v in (ds[0].values() if ds else []) if str(v).startswith("features=")), ""),
                         **label_structure(d, block, b["allele1"], b["allele2"])))
    write("evidence", rows)
    print(open(f"{out}/evidence.tsv").read())
    sys.exit(0)
