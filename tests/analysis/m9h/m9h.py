#!/usr/bin/env python3
"""M9-H analysis (frozen with DESIGN.md).

  m9h.py correction <out>  M9-G closure: unordered same-sequence-pair semantics with fixtures, re-applied to M9-G's saved
                           outputs (historical tables untouched)
  m9h.py models <out>      block-allele dumps (cyp2d6_block_alleles) for the LOO model and both new LZO models
  m9h.py freeze <out>      every prediction of both arms, before any truth-based evaluation
  m9h.py intervals <out>   sv:5..ctx:6 interval sequences (predicted, from the v2 representative labelling; truth, from the
                           block-sequences projection) and exact edit distances (genotype_probe edit-distances)
  m9h.py analyse <out>     paired outcomes, family classification and the predeclared batch reading
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
G9 = f"{R}/m9g_junction_integration_20261020/runs"
OUT = f"{R}/m9h_jp11_validation_20261020/runs"
ACOT_GFA = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
ACOT_REF = "GRCh38#0#chr14:73416151-73668669"
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l and not l.startswith("#")} if os.path.exists(p) else {}
labs = lambda s: [int(x) for x in s.split(",")]
CELLS = ["L1", "L2", "L3", "L4", "Z1", "Z2"]
LED = {r["unit"]: r for r in T(f"{OUT}/../ledger.tsv")}


def write(name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def same_sequence_pair(a, b):
    """Unordered haplotype-sequence pairs; UNKNOWN unless both are complete predictions with both haplotype digests."""
    def pair(x):
        if not x or x.get("prediction") != "PREDICTION" or not x.get("representative_hap1_sha256") or not x.get("representative_hap2_sha256"):
            return None
        return frozenset([x["representative_hap1_sha256"], x["representative_hap2_sha256"]]) if x["representative_hap1_sha256"] != x["representative_hap2_sha256"] \
            else ("hom", x["representative_hap1_sha256"])
    pa, pb = pair(a), pair(b)
    if pa is None or pb is None:
        return "UNKNOWN"
    return pa == pb


def prep():
    p = {}
    for l in open(f"{OUT}/prepare/prepare.tsv"):
        f = l.rstrip("\n").split("\t")
        if f[0].startswith("family_"):
            p[int(f[0].split("_")[1])] = dict(family=f[1], donor=f[2].split()[1])
    return p


def cell_info(cell):
    p = prep()
    i = int(cell[1])
    return dict(scheme="LOO" if cell[0] == "L" else "LZO", **p[i])


def dump(path):
    blocks = []
    for r in T(path):
        b = int(r["block_index"])
        while len(blocks) <= b:
            blocks.append([])
        blocks[b].append(r["sequence"])
    return blocks


if mode == "correction":
    fx = [("swapped haplotypes", same_sequence_pair({"prediction": "PREDICTION", "representative_hap1_sha256": "a", "representative_hap2_sha256": "b"},
                                                    {"prediction": "PREDICTION", "representative_hap1_sha256": "b", "representative_hap2_sha256": "a"}), True),
          ("different pair", same_sequence_pair({"prediction": "PREDICTION", "representative_hap1_sha256": "a", "representative_hap2_sha256": "b"},
                                                {"prediction": "PREDICTION", "representative_hap1_sha256": "a", "representative_hap2_sha256": "c"}), False),
          ("homozygous vs heterozygous", same_sequence_pair({"prediction": "PREDICTION", "representative_hap1_sha256": "a", "representative_hap2_sha256": "a"},
                                                            {"prediction": "PREDICTION", "representative_hap1_sha256": "a", "representative_hap2_sha256": "b"}), False),
          ("missing", same_sequence_pair({}, {"prediction": "PREDICTION", "representative_hap1_sha256": "a", "representative_hap2_sha256": "b"}), "UNKNOWN"),
          ("incomplete (provisional)", same_sequence_pair({"prediction": "PROVISIONAL", "representative_hap1_sha256": "a", "representative_hap2_sha256": "b"},
                                                          {"prediction": "PREDICTION", "representative_hap1_sha256": "a", "representative_hap2_sha256": "b"}), "UNKNOWN"),
          ("no prediction", same_sequence_pair({"prediction": "NO_PREDICTION"}, {"prediction": "NO_PREDICTION"}), "UNKNOWN")]
    rows = [dict(fixture=n, got=g, expected=e, result="PASS" if g == e else "FAIL") for n, g, e in fx]
    write("fixtures", rows)
    recon = []
    hist = {r["cell"]: r for r in T(f"{G9}/analyse/paired.tsv")}
    for cell in ["LOO.d1", "LOO.d3", "LOO.d6", "LZO.d1", "LZO.d2"]:
        a, b = KV(f"{G9}/call.off.{cell}/call.sequence_locus.tsv"), KV(f"{G9}/call.on.{cell}/call.sequence_locus.tsv")
        new = same_sequence_pair(a, b)
        recon.append(dict(cell=cell, historical_same_sequence_pair=hist[cell]["same_sequence_pair"], corrected=new,
                          conclusion_changed=str(new) != hist[cell]["same_sequence_pair"]))
    write("m9g_reapplied", recon)
    print(open(f"{out}/fixtures.tsv").read())
    print(open(f"{out}/m9g_reapplied.tsv").read())
    sys.exit(0 if all(r["result"] == "PASS" for r in rows) else 1)

if mode == "models":
    P = f"{OUT}/prepare"
    rows = []
    for name, cur in [("LOO", f"{P}/curation.LOO.txt")] + [(f"Z{i}", f"{P}/curation.LZO.{prep()[i]['donor']}.txt") for i in (1, 2)]:
        code = subprocess.run([f"{WT}/build/cyp2d6_block_alleles", ACOT_GFA, ACOT_REF, cur, f"{out}/{name}"], capture_output=True).returncode
        rows.append(dict(model=name, exit=code, fingerprint=KV(f"{out}/{name}.summary.tsv").get("model_fingerprint", "absent")))
    write("models", rows)
    print(open(f"{out}/models.tsv").read())
    sys.exit(0 if all(r["exit"] == 0 for r in rows) else 1)


def arm_dir(arm, cell):
    return f"{OUT}/call.{arm}.{cell}/call"


if mode == "freeze":
    rows = []
    for cell in CELLS:
        for arm in ("off", "on"):
            sl = KV(f"{arm_dir(arm, cell)}.sequence_locus.tsv")
            led = LED.get(f"call.{arm}.{cell}", {})
            rows.append(dict(cell=cell, arm=arm, **cell_info(cell), exit=led.get("exit", "not run"), timed_out=led.get("timed_out", ""),
                             prediction=sl.get("prediction", "absent"), sequence_status=sl.get("sequence_status", "NA"), provisional=sl.get("provisional", "NA"),
                             legacy_run_status=sl.get("legacy_run_status", "NA"), legacy_no_call_reason=sl.get("legacy_no_call_reason", "NA"),
                             hap1=sl.get("representative_hap1_sha256", "NA")[:16], hap2=sl.get("representative_hap2_sha256", "NA")[:16],
                             decision_set_groups=sl.get("decision_set_groups", "NA"),
                             sha_sequence_blocks=hashlib.sha256(open(f"{arm_dir(arm, cell)}.sequence_blocks.tsv", "rb").read()).hexdigest()[:16]
                             if os.path.exists(f"{arm_dir(arm, cell)}.sequence_blocks.tsv") else "absent"))
    write("freeze", rows)
    print(open(f"{out}/freeze.tsv").read())
    sys.exit(0)

if mode == "intervals":
    fasta, pairs, meta = [], [], []
    for cell in CELLS:
        info = cell_info(cell)
        model = "LOO" if info["scheme"] == "LOO" else cell
        blocks = dump(f"{OUT}/models/{model}.block_alleles.tsv")
        fp = KV(f"{OUT}/models/{model}.summary.tsv").get("model_fingerprint")
        tf = {}
        for blk in ("sv:5", "ctx:6"):
            recs, cur = {}, None
            for line in open(f"{OUT}/truth.{cell}/bs.{blk}.truth.fa"):
                if line.startswith(">"):
                    cur = line[1:].split()[0]
                    recs[cur] = ""
                else:
                    recs[cur] += line.strip()
            tf[blk] = recs
        truth = [tf["sv:5"][f"truth{k}"] + tf["ctx:6"][f"truth{k}"] for k in (1, 2)]
        fasta += [(f"{cell}_truth1", truth[0]), (f"{cell}_truth2", truth[1])]
        for arm in ("off", "on"):
            sl = KV(f"{arm_dir(arm, cell)}.sequence_locus.tsv")
            st = KV(f"{arm_dir(arm, cell)}.mosaic_status.tsv")
            if sl.get("prediction") not in ("PREDICTION", "PROVISIONAL"):
                meta.append(dict(cell=cell, arm=arm, status="no prediction"))
                continue
            ids = [r["block_id"] for r in T(f"{arm_dir(arm, cell)}.sequence_blocks.tsv")]
            b5, b6 = ids.index("sv:5"), ids.index("ctx:6")
            haps = [labs(sl["representative_labelling_hap1"]), labs(sl["representative_labelling_hap2"])]
            pred = [blocks[b5][h[b5]] + blocks[b6][h[b6]] for h in haps]
            fasta += [(f"{cell}_{arm}_h1", pred[0]), (f"{cell}_{arm}_h2", pred[1])]
            for h in (1, 2):
                for t in (1, 2):
                    pairs.append((f"{cell}|{arm}|h{h}|t{t}", f"{cell}_{arm}_h{h}", f"{cell}_truth{t}"))
            meta.append(dict(cell=cell, arm=arm, status="ok", model=model, model_fingerprint_matches_run=fp == st.get("model_fingerprint"),
                             sv5=f"{haps[0][b5]}/{haps[1][b5]}", ctx6=f"{haps[0][b6]}/{haps[1][b6]}"))
    with open(f"{out}/intervals.fa", "w") as f:
        for name, s in fasta:
            f.write(f">{name}\n{s}\n")
    with open(f"{out}/pairs.tsv", "w") as f:
        for lab, a, b in pairs:
            f.write(f"{lab}\t{out}/intervals.fa\t{a}\t{out}/intervals.fa\t{b}\n")
    code = subprocess.run([f"{WT}/build/genotype_probe", "edit-distances", "--pairs", f"{out}/pairs.tsv", "--out", f"{out}/distances.tsv", "--threads", "4"]).returncode
    dist = {}
    for r in T(f"{out}/distances.tsv"):
        dist[r["label"]] = int(r["distance"])
    for m in meta:
        if m["status"] != "ok":
            continue
        k = lambda h, t: dist[f"{m['cell']}|{m['arm']}|h{h}|t{t}"]
        m["interval_distance"] = min(k(1, 1) + k(2, 2), k(1, 2) + k(2, 1))
    write("intervals", meta)
    print(open(f"{out}/intervals.tsv").read())
    sys.exit(0 if code == 0 and all(m.get("model_fingerprint_matches_run", True) for m in meta) else 1)

if mode == "analyse":
    iv = {(r["cell"], r["arm"]): r for r in T(f"{OUT}/intervals/intervals.tsv")}
    res, changed_rows = [], []
    for cell in CELLS:
        info = cell_info(cell)
        arms = {}
        for arm in ("off", "on"):
            d = arm_dir(arm, cell)
            sl, st = KV(f"{d}.sequence_locus.tsv"), KV(f"{d}.mosaic_status.tsv")
            lg = (T(f"{d}.locus_genotype.tsv") or [{}])[0]
            ev = KV(f"{OUT}/eval.{arm}.{cell}/eval.summary.tsv")
            eb = {r["block_id"]: r for r in T(f"{OUT}/eval.{arm}.{cell}/eval.blocks.tsv")}
            sb = {r["block_id"]: r for r in T(f"{d}.sequence_blocks.tsv")}
            kinds = {r["block_id"]: r["kind"] for r in T(f"{d}.block_genotypes.tsv")} or {b: e.get("kind") for b, e in eb.items()}
            cnt = dict(exact=0, not_exact=0, unattainable=0, not_identified=0)
            for b, row in sb.items():
                if kinds.get(b) != "SV":
                    continue
                if row["block_status"] != "IDENTIFIED":
                    cnt["not_identified"] += 1
                    continue
                e = eb.get(b)
                if not e or e["truth1_closest_allele"] == "NA" or not (e["truth1_closest_distance"] == "0" and e["truth2_closest_distance"] == "0"):
                    cnt["unattainable"] += 1
                elif sorted(map(int, row["genotype"].split("/"))) == sorted((int(e["truth1_closest_allele"]), int(e["truth2_closest_allele"]))):
                    cnt["exact"] += 1
                else:
                    cnt["not_exact"] += 1
            moves = [f"{r['repaired_block']}:{r['allele_x']}/{r['allele_y']}" for r in T(f"{d}.repair_rounds.tsv") if r.get("moved") == "true"]
            arms[arm] = dict(sl=sl, st=st, lg=lg, ev=ev, sb=sb, cnt=cnt, moves=moves, iv=iv.get((cell, arm), {}))
        o, n = arms["off"], arms["on"]
        complete = all(a["sl"].get("prediction") == "PREDICTION" for a in (o, n))
        same = same_sequence_pair(o["sl"], n["sl"])
        changed = [b for b in o["sb"] if o["sb"][b].get("genotype") != n["sb"].get(b, {}).get("genotype")]
        for b in changed:
            changed_rows.append(dict(cell=cell, block=b, off=o["sb"][b].get("genotype"), on=n["sb"].get(b, {}).get("genotype"),
                                     off_error_mass_upper=o["sb"][b].get("error_mass_upper"), on_error_mass_upper=n["sb"].get(b, {}).get("error_mass_upper")))
        num = lambda a, k: int(a["ev"][k]) if a["ev"].get(k, "NA") not in ("NA", "") else None
        metrics = [("interval_distance", int(o["iv"]["interval_distance"]) if o["iv"].get("interval_distance") not in (None, "") else None,
                    int(n["iv"]["interval_distance"]) if n["iv"].get("interval_distance") not in (None, "") else None, "lower"),
                   ("exact_reportable_blocks", o["cnt"]["exact"], n["cnt"]["exact"], "higher"),
                   ("locus_excess", num(o, "excess_over_floor_total"), num(n, "excess_over_floor_total"), "lower")]
        better = [m for m, a, b, d in metrics if a is not None and b is not None and ((b < a) if d == "lower" else (b > a))]
        worse = [m for m, a, b, d in metrics if a is not None and b is not None and ((b > a) if d == "lower" else (b < a))]
        if not complete:
            fam = "INCOMPLETE_PAIR"
        elif same is True:
            fam = "UNCHANGED"
        elif better and worse:
            fam = "MIXED_WITHIN_FAMILY"
        elif better:
            fam = "IMPROVED"
        elif worse:
            fam = "DETERIORATED"
        else:
            fam = "CHANGED_NEUTRAL"
        g = lambda a, k: a["lg"].get(k, "NA")
        res.append(dict(cell=cell, scheme=info["scheme"], family=info["family"], donor=info["donor"], complete_pair=complete, same_sequence_pair=same,
                        family_outcome=fam, better=",".join(better) or "none", worse=",".join(worse) or "none", blocks_changed=",".join(changed) or "none",
                        off_prediction=f"{o['sl'].get('prediction')}/{o['sl'].get('sequence_status')}", on_prediction=f"{n['sl'].get('prediction')}/{n['sl'].get('sequence_status')}",
                        off_interval=o["iv"].get("interval_distance"), on_interval=n["iv"].get("interval_distance"),
                        off_sv5=o["iv"].get("sv5"), on_sv5=n["iv"].get("sv5"), off_ctx6=o["iv"].get("ctx6"), on_ctx6=n["iv"].get("ctx6"),
                        off_exact=o["cnt"]["exact"], on_exact=n["cnt"]["exact"], off_not_exact=o["cnt"]["not_exact"], on_not_exact=n["cnt"]["not_exact"],
                        unattainable=o["cnt"]["unattainable"], off_not_identified=o["cnt"]["not_identified"], on_not_identified=n["cnt"]["not_identified"],
                        floor=o["ev"].get("availability_floor_total"), off_distance=o["ev"].get("unordered_distance_total"), on_distance=n["ev"].get("unordered_distance_total"),
                        off_excess=o["ev"].get("excess_over_floor_total"), on_excess=n["ev"].get("excess_over_floor_total"),
                        off_total=g(o, "total_score"), on_total=g(n, "total_score"), off_fragment=g(o, "fragment_term"), on_fragment=g(n, "fragment_term"),
                        off_exposure=g(o, "exposure_term"), on_exposure=g(n, "exposure_term"), off_prior=g(o, "sequence_pair_prior"), on_prior=g(n, "sequence_pair_prior"),
                        off_moves=";".join(o["moves"]), on_moves=";".join(n["moves"]),
                        off_tie_branches=o["st"].get("tie_branches_run", "0"), on_tie_branches=n["st"].get("tie_branches_run", "0"),
                        off_phase=o["sl"].get("legacy_phase_support"), on_phase=n["sl"].get("legacy_phase_support"),
                        on_eligible=n["st"].get("junction_eligible_blocks"), on_invocations=n["st"].get("junction_invocations"),
                        on_cache_hits=n["st"].get("junction_cache_hits"), on_cutoff_ties=n["st"].get("junction_cutoff_tie_crossings"),
                        on_junction_seconds=n["st"].get("junction_seconds"),
                        off_seconds=LED.get(f"call.off.{cell}", {}).get("wall_seconds"), on_seconds=LED.get(f"call.on.{cell}", {}).get("wall_seconds"),
                        off_peak_gib=LED.get(f"call.off.{cell}", {}).get("peak_rss_gib"), on_peak_gib=LED.get(f"call.on.{cell}", {}).get("peak_rss_gib")))
    # predeclared batch reading
    loo = [r for r in res if r["scheme"] == "LOO"]
    lzo = [r for r in res if r["scheme"] == "LZO"]
    incomplete = [r["cell"] for r in res if not r["complete_pair"]]
    control_harm = [r["cell"] for r in lzo if r["family_outcome"] in ("DETERIORATED", "MIXED_WITHIN_FAMILY") or
                    (r["off_excess"] == "0" and r["family_outcome"] not in ("UNCHANGED", "INCOMPLETE_PAIR") and r["on_excess"] != "0")]
    changed_any = [r["cell"] for r in res if r["family_outcome"] not in ("UNCHANGED", "INCOMPLETE_PAIR")]
    improved = [r["cell"] for r in loo if r["family_outcome"] == "IMPROVED"]
    deteriorated = [r["cell"] for r in res if r["family_outcome"] in ("DETERIORATED", "MIXED_WITHIN_FAMILY")]
    if control_harm:
        verdict = "CONTROL_HARM"
    elif not changed_any:
        verdict = "NON_DISCRIMINATING"
    elif improved and deteriorated:
        verdict = "MIXED"
    elif improved and not deteriorated and all(r["family_outcome"] == "UNCHANGED" for r in lzo):
        verdict = "SUPPORTIVE"
    elif deteriorated:
        verdict = "MIXED" if improved else "DETERIORATION_ONLY"
    else:
        verdict = "CHANGED_WITHOUT_IMPROVEMENT"
    verdict_row = dict(verdict=verdict, improved=",".join(improved) or "none", deteriorated=",".join(deteriorated) or "none",
                       changed=",".join(changed_any) or "none", incomplete_pairs=",".join(incomplete) or "none",
                       limitation="incomplete pairs limit the conclusion and are never counted as unchanged" if incomplete else "none")
    write("paired", res)
    write("changed_blocks", changed_rows)
    write("verdict", [verdict_row])
    for r in res:
        print({k: r[k] for k in ("cell", "family", "donor", "family_outcome", "blocks_changed", "off_interval", "on_interval", "off_exact", "on_exact",
                                 "off_excess", "on_excess", "on_invocations")})
    print(verdict_row)
    sys.exit(0)
