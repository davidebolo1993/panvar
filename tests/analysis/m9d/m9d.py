#!/usr/bin/env python3
"""M9-D isolated attribution of the M9-C ACOT high-confidence errors (frozen with DESIGN.md). Read-only on M9-C.

  m9d.py models <out>     model dumps exist (cyp2d6_block_alleles, run by the runner); fingerprints equal each run's
                          model_fingerprint; degeneracy inventory result recorded (degeneracy_inventory.py, run by runner)
  m9d.py construct <out>  per case: context A (reported winner) and B (best evaluated carrier of the exact truth pair at the
                          block); within each, the called and the truth block pair in both homologue assignments, every other
                          block unchanged; LZO controls: the correct call vs the LOO-called pair mapped by exact sequence;
                          spelling gates (block bytes, unchanged flanks, offsets), sequence-level dedupe, pool/shortlist
                          membership; writes <case>.candidates.tsv and <case>.decompose for score-set
  m9d.py analyse <out>    reproduction gates, matched contrasts (truth - called) with fragment / exposure / prior / total and
                          numerical bounds, per-fragment localisation from the --decompose tables, conservative reading
"""
import csv
import hashlib
import math
import os
import sys

mode, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
R = "/Users/davide.bolognini/panvar/tests/results"
C = f"{R}/m9c_acot_reliability_20261020/runs"
OUT = f"{R}/m9d_attribution_20261020/runs"
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l and not l.startswith("#")} if os.path.exists(p) else {}
H = lambda s: hashlib.sha256(s.encode()).hexdigest()
# case, cell, block, kind, model dump name
CASES = [("E1", "LOO.d1", "sv:5", "error", "LOO"), ("E2", "LOO.d2", "sv:5", "error", "LOO"), ("E3", "LOO.d3", "sv:5", "error", "LOO"),
         ("E4", "LOO.d3", "sv:9", "error", "LOO"), ("E5", "LOO.d4", "sv:5", "error", "LOO"), ("E6", "LOO.d6", "sv:5", "error", "LOO"),
         ("C1", "LZO.d1", "sv:5", "control", "LZO.d1"), ("C2", "LZO.d2", "sv:5", "control", "LZO.d2")]
CONTROL_SOURCE = {"C1": "LOO.d1", "C2": "LOO.d2"}
EPS = 2.0 ** -52


def write(name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def dump(name):
    blocks = []
    for r in T(f"{OUT}/models/{name}.block_alleles.tsv"):
        b, a = int(r["block_index"]), int(r["allele_index"])
        while len(blocks) <= b:
            blocks.append([])
        assert a == len(blocks[b])
        blocks[b].append(r["sequence"])
    return blocks


def spell(blocks, labels):
    seq, spans, pos = [], [], 0
    for b, a in enumerate(labels):
        s = blocks[b][a]
        spans.append((pos, pos + len(s)))
        seq.append(s)
        pos += len(s)
    return "".join(seq), spans


labs = lambda s: [int(x) for x in s.split(",")]
key = lambda h1, h2: tuple(sorted((tuple(h1), tuple(h2))))

if mode == "models":
    import subprocess
    WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
    gfa = f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa"
    ref = "GRCh38#0#chr14:73416151-73668669"
    P = f"{C}/prepare"
    prep = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in open(f"{P}/prepare.tsv")}
    d1, d2 = prep["donor_1"][0], prep["donor_2"][0]
    for name, cur in (("LOO", f"{P}/curation.LOO.txt"), ("LZO.d1", f"{P}/curation.LZO.{d1}.txt"), ("LZO.d2", f"{P}/curation.LZO.{d2}.txt")):
        assert os.path.exists(cur), cur
        assert subprocess.run([f"{WT}/build/cyp2d6_block_alleles", gfa, ref, cur, f"{out}/{name}"]).returncode == 0
    rows = []
    for name, cell in (("LOO", "LOO.d1"), ("LZO.d1", "LZO.d1"), ("LZO.d2", "LZO.d2")):
        s = KV(f"{OUT}/models/{name}.summary.tsv")
        run_fps = {c: KV(f"{C}/call.{c}/call.mosaic_status.tsv")["model_fingerprint"] for c in
                   (["LOO.d1", "LOO.d2", "LOO.d3", "LOO.d4", "LOO.d6"] if name == "LOO" else [cell])}
        rows.append(dict(model=name, dump_fingerprint=s.get("fingerprint", s.get("model_fingerprint", "absent")),
                         run_fingerprints=";".join(f"{c}={v}" for c, v in run_fps.items()),
                         match=all(v == s.get("fingerprint", s.get("model_fingerprint")) for v in run_fps.values())))
    write("models", rows)
    print(open(f"{out}/models.tsv").read())
    sys.exit(0 if all(r["match"] for r in rows) else 1)

if mode == "degeneracy":   # model-only, truth-free; informs sequence-vs-label identity, never gates construction
    import subprocess
    WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
    M = f"{OUT}/models"
    code = subprocess.run([sys.executable, f"{WT}/tests/analysis/cyp2d6_pilot/degeneracy_inventory.py", f"{out}/degeneracy.tsv", f"{out}/degeneracy_summary.tsv",
                           f"LOO={M}/LOO", f"LZO.d1={M}/LZO.d1", f"LZO.d2={M}/LZO.d2"]).returncode
    print(open(f"{out}/degeneracy_summary.tsv").read() if os.path.exists(f"{out}/degeneracy_summary.tsv") else "no summary")
    sys.exit(code)

if mode == "construct":
    gates, members, manifest = [], [], []
    for case, cell, block, kind, model in CASES:
        blocks = dump(model)
        d = f"{C}/call.{cell}"
        ids = [r["block_id"] for r in T(f"{d}/call.sequence_blocks.tsv")]
        b = ids.index(block)
        lg = T(f"{d}/call.locus_genotype.tsv")[0]
        A = (labs(lg["hap1_block_path"]), labs(lg["hap2_block_path"]), float(lg["total_score"]))
        called = tuple(sorted((A[0][b], A[1][b])))
        cand = T(f"{d}/call.mosaic_candidates.tsv")
        sc = {r["sequence_candidate_id"]: r for r in T(f"{d}/call.mosaic_scores.tsv")}
        pool = {key(labs(r["hap1_alleles"]), labs(r["hap2_alleles"])): r for r in cand}
        contexts = {"A": A}
        if kind == "error":
            e = next(r for r in T(f"{C}/eval.{cell}/eval.blocks.tsv") if r["block_id"] == block)
            assert e["truth1_closest_distance"] == "0" and e["truth2_closest_distance"] == "0", "floor-zero case expected"
            other = tuple(sorted((int(e["truth1_closest_allele"]), int(e["truth2_closest_allele"]))))
            carriers = [r for r in cand if tuple(sorted((labs(r["hap1_alleles"])[b], labs(r["hap2_alleles"])[b]))) == other]
            best = max(carriers, key=lambda r: float(sc[r["sequence_candidate_id"]]["score"]))
            contexts["B"] = (labs(best["hap1_alleles"]), labs(best["hap2_alleles"]), float(sc[best["sequence_candidate_id"]]["score"]))
            names = {"called": called, "truth": other}
        else:
            # LZO control: the correct call vs the LOO-called pair of the same donor, mapped by exact sequence
            src = CONTROL_SOURCE[case]
            loo = dump("LOO")
            slg = T(f"{C}/call.{src}/call.locus_genotype.tsv")[0]
            lc = sorted((labs(slg["hap1_block_path"])[b], labs(slg["hap2_block_path"])[b]))
            mapped = []
            for a in lc:
                hits = [i for i, s in enumerate(blocks[b]) if s == loo[b][a]]
                gates.append(dict(case=case, gate="LOO allele mapped by exact sequence", detail=f"LOO {block}:{a} -> LZO {hits}", ok=len(hits) == 1))
                mapped.append(hits[0] if hits else None)
            other = tuple(sorted(mapped))
            names = {"called": called, "inserted_LOO_call": other}
        # constructed candidates: within each context, each named pair in both homologue assignments
        rows = []
        for cx, (h1, h2, expected) in contexts.items():
            orig = key(h1, h2)
            for pname, pair in names.items():
                for asg, (x, y) in (("ab", pair), ("ba", pair[::-1])):
                    n1, n2 = list(h1), list(h2)
                    n1[b], n2[b] = x, y
                    lab = f"{case}_{cx}_{pname}_{asg}"
                    rows.append(dict(label=lab, context=cx, pair=pname, assignment=asg, h1=n1, h2=n2, is_context=key(n1, n2) == orig,
                                     expected=expected if key(n1, n2) == orig else None))
                    # spelling gates: replacement bytes, unchanged flanks (offsets may shift)
                    for side, (old, new) in enumerate(((h1, n1), (h2, n2))):
                        so, spo = spell(blocks, old)
                        sn, spn = spell(blocks, new)
                        ok = (sn[:spn[b][0]] == so[:spo[b][0]] and sn[spn[b][1]:] == so[spo[b][1]:] and sn[spn[b][0]:spn[b][1]] == blocks[b][new[b]])
                        gates.append(dict(case=case, gate=f"{lab} hap{side + 1} flanks byte-identical, block = allele {new[b]}",
                                          detail=f"offset shift {len(sn) - len(so)} bp", ok=ok))
        # sequence-level dedupe and pool membership
        seen = {}
        for r in rows:
            s1, s2 = spell(blocks, r["h1"])[0], spell(blocks, r["h2"])[0]
            sk = tuple(sorted((H(s1), H(s2))))
            r["sequence_pair"] = f"{sk[0][:12]}|{sk[1][:12]}"
            r["duplicate_of"] = seen.get(sk, "")
            seen.setdefault(sk, r["label"])
            p = pool.get(key(r["h1"], r["h2"]))
            r["membership"] = ("proposal shortlist (mosaic beam)" if p["source"] == "mosaic" else f"scored pool ({p['source']}, {p['scored_in']})") if p else "absent from scored pool"
            r["pool_score"] = sc[p["sequence_candidate_id"]]["score"] if p else ""
            members.append({k: (",".join(map(str, v)) if isinstance(v, list) else v) for k, v in r.items() if k not in ("h1", "h2")} | dict(case=case, cell=cell, block=block))
        with open(f"{out}/{case}.candidates.tsv", "w") as f:
            for r in rows:
                if r["duplicate_of"]:
                    continue   # identical unordered sequence pair: scored once, the duplicate retained in membership.tsv
                f.write(f"{r['label']}\t{','.join(map(str, r['h1']))}\t{','.join(map(str, r['h2']))}\t{repr(r['expected']) if r['expected'] is not None else 'NA'}\n")
        # decompositions: truth (or inserted) minus called, within each context, every distinct assignment pairing
        uniq = [r for r in rows if not r["duplicate_of"]]
        pairs = []
        for cx in contexts:
            first = [r for r in uniq if r["context"] == cx and r["pair"] != "called"]
            second = [r for r in uniq if r["context"] == cx and r["pair"] == "called"]
            pairs += [f"{x['label']},{y['label']}" for x in first for y in second]
        open(f"{out}/{case}.decompose", "w").write(";".join(pairs))
        manifest.append(dict(case=case, cell=cell, block=block, kind=kind, model=model, block_index=b, called=f"{called[0]}/{called[1]}",
                             other=f"{other[0]}/{other[1]}" if None not in other else "unmapped", contexts=",".join(contexts),
                             candidates_scored=len(uniq), decompositions=len(pairs)))
    write("gates", gates)
    write("membership", members)
    write("cases", manifest)
    bad = [g for g in gates if not g["ok"]]
    print(open(f"{out}/cases.tsv").read())
    print(f"gates: {len(gates) - len(bad)}/{len(gates)} pass")
    sys.exit(0 if not bad else 1)

if mode == "intervals":   # Amendment 1: interval-level sequence comparison over the degenerate sv:5..ctx:6 (and sv:9 alone)
    blocks = dump("LOO")
    rows = []
    for case, cell, block, kind, model in CASES:
        if kind != "error":
            continue
        d = f"{C}/call.{cell}"
        ids = [r["block_id"] for r in T(f"{d}/call.sequence_blocks.tsv")]
        span = ["sv:5", "ctx:6"] if block == "sv:5" else [block]
        lg = T(f"{d}/call.locus_genotype.tsv")[0]
        h = (labs(lg["hap1_block_path"]), labs(lg["hap2_block_path"]))
        called = ["".join(blocks[ids.index(x)][hh[ids.index(x)]] for x in span) for hh in h]
        bs = f"{OUT}/interval.{cell}/bs"
        truth = []
        for k in (0, 1):
            seq = ""
            for x in span:
                recs, cur = {}, None
                for line in open(f"{bs}.{x}.truth.fa"):
                    if line.startswith(">"):
                        cur = line[1:].split()[0]
                        recs[cur] = ""
                    else:
                        recs[cur] += line.strip()
                seq += recs[f"truth{k + 1}"]   # Amendment 2: block-sequences names the records truth1, truth2
            truth.append(seq)
        def core(a, b):   # residual differing core after trimming the common prefix and suffix
            p = 0
            while p < min(len(a), len(b)) and a[p] == b[p]:
                p += 1
            q = 0
            while q < min(len(a), len(b)) - p and a[-1 - q] == b[-1 - q]:
                q += 1
            return len(a) - p - q, len(b) - p - q
        best = min(((0, 1), (1, 0)), key=lambda o: sum(sum(core(called[i], truth[o[i]])) for i in (0, 1)))
        for i in (0, 1):
            c, t = called[i], truth[best[i]]
            rows.append(dict(case=case, cell=cell, interval="..".join(span), haplotype=f"called hap{i + 1} vs truth{best[i] + 1}", identical=c == t,
                             called_bp=len(c), truth_bp=len(t), residual_core_called_bp=core(c, t)[0], residual_core_truth_bp=core(c, t)[1]))
    write("intervals", rows)
    print(open(f"{out}/intervals.tsv").read())
    sys.exit(0)

if mode == "analyse":
    rep, contrasts, frags = [], [], []
    mem = {r["label"]: r for r in T(f"{OUT}/construct/membership.tsv")}
    for case, cell, block, kind, model in CASES:
        d = f"{OUT}/score.{case}"
        sc = {r["label"]: r for r in T(f"{d}/score.scores_truth_free.tsv")}
        if not sc:
            contrasts.append(dict(case=case, status="no score-set output"))
            continue
        saved = {r["sequence_candidate_id"]: r for r in T(f"{C}/call.{cell}/call.mosaic_scores.tsv")}
        cand = {}
        for r in T(f"{C}/call.{cell}/call.mosaic_candidates.tsv"):
            cand[key(labs(r["hap1_alleles"]), labs(r["hap2_alleles"]))] = r["sequence_candidate_id"]
        # reproduction of the saved contexts (fragment, exposure, prior, total)
        for lab, m in mem.items():
            if m["case"] != case or m["is_context"] != "True" or m["duplicate_of"]:
                continue
            s = sc[lab]
            sid = cand.get(key(labs(s["hap1_alleles"]), labs(s["hap2_alleles"])))
            v = saved[sid]
            rep.append(dict(case=case, context=m["context"], label=lab, d_fragment=float(s["fragment_term"]) - float(v["fragment_term"]),
                            d_exposure=float(s["exposure_term"]) - float(v["exposure_term"]), d_prior=float(s["prior"]) - float(v["sequence_pair_prior"]),
                            d_total=float(s["total"]) - float(v["score"]), expected_total=s["expected_total"],
                            ok=all(abs(x) <= 1e-6 for x in (float(s["fragment_term"]) - float(v["fragment_term"]), float(s["exposure_term"]) - float(v["exposure_term"]),
                                                            float(s["prior"]) - float(v["sequence_pair_prior"]), float(s["total"]) - float(v["score"])))))
        # matched contrasts from the decompositions (truth/inserted minus called)
        for summ in sorted(f for f in os.listdir(d) if f.endswith(".summary.tsv") and ".decomp." in f):
            k = KV(f"{d}/{summ}")
            a, b2 = k["label_A"], k["label_B"]
            sa, sb = sc[a], sc[b2]
            fr = [(float(r["delta"]), r["class"], abs(float(r["contribution_A"])), abs(float(r["contribution_B"])))
                  for r in T(f"{d}/{summ.replace('.summary.tsv', '.fragments.tsv')}")]
            nfr = len(fr)
            abs_a = abs(float(sa["exposure_term"])) + sum(x[2] for x in fr)
            abs_b = abs(float(sb["exposure_term"])) + sum(x[3] for x in fr)
            prior_a, prior_b = float(sa["prior"]), float(sb["prior"])
            bound = lambda s_abs, pr: max((nfr + 16) * EPS * (s_abs + (abs(pr) if math.isfinite(pr) else 0.0)), 0.5e-9)
            bnd = bound(abs_a, prior_a) + bound(abs_b, prior_b)
            dt = (float(sa["total"]) - float(sb["total"])) if sa["total"] != "NA" and sb["total"] != "NA" else float("nan")
            df, de, dp = float(sa["fragment_term"]) - float(sb["fragment_term"]), float(sa["exposure_term"]) - float(sb["exposure_term"]), prior_a - prior_b
            deltas = sorted((x[0] for x in fr), key=lambda v: -abs(v))
            sabs = sum(abs(v) for v in deltas) or 1.0
            cls = {}
            for v, c, _, _ in fr:
                cls.setdefault(c, [0, 0.0])
                cls[c][0] += abs(v) > 1e-9
                cls[c][1] += v
            ctx = mem[a]["context"]
            absent = mem[a]["membership"].startswith("absent")
            if not math.isfinite(dt):
                reading = "no total (prior support undefined)"
            elif dt > bnd and absent:
                reading = "replacement wins total but was not evaluated: search/proposal gap in this context"
            elif dt > bnd:
                reading = "replacement already scored higher than the context: INCONSISTENCY, investigate first"
            elif df + de > 0 and dt < -bnd:
                reading = "replacement wins likelihood but loses total: prior effect (quantified, not retuned)"
            elif df + de < 0:
                reading = "replacement loses likelihood: inspect the contributing evidence (not by itself a bug)"
            else:
                reading = "within numerical bounds: tied"
            contrasts.append(dict(case=case, cell=cell, block=block, context=ctx, truth_or_inserted=a, called=b2,
                                  membership_truth=mem[a]["membership"], membership_called=mem[b2]["membership"],
                                  reading=reading, d_fragment=df, d_exposure=de, d_likelihood=df + de, d_prior=dp, d_total=dt, bound=bnd,
                                  certified=("yes" if abs(dt) > bnd else "no") if math.isfinite(dt) else "NA",
                                  fragments_nonzero=sum(1 for v in deltas if abs(v) > 1e-9), top5_share=sum(abs(v) for v in deltas[:5]) / sabs,
                                  top20_share=sum(abs(v) for v in deltas[:20]) / sabs, fragments_abs_gt_1nat=sum(1 for v in deltas if abs(v) > 1.0),
                                  by_class=";".join(f"{c}:{n} frags,{s:+.3f}" for c, (n, s) in sorted(cls.items())),
                                  placement_interpretable=k.get("placement_evidence_interpretable", "NA"),
                                  decomp_fragment_minus_scorer=float(k["fragment_term_A"]) - float(k["fragment_term_A_scorer"])))
            for r in T(f"{d}/{summ.replace('.summary.tsv', '.top_fragments.tsv')}")[:10]:
                frags.append(dict(case=case, contrast=f"{a} - {b2}", rank=r["rank"], fragment=r["fragment"], delta=r["delta"], cls=r["class"]))
    write("reproduction", rep)
    write("contrasts", contrasts)
    write("top_fragments", frags)
    print(open(f"{out}/reproduction.tsv").read())
    print(open(f"{out}/contrasts.tsv").read())
    sys.exit(0 if rep and all(r["ok"] for r in rep) else 1)
