#!/usr/bin/env python3
"""M9-F driver (frozen with DESIGN.md).

  m9f.py gates <out>            synthetic fixtures for rule JP1 against the independent oracle
  m9f.py propose <cell> <out>   JP1 shortlist at sv:5 in the authoritative v2 representative context; no truth read
  m9f.py freeze <out>           SHA-256 of every frozen shortlist, before any truth or full scoring
  m9f.py analyse <out>          after scoring: incumbent reproduction, improvements, M9-D pair ranks, truth agreement,
                                diversity / ties / overlap with the index-ordered list, LZO control behaviour, cost
"""
import csv
import hashlib
import math
import os
import random
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import jprop  # noqa: E402
from jprop import K  # noqa: E402

R = "/Users/davide.bolognini/panvar/tests/results"
C = f"{R}/m9c_acot_reliability_20261020/runs"
D = f"{R}/m9d_attribution_20261020/runs"
E9 = f"{R}/m9e_reachability_20261020/runs"
OUT = f"{R}/m9f_junction_proposals_20261020/runs"
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l and not l.startswith("#")} if os.path.exists(p) else {}
labs = lambda s: [int(x) for x in s.split(",")]
H = lambda s: hashlib.sha256(s.encode()).hexdigest()
CELLS = ["LOO.d1", "LOO.d3", "LOO.d6", "LZO.d1", "LZO.d2", "LOO.d2", "LOO.d4", "LOO.d5"]
MODEL = lambda cell: "LOO" if cell.startswith("LOO") else cell
SHORT = 32


def write(path, rows):
    with open(path, "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def dump(name):
    blocks = []
    for r in T(f"{D}/models/{name}.block_alleles.tsv"):
        b = int(r["block_index"])
        while len(blocks) <= b:
            blocks.append([])
        blocks[b].append(r["sequence"])
    return blocks


def fastq(path):
    with open(path) as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                yield line.strip()


def setup(blocks, labels, b):
    """(prefix, suffix) per homologue from the incumbent labels."""
    out = []
    for h in labels:
        out.append(("".join(blocks[i][h[i]] for i in range(b)), "".join(blocks[i][h[i]] for i in range(b + 1, len(blocks)))))
    return out


def propose(blocks, labels, b, reads, rate, bg):
    ctx = setup(blocks, labels, b)
    alleles = blocks[b]
    W1 = jprop.windows(ctx[0][0], ctx[0][1], alleles)
    W2 = jprop.windows(ctx[1][0], ctx[1][1], alleles)
    inventory = set().union(*W1, *W2)
    e = jprop.elsewhere([(ctx[0][0], alleles[labels[0][b]], ctx[0][1]), (ctx[1][0], alleles[labels[1][b]], ctx[1][1])], inventory)
    y = jprop.read_counts(reads, inventory)
    score = jprop.make_scorer(e, y, rate, bg)
    return jprop.rank_pairs(W1, W2, score), dict(inventory=inventory, e=e, y=y, ctx=ctx, W1=W1, W2=W2)


def oracle_check(blocks, labels, b, reads, rate, bg, rows, info, picks):
    """Score differences between picked candidates: fast delta vs the oracle on whole spelled haplotypes."""
    ctx, inv = info["ctx"], info["inventory"]
    vals = []
    for r in picks:
        x, y = r["x"], r["y"]
        h1 = ctx[0][0] + blocks[b][x] + ctx[0][1]
        h2 = ctx[1][0] + blocks[b][y] + ctx[1][1]
        vals.append((r["score_ab"], jprop.oracle_score([h1, h2], inv, reads, rate, bg)))
    return max(abs((vals[i][0] - vals[0][0]) - (vals[i][1] - vals[0][1])) for i in range(1, len(vals))) if len(vals) > 1 else 0.0


mode = sys.argv[1]

if mode == "gates":
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    rng = random.Random(20261020)
    rs = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    rc = lambda s: s.translate(jprop.COMP)[::-1]
    res = []
    ok = lambda name, cond, detail="": res.append(dict(gate=name, result="PASS" if cond else "FAIL", detail=detail))
    rate, bg = 1.6, 0.012

    def reads_from(haps, depth=6, L=60):
        rr = []
        for h in haps:
            for i in range(0, len(h) - L + 1, max(1, L // depth)):
                rr.append(h[i:i + L])
        return rr

    def model(alleles):
        return [[rs(120)], alleles, [rs(120)]]

    # F1 short alleles distinguishable only through junctions (5 bp: no internal 31-mer)
    al = [rs(5) for _ in range(6)]
    m = model(al)
    truth = (3, 4)
    reads = reads_from([m[0][0] + al[truth[0]] + m[2][0], m[0][0] + al[truth[1]] + m[2][0]])
    rows, info = propose(m, [[0, 0, 0], [0, 1, 0]], 1, reads, rate, bg)
    ok("F1 short allele via junctions", (rows[0]["x"], rows[0]["y"]) == truth and rows[0]["tie_class_size"] == 1,
       f"top {rows[0]['x']}/{rows[0]['y']} tie {rows[0]['tie_class_size']}")
    picks = [rows[0], rows[5], rows[-1]]
    d = oracle_check(m, [[0, 0, 0], [0, 1, 0]], 1, reads, rate, bg, rows, info, picks)
    ok("F1 oracle (delta vs whole-haplotype scan)", d < 1e-9, f"max |diff| {d:.2e}")
    # F2 empty allele
    al = ["", rs(4), rs(9), rs(12)]
    m = model(al)
    reads = reads_from([m[0][0] + "" + m[2][0], m[0][0] + al[2] + m[2][0]])
    rows, info = propose(m, [[0, 1, 0], [0, 3, 0]], 1, reads, rate, bg)
    ok("F2 empty allele", (rows[0]["x"], rows[0]["y"]) == (0, 2), f"top {rows[0]['x']}/{rows[0]['y']}")
    d = oracle_check(m, [[0, 1, 0], [0, 3, 0]], 1, reads, rate, bg, rows, info, [rows[0], rows[1], rows[-1]])
    ok("F2 oracle", d < 1e-9, f"{d:.2e}")
    # F3 a K-mer spanning both boundaries of a very short allele is counted once per occurrence
    al = [rs(3), rs(3)]
    m = model(al)
    W = jprop.windows(m[0][0], m[2][0], al)
    span_both = sum(1 for i in range(len(m[0][0][-(K - 1):] + al[0] + m[2][0][:K - 1]) - K + 1)
                    if i < K - 1 and i + K > K - 1 + 3)
    window = m[0][0][-(K - 1):] + al[0] + m[2][0][:K - 1]
    both = [jprop.canon(window[i:i + K]) for i in range(len(window) - K + 1) if i < K - 1 and i + K > K - 1 + 3]
    ok("F3 both-boundary K-mers counted", sum(W[0].values()) == 3 + K - 1 and span_both == (K - 1) - 3 and all(W[0][c] >= 1 for c in both),
       f"window K-mers {sum(W[0].values())} (expected {3 + K - 1}); spanning both boundaries {span_both}")
    # F4 repeated marker elsewhere in the locus: copy the junction context into the far flank
    left = rs(120)
    junction = left[-(K - 1):]
    right = rs(60) + junction + "A" + rs(60)
    al = ["A", "C"]
    m = [[left], al, [right]]
    reads = reads_from([left + "A" + right, left + "C" + right])
    rows, info = propose(m, [[0, 0, 0], [0, 1, 0]], 1, reads, rate, bg)
    km = jprop.canon(junction + "A")
    d = oracle_check(m, [[0, 0, 0], [0, 1, 0]], 1, reads, rate, bg, rows, info, rows)
    ok("F4 repeated marker elsewhere", info["e"].get(km, 0) >= 1 and d < 1e-9, f"elsewhere copies of the junction K-mer {info['e'].get(km, 0)}, oracle {d:.2e}")
    # F5 reverse-complement invariance (alleles, contexts and reads all reverse-complemented)
    al = [rs(rng.randint(0, 40)) for _ in range(5)]
    m = model(al)
    reads = reads_from([m[0][0] + al[1] + m[2][0], m[0][0] + al[4] + m[2][0]])
    rows, _ = propose(m, [[0, 0, 0], [0, 2, 0]], 1, reads, rate, bg)
    mr = [[rc(m[2][0])], [rc(a) for a in al], [rc(m[0][0])]]
    rows_rc, _ = propose(mr, [[0, 0, 0], [0, 2, 0]], 1, [rc(s) for s in reads], rate, bg)
    a = {(r["x"], r["y"]): r["best"] for r in rows}
    bb = {(r["x"], r["y"]): r["best"] for r in rows_rc}
    ok("F5 reverse-complement invariance", max(abs(a[k] - bb[k]) for k in a) < 1e-9, f"max |diff| {max(abs(a[k] - bb[k]) for k in a):.2e}")
    # F6 homologue swap of the incumbent
    rows_sw, _ = propose(m, [[0, 2, 0], [0, 0, 0]], 1, reads, rate, bg)
    c = {(r["x"], r["y"]): r["best"] for r in rows_sw}
    ok("F6 homologue-swap invariance", max(abs(a[k] - c[k]) for k in a) < 1e-9, f"{max(abs(a[k] - c[k]) for k in a):.2e}")
    # F7 allele-index permutation, sequences preserved
    perm = [3, 0, 4, 1, 2]
    mp = [m[0], [al[p] for p in perm], m[2]]
    inv = {p: i for i, p in enumerate(perm)}
    rows_p, _ = propose(mp, [[0, inv[0], 0], [0, inv[2], 0]], 1, reads, rate, bg)
    pm = {tuple(sorted((perm[r["x"]], perm[r["y"]]))): r["best"] for r in rows_p}
    ok("F7 allele-index permutation", max(abs(a[k] - pm[k]) for k in a) < 1e-9, f"{max(abs(a[k] - pm[k]) for k in a):.2e}")
    # F8 identical spelled sequence under different labels scores identically (and ties)
    al2 = al + [al[1]]
    m2 = model(al2)
    m2[0], m2[2] = m[0], m[2]
    rows8, _ = propose(m2, [[0, 0, 0], [0, 2, 0]], 1, reads, rate, bg)
    s8 = {(r["x"], r["y"]): r for r in rows8}
    ok("F8 same sequence, different labels", abs(s8[(1, 4)]["best"] - s8[(4, 5)]["best"]) < 1e-9 and s8[(1, 4)]["tie_class"] == s8[(4, 5)]["tie_class"],
       f"{s8[(1, 4)]['best']:.6f} vs {s8[(4, 5)]['best']:.6f}")
    # F9 genuinely uninformative: no read touches the block -> every pair tied
    # equal-length alleles: under JP1 a longer allele predicts more K-mers, so unequal lengths are informative even
    # without observations (expected-count term); the uninformative case is equal lengths and no read touching the block
    m9 = model([rs(20) for _ in range(5)])
    rows9, _ = propose(m9, [[0, 0, 0], [0, 2, 0]], 1, [rs(60) for _ in range(30)], rate, bg)
    ok("F9 uninformative stays tied", rows9[0]["tie_class_size"] == len(rows9), f"tie class {rows9[0]['tie_class_size']} of {len(rows9)}")
    write(f"{out}/gates.tsv", res)
    print(open(f"{out}/gates.tsv").read())
    sys.exit(0 if all(r["result"] == "PASS" for r in res) else 1)

if mode == "propose":
    cell, out = sys.argv[2], sys.argv[3]
    os.makedirs(out, exist_ok=True)
    blocks = dump(MODEL(cell))
    d = f"{C}/call.{cell}"
    sl = KV(f"{d}/call.sequence_locus.tsv")
    labels = [labs(sl["representative_labelling_hap1"]), labs(sl["representative_labelling_hap2"])]
    spelled = ["".join(blocks[i][h[i]] for i in range(len(blocks))) for h in labels]
    gate_context = [H(spelled[0]), H(spelled[1])] == [sl["representative_hap1_sha256"], sl["representative_hap2_sha256"]]
    ids = [r["block_id"] for r in T(f"{d}/call.sequence_blocks.tsv")]
    b = ids.index("sv:5")
    st = KV(f"{d}/call.mosaic_status.tsv")
    rate, bg = float(st["marker_depth_rate_per_copy"]), float(st["marker_background"])
    prep = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in open(f"{C}/prepare/prepare.tsv")}
    dn = prep[f"donor_{cell.split('.d')[1]}"][0]
    rd = f"{C}/prepare/reads/{dn}/cov5"
    reads = list(fastq(f"{rd}/reads_1.fq")) + list(fastq(f"{rd}/reads_2.fq"))
    import time
    t = time.time()
    rows, info = propose(blocks, labels, b, reads, rate, bg)
    secs = time.time() - t
    # empirical per-copy rate for all K-mers: K-mers present exactly once in each incumbent haplotype (invariant, 2 copies)
    c1, c2 = Counter(jprop.kmers(spelled[0])), Counter(jprop.kmers(spelled[1]))
    inv2 = {m for m in c1 if c1[m] == 1 and c2.get(m, 0) == 1}
    yy = jprop.read_counts(reads, inv2)
    emp = sum(yy.values()) / (2 * len(inv2))
    ratio = emp / rate
    # oracle on real data: three candidates at different ranks
    picks = [rows[0], rows[len(rows) // 2], rows[-1]]
    od = oracle_check(blocks, labels, b, reads, rate, bg, rows, info, picks)
    write(f"{out}/ranks.tsv", rows)
    short = rows[:SHORT]
    write(f"{out}/shortlist.tsv", [dict(rank=r["rank"], x=r["x"], y=r["y"], best=r["best"], score_ab=r["score_ab"], score_ba=r["score_ba"],
                                       tie_class=r["tie_class"], tie_class_size=r["tie_class_size"]) for r in short])
    cut = rows[SHORT - 1]
    summ = dict(cell=cell, model=MODEL(cell), block="sv:5", alleles=len(blocks[b]), pairs=len(rows), inventory=len(info["inventory"]),
                inventory_observed=sum(1 for m in info["inventory"] if info["y"].get(m, 0) > 0),
                inventory_elsewhere=sum(1 for m in info["inventory"] if info["e"].get(m, 0) > 0),
                context_sha_matches_v2=gate_context, sequence_status=sl["sequence_status"], rate=rate, background=bg,
                empirical_all_kmer_rate=emp, empirical_over_fitted=ratio, rate_gate=0.75 <= ratio <= 1.25,
                oracle_max_abs_diff=od, oracle_gate=od < 1e-6, distinct_scores=len({r["tie_class"] for r in rows}),
                top_tie_class_size=rows[0]["tie_class_size"], cutoff_tie_class_size=cut["tie_class_size"],
                cutoff_tie_spans_boundary=any(r["tie_class"] == cut["tie_class"] for r in rows[SHORT:]),
                shortlist_distinct_alleles=len({a for r in short for a in (r["x"], r["y"])}), proposal_seconds=secs)
    write(f"{out}/summary.tsv", [summ])
    print(open(f"{out}/summary.tsv").read())
    sys.exit(0 if gate_context and summ["rate_gate"] and summ["oracle_gate"] else 1)

if mode == "freeze":
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    rows = []
    for cell in CELLS:
        p = f"{OUT}/propose.{cell}/shortlist.tsv"
        rows.append(dict(cell=cell, shortlist_sha256=hashlib.sha256(open(p, "rb").read()).hexdigest() if os.path.exists(p) else "absent",
                         ranks_sha256=hashlib.sha256(open(p.replace("shortlist", "ranks"), "rb").read()).hexdigest()[:16] if os.path.exists(p) else "absent"))
        # score-set candidate file: the v2 representative (expected total) + both assignments of every shortlisted pair
        if not os.path.exists(p):
            continue
        sl = KV(f"{C}/call.{cell}/call.sequence_locus.tsv")
        h1, h2 = labs(sl["representative_labelling_hap1"]), labs(sl["representative_labelling_hap2"])
        b = [r["block_id"] for r in T(f"{C}/call.{cell}/call.sequence_blocks.tsv")].index("sv:5")
        seen = {tuple(sorted((tuple(h1), tuple(h2))))}
        lines = [f"incumbent\t{','.join(map(str, h1))}\t{','.join(map(str, h2))}\t{sl['representative_score']}"]
        for r in T(p):
            for asg, (x, y) in (("ab", (int(r["x"]), int(r["y"]))), ("ba", (int(r["y"]), int(r["x"])))):
                n1, n2 = list(h1), list(h2)
                n1[b], n2[b] = x, y
                k = tuple(sorted((tuple(n1), tuple(n2))))
                if k in seen:
                    continue
                seen.add(k)
                lines.append(f"P{r['rank']}_{asg}_{x}_{y}\t{','.join(map(str, n1))}\t{','.join(map(str, n2))}\tNA")
        open(f"{out}/{cell}.candidates.tsv", "w").write("\n".join(lines) + "\n")
    write(f"{out}/freeze.tsv", rows)
    print(open(f"{out}/freeze.tsv").read())
    sys.exit(0)

if mode == "analyse":
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    LED = {r["unit"]: r for r in T(f"{OUT}/../ledger.tsv")}
    improving = {"LOO.d1": [(6, 77)], "LOO.d3": [(34, 94)], "LOO.d6": [(13, 21)]}   # M9-D frozen improving pairs (LOO labels)
    res, cov = [], []
    for cell in CELLS:
        su = T(f"{OUT}/propose.{cell}/summary.tsv")
        if not su:
            res.append(dict(cell=cell, status="no proposal"))
            continue
        su = su[0]
        ranks = {(int(r["x"]), int(r["y"])): r for r in T(f"{OUT}/propose.{cell}/ranks.tsv")}
        short = T(f"{OUT}/propose.{cell}/shortlist.tsv")
        sl = KV(f"{C}/call.{cell}/call.sequence_locus.tsv")
        ids = [r["block_id"] for r in T(f"{C}/call.{cell}/call.sequence_blocks.tsv")]
        b = ids.index("sv:5")
        inc = tuple(sorted((labs(sl["representative_labelling_hap1"])[b], labs(sl["representative_labelling_hap2"])[b])))
        # original index-ordered list: the probe's ranking (M9-E for LOO; ranks.<LZO> here)
        orig_path = f"{E9}/ranks.{cell}/score.local_ranks.tsv" if cell in ("LOO.d1", "LOO.d3", "LOO.d6") else f"{OUT}/ranks.{cell}/score.local_ranks.tsv"
        orig = [r for r in T(orig_path) if r["block_id"] == "sv:5"]
        orig_top = {(int(r["allele_x"]), int(r["allele_y"])) for r in orig if r["in_local_list"] == "true"}
        orig_distinct = len({r["emission"] for r in orig}) if orig else "NA"
        new_top = {(int(r["x"]), int(r["y"])) for r in short}
        # truth (after freeze)
        ev = next((r for r in T(f"{C}/eval.{cell}/eval.blocks.tsv") if r["block_id"] == "sv:5"), None)
        tpair = tuple(sorted((int(ev["truth1_closest_allele"]), int(ev["truth2_closest_allele"])))) if ev else None
        attainable = ev is not None and ev["truth1_closest_distance"] == "0" and ev["truth2_closest_distance"] == "0"
        # full scoring
        sc = {r["label"]: r for r in T(f"{OUT}/score.{cell}/score.scores_truth_free.tsv")}
        row = dict(cell=cell, sequence_status=su["sequence_status"], incumbent_sv5=f"{inc[0]}/{inc[1]}", incumbent_correct=(inc == tpair) if attainable else "not attainable",
                   truth_sv5=f"{tpair[0]}/{tpair[1]}" if tpair else "NA", truth_attainable=attainable,
                   truth_pair_new_rank=ranks[tpair]["rank"] if tpair in ranks else "NA", truth_pair_in_new_32=tpair in new_top,
                   truth_pair_old_rank=next((r["rank"] for r in orig if tpair and (int(r["allele_x"]), int(r["allele_y"])) == tpair), "NA"),
                   m9d_improving_new_ranks=";".join(f"{x}/{y}:{ranks[(x, y)]['rank']}" for x, y in improving.get(cell, [])) or "n/a",
                   new_distinct_scores=su["distinct_scores"], new_top_tie=su["top_tie_class_size"], new_cutoff_tie=su["cutoff_tie_class_size"],
                   cutoff_tie_spans_boundary=su["cutoff_tie_spans_boundary"], shortlist_distinct_alleles=su["shortlist_distinct_alleles"],
                   overlap_with_index_list=len(new_top & orig_top), original_distinct_emissions=orig_distinct,
                   incumbent_pair_new_rank=ranks[inc]["rank"] if inc in ranks else "NA",
                   empirical_over_fitted_rate=su["empirical_over_fitted"], proposal_seconds=su["proposal_seconds"])
        if sc:
            inc_total = float(sc["incumbent"]["total"])
            best = max((r for k, r in sc.items() if k != "incumbent" and r["total"] != "NA"), key=lambda r: float(r["total"]), default=None)
            bl = best["label"].split("_") if best else None
            bpair = tuple(sorted((int(bl[2]), int(bl[3])))) if best else None
            imp = float(best["total"]) - inc_total if best else float("nan")
            row.update(incumbent_total_reproduced=abs(inc_total - float(sc["incumbent"]["expected_total"])) <= 1e-6 if sc["incumbent"]["expected_total"] != "NA" else "NA",
                       candidates_scored=len(sc) - 1, best_candidate=best["label"] if best else "NA", best_improvement=imp,
                       best_d_likelihood=float(best["likelihood"]) - float(sc["incumbent"]["likelihood"]) if best else "NA",
                       best_d_prior=float(best["prior"]) - float(sc["incumbent"]["prior"]) if best else "NA",
                       would_change=imp > 0, changed_sv5_to=f"{bpair[0]}/{bpair[1]}" if best and imp > 0 else "unchanged",
                       change_correct=(bpair == tpair) if (best and imp > 0 and attainable) else ("NA" if not (best and imp > 0) else "not attainable"),
                       change_is_real_sequence_change=({best["hap1_unique"], best["hap2_unique"]} != {sc["incumbent"]["hap1_unique"], sc["incumbent"]["hap2_unique"]}) if best and imp > 0 else "NA",
                       improvements_over_incumbent=sum(1 for k, r in sc.items() if k != "incumbent" and r["total"] != "NA" and float(r["total"]) > inc_total),
                       score_seconds=LED.get(f"score.{cell}", {}).get("wall_seconds", ""), score_peak_gib=LED.get(f"score.{cell}", {}).get("peak_rss_gib", ""))
        res.append(row)
    # Amendment 1 (analysis only): unordered sv:5 edit distance to the truth bytes (M9-D interval units, LOO cells only) for
    # the incumbent and the best replacement; overlap is NA where no original ranking was computed for the cell
    def ed(a, b):
        prev = list(range(len(b) + 1))
        for i, ca in enumerate(a, 1):
            cur = [i]
            for j, cb in enumerate(b, 1):
                cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
            prev = cur
        return prev[-1]
    for row in res:
        cell = row["cell"]
        if cell not in ("LOO.d1", "LOO.d3", "LOO.d6") and cell.startswith("LOO"):
            row["overlap_with_index_list"] = "NA (original ranking not computed for this cell)"
        tf = f"{D}/interval.{cell}/bs.sv:5.truth.fa"
        if not os.path.exists(tf):
            row["sv5_distance_incumbent"] = row["sv5_distance_best_change"] = "NA (no truth bytes)"
            continue
        recs, cur = {}, None
        for line in open(tf):
            if line.startswith(">"):
                cur = line[1:].split()[0]
                recs[cur] = ""
            else:
                recs[cur] += line.strip()
        t1, t2 = recs["truth1"], recs["truth2"]
        al = dump(MODEL(cell))[[r["block_id"] for r in T(f"{C}/call.{cell}/call.sequence_blocks.tsv")].index("sv:5")]
        dist = lambda x, y: min(ed(al[x], t1) + ed(al[y], t2), ed(al[y], t1) + ed(al[x], t2))
        ix, iy = (int(v) for v in row["incumbent_sv5"].split("/"))
        row["sv5_distance_incumbent"] = dist(ix, iy)
        if row.get("changed_sv5_to", "unchanged") != "unchanged":
            bx, by = (int(v) for v in row["changed_sv5_to"].split("/"))
            row["sv5_distance_best_change"] = dist(bx, by)
        else:
            row["sv5_distance_best_change"] = "unchanged"
    write(f"{out}/results.tsv", res)
    for r in res:
        print(r)
    sys.exit(0)
