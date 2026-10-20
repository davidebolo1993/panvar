#!/usr/bin/env python3
"""M10-D: COSIGT-style sequence-similarity evaluation with the existing exact edit distance (frozen with DESIGN.md).
Evaluation only: nothing here changes predictions, GQ, confidence, filters, priors or caller code, and no categorical cutoff,
cluster or posterior pooling is defined.

Contract, for predicted haplotypes A, B and truth C, D over one evaluated region with matched endpoints:
  D1 = d(A,C) + d(B,D);  D2 = d(A,D) + d(B,C);  D = min(D1, D2), one assignment for the whole region, chosen by summed edits.
  Tie (D1 == D2): assignment 1 (A-C, B-D), recorded as a tie.
  d = exact global (NW) edit distance, panvar::nw_edit_distance (edlib); L1, L2 = its alignment-column counts (the
  benchmark's S) for the selected assignment, from edlib's deterministic path (co-optimal paths may differ in column count:
  a mismatch is one column, an insertion/deletion pair two; recorded, not resolved).
  error_rate = D / (L1 + L2) (unfloored; exactly 0 for identical sequences; NA when L1 + L2 = 0)
  identity   = 1 - error_rate
  QV         = -10 log10(max(0.5, D) / (L1 + L2)); 60 when L1 + L2 = 0 -- the benchmark's qv_from (benchmark_command.cpp),
               combined length-weighted (sum delta / sum S) exactly as the benchmark combines haplotypes.

  m10d.py gate <out>      contract checks: fixtures (M10-C sequences reused + missing checks), tiny-sequence DP oracle,
                          normalization, empty and unsupported handling
  m10d.py truth <case> <out>   whole-locus truth units via the existing evaluation path (M10-C's truth(); block-sequences)
  m10d.py cases <out>     the frozen development cases
"""
import collections
import csv
import itertools
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m10c"))
import m10c  # noqa: E402  (frozen M10-C fixtures, case definitions, truth extraction)

WT, R = m10c.WT, m10c.R
OUT = f"{R}/m10d_sequence_similarity_20261020"
ALIGN = m10c.ALIGN
K = 31


def write(out, name, rows):
    m10c.write(out, name, rows)


def distances(pairs, d, tag):
    """pairs: (id, A, B). Exact edit distance and alignment columns per pair (m10c_align_eval --distance)."""
    os.makedirs(f"{d}/seq", exist_ok=True)
    man = f"{d}/{tag}.manifest.tsv"
    with open(man, "w") as f:
        for pid, a, b in pairs:
            pa, pb = f"{d}/seq/{pid}.A.txt", f"{d}/seq/{pid}.B.txt"
            open(pa, "w").write(a)
            open(pb, "w").write(b)
            f.write(f"{pid}\t{pa}\t{pb}\t-\n")
    subprocess.run([ALIGN, "--distance", man, f"{d}/{tag}.distance.tsv"], check=True)
    return {r["id"]: r for r in csv.DictReader(open(f"{d}/{tag}.distance.tsv"), delimiter="\t")}


def qv(delta, S):   # the benchmark's qv_from, verbatim
    if S == 0:
        return 60.0
    return -10.0 * math.log10(max(0.5, float(delta)) / float(S))


def diploid_pairs(pid, p, t):
    return [(f"{pid}.AC", p[0], t[0]), (f"{pid}.BD", p[1], t[1]), (f"{pid}.AD", p[0], t[1]), (f"{pid}.BC", p[1], t[0])]


def contract(rows, pid):
    e = {k: int(rows[f"{pid}.{k}"]["edits"]) for k in ("AC", "BD", "AD", "BC")}
    L = {k: int(rows[f"{pid}.{k}"]["aln_len"]) for k in ("AC", "BD", "AD", "BC")}
    D1, D2 = e["AC"] + e["BD"], e["AD"] + e["BC"]
    sel = ("AC", "BD") if D1 <= D2 else ("AD", "BC")
    D = min(D1, D2)
    S = L[sel[0]] + L[sel[1]]
    err = D / S if S else None
    return dict(D1=D1, D2=D2, selected="1 (A-C, B-D)" if sel == ("AC", "BD") else "2 (A-D, B-C)", tie=D1 == D2, D=D,
                h1_edits=e[sel[0]], h2_edits=e[sel[1]], L1=L[sel[0]], L2=L[sel[1]],
                error_rate="NA" if err is None else repr(err), identity="NA" if err is None else repr(1 - err), QV=f"{qv(D, S):.4f}",
                status="EMPTY" if S == 0 else "EVALUATED")


def lev(a, b):   # independent oracle: textbook DP
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
        prev = cur
    return prev[-1]


def kmers(s):
    return {s[i:i + K] for i in range(len(s) - K + 1)}


def orientation(pred, truth):
    """Whole-sequence orientation check: truth 31-mers found in the prediction as given vs in its reverse complement."""
    tk, pk, rk = kmers(truth), kmers(pred), kmers(m10c.rc(pred))
    f, r = len(tk & pk), len(tk & rk)
    return f, r, "AS_GIVEN" if f >= r else "REVERSE_COMPLEMENT_SHARES_MORE"


# ------------------------------------------------------------------------------------------------ gate
def gate(out):
    checks = []

    def check(name, ok, detail=""):
        checks.append(dict(check=name, ok=ok, detail=detail))
        print(("  PASS " if ok else "  FAIL ") + name + (f"  [{detail}]" if detail else ""))

    F = {fid: (p, t) for fid, _, p, t, _, _ in m10c.fixtures_def_a1()}
    extra = {}
    p6, t6 = F["F06"]
    extra["SWAP"] = (p6, (t6[1], t6[0]))                                # swapped truth homologues
    p2, t2 = F["F02"]
    extra["RC"] = (p2, (m10c.rc(t2[0]), m10c.rc(t2[1])))                 # whole-sequence reverse-complemented truth
    extra["E1"] = (("ACGTACGT", ""), ("", ""))                          # one empty truth pair
    extra["E2"] = (("", ""), ("", ""))                                  # everything empty
    extra["E3"] = (("", "ACGT"), ("ACGT", ""))                          # empty haplotypes matched by the assignment
    pairs = []
    for fid, (p, t) in list(F.items()) + list(extra.items()):
        pairs += diploid_pairs(fid, p, t)
    rows = distances(pairs, out, "fixtures")
    res = {fid: contract(rows, fid) for fid in list(F) + list(extra)}
    table = []
    for fid, v in res.items():
        p, t = (F.get(fid) or extra[fid])
        table.append(dict(fixture=fid, pred_bp=f"{len(p[0])}/{len(p[1])}", truth_bp=f"{len(t[0])}/{len(t[1])}", **v))
    write(out, "fixtures_contract", table)
    L0 = res["F01"]["L1"] + res["F01"]["L2"]
    check("identical DNA (node/block representation irrelevant): D = 0, error 0 exactly, QV floored",
          res["F01"]["D"] == 0 and res["F01"]["error_rate"] == "0.0" and res["F01"]["QV"] == f"{qv(0, L0):.4f}")
    check("swapped homologues: same D, the other assignment selected", res["SWAP"]["D"] == res["F06"]["D"] and
          res["SWAP"]["selected"] != res["F06"]["selected"], f"D {res['F06']['D']}")
    check("whole-sequence reverse complement: much larger D, and the orientation check prefers the given strand",
          res["RC"]["D"] > 10 * max(1, res["F02"]["D"]) and orientation(p2[0], t2[0])[2] == "AS_GIVEN" and
          orientation(p2[0], m10c.rc(t2[0]))[2] == "REVERSE_COMPLEMENT_SHARES_MORE", f"D {res['F02']['D']} vs rc {res['RC']['D']}")
    check("small substitutions and indels: 0 < D <= constructed 9 edits", 0 < res["F02"]["D"] <= 9, f"D {res['F02']['D']}")
    for fid, exact in (("F03", 55), ("F04", 600), ("F05", 1200), ("F11", 55), ("F14", 60)):
        p, t = F[fid]
        check(f"{fid} unequal lengths: D equals the length-difference lower bound {exact}", res[fid]["D"] == exact, f"D {res[fid]['D']}")
    check("F05 3+1 vs 2+2: D = 1200 although total length is equal (allocation counted by the whole-haplotype assignment)",
          res["F05"]["D"] == 1200 and len(F["F05"][0][0]) + len(F["F05"][0][1]) == len(F["F05"][1][0]) + len(F["F05"][1][1]))
    check("F06 order change: exact distance reported (0 < D <= two 600-bp moves)", 0 < res["F06"]["D"] <= 2400, f"D {res['F06']['D']}")
    check("F07 / F07b inversions: exact distance reported (0 < D <= inverted length)", 0 < res["F07"]["D"] <= 600 and 0 < res["F07b"]["D"] <= 60,
          f"D {res['F07']['D']} / {res['F07b']['D']}")
    check("assignment tie (F11): recorded, assignment 1 selected", res["F11"]["tie"] and res["F11"]["selected"].startswith("1"))
    check("empty truth pair (8 bp vs empty): D = 8, L = 8", res["E1"]["D"] == 8 and res["E1"]["L1"] + res["E1"]["L2"] == 8, f"D {res['E1']['D']}")
    check("all empty: status EMPTY, error NA, QV 60 (benchmark convention)", res["E2"]["status"] == "EMPTY" and res["E2"]["error_rate"] == "NA"
          and res["E2"]["QV"] == "60.0000")
    check("empty haplotypes matched by the assignment: D = 0 via assignment 2", res["E3"]["D"] == 0 and res["E3"]["selected"].startswith("2"))
    cols_ok = all(max(int(r["len_a"]), int(r["len_b"])) <= int(r["aln_len"]) <= int(r["len_a"]) + int(r["len_b"]) and
                  int(r["aln_len"]) >= int(r["edits"]) for r in rows.values())
    check("alignment columns within [max(len), len_a + len_b] and >= edits for every pair (no clipping)", cols_ok)
    # independent oracle on tiny sequences, including tiny order changes and inversions
    tiny, ids = [], []
    for i in range(400):
        a = m10c.rnd(f"ta{i}", i % 10)
        b = m10c.rnd(f"tb{i}", (i * 7) % 11)
        if i % 5 == 0:
            b = a[::-1].translate(str.maketrans("ACGT", "TGCA"))    # inversion of the whole tiny sequence
        if i % 5 == 1 and len(a) >= 6:
            b = a[3:] + a[:3]                                       # order change
        tiny.append((f"T{i}", a, b))
    trows = distances(tiny, out, "oracle")
    bad = [pid for pid, a, b in tiny if int(trows[pid]["edits"]) != lev(a, b)]
    check("tiny-sequence oracle: nw_edit_distance equals textbook DP on 400 pairs (incl. inversions, order changes, empty)",
          not bad, f"{len(bad)} mismatches")
    dip = []
    for i in range(60):
        hs = [m10c.rnd(f"d{i}{k}", 3 + (i + k) % 6) for k in range(4)]
        dip += diploid_pairs(f"P{i}", (hs[0], hs[1]), (hs[2], hs[3]))
    drows = distances(dip, out, "oracle_diploid")
    okd = True
    for i in range(60):
        hs = [m10c.rnd(f"d{i}{k}", 3 + (i + k) % 6) for k in range(4)]
        brute = min(lev(hs[0], hs[2]) + lev(hs[1], hs[3]), lev(hs[0], hs[3]) + lev(hs[1], hs[2]))
        okd &= contract(drows, f"P{i}")["D"] == brute
    check("diploid oracle: D equals the brute-force minimum over both whole-haplotype assignments (60 cases)", okd)
    check("normalization: QV floor 0.5, S = 0 -> 60, diploid QV = benchmark length-weighted combination",
          abs(qv(0, 1000) - 33.0103) < 1e-4 and qv(5, 0) == 60.0 and abs(qv(3 + 4, 100 + 200) - (-10 * math.log10(7 / 300))) < 1e-12)
    check("unsupported comparison (truth unavailable) yields no metric", True, "exercised in cases: status UNSUPPORTED, metrics NA")
    write(out, "gate_checks", checks)
    ok = all(c["ok"] for c in checks)
    print("m10d gate: all passed" if ok else "m10d gate: FAILURES")
    return 0 if ok else 1


# ------------------------------------------------------------------------------------------------ cases
HISTORICAL = {   # earlier per-block metrics, a DIFFERENT evaluation unit; kept separately, never subtracted
    ("C4_NA19650_allocation", "reported_3+1"): "M1-B sv:7 excess 65,444 bp (block unit; availability-floor based)",
    ("C4_NA19650_allocation", "saved_2+2"): "M1-B sv:7 excess 86 bp (block unit)",
    ("CYP_NA19650_dosage", "reported_5x_2+2"): "M1-B sv:9 excess 12,159 bp (block unit)",
    ("CYP_NA19650_dosage", "twin_30x_3+2"): "M1-B sv:9 excess 12 bp (block unit); also the M9-I sv:3..sv:5 synonymous case (compound target exact)",
    ("CYP_NA19650_dosage", "saved_2+3"): "M1-B sv:9 excess 281 bp (block unit)",
    ("ACOT_L4", "jp_on"): "M9-H sv:5..ctx:6 interval distance 0, excess 540 (interval unit)",
    ("ACOT_L4", "jp_off"): "M9-H sv:5..ctx:6 interval distance 14, excess 553 (interval unit)",
    ("ACOT_L2", "jp_on"): "M9-H sv:5..ctx:6 interval distance 3, excess 32 (interval unit)",
    ("ACOT_L2", "jp_off"): "M9-H sv:5..ctx:6 interval distance 5, excess 38 (interval unit)",
}


def cases(out):
    res, prov = [], []
    for cid, tcase, preds in m10c.CASES:
        dump = m10c.TRUTH[tcase][0]
        seq, ids = m10c.model_blocks(dump)
        tdir = f"{OUT}/runs/truth.{tcase}"
        th, ok = {"truth1": [], "truth2": []}, {"truth1": True, "truth2": True}
        for bid in ids:
            p = f"{tdir}/bs.{bid}.truth.fa"
            recs = m10c.fasta(p) if os.path.exists(p) else {}
            for k in th:
                th[k].append(recs.get(k))
                ok[k] &= k in recs
        n = len(ids)
        interior = tuple("".join(x for x in th[k][1:n - 1]) if all(x is not None for x in th[k][1:n - 1]) else None for k in ("truth1", "truth2"))
        prov.append(dict(evaluated_case=cid, **m10c.kv(f"{tdir}/provenance.tsv")))   # Amendment 1: the file has its own 'case' key
        for label, src in preds.items():
            tup = m10c.predictions(src)
            p = tuple("".join(seq[(b, tup[h][b])] for b in range(1, n - 1)) for h in (0, 1))
            base = dict(case=cid, prediction=label, scope=f"chain interior: blocks {ids[1]}..{ids[n - 2]} (matched chain endpoints; "
                        "the flank blocks carry each path's own assembly extent and are excluded identically)",
                        pred_bp=f"{len(p[0])}/{len(p[1])}", historical_block_metric=HISTORICAL.get((cid, label), "-"))
            if any(x is None for x in interior):
                res.append(dict(base, status="UNSUPPORTED", reason="truth projection unavailable at a block", truth_bp="-"))
                continue
            o = [orientation(p[h], interior[h]) for h in (0, 1)]
            rows = distances(diploid_pairs(f"{cid}.{label}", p, interior), out, f"{cid}.{label}")
            v = contract(rows, f"{cid}.{label}")
            res.append(dict(base, truth_bp=f"{len(interior[0])}/{len(interior[1])}", **v,
                            orientation_check=";".join(f"h{h + 1}:{o[h][2]}({o[h][0]}/{o[h][1]})" for h in (0, 1)), reason="-"))
            # flank-inclusive: endpoints are not matched (path extents); reported as unsupported, lengths only
            full_p = tuple("".join(seq[(b, tup[h][b])] for b in range(n)) for h in (0, 1))
            full_t = tuple("".join(x for x in th[k]) if all(x is not None for x in th[k]) else "" for k in ("truth1", "truth2"))
            res.append(dict(case=cid, prediction=label, scope="flank-inclusive (whole projected paths)", status="UNSUPPORTED",
                            reason="locus endpoints not matched: flank blocks are each path's own assembly extent",
                            pred_bp=f"{len(full_p[0])}/{len(full_p[1])}", truth_bp=f"{len(full_t[0])}/{len(full_t[1])}",
                            historical_block_metric="-"))
    write(out, "cases", res)
    write(out, "provenance", prov)
    for r in res:
        if r["status"] == "EVALUATED":
            print(f"{r['case']}\t{r['prediction']}\tD={r['D']} (D1 {r['D1']}, D2 {r['D2']}, sel {r['selected'][0]}{' tie' if r['tie'] else ''})"
                  f"\tL={r['L1']}+{r['L2']}\terr={float(r['error_rate']):.3e}\tQV={r['QV']}\t{r['orientation_check']}")
        else:
            print(f"{r['case']}\t{r['prediction']}\t{r['status']}\t{r['reason']}")
    return 0


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "truth":
        os.makedirs(sys.argv[3], exist_ok=True)
        sys.exit(m10c.truth(sys.argv[2], sys.argv[3]))
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"gate": gate, "cases": cases}[mode](out))
