#!/usr/bin/env python3
"""M10-C offline sequence-based structural evaluator study (frozen with DESIGN.md). Evaluation only: no caller, confidence
or posterior pooling reads anything here. A pairwise alignment-based match is NOT a transitive equivalence class.

  m10c.py gate <out>       aligner sanity on the frozen fixtures (determinism, reversal symmetry of the edit distance)
  m10c.py fixtures <out>   the frozen adversarial fixtures: measurements, the experimental verdict, required distinction met?
  m10c.py truth <case> <out>  whole-locus truth via the existing evaluation path (genotype_probe block-sequences, every block)
  m10c.py cases <out>      the inspected saved cases (development only; no accuracy claim)

Method (frozen): whole-locus global alignment (edlib NW via panvar::nw_edit_distance + one path walk; m10c_align_eval),
prediction = query, truth = target, both directions (forward, and plain-reversed sequences, which flips tie-breaking).
Diploid: the two WHOLE-haplotype assignments (p1-t1, p2-t2) and (p1-t2, p2-t1); the one with the smaller total forward edit
distance is used; a tie is reported.
Measurements (always reported): edits, aligned length, mismatches, inserted and deleted bases, indel runs >= 50 bp (E1),
contiguous non-match block bases >= 50 bp (E2, the benchmark's existing split), terminal indel runs, located-site results,
per-homologue net length difference.
Experimental verdict (EVENT = 50 bp: the existing configured event-size convention, used only as such):
  UNSUPPORTED  a truth haplotype is unavailable, or a terminal indel run >= EVENT (unmatched locus endpoints)
  AMBIGUOUS    the two assignments tie and their per-homologue structural flags differ, or forward and reverse disagree on
               whether a homologue has an E1 event
  DISCORDANT   some homologue has an E1 event (in both directions) or a located class-defining site differs
  CONCORDANT   otherwise
Flags: E1_E2_DIFFER (the two event detectors disagree), PLACEMENT_SENSITIVE (same E1 sizes, different positions between
directions), ASSIGNMENT_TIE, ENDPOINT_NOTE (a terminal indel run < EVENT), DOSAGE_BALANCED_ALLOCATION_DIFFERS (both
homologues carry E1 events of opposite net sign and the diploid net length difference is < EVENT).
"""
import collections
import csv
import hashlib
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
R = "/Users/davide.bolognini/panvar/tests/results"
OUT = f"{R}/m10c_sequence_evaluation_20261020"
EVENT = 50
ALIGN = f"{WT}/build/m10c_align_eval"
COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def rnd(tag, n):
    out, i = "", 0
    while len(out) < n:
        out += "".join("ACGT"[b % 4] for b in hashlib.sha256(f"{tag}:{i}".encode()).digest())
        i += 1
    return out[:n]


def mutate(s, positions, tag):
    s = list(s)
    for k, p in enumerate(positions):
        alt = [c for c in "ACGT" if c != s[p]]
        s[p] = alt[int(hashlib.sha256(f"{tag}:{k}".encode()).hexdigest(), 16) % 3]
    return "".join(s)


def write(out, name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


# ------------------------------------------------------------------------------------------------ alignment and verdict
def align(pairs, d, tag):
    """pairs: list of (id, A, B, sites). Returns {(id, direction): row}."""
    os.makedirs(f"{d}/seq", exist_ok=True)
    man = f"{d}/{tag}.manifest.tsv"
    with open(man, "w") as f:
        for pid, a, b, sites in pairs:
            pa, pb = f"{d}/seq/{pid}.A.txt", f"{d}/seq/{pid}.B.txt"
            open(pa, "w").write(a)
            open(pb, "w").write(b)
            f.write(f"{pid}\t{pa}\t{pb}\t{','.join(map(str, sites)) or '-'}\n")
    subprocess.run([ALIGN, man, f"{d}/{tag}.align.tsv"], check=True)
    return {(r["id"], r["direction"]): r for r in csv.DictReader(open(f"{d}/{tag}.align.tsv"), delimiter="\t")}


def runs(row):
    if row["indel_runs_ge50"] == "-":
        return []
    out = []
    for x in row["indel_runs_ge50"].split(";"):
        op_len, pos = x.split("@")
        op, n = op_len.split(":")
        out.append((op, int(n), int(pos)))
    return out


def term(t):
    return 0 if t == "0" else int(t.split(":")[1])


def homologue(fw, rv):
    """Per-homologue measurements and flags from the forward and reverse rows."""
    e1f, e1r = runs(fw), runs(rv)
    endpoint = max(term(fw["terminal_left"]), term(fw["terminal_right"]), term(rv["terminal_left"]), term(rv["terminal_right"]))
    e1f_core = [x for x in e1f]
    sites_bad = [s for s in fw["site_results"].split(";") if s != "-" and s.split(":")[1] != s.split(":")[2]]
    return dict(edits=int(fw["edits"]), edits_reverse=int(rv["edits"]), aln_len=int(fw["aln_len"]), mismatches=int(fw["mismatches"]),
                ins_bp=int(fw["ins_bp"]), del_bp=int(fw["del_bp"]), e1_forward=sorted((op, n) for op, n, _ in e1f_core),
                e1_reverse=sorted((op, n) for op, n, _ in e1r), e1_positions_forward=[p for _, _, p in e1f], e1_positions_reverse=[p for _, _, p in e1r],
                e2_forward=int(fw["over_threshold_bp"]), e2_reverse=int(rv["over_threshold_bp"]), endpoint=endpoint,
                terminal=f"{fw['terminal_left']}|{fw['terminal_right']}", sites_bad=sites_bad,
                net=int(fw["len_a"]) - int(fw["len_b"]), max_run=int(fw["max_indel_run"]))


def diploid(rows, pid, truth_ok=(True, True)):
    """Verdict for prediction pid over the two whole-haplotype assignments (pairs pid.11, .22, .12, .21)."""
    if not all(truth_ok):
        return dict(verdict="UNSUPPORTED", flags="TRUTH_UNAVAILABLE", assignment="-")
    h = {k: homologue(rows[(f"{pid}.{k}", "forward")], rows[(f"{pid}.{k}", "reverse")]) for k in ("11", "22", "12", "21")}
    t1, t2 = h["11"]["edits"] + h["22"]["edits"], h["12"]["edits"] + h["21"]["edits"]
    best = ("11", "22") if t1 <= t2 else ("12", "21")
    other = ("12", "21") if best == ("11", "22") else ("11", "22")
    flags = []
    tie = t1 == t2

    def struct(k):
        x = h[k]
        return (bool(x["e1_forward"]), bool(x["e1_reverse"]), bool(x["sites_bad"]))
    if tie:
        flags.append("ASSIGNMENT_TIE")
    hb = [h[k] for k in best]
    verdict = "CONCORDANT"
    if any(x["endpoint"] >= EVENT for x in hb):
        verdict = "UNSUPPORTED"
        flags.append("ENDPOINT_UNMATCHED")
    elif tie and sorted(struct(k) for k in best) != sorted(struct(k) for k in other):
        verdict = "AMBIGUOUS"
    elif any(bool(x["e1_forward"]) != bool(x["e1_reverse"]) for x in hb):
        verdict = "AMBIGUOUS"
        flags.append("DIRECTION_CHANGES_EVENTS")
    elif any(x["e1_forward"] for x in hb) or any(x["sites_bad"] for x in hb):
        verdict = "DISCORDANT"
    if any(x["endpoint"] and x["endpoint"] < EVENT for x in hb):
        flags.append("ENDPOINT_NOTE")
    if any(bool(x["e1_forward"]) != (x["e2_forward"] > 0) for x in hb):
        flags.append("E1_E2_DIFFER")
    if any(x["e1_forward"] == x["e1_reverse"] and x["e1_positions_forward"] != x["e1_positions_reverse"] for x in hb):
        flags.append("PLACEMENT_SENSITIVE")
    if any(x["e1_forward"] != x["e1_reverse"] for x in hb):
        flags.append("E1_SIZES_DIRECTION_DEPENDENT")
    signs = [sum((n if op == "I" else -n) for op, n in x["e1_forward"]) for x in hb]
    if all(signs) and signs[0] * signs[1] < 0 and abs(hb[0]["net"] + hb[1]["net"]) < EVENT:
        flags.append("DOSAGE_BALANCED_ALLOCATION_DIFFERS")
    if any(x["sites_bad"] for x in hb):
        flags.append("CLASS_SITE_DISCORDANT")
    return dict(verdict=verdict, flags=",".join(flags) or "-", assignment="p1-t1,p2-t2" if best == ("11", "22") else "p1-t2,p2-t1",
                total_edits_best=min(t1, t2), total_edits_other=max(t1, t2),
                **{f"h{i + 1}_{k}": (";".join(map(str, v)) if isinstance(v, list) else v) for i, x in enumerate(hb) for k, v in x.items()
                   if k in ("edits", "edits_reverse", "mismatches", "ins_bp", "del_bp", "e1_forward", "e1_reverse", "e1_positions_forward",
                            "e1_positions_reverse", "e2_forward", "e2_reverse", "terminal", "net", "sites_bad", "max_run")})


def diploid_pairs(pid, p, t, sites=((), ())):
    return [(f"{pid}.11", p[0], t[0], sites[0]), (f"{pid}.22", p[1], t[1], sites[1]),
            (f"{pid}.12", p[0], t[1], sites[1]), (f"{pid}.21", p[1], t[0], sites[0])]


# ------------------------------------------------------------------------------------------------ Amendment 1: method B
# Method A's path cannot locate events (unit-cost gaps fragment: F03's 55-bp deletion has longest run 8). Method B uses the
# existing panvar::minimap2_hits (asm20; affine gaps; all regions, unflagged) in both roles. Its frozen event rule, per role:
#   no hit -> NO_ALIGNMENT (UNSUPPORTED); a reverse-strand hit spanning >= EVENT query bases -> INVERSION;
#   chain = forward hits in decreasing matches, dropping any that overlaps a kept hit by >= EVENT query bases (those are
#   ALTERNATIVE_PLACEMENT, a flag); in the chain ordered by query start: a target position decrease -> ORDER, a target
#   overlap >= EVENT -> DUPLICATION, a query or target gap >= EVENT between neighbours -> GAP; any CIGAR I/D op >= EVENT in a
#   chain hit -> INDEL; uncovered query or target prefix/suffix >= EVENT -> ENDPOINT_UNMATCHED (< EVENT: ENDPOINT_NOTE).
# Per homologue: UNSUPPORTED (no alignment / endpoint) > AMBIGUOUS (the two roles disagree on whether an event exists) >
# DISCORDANT (an event, or a located-site difference from method A's path) > CONCORDANT. The assignment is method A's.
def align_b(pairs, d, tag):
    man = f"{d}/{tag}.manifest.tsv"
    subprocess.run([ALIGN, "--minimap2", man, f"{d}/{tag}.minimap2.tsv"], check=True)
    out = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{d}/{tag}.minimap2.tsv"), delimiter="\t"):
        out[(r["id"], r["role"])].append(r)
    return out


def events_b(hits, qlen, tlen):
    if not hits or hits[0]["hit"] == "none":
        return {"NO_ALIGNMENT"}, set(), 0
    ev, flags = set(), set()
    H = [dict(strand=h["strand"], qs=int(h["qs"]), qe=int(h["qe"]), ts=int(h["ts"]), te=int(h["te"]), m=int(h["matches"]), big=h["big_ops"]) for h in hits]
    if any(h["strand"] == "-" and h["qe"] - h["qs"] >= EVENT for h in H):
        ev.add("INVERSION")
    chain = []
    for h in sorted((h for h in H if h["strand"] == "+"), key=lambda h: -h["m"]):
        if any(min(h["qe"], c["qe"]) - max(h["qs"], c["qs"]) >= EVENT for c in chain):
            flags.add("ALTERNATIVE_PLACEMENT")
            continue
        chain.append(h)
    if not chain:
        return ev | {"NO_FORWARD_ALIGNMENT"}, flags, 0
    chain.sort(key=lambda h: h["qs"])
    for x, y in zip(chain, chain[1:]):
        if y["ts"] < x["ts"]:
            ev.add("ORDER")
        if x["te"] - y["ts"] >= EVENT:
            ev.add("DUPLICATION")
        if y["qs"] - x["qe"] >= EVENT or y["ts"] - x["te"] >= EVENT:
            ev.add("GAP")
    if any(h["big"] != "-" for h in chain):
        ev.add("INDEL")
    ends = [chain[0]["qs"], qlen - chain[-1]["qe"], min(h["ts"] for h in chain), tlen - max(h["te"] for h in chain)]
    if max(ends) >= EVENT:
        ev.add("ENDPOINT_UNMATCHED")
    elif max(ends) > 0:
        flags.add("ENDPOINT_NOTE")
    return ev, flags, max(ends)


def diploid_b(rows_a, hits, pid, truth_ok=(True, True)):
    if not all(truth_ok):
        return dict(verdict="UNSUPPORTED", flags="TRUTH_UNAVAILABLE")
    tot = {k: int(rows_a[(f"{pid}.{k}", "forward")]["edits"]) for k in ("11", "22", "12", "21")}
    t1, t2 = tot["11"] + tot["22"], tot["12"] + tot["21"]
    best, other = (("11", "22"), ("12", "21")) if t1 <= t2 else (("12", "21"), ("11", "22"))

    def hom(k):
        fw = rows_a[(f"{pid}.{k}", "forward")]
        qlen, tlen = int(fw["len_a"]), int(fw["len_b"])
        e1, f1, _ = events_b(hits.get((f"{pid}.{k}", "pred_on_truth"), []), qlen, tlen)
        e2, f2, _ = events_b(hits.get((f"{pid}.{k}", "truth_on_pred"), []), tlen, qlen)
        sites = [x for x in fw["site_results"].split(";") if x != "-" and x.split(":")[1] != x.split(":")[2]]
        return dict(ev_pred_on_truth=sorted(e1), ev_truth_on_pred=sorted(e2), flags=sorted(f1 | f2), sites=sites,
                    struct=(bool(e1 - {"ENDPOINT_UNMATCHED"}), bool(e2 - {"ENDPOINT_UNMATCHED"}), bool(sites)), net=qlen - tlen)
    hb, ho = [hom(k) for k in best], [hom(k) for k in other]
    flags = set(f for x in hb for f in x["flags"])
    if t1 == t2:
        flags.add("ASSIGNMENT_TIE")
    allev = [set(x["ev_pred_on_truth"]) | set(x["ev_truth_on_pred"]) for x in hb]
    if any({"NO_ALIGNMENT", "NO_FORWARD_ALIGNMENT"} & e for e in allev) and not any({"INVERSION"} & e for e in allev):
        verdict = "UNSUPPORTED"
    elif any("ENDPOINT_UNMATCHED" in e for e in allev):
        verdict = "UNSUPPORTED"
        flags.add("ENDPOINT_UNMATCHED")
    elif t1 == t2 and sorted(x["struct"] for x in hb) != sorted(x["struct"] for x in ho):
        verdict = "AMBIGUOUS"
    elif any(x["struct"][0] != x["struct"][1] for x in hb):
        verdict = "AMBIGUOUS"
        flags.add("ROLES_DISAGREE")
    elif any(x["struct"][0] or x["sites"] for x in hb):
        verdict = "DISCORDANT"
    else:
        verdict = "CONCORDANT"
    if any(x["sites"] for x in hb):
        flags.add("CLASS_SITE_DISCORDANT")
    if all(x["struct"][0] for x in hb) and hb[0]["net"] * hb[1]["net"] < 0 and abs(hb[0]["net"] + hb[1]["net"]) < EVENT:
        flags.add("DOSAGE_BALANCED_ALLOCATION_DIFFERS")
    return dict(verdict=verdict, flags=",".join(sorted(flags)) or "-", assignment="p1-t1,p2-t2" if best == ("11", "22") else "p1-t2,p2-t1",
                total_edits_best=min(t1, t2), total_edits_other=max(t1, t2),
                h1_events=";".join(sorted(allev[0])) or "-", h2_events=";".join(sorted(allev[1])) or "-",
                h1_edits=tot[best[0]], h2_edits=tot[best[1]], h1_net=hb[0]["net"], h2_net=hb[1]["net"])


B_REQUIRED = ("F01", "F02", "F03", "F03b", "F04", "F05", "F06", "F07", "F07b", "F08", "F08b", "F09", "F10b", "F10c", "F11", "F11b", "F14")


def fixtures_b(out):
    F = fixtures_def_a1()
    pairs = []
    for fid, _, p, t, sites, _ in F:
        pairs += diploid_pairs(fid, p, t, sites or ((), ()))
    rows = align(pairs, out, "fixtures")
    hits = align_b(pairs, out, "fixtures")
    res = []
    for fid, desc, p, t, sites, req in F:
        v = diploid_b(rows, hits, fid)
        met = "n/a (definition required)" if req == "DEFINITION_REQUIRED" else ("yes" if v["verdict"] == req else "NO")
        res.append(dict(fixture=fid, description=desc, required=req, verdict=v["verdict"], required_met=met, **{k: x for k, x in v.items() if k != "verdict"}))
    write(out, "fixtures_b", res)
    for r in res:
        print(f"{r['fixture']}\t{r['required']}\t{r['verdict']}\t{r['required_met']}\t{r['flags']}\t{r['h1_events']}|{r['h2_events']}")
    supported = all(r["required_met"] == "yes" for r in res if r["fixture"] in B_REQUIRED)
    open(f"{out}/method_b_supported.txt", "w").write(f"{supported}\n")
    print(f"method B meets every structural and concordance fixture ({','.join(B_REQUIRED)}): {supported}")
    return 0


def fixtures_def_a1():
    F = fixtures_def()
    L, Rf, U = rnd("L", 800), rnd("R", 800), rnd("U", 600)
    copy = [mutate(U, [100 + 37 * k, 400 + 11 * k], f"copy{k}") for k in range(4)]
    base = L + copy[0] + copy[1] + Rf
    site = len(L) + 250
    F.append(("F10c", "1-bp change at a class site declared for the locus (on both truth haplotypes)", (mutate(base, [site], "class"), base),
              (base, base), ([site], [site]), "DISCORDANT"))
    return F


# ------------------------------------------------------------------------------------------------ Amendment 2: method C
# Method B fails F04/F05/F06 (minimap2 aligns only part of a sequence with repeated units) and F07b. Method C uses the
# existing panvar::poa_msa (abPOA, global, affine gaps) as a two-row alignment, forward and on plain-reversed sequences;
# events are gap runs >= EVENT (contiguous by the affine model); the verdict rule is method A's with these runs; the
# assignment remains method A's exact edit distance. Inversions are not representable in a single-strand global alignment.
def align_c(d, tag):
    subprocess.run([ALIGN, "--poa", f"{d}/{tag}.manifest.tsv", f"{d}/{tag}.poa.tsv"], check=True)
    return {(r["id"], r["direction"]): r for r in csv.DictReader(open(f"{d}/{tag}.poa.tsv"), delimiter="\t")}


def diploid_c(rows_a, rows_c, pid, truth_ok=(True, True)):
    if not all(truth_ok):
        return dict(verdict="UNSUPPORTED", flags="TRUTH_UNAVAILABLE")
    conv = {}
    for k, r in rows_c.items():
        conv[k] = dict(r, indel_runs_ge50=r["gap_runs_ge50"], edits=rows_a[k]["edits"], aln_len=r["columns"],
                       over_threshold_bp=str(sum(int(x.split(":")[1].split("@")[0]) for x in r["gap_runs_ge50"].split(";")) if r["gap_runs_ge50"] != "-" else 0),
                       max_indel_run=r["max_gap_run"])
    return diploid(conv, pid, truth_ok)


C_REQUIRED = B_REQUIRED


def fixtures_c(out):
    F = fixtures_def_a1()
    pairs = []
    for fid, _, p, t, sites, _ in F:
        pairs += diploid_pairs(fid, p, t, sites or ((), ()))
    rows = align(pairs, out, "fixtures")
    rc_ = align_c(out, "fixtures")
    res = []
    for fid, desc, p, t, sites, req in F:
        v = diploid_c(rows, rc_, fid)
        met = "n/a (definition required)" if req == "DEFINITION_REQUIRED" else ("yes" if v["verdict"] == req else "NO")
        res.append(dict(fixture=fid, description=desc, required=req, verdict=v["verdict"], required_met=met, **{k: x for k, x in v.items() if k != "verdict"}))
    write(out, "fixtures_c", res)
    for r in res:
        print(f"{r['fixture']}\t{r['required']}\t{r['verdict']}\t{r['required_met']}\t{r['flags']}\t{r.get('h1_e1_forward', '-')}|{r.get('h2_e1_forward', '-')}")
    full = all(r["required_met"] == "yes" for r in res if r["fixture"] in C_REQUIRED)
    narrow = all(r["required_met"] == "yes" for r in res if r["fixture"] in C_REQUIRED and r["fixture"] not in ("F07", "F07b"))
    status = "FULL" if full else "INVERSION_BLIND" if narrow else "NONE"
    open(f"{out}/method_c_supported.txt", "w").write(status + "\n")
    print(f"method C support: {status}")
    return 0


# ------------------------------------------------------------------------------------------------ frozen fixtures
def fixtures_def():
    L, Rf, U = rnd("L", 800), rnd("R", 800), rnd("U", 600)
    copy = [mutate(U, [100 + 37 * k, 400 + 11 * k], f"copy{k}") for k in range(4)]   # copy-specific variants
    A, B = rnd("A", 600), rnd("B", 600)
    base = L + copy[0] + copy[1] + Rf
    X55, Y55 = rnd("X55", 55), rnd("Y55", 55)
    mid = len(L) + 300
    F = []   # id, description, prediction (p1, p2), truth (t1, t2), sites (per truth haplotype), required verdict
    F.append(("F01", "identical DNA (different node splits / block boundaries are representation only)", (base, base), (base, base), None, "CONCORDANT"))
    neutral = mutate(base, [150, 900, 1500, 2100], "neutral")
    neutral = neutral[:1000] + neutral[1002:1700] + "GTC" + neutral[1700:]
    F.append(("F02", "neutral substitutions (4) and small indels (2-bp deletion, 3-bp insertion)", (neutral, base), (base, base), None, "CONCORDANT"))
    F.append(("F03", "55-bp deletion (in a graph: many short nodes; sequence-level here)", (base[:mid] + base[mid + 55:], base), (base, base), None, "DISCORDANT"))
    rep = L + "CA" * 60 + copy[0] + Rf
    rep_del = L + "CA" * 32 + copy[0] + Rf   # 56 bp of a dinucleotide repeat removed: placement ambiguous, size fixed
    F.append(("F03b", "56-bp deletion inside a (CA)n repeat (placement ambiguous)", (rep_del, rep), (rep, rep), None, "DISCORDANT"))
    cn3 = L + copy[0] + copy[1] + copy[2] + Rf
    F.append(("F04", "repeat copy number 3 vs 2 on one homologue", (cn3, base), (base, base), None, "DISCORDANT"))
    cn1 = L + copy[0] + Rf
    F.append(("F05", "equal total dosage, 3+1 vs 2+2", (cn3, cn1), (base, base), None, "DISCORDANT"))
    F.append(("F06", "order change A,B vs B,A (600 bp each)", (L + B + A + Rf, base), (L + A + B + Rf, base), None, "DISCORDANT"))
    F.append(("F07", "600-bp inversion", (L + rc(A) + Rf, base), (L + A + Rf, base), None, "DISCORDANT"))
    small = rnd("inv60", 60)
    F.append(("F07b", "60-bp inversion", (L + rc(small) + Rf, base), (L + small + Rf, base), None, "DISCORDANT"))
    same3 = L + U + U + U + Rf
    F.append(("F08", "shifted representation of one tandem expansion (identical copies; same DNA)", (same3, base), (same3, base), None, "CONCORDANT"))
    msat_l, msat_r = L + "AC" * 34 + Rf, L + "AC" * 34 + Rf
    F.append(("F08b", "microsatellite expansion placed left vs right (same DNA)", (msat_l, base), (msat_r, base), None, "CONCORDANT"))
    F.append(("F09", "novel small differences absent from any panel (6 SNPs)", (mutate(base, [200, 777, 1234, 1600, 1999, 2500], "novel"), base),
              (base, base), None, "CONCORDANT"))
    site = len(L) + 250
    cls = mutate(base, [site], "class")
    F.append(("F10", "short located class-defining difference (1 bp at a declared site)", (cls, base), (base, base), ([site], []), "DISCORDANT"))
    F.append(("F10b", "the same 1-bp change at an undeclared position", (mutate(base, [site + 3], "class"), base), (base, base), None, "CONCORDANT"))
    withX = base[:mid] + X55 + base[mid:]
    F.append(("F11", "assignment tie: prediction (S, S) vs truth (S+55 bp, S)", (base, base), (withX, base), None, "DISCORDANT"))
    cn3b = L + copy[0] + copy[1] + copy[2] + Rf
    two = L + copy[0] + copy[2] + Rf   # the middle copy missing: which copy is deleted depends on placement
    F.append(("F11b", "one of three near-identical copies deleted (allocation of the deleted copy ambiguous)", (two, base), (cn3b, base), None, "DISCORDANT"))
    F.append(("F12", "unmatched locus endpoint: truth extends 200 bp beyond the right anchor", (base, base), (base + rnd("ext", 200), base), None, "UNSUPPORTED"))
    F.append(("F12b", "unmatched locus endpoint: truth lacks the first 20 bp", (base, base), (base[20:], base), None, "UNSUPPORTED"))
    F.append(("F13", "equal-length 55-bp divergent replacement (definition required)", (base[:mid] + Y55 + base[mid + 55:], base),
              (base[:mid] + X55 + base[mid + 55:], base), None, "DEFINITION_REQUIRED"))
    near = mutate(base[mid:mid + 60], [10, 30, 50], "near")
    tandem = base[:mid] + near + base[mid:]   # a 60-bp near-copy inserted next to its source
    F.append(("F14", "60-bp near-tandem duplication (3 mismatches to its source) present in truth only", (base, base), (tandem, base), None, "DISCORDANT"))
    return F


def fixtures(out, gate=False):
    F = fixtures_def()
    pairs = []
    for fid, _, p, t, sites, _ in F:
        pairs += diploid_pairs(fid, p, t, sites or ((), ()))
    rows = align(pairs, out, "fixtures")
    if gate:
        rows2 = align(pairs, out, "fixtures_repeat")
        det = all(rows[k] == rows2[k] for k in rows)
        sym = all(rows[(i, "forward")]["edits"] == rows[(i, "reverse")]["edits"] for i, _ in rows)
        zero = rows[("F01.11", "forward")]["edits"] == "0"
        print(f"gate: deterministic {det}; reversal-symmetric edit distance {sym}; identical -> 0 edits {zero}")
        return 0 if det and sym and zero else 1
    res = []
    for fid, desc, p, t, sites, req in F:
        v = diploid(rows, fid)
        met = "n/a (definition required)" if req == "DEFINITION_REQUIRED" else ("yes" if v["verdict"] == req else "NO")
        res.append(dict(fixture=fid, description=desc, required=req, verdict=v["verdict"], required_met=met, **{k: x for k, x in v.items() if k != "verdict"}))
    write(out, "fixtures", res)
    for r in res:
        print(f"{r['fixture']}\t{r['required']}\t{r['verdict']}\t{r['required_met']}\t{r['flags']}")
    return 0


# ------------------------------------------------------------------------------------------------ saved cases
M1B, M9H, M10A = f"{R}/m1b_structural_20261019/runs", f"{R}/m9h_jp11_validation_20261020/runs", f"{R}/m10a_dev_configuration_20261020/runs"
DUMPS = f"{R}/m10b1_projection_hardening_20261020/runs"
C4REF, CYPREF, ACOTREF = "HG00096#1#haplotype1-0000024:31848049-32074186", "grch38#1#chr22:42031864-42245566", "GRCh38#0#chr14:73416151-73668669"
CYPCUR = "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"
ACOTCUR = f"{M9H}/prepare/curation.LOO.txt"
TRUTH = {   # case: model dump, model gfa, reference, curation, truth gfa (verified input), donor
    "c4_na19650": ("c4loo", f"{R}/genotype_c4_production_11fe681/graphs/NA19650.loo.gfa", C4REF, "", f"{R}/genotype_c4_production_11fe681/graphs/c4.gfa", "NA19650"),
    "c4_hg00096": ("c4dev", f"{R}/genotype_c4_production_11fe681/graphs/c4.gfa", C4REF, "", f"{R}/genotype_c4_production_11fe681/graphs/c4.gfa", "HG00096"),
    "cyp_na19650": ("cyploo", f"{R}/cyp2d6_shortread_pilot_v1/data/graphs/loo.NA19650.gfa", CYPREF, CYPCUR, f"{R}/cyp2d6_shortread_pilot_v1/data/graphs/cyp2d6.gfa", "NA19650"),
    "acot_hg02083": ("acotloo", f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", ACOTREF, ACOTCUR, f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", "HG02083"),
    "acot_hg02074": ("acotloo", f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", ACOTREF, ACOTCUR, f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", "HG02074"),
}
CASES = [   # id, truth case, predictions (label -> source)
    ("C4_NA19650_allocation", "c4_na19650", {"reported_3+1": ("alts", f"{M1B}/alts.c4_10/c4_10.candidates.tsv", "reported"),
                                            "saved_2+2": ("alts", f"{M1B}/alts.c4_10/c4_10.candidates.tsv", "saved_2+2")}),
    ("C4_HG00096_control", "c4_hg00096", {"call": ("call", f"{M10A}/D.c4_5")}),
    ("CYP_NA19650_dosage", "cyp_na19650", {"reported_5x_2+2": ("call", f"{M1B}/call.cyp_5"), "twin_30x_3+2": ("call", f"{M1B}/call.cyp_30"),
                                          "saved_2+3": ("alts", f"{M1B}/alts.cyp_5/cyp_5.candidates.tsv", "saved_2+3")}),
    ("ACOT_L4", "acot_hg02083", {"jp_on": ("call", f"{M9H}/call.on.L4"), "jp_off": ("call", f"{M9H}/call.off.L4")}),
    ("ACOT_L2", "acot_hg02074", {"jp_on": ("call", f"{M9H}/call.on.L2"), "jp_off": ("call", f"{M9H}/call.off.L2")}),
]


def kv(p):
    return dict(l.rstrip("\n").split("\t", 1) for l in open(p) if "\t" in l and not l.startswith("#"))


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def truth(case, out):
    dump, mgfa, ref, cur, tgfa, donor = TRUTH[case]
    blocks = []
    for r in csv.DictReader(open(f"{DUMPS}/dump.{dump}/{dump}.block_alleles.tsv"), delimiter="\t"):
        if not blocks or blocks[-1] != r["block_id"]:
            blocks.append(r["block_id"])
    cmd = [f"{WT}/build/genotype_probe", "block-sequences", "--out", f"{out}/bs", "--reference", ref, "--model-gfa", mgfa,
           "--blocks", ",".join(blocks), "--truth-gfa", tgfa, "--donor", donor] + (["--panel-curation", cur] if cur else [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    print(r.stdout + r.stderr)
    fp = kv(f"{DUMPS}/dump.{dump}/{dump}.summary.tsv")["model_fingerprint"]
    ok = r.returncode == 0 and f"model fingerprint {fp}" in r.stdout
    open(f"{out}/provenance.tsv", "w").write(f"key\tvalue\ncase\t{case}\ndonor\t{donor}\ntruth_gfa\t{tgfa}\ntruth_gfa_sha256\t{sha(tgfa)}\n"
                                             f"model_gfa\t{mgfa}\nmodel_fingerprint\t{fp}\nfingerprint_matches_dump\t{ok}\n")
    return 0 if ok else 1


def model_blocks(dump):
    seq, ids = {}, []
    for r in csv.DictReader(open(f"{DUMPS}/dump.{dump}/{dump}.block_alleles.tsv"), delimiter="\t"):
        seq[(int(r["block_index"]), int(r["allele_index"]))] = r["sequence"]
        if len(ids) <= int(r["block_index"]):
            ids.append(r["block_id"])
    return seq, ids


def fasta(p):
    out, name = {}, None
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line)
    return {k: "".join(v) for k, v in out.items()}


def predictions(src):
    if src[0] == "call":
        sl = kv(f"{src[1]}/call.sequence_locus.tsv")
        return [list(map(int, sl[f"representative_labelling_hap{h}"].split(","))) for h in (1, 2)]
    for line in open(src[1]):
        f = line.rstrip("\n").split("\t")
        if f[0] == src[2]:
            return [list(map(int, f[1].split(","))), list(map(int, f[2].split(",")))]
    raise KeyError(src)


def cases(out):
    res, prov = [], []
    for cid, tcase, preds in CASES:
        dump = TRUTH[tcase][0]
        seq, ids = model_blocks(dump)
        tdir = f"{OUT}/runs/truth.{tcase}.a2"   # Amendment 2: truth extraction runs as .a2
        th = {"truth1": [], "truth2": []}
        avail = {"truth1": True, "truth2": True}
        for b, bid in enumerate(ids):
            recs = fasta(f"{tdir}/bs.{bid}.truth.fa") if os.path.exists(f"{tdir}/bs.{bid}.truth.fa") else {}
            for k in th:
                if k in recs:
                    th[k].append(recs[k])
                else:
                    avail[k] = False
                    th[k].append(None)
        for scope, lo, hi in (("interior (matched chain endpoints)", 1, len(ids) - 1), ("with flank blocks (path extents)", 0, len(ids))):
            t = tuple("".join(x for x in th[k][lo:hi]) if all(x is not None for x in th[k][lo:hi]) else "" for k in ("truth1", "truth2"))
            ok = tuple(all(x is not None for x in th[k][lo:hi]) for k in ("truth1", "truth2"))
            for label, src in preds.items():
                tup = predictions(src)
                p = tuple("".join(seq[(b, tup[h][b])] for b in range(lo, hi)) for h in (0, 1))
                pid = f"{cid}.{label}.{'int' if lo else 'full'}"
                rows = align(diploid_pairs(pid, p, t), out, pid) if all(ok) else {}
                v = diploid(rows, pid, ok)
                vb = diploid_c(rows, align_c(out, pid), pid, ok) if all(ok) else dict(verdict="UNSUPPORTED", flags="TRUTH_UNAVAILABLE")   # A2: method C
                res.append(dict(case=cid, prediction=label, scope=scope, truth_available=f"{ok[0]}/{ok[1]}",
                                pred_bp=f"{len(p[0])}/{len(p[1])}", truth_bp=f"{len(t[0])}/{len(t[1])}",
                                method_c_verdict=vb["verdict"], method_c_flags=vb["flags"], method_c_h1_runs=vb.get("h1_e1_forward", "-"),
                                method_c_h2_runs=vb.get("h2_e1_forward", "-"), method_c_h1_runs_reverse=vb.get("h1_e1_reverse", "-"),
                                method_c_h2_runs_reverse=vb.get("h2_e1_reverse", "-"), method_a_verdict=v["verdict"], method_a_flags=v["flags"],
                                **{k: x for k, x in v.items() if k not in ("verdict", "flags")}))
        prov.append(dict(case=cid, **kv(f"{tdir}/provenance.tsv")))
    write(out, "cases", res)
    write(out, "provenance", prov)
    for r in res:
        print(f"{r['case']}\t{r['prediction']}\t{r['scope'][:8]}\tC={r['method_c_verdict']}\t{r['method_c_flags']}\t{r['method_c_h1_runs']}|{r['method_c_h2_runs']}\tA={r['method_a_verdict']}\t{r.get('total_edits_best', '-')}")
    return 0


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "truth":
        os.makedirs(sys.argv[3], exist_ok=True)
        sys.exit(truth(sys.argv[2], sys.argv[3]))
    out = sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"gate": lambda o: fixtures(o, gate=True), "fixtures": fixtures, "fixtures_b": fixtures_b, "fixtures_c": fixtures_c, "cases": cases}[mode](out))
