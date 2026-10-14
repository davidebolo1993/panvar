#!/usr/bin/env python3
"""CYP2D6/D7 pilot: prospective mosaic-donor selection and the six off-panel mosaic plans (truth construction only).

No reads are simulated and no call is made. Selection uses only structure (orientation, completeness, respelling, BED
census, CYP2D6 copy allocation) and donor names; never genotyping performance or reads.

Selection (AMENDMENT_2.md):
  rules 1-4  eligibility from preflight (two complete assigned-sample rows, exact respelling, one CYP2D7 per homologue),
             forward/forward orientation; reference assemblies and unassigned paths excluded;
  rule 5     strata = CYP2D6 H1/H2 copy allocation;
  rule 6     NA19705 first;
  rule 7     structural diversity: the remaining donors are taken round-robin over strata, rarest stratum first
             (ties by stratum label), so every stratum is represented before the majority stratum repeats;
  rule 8     within a stratum, lexicographic donor order.
  The ordered list L gives donor A of plan k as L[k] and donor B as L[k + 6] (k = 0..5); L[12:] are reserves.

Plans, in this order: one_1, both_1, three_1, one_2, both_2, three_2 (one homologue switched once; both homologues
switched once at different boundaries; three switches on each homologue at the same three boundaries). Boundary i
(1..18) means block i is the first block from the other donor. Candidate boundaries are ranked by
sha256("panvar-cyp2d6-pilot-v1:<plan>:<i>"). A plan's switch set is valid when every switched homologue is in no panel
row, every donor segment of a switched homologue contains a block where the two source rows differ (no switch is a
no-op), and the two mosaic homologues differ. one: the first valid boundary on hap1, hap2 donor A's row. both: hap1 at
the first valid boundary, hap2 at the next valid boundary different from it. three: the first valid 3-combination of
the ranked list (itertools order), used on both homologues. If a pair admits no valid plan, donor B is replaced by
the next unused reserve; if the reserves run out, the construction stops.

Usage: select_and_plan.py <genotype_probe> <truth-construction gfa> <original cyp2d6.gfa.gz> <eligibility.tsv>
                          <panel_alleles.tsv> <out dir>
"""
import csv
import gzip
import hashlib
import itertools
import os
import subprocess
import sys

REF = "grch38#1#chr22:42031864-42245566"
SALT = "panvar-cyp2d6-pilot-v1"
PLANS = [("one_1", "one"), ("both_1", "both"), ("three_1", "three"), ("one_2", "one"), ("both_2", "both"),
         ("three_2", "three")]
REQUIRED = "NA19705"
HERE = os.path.dirname(os.path.abspath(__file__))
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


class Stop(Exception):
    pass


def read_gfa(path):
    op = gzip.open if path.endswith(".gz") else open
    seq, paths = {}, {}
    with op(path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] == "S":
                seq[p[1]] = p[2]
            elif p[0] == "P":
                paths[p[1]] = p[2].split(",")
    return seq, paths


def spell(seq, steps):
    return "".join(seq[s[:-1]].translate(COMP)[::-1] if s.endswith("-") else seq[s[:-1]] for s in steps)


def panel_alleles(path):
    out = {}
    for line in open(path):
        if line.startswith(("#", "template")) or not line.strip():
            continue
        name, alleles = line.rstrip("\n").split("\t")[:2]
        out[name] = [int(x) for x in alleles.split(",")]
    return out


def rank(plan, nblocks):
    return sorted(range(1, nblocks), key=lambda i: hashlib.sha256(f"{SALT}:{plan}:{i}".encode()).hexdigest())


def donor_of_block(switches, nblocks):
    cur, out = 0, []
    for b in range(nblocks):
        if b in switches:
            cur = 1 - cur
        out.append(cur)
    return out


def haplotype(rows, switches, nblocks):
    src = donor_of_block(set(switches), nblocks)
    return [rows[src[b]][b] for b in range(nblocks)], src


def segments_effective(rows, switches, nblocks):
    """Every maximal single-donor segment contains a block where the two source rows differ."""
    cuts = [0] + sorted(switches) + [nblocks]
    return all(any(rows[0][b] != rows[1][b] for b in range(lo, hi)) for lo, hi in zip(cuts, cuts[1:]))


def valid_hap(rows, switches, nblocks, panel_rows):
    vec, _ = haplotype(rows, switches, nblocks)
    return vec not in panel_rows and segments_effective(rows, switches, nblocks), vec


def choose(kind, plan, pa, a, b, nblocks, panel_rows):
    """Returns (hap1 switches, hap2 switches, trace lines) or None."""
    ha = sorted(n for n in pa if n.startswith(a + "#"))
    hb = sorted(n for n in pa if n.startswith(b + "#"))
    rows = [[pa[ha[h]], pa[hb[h]]] for h in range(2)]   # rows[h] = (donor A row h, donor B row h); swap_hap_b 0
    order = rank(plan, nblocks)
    trace = [f"{plan}: boundary rank order {','.join(map(str, order))}"]
    if kind == "one":
        for i in order:
            ok1, v1 = valid_hap(rows[0], [i], nblocks, panel_rows)
            ok = ok1 and v1 != rows[1][0]
            trace.append(f"{plan}: hap1 boundary {i}: {'valid' if ok else 'invalid'}")
            if ok:
                return [i], [], trace
        return None
    if kind == "both":
        first = None
        for i in order:
            ok1, _ = valid_hap(rows[0], [i], nblocks, panel_rows)
            trace.append(f"{plan}: hap1 boundary {i}: {'valid' if ok1 else 'invalid'}")
            if ok1:
                first = i
                break
        if first is None:
            return None
        v1 = haplotype(rows[0], [first], nblocks)[0]
        for j in order:
            if j == first:
                continue
            ok2, v2 = valid_hap(rows[1], [j], nblocks, panel_rows)
            ok2 = ok2 and v2 != v1
            trace.append(f"{plan}: hap2 boundary {j}: {'valid' if ok2 else 'invalid'}")
            if ok2:
                return [first], [j], trace
        return None
    for combo in itertools.combinations(order, 3):
        s = sorted(combo)
        ok1, v1 = valid_hap(rows[0], s, nblocks, panel_rows)
        ok2, v2 = valid_hap(rows[1], s, nblocks, panel_rows)
        ok = ok1 and ok2 and v1 != v2
        trace.append(f"{plan}: boundaries {','.join(map(str, s))}: {'valid' if ok else 'invalid'}")
        if ok:
            return s, s, trace
    return None


def run(cmd, log):
    p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    with open(log, "w") as f:
        f.write("$ " + " ".join(map(str, cmd)) + "\n" + p.stdout + p.stderr)
    return p.returncode, (p.stdout + p.stderr).strip()


trace, validation = [], []


def write_records(out):
    with open(os.path.join(out, "selection_trace.txt"), "w") as f:
        f.write("\n".join(trace) + "\n")
    with open(os.path.join(out, "plan_validation.tsv"), "w") as f:
        f.write("plan\tstatus\tcheck\tdetail\n" + "".join("\t".join(map(str, r)) + "\n" for r in validation))


def main(probe, truth_gfa, original, elig_path, pa_path, out):

    # ---- selection
    with open(elig_path) as f:
        elig = list(csv.DictReader(f, delimiter="\t"))
    pa = panel_alleles(pa_path)
    nblocks = len(next(iter(pa.values())))
    panel_rows = [v for v in pa.values()]
    donors = []
    for r in elig:
        reasons = []
        if r["orientation"] != "forward/forward":
            reasons.append("not forward/forward")
        for col, rule in (("complete_two_rows", "rule 1"), ("respell_exact", "rule 2"),
                          ("assigned_sample_path", "rule 3"), ("bed_defined", "rule 4")):
            if r[col] != "True":
                reasons.append(f"{rule} ({col})")
        if r["donor"].lower().startswith(("grch38", "chm13")):
            reasons.append("rule 3 (reference assembly)")
        if r["CYP2D7_H1H2"] != "1/1":
            reasons.append("rule 4 (CYP2D7 copies per homologue != 1)")
        trace.append(f"eligibility {r['donor']}: {'eligible' if not reasons else 'excluded: ' + '; '.join(reasons)}")
        if not reasons:
            donors.append(r)
    strata = {}
    for r in donors:
        strata.setdefault(r["CYP2D6_H1H2"], []).append(r["donor"])
    for k in strata:
        strata[k].sort()
    if REQUIRED not in {r["donor"] for r in donors}:
        raise Stop(f"{REQUIRED} is not eligible")
    ordered = [REQUIRED]
    trace.append(f"rule 6: {REQUIRED} first")
    queues = {k: [d for d in v if d != REQUIRED] for k, v in strata.items()}
    stratum_order = sorted(queues, key=lambda k: (len(strata[k]), k))
    trace.append("rule 5/7: strata (rarest first) " + "; ".join(f"{k}: {len(strata[k])} ({','.join(strata[k])})"
                                                                  for k in stratum_order))
    while any(queues.values()):
        for k in stratum_order:
            if queues[k]:
                d = queues[k].pop(0)
                ordered.append(d)
                trace.append(f"rule 7/8: next from stratum {k}: {d}")
    if len(ordered) < 12:
        raise Stop(f"only {len(ordered)} eligible donors; 12 are needed")
    a_list, b_list, reserves = ordered[:6], ordered[6:12], ordered[12:]
    stratum_of = {r["donor"]: r["CYP2D6_H1H2"] for r in donors}
    with open(os.path.join(out, "eligible_donors.tsv"), "w") as f:
        f.write("order\tdonor\tH1\tH2\tCYP2D6_H1H2\tCYP2D7_H1H2\trole\n")
        byname = {r["donor"]: r for r in donors}
        for i, d in enumerate(ordered):
            role = f"A{i + 1}" if i < 6 else f"B{i - 5}" if i < 12 else "reserve"
            r = byname[d]
            f.write(f"{i + 1}\t{d}\t{r['H1']}\t{r['H2']}\t{r['CYP2D6_H1H2']}\t{r['CYP2D7_H1H2']}\t{role}\n")

    # ---- plans
    oseq, opaths = read_gfa(original)
    template_spellings = {spell(oseq, s): n for n, s in opaths.items()}
    rows_out = []

    def check(plan, ok, what, detail=""):
        validation.append((plan, "PASS" if ok else "FAIL", what, detail))
        if not ok:
            raise Stop(f"{plan}: {what} {detail}")

    used = set(a_list + b_list)
    for k, (plan, kind) in enumerate(PLANS):
        a, b = a_list[k], b_list[k]
        chosen = choose(kind, plan, pa, a, b, nblocks, panel_rows)
        while chosen is None:
            trace.append(f"{plan}: no valid {kind} plan for {a} x {b}")
            nxt = next((d for d in reserves if d not in used), None)
            if nxt is None:
                raise Stop(f"{plan}: no valid plan and no reserve donor left")
            used.add(nxt)
            trace.append(f"{plan}: donor B {b} replaced by reserve {nxt}")
            b = nxt
            chosen = choose(kind, plan, pa, a, b, nblocks, panel_rows)
        s1, s2, t = chosen
        trace.extend(t)
        require = "one" if kind == "one" else "both"
        d = os.path.join(out, "plans", plan)
        os.makedirs(d, exist_ok=True)
        cmd = [probe, "mosaic-truth", "--gfa", truth_gfa, "--reference", REF, "--donor-a", a, "--donor-b", b,
               "--switches-hap1", ",".join(map(str, s1)), "--switches-hap2", ",".join(map(str, s2)),
               "--require-off-panel", require]
        rc, msg = run(cmd + ["--out-gfa", os.path.join(d, "staged.gfa"), "--out", os.path.join(d, "truth.tsv")],
                      os.path.join(d, "mosaic_truth.log"))
        check(plan, rc == 0, "mosaic-truth accepted (two-route self-check)", msg.splitlines()[-1] if msg else "")
        rc2, _ = run(cmd + ["--out-gfa", os.path.join(d, "staged.again.gfa"), "--out", os.path.join(d, "truth.again.tsv")],
                     os.path.join(d, "mosaic_truth.again.log"))
        same = rc2 == 0 and all(open(os.path.join(d, x), "rb").read() == open(os.path.join(d, y), "rb").read()
                                for x, y in (("staged.gfa", "staged.again.gfa"), ("truth.tsv", "truth.again.tsv")))
        check(plan, same, "mosaic-truth rebuilt byte-identically")
        os.remove(os.path.join(d, "staged.again.gfa"))
        os.remove(os.path.join(d, "truth.again.tsv"))
        rc, msg = run([sys.executable, "-I", os.path.join(HERE, "translate_mosaic_frame.py"), original,
                       os.path.join(d, "staged.gfa"), os.path.join(d, "truth.gfa"), os.path.join(d, "translate.tsv")],
                      os.path.join(d, "translate.log"))
        check(plan, rc == 0, "MOSAIC paths translated to the original node frame, spelling identical", msg[-200:])
        # independent checks on the staged graph's panel alleles
        rc, msg = run([probe, "panel-alleles", "--gfa", os.path.join(d, "staged.gfa"), "--reference", REF,
                       "--out", os.path.join(d, "panel_alleles.tsv")], os.path.join(d, "panel_alleles.log"))
        check(plan, rc == 0, "panel-alleles over the staged graph")
        spa = panel_alleles(os.path.join(d, "panel_alleles.tsv"))
        check(plan, len(next(iter(spa.values()))) == nblocks, "staged model keeps the 19-block chain", str(nblocks))
        donor_rows = {n: v for n, v in spa.items() if not n.startswith("MOSAIC#")}
        check(plan, donor_rows == pa, "staged panel rows equal the calling panel's rows (126 templates, same alleles)")
        ha = sorted(n for n in pa if n.startswith(a + "#"))
        hb = sorted(n for n in pa if n.startswith(b + "#"))
        mos = [spa["MOSAIC#1#synthetic:0-0"], spa["MOSAIC#2#synthetic:0-0"]]
        srcs = []
        for h, sw in ((0, s1), (1, s2)):
            intended, src = haplotype([pa[ha[h]], pa[hb[h]]], sw, nblocks)
            srcs.append(src)
            check(plan, mos[h] == intended, f"MOSAIC#{h + 1}: every block carries its intended donor's allele")
            check(plan, all(0 < i < nblocks for i in sw), f"MOSAIC#{h + 1}: switches at internal boundaries", str(sw))
            if sw:
                check(plan, mos[h] not in panel_rows, f"MOSAIC#{h + 1} (switched) is in no panel row")
                check(plan, segments_effective([pa[ha[h]], pa[hb[h]]], sw, nblocks),
                      f"MOSAIC#{h + 1}: every donor segment is distinguishable (no no-op switch)")
            else:
                check(plan, mos[h] == pa[ha[h]], f"MOSAIC#{h + 1} (unswitched) is donor A's row")
        check(plan, mos[0] != mos[1], "the two mosaic homologues differ")
        diplotype_off = not any(sorted([mos[0], mos[1]]) == sorted([pa[x], pa[y]])
                                for x in pa for y in pa)
        check(plan, diplotype_off, "the mosaic diplotype is no pair of panel rows")
        # sequence level, on the original frame
        tseq, tpaths = read_gfa(os.path.join(d, "truth.gfa"))
        sseq, spaths = read_gfa(os.path.join(d, "staged.gfa"))
        for h in (1, 2):
            name = f"MOSAIC#{h}#synthetic:0-0"
            t = spell(tseq, tpaths[name])
            check(plan, t == spell(sseq, spaths[name]), f"{name}: original-frame spelling == staged spelling", f"{len(t)} bp")
            sw = s1 if h == 1 else s2
            if sw:
                check(plan, t not in template_spellings, f"{name}: sequence equals no template of the original graph")
        others = {n: p for n, p in tpaths.items() if not n.startswith("MOSAIC#")}
        check(plan, others == opaths and tseq == oseq, "truth graph = original graph + the two MOSAIC paths")
        rc, msg = run([probe, "respell", "--gfa", os.path.join(d, "staged.gfa"), "--reference", REF,
                       "--out", os.path.join(d, "respell.tsv")], os.path.join(d, "respell.log"))
        failing = next((l.split("\t")[1].strip() for l in open(os.path.join(d, "respell.tsv")) if l.startswith("#failing")), "?")
        check(plan, rc == 0 and failing == "0", "respell over the staged graph: 0 failing", failing)
        rows_out.append({
            "plan": plan, "kind": {"one": "one_homologue_one_switch", "both": "both_homologues_one_switch",
                                   "three": "three_switches_per_homologue"}[kind],
            "donor_a": a, "donor_b": b, "donor_a_stratum": stratum_of[a], "donor_b_stratum": stratum_of[b],
            "swap_hap_b": 0, "switches_hap1": ",".join(map(str, s1)) or "none",
            "switches_hap2": ",".join(map(str, s2)) or "none", "require_off_panel": require,
            "hap1_sources": "".join("AB"[x] for x in srcs[0]), "hap2_sources": "".join("AB"[x] for x in srcs[1]),
            "hap1_bp": len(spell(tseq, tpaths["MOSAIC#1#synthetic:0-0"])),
            "hap2_bp": len(spell(tseq, tpaths["MOSAIC#2#synthetic:0-0"])),
            "staged_gfa_sha256": hashlib.sha256(open(os.path.join(d, "staged.gfa"), "rb").read()).hexdigest(),
            "truth_gfa_sha256": hashlib.sha256(open(os.path.join(d, "truth.gfa"), "rb").read()).hexdigest(),
            "truth_tsv_sha256": hashlib.sha256(open(os.path.join(d, "truth.tsv"), "rb").read()).hexdigest()})
    donors_used = [x for r in rows_out for x in (r["donor_a"], r["donor_b"])]
    check("all", len(set(donors_used)) == 12, "twelve distinct donors across the six plans")
    with open(os.path.join(out, "mosaic_plans.tsv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows_out)
    print(f"six plans constructed; {len(validation)} checks pass")


if __name__ == "__main__":
    os.makedirs(sys.argv[6], exist_ok=True)
    try:
        main(*sys.argv[1:7])
    except Stop as e:
        trace.append(f"STOP: {e}")
        print(f"STOP: {e}", file=sys.stderr)
        sys.exit(3)
    finally:
        write_records(sys.argv[6])
