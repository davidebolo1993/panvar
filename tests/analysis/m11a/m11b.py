#!/usr/bin/env python3
"""M11-B: offline development audit of per-block truth-free confidence against the established floor-adjusted block
outcome (descriptive only; no fitting, thresholds, calibration, validation, or prediction changes).

  m11b.py manifest <out>   freeze the input manifest from saved evaluations (existing records only)
  m11b.py audit <out>      join each evaluated block to its run's per-block truth-free fields and describe the relation

Outcome contract (reused verbatim from scripts/genotype_report/load.py and aggregate.py):
  truth_defined          both truth1_bp and truth2_bp present
  excess                 unordered_distance - availability_floor (the evaluator's recorded field), bp
  normalized_excess      excess / max(truth1_bp + truth2_bp, called1_bp + called2_bp), truth-defined and denominator > 0
  floor_attained         excess == 0 (truth-defined)
  catastrophic           excess >= 1000 bp (the report's existing catastrophic_bp default; descriptive only)
Per-block truth-free quantity: sequence_first_v2 runs -- error_mass_upper (projection-aware, pool-conditional, from
<run>.sequence_blocks.tsv; error_mass_lower beside it); every run -- conditional_GQ / conditional_genotype_posterior and
genotype_support from <run>.block_genotypes.tsv. Locus-level support (decision_set_probability; legacy
winner_candidate_set_probability) is tabulated separately and never broadcast onto blocks.
Strata are configurations (never pooled); donors, not blocks, are the independent biological units.
"""
import collections
import csv
import hashlib
import math
import os
import re
import sys

R = "/Users/davide.bolognini/panvar/tests/results"
OUT = f"{R}/m11a_confidence_design_20261020/m11b"
CATASTROPHIC_BP = 1000.0
M9I = f"{R}/m9i1_compound_contract_20261020/runs/reapply/reapply.tsv"
INTERVALS = {"ACOT": ("sv:5", "ctx:6"), "CYP2D6": ("sv:3", "ctx:4", "sv:5")}   # model degeneracy intervals (M9-I inventory)
CYP_DEGENERATE_FP = {"2a5b720a1c5bf669", "82f7142952e845a3", "a12cabc3b8dd6766", "a2e2cc7795e47893", "e41f77fda080c143",
                     "8c7141b7b464601b", "0c24de16db70e642"}   # every CYP2D6 pilot model carries sv:3..sv:5 (41-42 equivalences)
ACOT_DEGENERATE_FP = {"12cbe6beb1ad28d1", "4606e1bd4b5e8d91", "545ff712f4e13b1c", "e560b1a604376290", "5bd687ecea4af414"}
SOURCES = [   # stratum, locus, package glob root, eval summary pattern, case-label function
    ("v2:M1-B", None, f"{R}/m1b_structural_20261019/runs", r"^eval\.(c4|cyp)_\d+$"),
    ("v2:M9-C", "ACOT", f"{R}/m9c_acot_reliability_20261020/runs", r"^eval\."),
    ("v2:M9-G", "ACOT", f"{R}/m9g_junction_integration_20261020/runs", r"^eval\."),
    ("v2:M9-H", "ACOT", f"{R}/m9h_jp11_validation_20261020/runs", r"^eval\."),
    ("legacy:C4-production-11fe681", "C4", f"{R}/genotype_c4_production_11fe681", None),
    ("legacy:CYP2D6-pilot-a5", "CYP2D6", f"{R}/cyp2d6_shortread_pilot_v1/data/eval_a5", r"\.reported\.summary\.tsv$"),
]


def T(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []


def KV(p):
    return dict(l.rstrip("\n").split("\t", 1) for l in open(p) if "\t" in l and not l.startswith("#")) if os.path.exists(p) else {}


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def write(out, name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def summaries(root, pattern):
    out = []
    if pattern is None:   # the C4 production campaign: every evaluator summary under the package
        for dp, _, fs in os.walk(root):
            for f in fs:
                if f.endswith(".summary.tsv") and "eval" in f and "/sentinel" not in dp:
                    out.append(os.path.join(dp, f))
        return sorted(out)
    if pattern.startswith("^eval"):
        for d in sorted(os.listdir(root)):
            if re.match(pattern, d) and os.path.exists(f"{root}/{d}/eval.summary.tsv"):
                out.append(f"{root}/{d}/eval.summary.tsv")
        return out
    return sorted(f"{root}/{f}" for f in os.listdir(root) if re.search(pattern, f)) if os.path.isdir(root) else []


def case_label(stratum, s):
    run = s["run"]
    m = re.search(r"/call\.(on|off)\.(L\d|Z\d)/call$", run)
    if m and stratum == "v2:M9-H":
        return f"M9-H {m.group(1)} {m.group(2)}"
    m = re.search(r"/call\.(on|off)\.(LOO|LZO)\.(d\d)/call$", run)
    if m:
        return f"M9-G {m.group(1)} {m.group(2)}.{m.group(3)}"
    m = re.search(r"/call\.(LOO|LZO)\.(d\d)/call$", run)
    if m:
        return f"M9-C {m.group(1)}.{m.group(2)}"
    m = re.search(r"/call\.(cyp_\d+)/call$", run)
    if m:
        return {"cyp_5": "CYP M1-B cyp_5 (NA19650 LOO 5x)", "cyp_30": "CYP M1-B cyp_30 (NA19650 LOO 30x)"}[m.group(1)]
    return "-"


def manifest(out):
    rows = []
    for stratum, locus, root, pattern in SOURCES:
        for sp in summaries(root, pattern):
            s = KV(sp)
            run = s.get("run", "")
            bp = sp.replace("summary.tsv", "blocks.tsv")
            if s.get("run_no_call") == "true":   # Amendment 1: a no-call run stays visible (incomplete stratum), never "missing"
                rows.append(dict(stratum=stratum, locus=locus or "-", summary=sp, run=run, status="NO_CALL",
                                 reason=f"run emitted no block genotype ({s.get('status', '-')})", donor=s.get("donor", "-"),
                                 scheme=s.get("scheme", "-"), model_fingerprint=s.get("model_fingerprint", "-"), run_status=s.get("status", "-")))
                continue
            if not run or not os.path.exists(bp) or not os.path.exists(run + ".block_genotypes.tsv"):
                why = "no run prefix in the summary" if not run else "no evaluator block table" if not os.path.exists(bp) else "no run block table"
                rows.append(dict(stratum=stratum, summary=sp, status="EXCLUDED", reason=why))
                continue
            loc = locus or ("C4" if "/call.c4_" in run else "CYP2D6")
            st = KV(run + ".mosaic_status.tsv")
            rows.append(dict(stratum=stratum, locus=loc, summary=sp, blocks=bp, run=run, status="INCLUDED", reason="-",
                             donor=s.get("donor", "-"), scheme=s.get("scheme", "-"), model_fingerprint=s.get("model_fingerprint", "-"),
                             profile_sha256=st.get("profile_sha256", "-"), reads1_sha256=st.get("sha256_reads1", "-"), run_status=s.get("status", "-"),
                             v2=os.path.exists(run + ".sequence_blocks.tsv"), case_label=case_label(stratum, s),
                             eval_blocks_sha256=sha(bp)))
    # duplicates: the same evaluation of the same reads, model and profile (M9-G off = M9-C). Amendment 1: the reads' SHA-256
    # is part of the key, so coverage replicates (different reads) with identical evaluator tables stay distinct runs
    seen = {}
    for r in rows:
        if r["status"] != "INCLUDED":
            continue
        key = (r["model_fingerprint"], r["profile_sha256"], r["reads1_sha256"], r["eval_blocks_sha256"])
        if key in seen:
            r["status"], r["reason"] = "DUPLICATE", f"identical evaluator block table, model and profile as {seen[key]}"
        else:
            seen[key] = r["run"]
    write(out, "manifest", rows)
    c = collections.Counter((r["stratum"], r["status"]) for r in rows)
    for k in sorted(c):
        print(f"{k[0]}\t{k[1]}\t{c[k]}")
    return 0


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def mass_bin(m):
    if m is None:
        return "NA"
    if m == 0.0:
        return "printed 0 (rounding/underflow; not certainty)"
    e = math.floor(math.log10(m))
    return "1e-1..1" if e >= -1 else f"1e{e}..1e{e + 1}" if e >= -9 else "<1e-9"


def gq_bin(g):
    if g is None:
        return "NA"
    if math.isinf(g):
        return "inf (posterior printed as 1)"
    return f"{int(g // 10) * 10}-{int(g // 10) * 10 + 10}" if g < 100 else ">=100"


def audit(out):
    allm = T(f"{OUT}/runs/manifest.a1/manifest.tsv")   # Amendment 1
    man = [r for r in allm if r["status"] == "INCLUDED"]
    comp = {r["case"]: num(r["compound_distance_m9i"]) for r in T(M9I)}
    obs, runs = [], []
    for m in man:
        ev = {r["block_id"]: r for r in T(m["blocks"])}
        bg = {r["block_id"]: r for r in T(m["run"] + ".block_genotypes.tsv")}
        sb = {r["block_id"]: r for r in T(m["run"] + ".sequence_blocks.tsv")} if m["v2"] == "True" else {}
        sl = KV(m["run"] + ".sequence_locus.tsv")
        s = KV(m["summary"])
        loc = m["locus"]
        degen = INTERVALS.get(loc) if ((loc == "ACOT" and m["model_fingerprint"] in ACOT_DEGENERATE_FP) or
                                       (loc == "CYP2D6" and m["model_fingerprint"] in CYP_DEGENERATE_FP)) else None
        dsum = sum(num(ev[b]["unordered_distance"]) or 0 for b in (degen or ()) if b in ev and num(ev[b]["unordered_distance"]) is not None)
        cd = comp.get(m["case_label"])
        artefact_bp = max(0.0, dsum - cd) if (degen and cd is not None) else None
        for bid, e in ev.items():
            t1, t2, c1, c2 = (num(e.get(k)) for k in ("truth1_bp", "truth2_bp", "called1_bp", "called2_bp"))
            td = t1 is not None and t2 is not None
            ex = num(e["excess"]) if td else None
            den = max((t1 or 0) + (t2 or 0), (c1 or 0) + (c2 or 0)) if td else None
            g = bg.get(bid, {})
            gq_s = g.get("conditional_GQ", "NA")
            gq = float("inf") if gq_s in ("inf", "Infinity") else num(gq_s)
            mu = num(sb.get(bid, {}).get("error_mass_upper")) if sb else None
            in_deg = bool(degen and bid in degen)
            obs.append(dict(stratum=m["stratum"], locus=loc, donor=m["donor"], scheme=m["scheme"], run=m["run"], run_status=m["run_status"],
                            model_fingerprint=m["model_fingerprint"], block_id=bid, kind=e.get("kind", "-"), truth_defined=td,
                            excess=ex, normalized_excess=(ex / den) if (td and den) else None, floor_attained=(ex == 0) if td else None,
                            catastrophic=(ex >= CATASTROPHIC_BP) if td else None,
                            error_mass_upper=mu, error_mass_lower=num(sb.get(bid, {}).get("error_mass_lower")) if sb else None,
                            error_mass_printed=sb.get(bid, {}).get("error_mass_upper", "NA") if sb else "NA",
                            block_status=sb.get(bid, {}).get("block_status", "NA") if sb else "NA",
                            conditional_GQ=gq, gq_printed=gq_s, posterior_printed=g.get("conditional_genotype_posterior", "NA"),
                            genotype_support=g.get("genotype_support", "NA"), quality_filter=g.get("quality_filter", "NA"),
                            in_degeneracy_interval=in_deg,
                            synonymous_boundary_artefact=("YES" if in_deg and artefact_bp and (num(e["unordered_distance"]) or 0) > 0 else
                                                          "UNQUANTIFIED" if in_deg and artefact_bp is None and (ex or 0) > 0 else "no"),
                            interval_artefact_bp=artefact_bp if in_deg else None))
        runs.append(dict(stratum=m["stratum"], locus=loc, donor=m["donor"], scheme=m["scheme"], run_status=m["run_status"],
                         completion=sl.get("completion", "legacy"), decision_set_probability=sl.get("decision_set_probability", "NA"),
                         winner_candidate_set_probability=s.get("winner_candidate_set_probability", "NA"),
                         excess_over_floor_total=s.get("excess_over_floor_total", "NA"), undefined_truth_blocks=s.get("undefined_truth_blocks", "NA"),
                         case_label=m["case_label"], interval_artefact_bp=artefact_bp))
    for m in allm:   # Amendment 1: no-call runs stay in the run-level denominators
        if m["status"] == "NO_CALL":
            runs.append(dict(stratum=m["stratum"], locus=m["locus"], donor=m["donor"], scheme=m["scheme"], run_status=m["run_status"],
                             completion="NO_CALL", decision_set_probability="NA", winner_candidate_set_probability="NA",
                             excess_over_floor_total="NA", undefined_truth_blocks="NA", case_label="-", interval_artefact_bp=None))
    write(out, "block_observations", obs)
    write(out, "locus_level", runs)
    summary, extremes = [], []
    for stratum in dict.fromkeys(o["stratum"] for o in obs):
        S = [o for o in obs if o["stratum"] == stratum and o["truth_defined"]]
        for qname, binf, key in (("error_mass_upper", mass_bin, "error_mass_upper"), ("conditional_GQ", gq_bin, "conditional_GQ")):
            Q = [o for o in S if o[key] is not None]
            if not Q:
                continue
            bins = collections.defaultdict(list)
            for o in Q:
                bins[binf(o[key])].append(o)
            for b, L in sorted(bins.items()):
                xs = sorted(o["excess"] for o in L)
                summary.append(dict(stratum=stratum, quantity=qname, bin=b, blocks=len(L), donors=len({o["donor"] for o in L}),
                                    floor_attained=sum(o["floor_attained"] for o in L), floor_attained_rate=f"{sum(o['floor_attained'] for o in L) / len(L):.3f}",
                                    excess_median=xs[len(xs) // 2], excess_max=xs[-1], catastrophic=sum(o["catastrophic"] for o in L),
                                    normalized_excess_median=sorted(o["normalized_excess"] or 0 for o in L)[len(L) // 2],
                                    synonymous_artefact_blocks=sum(o["synonymous_boundary_artefact"] == "YES" for o in L)))
            # descriptive rank separation: P(error-mass of a non-floor block > that of a floor block) (ties count 1/2);
            # for GQ the direction is reversed. Blocks are not independent; no test is made.
            pos = [o[key] for o in Q if not o["floor_attained"]]
            neg = [o[key] for o in Q if o["floor_attained"]]
            if pos and neg:
                sgn = 1 if qname == "error_mass_upper" else -1
                wins = sum((1.0 if sgn * (p - n) > 0 else 0.5 if p == n else 0.0) for p in pos for n in neg) / (len(pos) * len(neg))
            else:
                wins = None
            summary.append(dict(stratum=stratum, quantity=qname, bin="ALL (rank separation, descriptive)", blocks=len(Q),
                                donors=len({o["donor"] for o in Q}), floor_attained=len(neg),
                                floor_attained_rate=f"{len(neg) / len(Q):.3f}", excess_median="-", excess_max=max(o["excess"] for o in Q),
                                catastrophic=sum(o["catastrophic"] for o in Q), normalized_excess_median="-",
                                synonymous_artefact_blocks=sum(o["synonymous_boundary_artefact"] == "YES" for o in Q),
                                rank_separation_nonfloor_vs_floor="NA" if wins is None else f"{wins:.3f}"))
            hi_err = sorted((o for o in Q if not o["floor_attained"]), key=lambda o: (o[key] if qname == "error_mass_upper" else -o[key], -o["excess"]))[:12]
            lo_ok = sorted((o for o in Q if o["floor_attained"]), key=lambda o: (-o[key] if qname == "error_mass_upper" else o[key]))[:8]
            for kind, L in (("HIGH_SUPPORT_ERROR", hi_err), ("LOW_SUPPORT_AT_FLOOR", lo_ok)):
                for o in L:
                    extremes.append(dict(stratum=stratum, quantity=qname, kind=kind, donor=o["donor"], scheme=o["scheme"], block_id=o["block_id"],
                                         value=o["error_mass_printed"] if qname == "error_mass_upper" else o["gq_printed"], excess=o["excess"],
                                         normalized_excess=o["normalized_excess"], genotype_support=o["genotype_support"],
                                         synonymous_boundary_artefact=o["synonymous_boundary_artefact"], run=o["run"][-70:]))
    write(out, "binned_summary", summary)
    write(out, "extremes", extremes)
    rounding = []
    for stratum in dict.fromkeys(o["stratum"] for o in obs):
        S = [o for o in obs if o["stratum"] == stratum]
        rounding.append(dict(stratum=stratum, blocks=len(S), truth_undefined=sum(not o["truth_defined"] for o in S),
                             error_mass_printed_zero=sum(o["error_mass_printed"] in ("0", "0.0") for o in S),
                             posterior_printed_one=sum(o["posterior_printed"] in ("1", "1.0") for o in S),
                             gq_inf_or_missing=sum(o["conditional_GQ"] is None or (o["conditional_GQ"] is not None and math.isinf(o["conditional_GQ"])) for o in S),
                             in_degeneracy_interval=sum(o["in_degeneracy_interval"] for o in S),
                             synonymous_artefact_yes=sum(o["synonymous_boundary_artefact"] == "YES" for o in S),
                             synonymous_artefact_unquantified=sum(o["synonymous_boundary_artefact"] == "UNQUANTIFIED" for o in S),
                             donors=len({o["donor"] for o in S}), runs=len({o["run"] for o in S})))
    write(out, "rounding_and_artefacts", rounding)
    for n in ("rounding_and_artefacts", "binned_summary"):
        print(f"## {n}\n" + open(f"{out}/{n}.tsv").read())
    return 0


if __name__ == "__main__":
    mode, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    sys.exit({"manifest": manifest, "audit": audit}[mode](out))
