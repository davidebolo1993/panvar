#!/usr/bin/env python3
"""M9-C freeze and descriptive analysis (frozen with DESIGN.md).

  m9c.py freeze <out>    after every call, before any evaluation: per cell the outcome (resource failure / NO_PREDICTION /
                         PROVISIONAL / PREDICTION with its sequence_status), the v2 prediction (representative SHA-256s,
                         origins, FASTA decision), SHA-256 of the authoritative sequence outputs, and the effective-settings
                         check; reads no truth
  m9c.py analyse <out>   reliability of the v2 projection-aware block error masses at reportable blocks with IDENTIFIED
                         genotypes, against exact block-allele agreement from the archived evaluator's truth columns;
                         predeclared bins; floor-zero vs truth-unavailable strata; overall accounting; family-level outcomes;
                         legacy GQ comparison where assumptions permit; high-confidence error inspection; phase separately
"""
import csv
import hashlib
import math
import os
import sys

mode, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
OUT = "/Users/davide.bolognini/panvar/tests/results/m9c_acot_reliability_20261020"
RUNS = f"{OUT}/runs"
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l and not l.startswith("#")} if os.path.exists(p) else {}
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:16] if os.path.exists(p) else "absent"
prep = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in open(f"{RUNS}/prepare/prepare.tsv")}
DONORS = [prep[f"donor_{i}"][0] for i in range(1, 7)]
FAMILY = {prep[f"donor_{i}"][0]: prep[f"donor_{i}"][1].split()[1] for i in range(1, 7)}
CELLS = [(f"LOO.d{i}", DONORS[i - 1], "LOO") for i in range(1, 7)] + [(f"LZO.d{i}", DONORS[i - 1], "LZO") for i in (1, 2)]
LED = {r["unit"]: r for r in T(f"{OUT}/ledger.tsv")}
EXPECT = {"output_schema": "sequence_first_v2", "switch_probability": "0.01", "copying_transition": "constant_v1", "repair_joint_search": "off",
          "threads": "4", "time_limit_seconds": "900", "footprint_limit_bytes": "6442450944"}
BINS = [("0 (printed zero)", lambda u: u == 0.0), ("(0, 1e-6)", lambda u: 0 < u < 1e-6), ("[1e-6, 1e-3)", lambda u: 1e-6 <= u < 1e-3),
        ("[1e-3, 1e-2)", lambda u: 1e-3 <= u < 1e-2), ("[1e-2, 1e-1)", lambda u: 1e-2 <= u < 1e-1), ("[0.1, 1]", lambda u: u >= 0.1)]
HIGH = 1e-2   # predeclared: a high-confidence block has error_mass_upper < 0.01


def write(name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def outcome(cell):
    d = f"{RUNS}/call.{cell}"
    led = LED.get(f"call.{cell}", {})
    sl = KV(f"{d}/call.sequence_locus.tsv")
    if not led:
        return "NOT_RUN", sl
    if led["timed_out"] == "1":
        return "RESOURCE_FAILURE (watchdog)", sl
    if not sl:
        st = KV(f"{d}/call.mosaic_status.tsv")
        return f"NO_SEQUENCE_REPORT (exit {led['exit']}, status {st.get('status', 'absent')})", sl
    return sl.get("prediction", "absent"), sl


if mode == "freeze":
    rows = []
    for cell, dn, scheme in CELLS:
        d = f"{RUNS}/call.{cell}"
        oc, sl = outcome(cell)
        settings = {r["key"]: r["value"] for r in T(f"{d}/call.settings.tsv")}
        st = KV(f"{d}/call.mosaic_status.tsv")
        seen = {k: settings.get(k, st.get(k, sl.get(k, "absent"))) for k in EXPECT}
        same = lambda k, v: (float(seen[k]) == float(v)) if k == "switch_probability" else seen[k] == v
        bad = [k for k, v in EXPECT.items() if seen[k] != "absent" and not same(k, v)]
        unverified = [k for k in EXPECT if seen[k] == "absent"]
        haps = T(f"{d}/call.sequence_haplotypes.tsv")
        rows.append(dict(cell=cell, family=FAMILY[dn], outcome=oc, sequence_status=sl.get("sequence_status", "NA"), provisional=sl.get("provisional", "NA"),
                         legacy_run_status=sl.get("legacy_run_status", st.get("status", "NA")), legacy_no_call_reason=sl.get("legacy_no_call_reason", "NA"),
                         dosage_mode=sl.get("dosage_mode", "NA"), decision_set_groups=sl.get("decision_set_groups", "NA"),
                         representative_probability=sl.get("representative_probability", "NA"), decision_set_probability=sl.get("decision_set_probability", "NA"),
                         identified_blocks=sl.get("identified_blocks", "NA"), representation_ambiguous_blocks=sl.get("representation_ambiguous_blocks", "NA"),
                         genotype_unresolved_blocks=sl.get("genotype_unresolved_blocks", "NA"), provisional_label_blocks=sl.get("provisional_label_blocks", "NA"),
                         representation_intervals=sl.get("representation_intervals", "NA"),
                         hap1_sha256=sl.get("representative_hap1_sha256", "NA")[:16], hap2_sha256=sl.get("representative_hap2_sha256", "NA")[:16],
                         hap_origins=",".join(h["origin"] for h in haps) or "NA", sequence_fasta=st.get("sequence_fasta", "NA"),
                         sequence_fasta_reason=st.get("sequence_fasta_reason", "NA"), fasta_present=os.path.exists(f"{d}/call.sequence_haplotypes.fa"),
                         legacy_fasta_present=os.path.exists(f"{d}/call.mosaic_haplotypes.fa"), legacy_phase_support=sl.get("legacy_phase_support", "NA"),
                         sha_sequence_locus=sha(f"{d}/call.sequence_locus.tsv"), sha_sequence_blocks=sha(f"{d}/call.sequence_blocks.tsv"),
                         sha_sequence_haplotypes=sha(f"{d}/call.sequence_haplotypes.tsv"),
                         settings_seen=";".join(f"{k}={v}" for k, v in seen.items()), settings_mismatch=",".join(bad) or "none", settings_unverified=",".join(unverified) or "none",
                         wall_seconds=LED.get(f"call.{cell}", {}).get("wall_seconds", ""), peak_rss_gib=LED.get(f"call.{cell}", {}).get("peak_rss_gib", ""),
                         exit=LED.get(f"call.{cell}", {}).get("exit", "")))
    write("freeze", rows)
    print(open(f"{out}/freeze.tsv").read())
    sys.exit(0)

if mode == "analyse":
    fz = {r["cell"]: r for r in T(f"{RUNS}/freeze/freeze.tsv")}
    blocks, fam, hc, acct, gq = [], [], [], {}, []
    for cell, dn, scheme in CELLS:
        d = f"{RUNS}/call.{cell}"
        sb = T(f"{d}/call.sequence_blocks.tsv")
        ev = {r["block_id"]: r for r in T(f"{RUNS}/eval.{cell}/eval.blocks.tsv")}
        lg = {r["block_id"]: r for r in T(f"{d}/call.block_genotypes.tsv")}
        sl = KV(f"{d}/call.sequence_locus.tsv")
        f = fz.get(cell, {})
        kinds = {b: (ev[b]["kind"] if b in ev else lg[b]["kind"] if b in lg else "unknown") for b in [r["block_id"] for r in sb]}
        # candidate genotypes at each block, for the high-confidence error inspection (evaluated label pairs and their groups)
        cand, scores = None, None
        n_id = n_exact = n_err = n_hc_err = n_rep = 0
        for r in sb:
            b = r["block_id"]
            kind = kinds.get(b, "unknown")
            key = (scheme, kind, r["block_status"])
            acct[key] = acct.get(key, 0) + 1
            if kind != "SV":
                continue
            n_rep += 1
            if r["block_status"] != "IDENTIFIED":
                continue
            n_id += 1
            a1, a2 = sorted(int(x) for x in r["genotype"].split("/"))
            lo, up = float(r["error_mass_lower"]), float(r["error_mass_upper"])
            e = ev.get(b)
            if e is None or e["truth1_closest_allele"] == "NA":
                stratum, exact = ("not evaluated" if e is None else "truth undefined"), None
            else:
                avail = e["truth1_closest_distance"] == "0" and e["truth2_closest_distance"] == "0"
                t = sorted((int(e["truth1_closest_allele"]), int(e["truth2_closest_allele"])))
                stratum = "floor-zero (both truth alleles in panel)" if avail else "truth allele unavailable in calling panel"
                exact = avail and [a1, a2] == t
            row = dict(cell=cell, family=FAMILY[dn], scheme=scheme, block=b, genotype=f"{a1}/{a2}", error_mass_lower=r["error_mass_lower"],
                       error_mass_upper=r["error_mass_upper"], bin=next(n for n, fn in BINS if fn(up)), legacy_GQ=r["legacy_representative_GQ"],
                       representation_interval=r["representation_interval"], stratum=stratum, exact=exact,
                       truth_closest=f"{e['truth1_closest_allele']}/{e['truth2_closest_allele']}" if e else "NA",
                       truth_closest_distances=f"{e['truth1_closest_distance']}/{e['truth2_closest_distance']}" if e else "NA",
                       legacy_unordered_distance=e["unordered_distance"] if e else "NA", floor=e["availability_floor"] if e else "NA")
            blocks.append(row)
            if exact is True:
                n_exact += 1
            elif exact is False:
                n_err += 1
            # legacy GQ comparison where the assumptions coincide: no mixed groups at this block (lower == upper) and a
            # single decision-set group (the legacy representative is then the v2 representative's group)
            g = r["legacy_representative_GQ"]
            if g not in ("NA", "inf") and lo == up and sl.get("decision_set_groups") == "1":
                eg = 10 ** (-float(g) / 10)
                gq.append(dict(cell=cell, block=b, error_mass_upper=up, error_from_legacy_GQ=eg,
                               log10_ratio=(math.log10(up) - math.log10(eg)) if up > 0 and eg > 0 else "NA (zero)"))
            elif g not in ("NA",):
                gq.append(dict(cell=cell, block=b, error_mass_upper=up, error_from_legacy_GQ="not compared",
                               log10_ratio="assumptions differ (mixed groups or decision set > 1)"))
            if up < HIGH and exact is False:
                n_hc_err += 1
                if cand is None:
                    sc = {x["sequence_candidate_id"]: x for x in T(f"{d}/call.mosaic_scores.tsv")}
                    cand = [(x["sequence_candidate_id"], x["hap1_alleles"].split(","), x["hap2_alleles"].split(",")) for x in T(f"{d}/call.mosaic_candidates.tsv")]
                    scores = sc
                    idx = {x["block_id"]: i for i, x in enumerate(sb)}
                bi = idx[b]
                carriers = {gid for gid, h1, h2 in cand if sorted((int(h1[bi]), int(h2[bi]))) == t}
                rep = scores.get(sl.get("representative_group", ""), max(scores.values(), key=lambda x: float(x["score"])))
                lik = lambda x: float(x["fragment_term"]) + float(x["exposure_term"])
                if r["representation_interval"] not in ("", "NA") or lo != up:
                    cat = "representation ambiguity"
                elif not carriers:
                    cat = "correct block pair absent from evaluated candidates"
                else:
                    cat = "correct pair present but losing as whole hypotheses"
                best = max((scores[g2] for g2 in carriers if g2 in scores), key=lambda x: float(x["score"]), default=None)
                hc.append(dict(cell=cell, family=FAMILY[dn], block=b, called=f"{a1}/{a2}", truth=f"{t[0]}/{t[1]}", error_mass_upper=r["error_mass_upper"],
                               legacy_GQ=g, category=cat, evaluated_hypotheses_with_truth_pair=len(carriers),
                               best_truth_pair_gap_total=float(rep["score"]) - float(best["score"]) if best else "",
                               best_truth_pair_gap_likelihood=lik(rep) - lik(best) if best else "",
                               best_truth_pair_gap_prior=float(rep["sequence_pair_prior"]) - float(best["sequence_pair_prior"]) if best else "",
                               attribution="whole-hypothesis contrast only; not block-isolated (no block-local likelihood claim)"))
        fam.append(dict(cell=cell, family=FAMILY[dn], scheme=scheme, outcome=f.get("outcome", "NA"), sequence_status=f.get("sequence_status", "NA"),
                        legacy_run_status=f.get("legacy_run_status", "NA"), reportable_blocks=n_rep, identified_reportable=n_id, exact=n_exact, errors=n_err,
                        high_confidence_errors=n_hc_err, phase_support=f.get("legacy_phase_support", "NA"),
                        wall_seconds=f.get("wall_seconds", ""), peak_rss_gib=f.get("peak_rss_gib", "")))
    # reliability table: predeclared bins x stratum, LOO and LZO separately
    rel = []
    for scheme in ("LOO", "LZO"):
        for stratum in ("floor-zero (both truth alleles in panel)", "truth allele unavailable in calling panel", "truth undefined", "not evaluated"):
            for name, _ in BINS:
                sel = [x for x in blocks if x["scheme"] == scheme and x["stratum"] == stratum and x["bin"] == name]
                if sel:
                    rel.append(dict(scheme=scheme, stratum=stratum, bin=name, blocks=len(sel), exact=sum(x["exact"] is True for x in sel),
                                    not_exact=sum(x["exact"] is False for x in sel), families=len({x["family"] for x in sel}),
                                    mean_error_mass_upper=sum(float(x["error_mass_upper"]) for x in sel) / len(sel)))
    write("blocks", blocks)
    write("reliability", rel)
    write("families", fam)
    write("high_confidence_errors", hc)
    write("legacy_gq_comparison", gq)
    write("accounting", [dict(scheme=k[0], kind=k[1], block_status=k[2], blocks=v) for k, v in sorted(acct.items())])
    for n in ("families", "reliability", "high_confidence_errors", "accounting"):
        print(f"## {n}\n" + open(f"{out}/{n}.tsv").read())
    sys.exit(0)
