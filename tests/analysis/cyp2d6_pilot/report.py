#!/usr/bin/env python3
"""CYP2D6/D7 short-read pilot: tables, Figures A and B and caption numbers from the sealed, evaluated cells (read-only).

Usage: report.py <pilot package> <scripts dir> <out dir> <tag>
Reads jobs.tsv, call_classification_<tag>.tsv (classify_ties.py), degeneracy/ (model dumps and inventory),
data/runs_<tag>/<job>.* (call outputs), data/eval_<tag>/<job>.{raw,reported}.{blocks,summary}.tsv and the simulated truth
haplotypes data/reads/<sample>/cov30/truth.fa. Writes <out dir>/figures/ and <out dir>/data/*.tsv.

Kept separate, as AMENDMENT_5 requires:
  callability   every cell by class: strict, genotype-equivalent, SEQUENCE_IDENTICAL_BLOCK_LABEL_AMBIGUITY,
                DISTINCT_SEQUENCE_TIE, other no-call, resource limit. Every no-call stays in its denominator.
  block accuracy only where a block genotype was reported (calls); raw = the pre-dosage-repair outputs (.original),
                reported = the final outputs, retained = reported blocks whose caller genotype_filter is PASS (read
                from the caller's block table before any truth column is joined). No block value is made for a no-call.
  whole sequence strict calls: the called unordered haplotype pair against the truth pair (orientation-agnostic);
                label-ambiguity no-calls: the identified sequence (all tied members spell it) against truth, reported
                as the sequence only, no block labelling chosen; genotype-equivalent calls and distinct-sequence ties:
                only whether the truth pair is among the tied members, never counted as a correct sequence call.
  sv:3..sv:5    cells and samples with label ambiguity; truth haplotypes whose labelling has a same-sequence
                alternative; calls whose block excess at sv:3/sv:5 is purely a labelling difference.
"""
import csv
import os
import sys

import pandas as pd

PKG, SCRIPTS, OUT, TAG = sys.argv[1:5]
RUNS, EVAL = f"runs_{TAG}", f"eval_{TAG}"
sys.path.insert(0, os.path.abspath(SCRIPTS))
from genotype_report import plots as P  # noqa: E402


def compact_range(lo, hi):
    """Truth-length range for a block tick, compact enough for 19 blocks on one slide: 13–308bp, 64.3–64.4kb."""
    unit = "bp" if hi < 1000 else "kb"
    if unit == "kb" and lo < 1000:
        return f"{lo}bp–{hi / 1000:.1f}kb"
    if unit == "bp":
        return f"{lo}bp" if lo == hi else f"{lo}–{hi}bp"
    # whole kb from 10 kb up when the range stays distinct; one decimal otherwise
    places = 0 if lo >= 10000 and round(lo / 1000) != round(hi / 1000) else 1
    a, b = f"{lo / 1000:.{places}f}", f"{hi / 1000:.{places}f}"
    return f"{a}kb" if a == b else f"{a}–{b}kb"

ROWS = [("lzo", "LZO"), ("mosaic", "Mosaic LZO"), ("loo", "LOO")]
COVS = [5, 30]
BIG = 1000
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
LABEL = "SEQUENCE_IDENTICAL_BLOCK_LABEL_AMBIGUITY"
arm_order = {k: i for i, (k, _) in enumerate(ROWS)}


def tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def kv(path):
    out = {}
    for line in open(path).read().splitlines():
        p = line.split("\t")
        if len(p) >= 2:
            out.setdefault(p[0], p[1])
    return out


def canon(s):
    r = s.translate(COMP)[::-1]
    return min(s, r)


def fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line)
    return {k: "".join(v) for k, v in out.items()}


def sort_arm(df, extra=()):
    return df.sort_values(["arm", "coverage", *extra], key=lambda s: s.map(arm_order) if s.name == "arm" else s,
                          kind="stable")


model_cache = {}


def model(name):
    if name not in model_cache:
        blocks, ids = [], []
        with open(os.path.join(PKG, "degeneracy", "models", f"{name}.block_alleles.tsv")) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                b = int(r["block_index"])
                while len(blocks) <= b:
                    blocks.append([])
                    ids.append(r["block_id"])
                blocks[b].append(r["sequence"])
        model_cache[name] = (blocks, ids)
    return model_cache[name]


def spell(model_name, alleles):
    blocks, _ = model(model_name)
    return "".join(blocks[b][int(a)] for b, a in enumerate(alleles.split(",")))


jobs = tsv(os.path.join(PKG, "jobs.tsv"))
jobs["coverage"] = jobs.coverage.astype(int)
jobs["out_prefix"] = jobs.out_prefix.str.replace("data/runs/", f"data/{RUNS}/", n=1, regex=False)
cls = tsv(os.path.join(PKG, f"call_classification_{TAG}.tsv")).set_index("job")
inventory = tsv(os.path.join(PKG, "degeneracy", "block_representation_degeneracies.tsv"))
D = os.path.join(OUT, "data")
os.makedirs(D, exist_ok=True)
os.makedirs(os.path.join(OUT, "figures"), exist_ok=True)

# ---- cells, block observations (calls only), whole-sequence outcomes -------------------------------------
cells, obs, wseq, undefined = [], [], [], []
for j in jobs.itertuples():
    c = cls.loc[j.job]
    run = tsv(os.path.join(PKG, j.out_prefix + ".run.tsv")).iloc[0]
    st = kv(os.path.join(PKG, j.out_prefix + ".mosaic_status.tsv"))
    klass = c["classification"]
    cells.append({"job": j.job, "arm": j.arm, "sample": j.sample, "coverage": j.coverage, "exit": run["exit"],
                  "classification": klass, "caller_status": st.get("status"), "no_call_reason": st.get("no_call_reason"),
                  "wall_seconds": float(run.wall_seconds), "peak_rss_GiB": round(int(run.peak_rss_bytes) / 2 ** 30, 3),
                  "statistical_state_classes": st.get("statistical_state_classes"),
                  "physical_origins": st.get("physical_origins"), "repair_rounds_run": st.get("repair_rounds_run"),
                  "tie_branches_run": st.get("tie_branches_run")})
    mname = c["model"]
    truth = fasta(os.path.join(PKG, "data", "reads", j.sample, "cov30", "truth.fa"))
    tpair = sorted(canon(s) for s in truth.values())
    eq = list(csv.DictReader(open(os.path.join(PKG, j.out_prefix + ".equivalence_set.tsv")), delimiter="\t")) \
        if os.path.exists(os.path.join(PKG, j.out_prefix + ".equivalence_set.tsv")) else []
    members = [sorted(canon(spell(mname, r[h])) for h in ("hap1_alleles", "hap2_alleles")) for r in eq]
    row = {"job": j.job, "arm": j.arm, "sample": j.sample, "coverage": j.coverage, "classification": klass,
           "whole_sequence": "NA", "truth_pair_among_tied_members": "NA", "fasta_matches_model_spelling": "NA"}
    if klass == "CALL_strict_call":
        fa = fasta(os.path.join(PKG, j.out_prefix + ".mosaic_haplotypes.fa"))
        called = sorted(canon(s) for s in fa.values())
        row["fasta_matches_model_spelling"] = str(len(members) >= 1 and called == members[0]).lower()
        row["whole_sequence"] = "exact" if called == tpair else "different"
    elif klass == LABEL:
        row["whole_sequence"] = "exact (identified sequence)" if members[0] == tpair else "different (identified sequence)"
    elif klass in ("CALL_genotype_equivalent_call", "DISTINCT_SEQUENCE_TIE"):
        row["truth_pair_among_tied_members"] = str(tpair in members).lower()
        row["whole_sequence"] = "not a sequence call"
    wseq.append(row)
    if run["exit"] != "0":
        continue
    filt = tsv(os.path.join(PKG, j.out_prefix + ".block_genotypes.tsv"))
    member = dict(zip(filt.block_id, filt.genotype_filter))       # decided here, before any truth column
    for layer in ("raw", "reported"):
        b = tsv(os.path.join(PKG, "data", EVAL, f"{j.job}.{layer}.blocks.tsv"))
        for i, r in enumerate(b.itertuples()):
            if r.truth1_bp == "NA" or r.truth2_bp == "NA":
                # truth cannot be projected onto this block of the call's model: the block genotype is reported but
                # unscorable; tabulated, never given a value
                ctab = filt if layer == "reported" else tsv(os.path.join(PKG, j.out_prefix + ".original.block_genotypes.tsv"))
                crow = ctab[ctab.block_id == r.block_id].iloc[0]
                undefined.append({"job": j.job, "layer": layer, "arm": j.arm, "sample": j.sample, "coverage": j.coverage,
                                  "block_id": r.block_id, "called1_bp": crow.allele1_bp, "called2_bp": crow.allele2_bp,
                                  "genotype_filter": member[r.block_id]})
                continue
            if r.called1_bp == "NA":
                sys.exit(f"{j.job} {layer}: a call without a block genotype at {r.block_id}")
            t = int(r.truth1_bp) + int(r.truth2_bp)
            cs = int(r.called1_bp) + int(r.called2_bp)
            e = float(r.excess)
            obs.append({"run_id": j.job, "layer": layer, "row": j.arm, "arm": j.arm, "sample": j.sample,
                        "coverage": j.coverage, "block_index": i, "block_id": r.block_id, "kind": r.kind,
                        "truth1_bp": int(r.truth1_bp), "truth2_bp": int(r.truth2_bp),
                        "called1_bp": int(r.called1_bp), "called2_bp": int(r.called2_bp),
                        "called_alleles": f"{r.called_allele1},{r.called_allele2}",
                        "truth_alleles": f"{r.truth1_closest_allele},{r.truth2_closest_allele}",
                        "distance": float(r.unordered_distance), "floor": float(r.availability_floor),
                        "excess": e, "excess_norm": e / max(t, cs) if max(t, cs) > 0 else 0.0,
                        "genotype_filter": member[r.block_id] if layer == "reported" else "NA"})
cells, obs, wseq = pd.DataFrame(cells), pd.DataFrame(obs), pd.DataFrame(wseq)
pd.DataFrame(undefined, columns=["job", "layer", "arm", "sample", "coverage", "block_id", "called1_bp", "called2_bp",
                                 "genotype_filter"]).to_csv(os.path.join(D, "truth_undefined_blocks.tsv"), sep="\t", index=False)
cells.to_csv(os.path.join(D, "cells.tsv"), sep="\t", index=False)
obs.to_csv(os.path.join(D, "block_observations.tsv"), sep="\t", index=False)
wseq.to_csv(os.path.join(D, "whole_sequence.tsv"), sep="\t", index=False)

# ---- callability ---------------------------------------------------------------------------------------
classes = [("strict", "CALL_strict_call"), ("genotype_equivalent", "CALL_genotype_equivalent_call"),
           ("label_ambiguity_no_call", LABEL), ("distinct_sequence_tie_no_call", "DISTINCT_SEQUENCE_TIE")]
call_rows = []
for (arm, cov), x in cells.groupby(["arm", "coverage"], sort=False):
    r = {"arm": arm, "coverage": cov, "cells": len(x)}
    for name, k in classes:
        r[name] = int((x.classification == k).sum())
    r["other_no_call"] = int(x.classification.str.startswith("NO_CALL_").sum())
    r["not_completed"] = int(x.classification.str.startswith("NOT_COMPLETED").sum())
    r["block_genotype_reported"] = r["strict"] + r["genotype_equivalent"]
    r["callability"] = f"{r['block_genotype_reported']}/{r['cells']}"
    call_rows.append(r)
callab = sort_arm(pd.DataFrame(call_rows))
callab.to_csv(os.path.join(D, "callability.tsv"), sep="\t", index=False)

# ---- block accuracy (calls only) -----------------------------------------------------------------------
raw = obs[obs.layer == "raw"].copy()
rep = obs[obs.layer == "reported"].copy()
kept = rep[rep.genotype_filter == "PASS"].copy()


def block_summary(x, keys):
    g = x.groupby(keys, sort=False)
    return g.agg(block_calls=("excess", "size"), at_floor=("excess", lambda s: int((s == 0).sum())),
                 positive_excess=("excess", lambda s: int((s > 0).sum())),
                 above_1kb=("excess", lambda s: int((s > BIG).sum())),
                 excess_total_bp=("excess", "sum"), excess_max_bp=("excess", "max"),
                 excess_norm_median=("excess_norm", "median")).reset_index()


for keys, name in ((["arm", "coverage"], "by_arm_coverage"), (["arm", "coverage", "block_index", "block_id"], "by_block")):
    parts = []
    for layer, x in (("raw", raw), ("reported", rep), ("retained", kept)):
        t = block_summary(x, keys)
        t.insert(0, "layer", layer)
        parts.append(t)
    t = pd.concat(parts, ignore_index=True)
    t = t.sort_values(["layer", "arm", "coverage"] + (["block_index"] if "block_index" in keys else []),
                      key=lambda s: s.map(arm_order) if s.name == "arm" else s, kind="stable")
    t.to_csv(os.path.join(D, f"block_accuracy_{name}.tsv"), sep="\t", index=False)

f = []
for (arm, cov), x in rep.groupby(["arm", "coverage"], sort=False):
    ret, ok = x.genotype_filter == "PASS", x.excess == 0
    f.append({"arm": arm, "coverage": cov, "block_calls": len(x), "retained": int(ret.sum()), "excluded": int((~ret).sum()),
              "retained_at_floor": int((ret & ok).sum()), "retained_positive_excess": int((ret & ~ok).sum()),
              "excluded_at_floor": int((~ret & ok).sum()), "excluded_positive_excess": int((~ret & ~ok).sum()),
              "filter_precision": f"{int((ret & ok).sum())}/{int(ret.sum())}",
              "filter_recall": f"{int((ret & ok).sum())}/{int(ok.sum())}",
              "errors_removed": f"{int((~ret & ~ok).sum())}/{int((~ok).sum())}",
              "above_1kb_removed": f"{int((~ret & (x.excess > BIG)).sum())}/{int((x.excess > BIG).sum())}"})
sort_arm(pd.DataFrame(f)).to_csv(os.path.join(D, "filter_evaluation.tsv"), sep="\t", index=False)

# locus-wide outcome among calls, separate from block accuracy
loc = []
for (arm, cov), x in cells.groupby(["arm", "coverage"], sort=False):
    calls = x[x.exit == "0"].job
    r = {"arm": arm, "coverage": cov, "calls": len(calls)}
    for layer, o in (("raw", raw), ("reported", rep)):
        g = o[o.run_id.isin(calls)].groupby("run_id")
        r[f"{layer}_all_blocks_at_floor"] = int(g.excess.apply(lambda s: (s == 0).all()).sum()) if len(calls) else 0
        r[f"{layer}_all_blocks_exact"] = int(g.distance.apply(lambda s: (s == 0).all()).sum()) if len(calls) else 0
    g = rep[rep.run_id.isin(calls)].groupby("run_id")
    r["reported_all_blocks_retained"] = int(g.genotype_filter.apply(lambda s: (s == "PASS").all()).sum()) if len(calls) else 0
    loc.append(r)
sort_arm(pd.DataFrame(loc)).to_csv(os.path.join(D, "locus_outcomes_among_calls.tsv"), sep="\t", index=False)

# ---- whole-sequence summary ----------------------------------------------------------------------------
ws = []
for (arm, cov), x in wseq.groupby(["arm", "coverage"], sort=False):
    s = x[x.classification == "CALL_strict_call"]
    a = x[x.classification == LABEL]
    o = x[x.classification.isin(["CALL_genotype_equivalent_call", "DISTINCT_SEQUENCE_TIE"])]
    ws.append({"arm": arm, "coverage": cov,
               "strict_calls_whole_sequence_exact": f"{int((s.whole_sequence == 'exact').sum())}/{len(s)}",
               "label_ambiguity_identified_sequence_exact": f"{int(a.whole_sequence.str.startswith('exact').sum())}/{len(a)}",
               "equivalent_or_tie_truth_among_tied_members": f"{int((o.truth_pair_among_tied_members == 'true').sum())}/{len(o)}"})
sort_arm(pd.DataFrame(ws)).to_csv(os.path.join(D, "whole_sequence_summary.tsv"), sep="\t", index=False)

# ---- sv:3..sv:5 ----------------------------------------------------------------------------------------
deg = {}
for r in inventory.itertuples():
    lo, hi = int(r.first_block_index), int(r.last_block_index)
    deg.setdefault((r.model, lo, hi), set()).update({r.tuple_a, r.tuple_b})
iv = {}
for (m, lo, hi) in deg:
    iv.setdefault(m, set()).add((lo, hi))
run_exit = dict(zip(cells.job, cells.exit))
prefix_of = dict(zip(jobs.job, jobs.out_prefix))
sv = []
for j in jobs.itertuples():
    c = cls.loc[j.job]
    mname = c["model"]
    b = tsv(os.path.join(PKG, "data", EVAL, f"{j.job}.reported.blocks.tsv"))
    exposed = 0
    for lo, hi in iv.get(mname, ()):
        for h in (1, 2):
            tup = ",".join(b[f"truth{h}_closest_allele"].iloc[lo:hi + 1])
            dist = b[f"truth{h}_closest_distance"].iloc[lo:hi + 1]
            if all(d == "0" for d in dist) and tup in deg[(mname, lo, hi)]:
                exposed += 1
    label_only, segment, seg_excess = "NA", "NA", "NA"
    if c["classification"] == "CALL_strict_call":
        w = wseq[wseq.job == j.job].iloc[0]
        x = rep[(rep.run_id == j.job) & rep.block_id.isin(["sv:3", "ctx:4", "sv:5"])]
        if w.whole_sequence == "exact":
            label_only = str(bool((x.excess > 0).any())).lower()
    if run_exit[j.job] == "0":
        # the sequence over each degeneracy interval, called (phased tied members) against truth (the truth's own
        # projected alleles, only where all are model alleles at distance 0); no labelling is chosen
        blocks_m, _ = model(mname)
        x = rep[(rep.run_id == j.job)]
        eqrows = list(csv.DictReader(open(os.path.join(PKG, prefix_of[j.job] + ".equivalence_set.tsv")), delimiter="\t"))
        verdicts = []
        for lo, hi in sorted(iv.get(mname, ())):
            dist = [b[f"truth{h}_closest_distance"].iloc[k] for h in (1, 2) for k in range(lo, hi + 1)]
            if any(d != "0" for d in dist):
                verdicts.append("truth_not_spellable")
                continue
            tseg = sorted("".join(blocks_m[k][int(b[f"truth{h}_closest_allele"].iloc[k])] for k in range(lo, hi + 1))
                          for h in (1, 2))
            members = {tuple(sorted("".join(blocks_m[k][int(r[f"hap{h}_alleles"].split(",")[k])] for k in range(lo, hi + 1))
                                    for h in (1, 2))) for r in eqrows}
            verdicts.append("exact" if members == {tuple(tseg)} else "different" if tuple(tseg) not in members else "mixed")
        segment = ";".join(verdicts) or "NA"
        ex = x[x.block_id.isin(["sv:3", "sv:5"])].excess
        seg_excess = str(bool((ex > 0).any())).lower()
    sv.append({"job": j.job, "arm": j.arm, "sample": j.sample, "coverage": j.coverage,
               "classification": c["classification"], "truth_haplotypes_with_alternative_labelling": exposed,
               "exact_sequence_call_with_label_only_block_excess": label_only,
               "sv3_sv5_segment_sequence": segment, "sv3_or_sv5_positive_block_excess": seg_excess})
sv = pd.DataFrame(sv)
sv.to_csv(os.path.join(D, "sv3_sv5_cells.tsv"), sep="\t", index=False)
svs = []
for (arm, cov), x in sv.groupby(["arm", "coverage"], sort=False):
    svs.append({"arm": arm, "coverage": cov, "cells": len(x),
                "label_ambiguity_no_calls": int((x.classification == LABEL).sum()),
                "samples_with_truth_alternative_labelling": int((x.truth_haplotypes_with_alternative_labelling > 0).sum()),
                "exact_sequence_calls_with_label_only_excess": int((x.exact_sequence_call_with_label_only_block_excess == "true").sum()),
                "calls_with_sv3_or_sv5_excess": int((x.sv3_or_sv5_positive_block_excess == "true").sum()),
                "of_which_segment_sequence_exact": int(((x.sv3_or_sv5_positive_block_excess == "true")
                                                       & (x.sv3_sv5_segment_sequence == "exact")).sum()),
                "of_which_segment_truth_not_spellable": int(((x.sv3_or_sv5_positive_block_excess == "true")
                                                            & (x.sv3_sv5_segment_sequence == "truth_not_spellable")).sum())})
sort_arm(pd.DataFrame(svs)).to_csv(os.path.join(D, "sv3_sv5_summary.tsv"), sep="\t", index=False)
affected_samples = sorted(set(sv[sv.classification == LABEL]["sample"]))

# ---- figures -------------------------------------------------------------------------------------------
order = list(rep.drop_duplicates("block_id").sort_values("block_index").block_id)
alltruth = []
for j in jobs.itertuples():
    b = tsv(os.path.join(PKG, "data", EVAL, f"{j.job}.reported.blocks.tsv"))
    alltruth.append(b[["block_id", "truth1_bp", "truth2_bp"]])
alltruth = pd.concat(alltruth)
alltruth = alltruth[(alltruth.truth1_bp != "NA") & (alltruth.truth2_bp != "NA")]
lengths = pd.DataFrame([{"block_id": bid, "displayed_length_range": compact_range(
    pd.concat([alltruth[alltruth.block_id == bid].truth1_bp, alltruth[alltruth.block_id == bid].truth2_bp]).astype(int).min(),
    pd.concat([alltruth[alltruth.block_id == bid].truth1_bp, alltruth[alltruth.block_id == bid].truth2_bp]).astype(int).max())}
    for bid in order])
lengths.to_csv(os.path.join(D, "block_lengths.tsv"), sep="\t", index=False)
ticks = P.block_labels(order, lengths)
ylims = {}
for rv, _ in ROWS:
    top = max([10.0] + [float(x[x.row == rv].excess.max()) for x in (raw, kept) if (x.row == rv).any()])
    ylims[rv] = (-0.4, top * 1.8)
pd.DataFrame([(k, *v) for k, v in ylims.items()], columns=["row", "ymin", "ymax"]).to_csv(
    os.path.join(D, "panel_axes.tsv"), sep="\t", index=False)


def draw(points, stem, title, footer):
    fig, drawn = P.fig_block_raincloud_rows(points, order, ROWS, COVS, ticks, title, ylims, footer=footer,
                                            coverage_markers=True)
    clipped, overlapping = P.text_layout_problems(fig)
    P.save(fig, os.path.join(OUT, "figures", stem), dpi=P.SLIDE["dpi"])
    with open(os.path.join(D, f"{stem}_layout.txt"), "w") as fh:
        fh.write("clipped:\n" + "\n".join(clipped) + "\noverlapping:\n" + "\n".join(overlapping) + "\n")
    pd.DataFrame(drawn["points"], columns=["row", "metric", "run_id", "block_id", "coverage", "x", "excess"]).to_csv(
        os.path.join(D, f"{stem}_points.tsv"), sep="\t", index=False)
    pd.DataFrame(drawn["violins"], columns=["row", "metric", "block_id", "coverage", "kind", "n", "q1", "median", "q3"]).to_csv(
        os.path.join(D, f"{stem}_violins.tsv"), sep="\t", index=False)
    return len(drawn["points"]), len(clipped), len(overlapping)


def per_cell_n(x):
    g = x.groupby(["row", "coverage"]).run_id.nunique()
    return ", ".join(f"{dict(ROWS)[r]} {c}x {n}" for (r, c), n in sorted(g.items(), key=lambda t: (arm_order[t[0][0]], t[0][1])))


nocalls = int((cells.exit != "0").sum())
fa = draw(raw, "figure_A_raw_block_excess", "Raw block excess at CYP2D6/D7",
          f"Calls only (raw, pre-dosage-repair); {nocalls} no-calls not drawn; violins descriptive only. Pilot profile, "
          f"not built-in defaults.")
captions = {"A_calls_per_panel": per_cell_n(raw), "B_calls_per_panel": per_cell_n(kept)}
excluded = rep[rep.genotype_filter != "PASS"]
fb = draw(kept, "figure_B_retained_block_excess", "Retained block excess at CYP2D6/D7",
          f"GENOTYPE_FILTER == PASS: {len(kept):,}/{len(rep):,} block genotypes retained; {nocalls} no-calls not drawn; "
          f"axes fixed to Figure A's.")
pd.DataFrame([captions]).to_csv(os.path.join(D, "caption_numbers.tsv"), sep="\t", index=False)
excluded.to_csv(os.path.join(D, "excluded_block_calls.tsv"), sep="\t", index=False)
pd.DataFrame([{"figure": "A", "points": fa[0], "expected": len(raw), "clipped_texts": fa[1], "overlapping_text_pairs": fa[2]},
              {"figure": "B", "points": fb[0], "expected": len(kept), "clipped_texts": fb[1], "overlapping_text_pairs": fb[2]}]
             ).to_csv(os.path.join(D, "figure_check.tsv"), sep="\t", index=False)

# ---- built-in defaults against the amended profile, the three sentinels (truth-free) -------------------
cmp_rows = []
for job in ("lzo.HG02769.5", "mosaic.mosaic_one_1.5", "loo.NA19983.5", "lzo.HG02769.30", "mosaic.mosaic_one_1.30",
            "loo.NA19983.30"):
    for label, d in (("default_neighbour_cap_670", "runs"), ("neighbours_958_state_classes_200M", "runs_a3"),
                     (f"pilot_profile_958_600M", RUNS)):
        p = os.path.join(PKG, "data", d, job + ".run.tsv")
        if not os.path.exists(p):
            continue
        rec = tsv(p).iloc[0]
        st = kv(os.path.join(PKG, "data", d, job + ".mosaic_status.tsv"))
        cmp_rows.append({"job": job, "profile": label, "exit": rec["exit"], "wall_seconds": rec.wall_seconds,
                         "peak_rss_GiB": round(int(rec.peak_rss_bytes) / 2 ** 30, 3), "status": st.get("status"),
                         "call_class": st.get("call_class"), "no_call_reason": st.get("no_call_reason", "NA"),
                         "note": st.get("note", "")[:120], "repair_rounds_run": st.get("repair_rounds_run"),
                         "pool_repair_scored": st.get("pool_repair_scored", "NA"),
                         "statistical_state_classes": st.get("statistical_state_classes", "NA")})
pd.DataFrame(cmp_rows).to_csv(os.path.join(D, "default_vs_pilot_profile.tsv"), sep="\t", index=False)
print(f"{len(cells)} cells; {len(raw)} raw and {len(rep)} reported block genotypes from "
      f"{int((cells.exit == '0').sum())} calls; {len(kept)} retained; label-ambiguity samples {affected_samples}")
