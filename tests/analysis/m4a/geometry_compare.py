#!/usr/bin/env python3
"""M4-A geometry comparison (diagnostic; reads prior-geometry outputs only; no prior is evaluated on reads).

Per model: block/boundary counts, total distance under each candidate geometry, the share of total carried by the largest
boundary, non-positive core-midpoint distances, the ratio of reference core-midpoint distance to M3-A's panel-mean span per
boundary, and the hazard-matched copying length a diagnostic control would use (total / ((B-1) * -log(1-0.01))), which
is block-count dependent BY CONSTRUCTION. Per declared pair (same graph curated vs uncurated; full graph vs a rebuilt
leave-one-out graph): whether the block frame is the same, and per geometry the largest per-boundary change.

  geometry_compare.py <out prefix> <label>=<probe prefix> ... -- <labelA>,<labelB> ...
"""
import csv
import math
import sys

prefix = sys.argv[1]
sep = sys.argv.index("--")
models = dict(a.split("=", 1) for a in sys.argv[2:sep])
pairs = [p.split(",") for p in sys.argv[sep + 1:]]
G = ["d_panel_mean_m3a", "d_ref_start", "d_ref_block_mid", "d_ref_core_mid", "d_invariant_gap"]
T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l}
data = {k: (KV(v + ".summary.tsv"), T(v + ".boundaries.tsv"), T(v + ".blocks.tsv")) for k, v in models.items()}
out = ["model\tblocks\ttemplates\treference_tiles\tref_template_agrees\tgeometry\ttotal_bp\tlargest_boundary_share\tnonpositive\tratio_to_panel_mean_min\t"
       "ratio_to_panel_mean_median\tratio_to_panel_mean_max\thazard_matched_length_bp(control)"]
for k, (s, bd, bl) in data.items():
    n = len(bd)
    for g in G:
        vals = [float(r[g]) for r in bd]
        tot = sum(vals)
        ratios = sorted(float(r[g]) / float(r["d_panel_mean_m3a"]) for r in bd if float(r["d_panel_mean_m3a"]) > 0)
        out.append("\t".join(map(str, [k, s["blocks"], s["templates"], s["reference_tiles_chain"], s["reference_template_alleles_agree"], g, f"{tot:.1f}",
                                       f"{max(vals) / tot:.3f}" if tot > 0 else "NA", sum(v <= 0 for v in vals), f"{ratios[0]:.3f}" if ratios else "NA",
                                       f"{ratios[len(ratios) // 2]:.3f}" if ratios else "NA", f"{ratios[-1]:.3f}" if ratios else "NA",
                                       f"{tot / (n * -math.log1p(-0.01)):.1f}" if n else "NA"])))
out.append("")
out.append("pair\tsame_block_frame\tgeometry\tmax_abs_boundary_change_bp\ttotal_change_bp")
for a, b in pairs:
    sa, ba, la = data[a]
    sb, bb, lb = data[b]
    same = [r["block"] for r in la] == [r["block"] for r in lb] and [r["ref_bp"] for r in la] == [r["ref_bp"] for r in lb]
    for g in G:
        if same:
            ch = [abs(float(x[g]) - float(y[g])) for x, y in zip(ba, bb)]
            out.append(f"{a} vs {b}\t{same}\t{g}\t{max(ch):.2f}\t{sum(float(y[g]) for y in bb) - sum(float(x[g]) for x in ba):.2f}")
        else:
            out.append(f"{a} vs {b}\t{same}\t{g}\tNA (frame differs: {len(la)} vs {len(lb)} blocks)\t"
                       f"{sum(float(y[g]) for y in bb) - sum(float(x[g]) for x in ba):.2f}")
open(prefix + ".tsv", "w").write("\n".join(out) + "\n")
print("\n".join(out))
