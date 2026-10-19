#!/usr/bin/env python3
"""M4-B: freeze the hazard-matched copying length per cell from the IMPLEMENTED adapter geometry (prior-geometry outputs),
bound to the model fingerprint and the geometry/frame digests, before any matched run. lambda_HM = sum_b d_b / ((B-1) *
-log(1 - 0.01)), r = 0.01 being every cell's constant switch probability. A DIAGNOSTIC CONTROL: block-count dependent by
construction, never fitted to accuracy. The decimal written (17 significant digits) must round-trip to the probe's exact double.

  freeze_lambda.py <runs dir> <out LAMBDA.tsv> cell...
"""
import sys

runs, out, cells = sys.argv[1], sys.argv[2], sys.argv[3:]
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l}
rows = ["cell\tmodel_fingerprint\tgeometry_digest\tframe_digest\tgeometry_reference\tconstant_r\tlambda_bp\tlambda_hex\tboundaries\ttotal_distance_bp\t"
        "zero_distance_boundaries\tzero_distance_informative_boundaries"]
for c in cells:
    s = KV(f"{runs}/geometry.{c}/geom.summary.tsv")
    assert s["adapter_geometry_status"] == "available", (c, s["adapter_geometry_status"])
    lam, hx = s["hazard_matched_copying_length_bp"], s["hazard_matched_copying_length_hex"]
    assert float(lam) == float.fromhex(hx), (c, lam, hx)
    rows.append("\t".join([c, s["model_fingerprint"], s["adapter_geometry_digest"], s["adapter_frame_digest"], s["adapter_geometry_reference"],
                           s["hazard_matched_constant_r"], lam, hx, str(int(s["blocks"]) - 1), s["adapter_total_distance_bp"],
                           s["adapter_zero_distance_boundaries"], s["adapter_zero_distance_informative_boundaries"]]))
open(out, "w").write("\n".join(rows) + "\n")
print("\n".join(rows))
