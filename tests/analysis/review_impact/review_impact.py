#!/usr/bin/env python3
"""Impact screens for review findings C2 and C3 on existing sealed outputs (read-only; no truth is read).

These screens bound exposure; they do not show that a finding changed any historical call.

C2 (representative-only block support/GQ for lower-ranked sequence groups): for every call, counts lower-ranked sequence
groups (not in the winner set) that are projection-ambiguous (several block projections) within W nats of the best score.
A group whose representative CARRIES the called genotype is "susceptible": another of its projections may carry a
different block genotype that the current annotation never counts. A group whose representative does not carry it is
already counted as error mass.
C3 (dosage route depends on node representation): inventories, per block, the dosage census route and its graph inputs
(max node multiplicity, max repeated bp) from the sealed C4 dosage replay; routes that rest on node multiplicity are
representation-sensitive by construction. Also tallies the dosage audit status of the CYP2D6 pilot calls.

Usage: review_impact.py <out dir> <results root>
"""
import csv
import glob
import os
import sys
from collections import Counter

OUT, ROOT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
WINDOWS = (10, 20, 50)


def kv(path):
    out = {}
    for line in open(path).read().splitlines():
        p = line.split("\t")
        if len(p) >= 2:
            out.setdefault(p[0], p[1])
    return out


def c2_rows(label, files):
    rows = []
    for f in sorted(files):
        st = kv(f.replace(".mosaic_scores.tsv", ".mosaic_status.tsv"))
        if st.get("call_class") not in ("strict_call", "genotype_equivalent_call"):
            continue
        with open(f) as h:
            scores = list(csv.DictReader(h, delimiter="\t"))
        top = max(float(r["score"]) for r in scores)
        r = {"set": label, "run": os.path.relpath(f, ROOT), "call_class": st["call_class"], "groups": len(scores),
             "projection_ambiguous_groups": sum(x["projection_ambiguous"] == "true" for x in scores)}
        for w in WINDOWS:
            lower = [x for x in scores if x["in_winner_set"] != "true" and x["projection_ambiguous"] == "true"
                     and top - float(x["score"]) <= w]
            r[f"lower_ambiguous_within_{w}"] = len(lower)
            r[f"susceptible_within_{w}"] = sum(x["carries_called_genotype"] == "true" for x in lower)
        rows.append(r)
    return rows


c2 = c2_rows("cyp2d6_pilot_a5", glob.glob(os.path.join(ROOT, "cyp2d6_shortread_pilot_v1/data/runs_a5/*.mosaic_scores.tsv")))
c2 += c2_rows("c4_production_11fe681", [f for f in glob.glob(os.path.join(ROOT, "genotype_c4_production_11fe681/*/*/call.*.mosaic_scores.tsv"))
                                         if ".original." not in f])
c2 = [r for r in c2 if ".original." not in r["run"]]
with open(os.path.join(OUT, "c2_lower_group_exposure.tsv"), "w") as f:
    w = csv.DictWriter(f, fieldnames=list(c2[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(c2)
summary = []
for label in sorted({r["set"] for r in c2}):
    x = [r for r in c2 if r["set"] == label]
    s = {"set": label, "calls": len(x)}
    for wd in WINDOWS:
        s[f"calls_with_lower_ambiguous_within_{wd}"] = sum(r[f"lower_ambiguous_within_{wd}"] > 0 for r in x)
        s[f"calls_with_susceptible_within_{wd}"] = sum(r[f"susceptible_within_{wd}"] > 0 for r in x)
    summary.append(s)
with open(os.path.join(OUT, "c2_summary.tsv"), "w") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(summary)

# C3
census = {}
for f in glob.glob(os.path.join(ROOT, "genotype_c4_production_11fe681_dosage_repair/final_26401b4/replay/campaign/*.dosage_census.tsv")):
    model = "loo" if os.path.basename(f).startswith("loo.") else "full_panel"
    with open(f) as h:
        for r in csv.DictReader(h, delimiter="\t"):
            census.setdefault((model, r["block_id"]), Counter())[(r["census_route"], r["max_node_multiplicity"], r["max_repeated_bp"])] += 1
c3 = [{"model_kind": m, "block_id": b, "route_inputs_observed": "; ".join(f"{k[0]} (node multiplicity {k[1]}, repeated {k[2]} bp) x{n}"
                                                                    for k, n in sorted(c.items()))}
      for (m, b), c in sorted(census.items())]
with open(os.path.join(OUT, "c3_c4_census_routes.tsv"), "w") as f:
    w = csv.DictWriter(f, fieldnames=list(c3[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(c3)
cyp = Counter()
for f in glob.glob(os.path.join(ROOT, "cyp2d6_shortread_pilot_v1/data/runs_a5/*.dosage_status.tsv")):
    s = kv(f)
    cyp[(s.get("dosage_audit_status"), s.get("dosage_audit_note", "")[:90])] += 1
with open(os.path.join(OUT, "c3_cyp2d6_dosage_audit.tsv"), "w") as f:
    f.write("dosage_audit_status\tnote\tcalls\n" + "".join(f"{k[0]}\t{k[1]}\t{n}\n" for k, n in sorted(cyp.items())))
for s in summary:
    print(s)
print("C3 C4 census blocks:", len(c3), "| CYP2D6 dosage audit:", dict(cyp))
