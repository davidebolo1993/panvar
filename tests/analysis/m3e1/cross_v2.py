#!/usr/bin/env python3
"""M3-E.1: the C++ block_data_support_v2 (support-audit probe) against the independent Python v2 audit, label by label,
on the same saved tables. Usage: cross_v2.py <python blocks.tsv> <probe dir> <out.tsv>"""
import csv
import os
import sys

py, probe_dir, out = sys.argv[1:4]
ref = {(r["run"], r["block"]): r for r in csv.DictReader(open(py), delimiter="\t")}
lines, failed = ["run\tblock\tpython_v2\tcpp_v2\tagree"], 0
for f in sorted(os.listdir(probe_dir)):
    if not f.endswith(".support_v2.tsv"):
        continue
    run = f[: -len(".support_v2.tsv")]
    for r in csv.DictReader(open(os.path.join(probe_dir, f)), delimiter="\t"):
        p = ref[(run, r["block"])]["v2_label"]
        ok = p == r["label"]
        failed += not ok
        lines.append(f"{run}\t{r['block']}\t{p}\t{r['label']}\t{ok}")
open(out, "w").write("\n".join(lines) + "\n")
print(f"{len(lines) - 1} blocks compared, {failed} disagree")
sys.exit(1 if failed else 0)
