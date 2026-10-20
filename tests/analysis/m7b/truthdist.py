#!/usr/bin/env python3
"""M7-B evaluation-only target-block truth distance (best homologue assignment), with the availability floor implied by
M6-A's recorded block excesses (every M6-A row must imply the same floor).

  truthdist.py <unit dir>
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m7b_common import *  # noqa: E402,F401,F403

d = sys.argv[1]
log = f"{d}/commands.log"
P = f"{OUT}/runs/pools"
mapping = table(f"{P}/mapping.tsv")
excess = {(r["cell"], r["label"]): int(r["block_excess"]) for r in table(f"{R}/m6a_count_model_20261019/runs/summary.a1/m6a.tsv") if r.get("block_excess", "") not in ("", "NA")}
out, ok_all = [], True
for cell, c in CELLS.items():
    pairs = f"{d}/{cell}.pairs.tsv"
    rows = [m for m in mapping if m["cell"] == cell and m["block_interval_1"] != "-"]
    with open(pairs, "w") as f:
        for m in rows:
            for k in (1, 2):
                for t in ("truth1", "truth2"):
                    f.write(f"row{m['row']}.{k}:{t}\t{P}/{cell}.blocks.fa\trow{m['row']}.{k}\t{c['truth_block']}\t{t}\n")
    code = run([SCORING_PROBE, "edit-distances", "--pairs", pairs, "--out", f"{d}/{cell}.dist.tsv", "--threads", 4], log)
    assert code == 0
    dist = {r["label"]: int(r["distance"]) for r in table(f"{d}/{cell}.dist.tsv")}
    floors = set()
    for m in rows:
        i = m["row"]
        dd = min(dist[f"row{i}.1:truth1"] + dist[f"row{i}.2:truth2"], dist[f"row{i}.1:truth2"] + dist[f"row{i}.2:truth1"])
        m6 = [l[4:] for l in m["labels"].split("|") if l.startswith("m6a_")]
        rec = excess.get((cell, m6[0])) if m6 else None
        if rec is not None:
            floors.add(dd - rec)
        out.append(dict(cell=cell, row=i, canonical_diplotype=m["canonical_diplotype"], tier=m["tier"], unit_class=m["unit_class"], truth_distance=dd,
                        m6a_label=m6[0] if m6 else "", m6a_recorded_excess="" if rec is None else rec))
    floor = floors.pop() if len(floors) == 1 else None
    ok_all &= floor is not None
    for r in out:
        if r["cell"] == cell:
            r["floor_implied_by_m6a"] = floor if floor is not None else f"INCONSISTENT {sorted(floors)}"
            r["block_excess"] = r["truth_distance"] - floor if floor is not None else ""
write_tsv(f"{d}/truthdist.tsv", out)
print(open(f"{d}/truthdist.tsv").read()[:4000])
sys.exit(0 if ok_all else 1)
