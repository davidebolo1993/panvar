#!/usr/bin/env python3
"""M7-B.1 reuse check: bind the stopped M7-B artefacts that are reused (provenance verdicts, canonical pools, Panvar
scores, truth distances) by SHA-256, confirm their recorded gates passed, and confirm every executable is unchanged.

  reuse.py <unit dir>
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m7b_common import *  # noqa: E402,F401,F403

d = sys.argv[1]
rows, ok = [], True
led = {r["unit"]: r for r in table(f"{OLD}/ledger.tsv")}
need = ["provenance", "pools", "truthdist.a2"] + [f"panvar.D.{c}" for c in CELLS] + ["panvar.X.cyp_loo_na19650_5"]
for u in need:
    good = led.get(u, {}).get("exit") == "0"
    rows.append(dict(item=f"M7-B unit {u}", value=f"exit {led.get(u, {}).get('exit')}", ok=good))
    ok &= good
for r in table(f"{OLD}/runs/provenance/lineage.tsv"):
    rows.append(dict(item=f"lineage {r['lineage']}", value=f"{r['verdict']} ({r['simulator']})", ok=r["verdict"] == "PASS"))
    ok &= r["verdict"] == "PASS"
for r in table(f"{OLD}/runs/pools/pools_checks.tsv"):
    rows.append(dict(item=f"pools check {r['cell']}: {r['check'][:60]}", value=r.get("value", ""), ok=r["ok"] == "True"))
    ok &= r["ok"] == "True"
for c in list(CELLS) + ["cyp_loo_na19650_5:X"]:
    cell, tier = (c.split(":") + ["D"])[:2]
    lg = open(f"{OLD}/runs/panvar.{tier}.{cell}/run.log").read()
    m = [l for l in lg.splitlines() if "declared" in l and "expectations reproduced" in l]
    rows.append(dict(item=f"Panvar score-set {tier} {cell}", value=m[0].split(": ", 1)[1] if m else "missing", ok=bool(m)))
    ok &= bool(m)
for r in table(f"{OLD}/runs/provenance/executables.tsv"):
    if r["path"].startswith("/") and os.path.isfile(r["path"]):
        same = sha(r["path"]) == r["sha256"]
        rows.append(dict(item=f"executable unchanged: {r['tool']}", value=r["sha256"][:16], ok=same))
        ok &= same
bind = []
for root, _, files in os.walk(f"{OLD}/runs"):
    for f in sorted(files):
        p = os.path.join(root, f)
        rel = os.path.relpath(p, OLD)
        if rel.split("/")[1] in ("provenance", "pools", "truthdist.a2", "compare.a2") or rel.split("/")[1].startswith("panvar."):
            if not rel.endswith((".fq", ".gfa")) and os.path.getsize(p) < 200_000_000:
                bind.append(dict(file=rel, sha256=sha(p)))
write_tsv(f"{d}/reused_bindings.tsv", bind)
write_tsv(f"{d}/reuse_checks.tsv", rows)
print(open(f"{d}/reuse_checks.tsv").read())
print(f"bound {len(bind)} reused files")
sys.exit(0 if ok else 1)
