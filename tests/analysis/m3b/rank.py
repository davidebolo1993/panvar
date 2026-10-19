#!/usr/bin/env python3
"""M3-B selection (TRUTH-FREE; frozen before evaluation). For each cell: rank every scored catalogue member by total
(fragment + exposure + sequence-pair prior at r = 0.01, from score-set), report the winner, the incumbent's rank, the
top 10, and -- when the winner came from a two-block move -- the incumbent, both single-block intermediates and the joint
candidate decomposed into fragment, exposure and prior. Ties are reported (equal totals), never broken silently: the
order among exact ties is label order. Writes <out>.selection.tsv, <out>.top.tsv, <out>.decomposition.tsv and a digest.

  rank.py <runs dir> <out prefix> cell...
"""
import csv
import hashlib
import re
import sys

runs, prefix, cells = sys.argv[1], sys.argv[2], sys.argv[3:]


def rows(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


sel, top, dec = [], [], []
for cell in cells:
    # Amendment 1: a cell may be scored whole or in two parts; every part carries the incumbent, which must agree exactly
    import os
    whole = f"{runs}/score.{cell}/score.scores_truth_free.tsv"
    files = [whole] if os.path.exists(whole) else [f"{runs}/score.{cell}.part{i}/score.scores_truth_free.tsv" for i in (1, 2)]
    scores, seen = [], {}
    for f in files:
        for r in rows(f):
            if r["label"] in seen:
                assert r["label"] == "incumbent" and r["total"] == seen[r["label"]]["total"], f"{cell}: part disagreement on {r['label']}"
                continue
            seen[r["label"]] = r
            scores.append(r)
    with open(f"{os.path.dirname(prefix)}/merged.{cell}.tsv", "w") as out:
        keys = list(scores[0].keys())
        out.write("\t".join(keys) + "\n")
        for r in scores:
            out.write("\t".join(r[k] for k in keys) + "\n")
    prov = {r["label"]: r for r in rows(f"{runs}/cat.{cell}/cat.provenance.tsv")}
    pre = dict(r for r in (l.rstrip("\n").split("\t", 1) for l in open(f"{runs}/cat.{cell}/cat.preflight.tsv")) if len(r) == 2)
    for r in scores:
        r["t"] = float(r["total"])
    ranked = sorted(scores, key=lambda r: (-r["t"], r["label"]))
    by = {r["label"]: r for r in scores}
    w = ranked[0]
    inc = by["incumbent"]
    inc_rank = next(i + 1 for i, r in enumerate(ranked) if r["label"] == "incumbent")
    ties = sum(1 for r in ranked if r["t"] == w["t"])
    better = sum(1 for r in ranked if r["t"] > inc["t"])
    kind = "incumbent" if w["label"] == "incumbent" else prov[w["label"]]["first_move_kind"]
    sel.append({"cell": cell, "catalogue_scored": len(scores), "raw_candidates": pre.get("raw_candidates"),
                "unique_candidates": pre.get("unique_candidates_excluding_incumbent"), "distinct_haplotypes": pre.get("distinct_haplotypes_including_incumbent"),
                "winner": w["label"], "winner_kind": kind, "winner_moves": "-" if w["label"] == "incumbent" else prov[w["label"]]["moves"],
                "winner_minus_incumbent": w["t"] - inc["t"], "exact_ties_at_top": ties, "incumbent_rank": inc_rank,
                "candidates_above_incumbent": better, "incumbent_total_reproduced": inc["total"]})
    for i, r in enumerate(ranked[:10]):
        top.append({"cell": cell, "rank": i + 1, "label": r["label"], "moves": "-" if r["label"] == "incumbent" else prov[r["label"]]["moves"],
                    "fragment": r["fragment_term"], "exposure": r["exposure_term"], "prior": r["prior"], "total": r["total"],
                    "minus_incumbent": r["t"] - inc["t"]})
    # decomposition for the best two-block (joint) candidate, winner or not
    best_pair = next((r for r in ranked if r["label"] != "incumbent" and prov[r["label"]]["first_move_kind"] == "pair"), None)
    if best_pair:
        move = next(m for m in prov[best_pair["label"]]["moves"].split("|") if m.startswith("pair:"))
        h, rest = move.split(":", 2)[1], move.split(":", 2)[2]
        parts = rest.split("+")
        inter = []
        for part in parts:
            # Amendment 2: block ids contain ':' (e.g. ctx:0), so split from the right: <block id>:<from>><to>:r<rank>
            blk, change, _ = part.rsplit(":", 2)
            want = f"single:{h}:{blk}:{change}"
            lab = next((l for l, p in prov.items() if any(m.startswith(want + ":") for m in p["moves"].split("|"))), None)
            inter.append(lab)
        for role, lab in [("incumbent", "incumbent"), ("intermediate_1", inter[0]), ("intermediate_2", inter[1]), ("joint", best_pair["label"])]:
            r = by.get(lab) if lab else None
            dec.append({"cell": cell, "role": role, "label": lab or "absent", "move": move if role == "joint" else (prov[lab]["moves"] if lab and lab != "incumbent" else "-"),
                        "fragment": r["fragment_term"] if r else "", "exposure": r["exposure_term"] if r else "", "prior": r["prior"] if r else "",
                        "total": r["total"] if r else "", "minus_incumbent": (r["t"] - inc["t"]) if r else "",
                        "joint_is_winner": best_pair["label"] == w["label"]})


def write(name, items):
    keys = list(items[0].keys()) if items else []
    with open(f"{prefix}.{name}.tsv", "w") as f:
        f.write("\t".join(keys) + "\n")
        for it in items:
            f.write("\t".join(str(it[k]) for k in keys) + "\n")


write("selection", sel)
write("top", top)
write("decomposition", dec)
digest = hashlib.sha256(open(f"{prefix}.selection.tsv", "rb").read() + open(f"{prefix}.decomposition.tsv", "rb").read()).hexdigest()
open(f"{prefix}.selection.sha256", "w").write(digest + "\n")
for s in sel:
    print("\t".join(f"{k}={v}" for k, v in s.items()))
print("selection digest", digest)
