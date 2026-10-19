#!/usr/bin/env python3
"""M3-E saved-evidence audit (independent of the C++; DESIGN.md section 1). TRUTH-FREE.

Per run (a saved M3-D caller output), per block: the legacy fragment-only support (from .block_genotypes.tsv), the
experimental full-data-likelihood genotype support recomputed from .mosaic_scores.tsv + .mosaic_candidates.tsv, the
phase support (unchanged), and the strongest disagreeing group with call-minus-competitor F, E, L and P differences.
The saved tables are rounded (scores ~1e-5, terms ~1e-7..1e-11): a comparison whose gap is within ROUNDING of the
certification bounds is flagged NEAR_TIE_UNCERTIFIABLE_FROM_ROUNDED rather than certified either way.

  audit_saved.py <out prefix> <label>=<run prefix> ...
"""
import csv
import sys

ROUNDING = 1e-4


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def audit(run):
    cands = {r["candidate_id"]: r for r in rows(run + ".mosaic_candidates.tsv")}
    groups = {}
    for r in rows(run + ".mosaic_scores.tsv"):
        groups[r["sequence_candidate_id"]] = {"F": float(r["fragment_term"]), "E": float(r["exposure_term"]),
                                              "P": float(r["sequence_pair_prior"]), "u": float(r["numeric_error_bound"]),
                                              "S": float(r["score"]), "winner": r["in_winner_set"] == "true", "rows": set()}
    for c in cands.values():
        g = groups.get(c["sequence_candidate_id"])
        if g is None:
            continue
        h1 = [int(x) for x in c["hap1_alleles"].split(",")]
        h2 = [int(x) for x in c["hap2_alleles"].split(",")]
        g["rows"].add(tuple(tuple(sorted(p)) for p in zip(h1, h2)))
    winners = [k for k, g in groups.items() if g["winner"]]
    assert len(winners) == 1, f"{run}: expected a unique winner group, got {len(winners)}"
    w = groups[winners[0]]
    blocks = len(next(iter(w["rows"])))
    legacy = {r["block_id"]: r for r in rows(run + ".block_genotypes.tsv")}
    block_ids = [r["block_id"] for r in rows(run + ".block_genotypes.tsv")]
    out = []
    for b in range(blocks):
        L = lambda g: g["F"] + g["E"]
        wg = {row[b] for row in w["rows"]}
        rec = {"block": block_ids[b], "legacy": legacy[block_ids[b]]["genotype_support"], "phase": legacy[block_ids[b]]["phase_support"]}
        if len(wg) > 1:
            rec.update(label="AMBIGUOUS_PROJECTION")
            out.append(rec)
            continue
        gb = next(iter(wg))
        same, diff, mixed = [], [], []
        for k, g in groups.items():
            at = {row[b] for row in g["rows"]}
            if at == {gb}:
                same.append((k, g))
            elif gb not in at:
                diff.append((k, g))
            else:
                mixed.append((k, g))
        rec["disagreeing"], rec["mixed"] = len(diff), len(mixed)
        if not diff:
            rec.update(label="NO_EVALUATED_ALTERNATIVE")
            out.append(rec)
            continue
        sk, s = max(same, key=lambda kg: (L(kg[1]), kg[0]))
        dk, d = max(diff, key=lambda kg: (L(kg[1]), kg[0]))
        margin = L(s) - L(d)
        above = L(s) - s["u"] > L(d) + d["u"]
        below = L(d) - d["u"] > L(s) + s["u"]
        mixed_ok = all(L(s) - s["u"] > L(g) + g["u"] for _, g in mixed)
        near = abs(margin) <= s["u"] + d["u"] + ROUNDING
        if above and mixed_ok:
            label = "DATA_SUPPORTED"
        elif below:
            label = "PRIOR_DEPENDENT"
        elif not mixed_ok:
            label = "AMBIGUOUS_COMPETITOR"
        else:
            label = "DATA_TIED"
        rec.update(label=label, margin=margin, certified=above or below, near_tie=near, competitor=dk,
                   dF=w["F"] - d["F"], dE=w["E"] - d["E"], dL=L(w) - L(d), dP=w["P"] - d["P"])
        out.append(rec)
    # phase competitors (legacy fragment contenders), classified
    wrow = next(iter(w["rows"])) if len(w["rows"]) == 1 else None
    changing = preserving = mix = 0
    for k, g in groups.items():
        if g is w or g["F"] + g["u"] < w["F"] - w["u"]:
            continue
        kinds = {("same" if wrow is not None and row == wrow else "diff") for row in g["rows"]}
        if kinds == {"diff"}:
            changing += 1
        elif kinds == {"same"}:
            preserving += 1
        else:
            mix += 1
    return out, (changing, preserving, mix), len(groups)


if __name__ == "__main__":
    prefix = sys.argv[1]
    lines = ["run\tblock\tlegacy_genotype_support\tgenotype_support_v1x\tphase_support\tmargin_L\tcertified\tnear_tie\tcompetitor_group\t"
             "delta_F\tdelta_E\tdelta_L\tdelta_P\tdisagreeing_groups\tmixed_groups"]
    summary = ["run\tgroups\tlabels_v1x\tlegacy_labels\tphase_competitors_changing\tphase_competitors_preserving\tphase_competitors_mixed"]
    for arg in sys.argv[2:]:
        label, run = arg.split("=", 1)
        recs, phase, n = audit(run)
        for r in recs:
            f = lambda k, fmt="{:.6f}": (fmt.format(r[k]) if isinstance(r.get(k), float) else str(r.get(k, "-")))
            lines.append("\t".join([label, r["block"], r["legacy"], r["label"], r["phase"], f("margin"), f("certified"), f("near_tie"),
                                    f("competitor"), f("dF"), f("dE"), f("dL"), f("dP"), f("disagreeing"), f("mixed")]))
        from collections import Counter
        summary.append("\t".join([label, str(n), ",".join(f"{k}:{v}" for k, v in sorted(Counter(r["label"] for r in recs).items())),
                                  ",".join(f"{k}:{v}" for k, v in sorted(Counter(r["legacy"] for r in recs).items())), *map(str, phase)]))
    open(prefix + ".blocks.tsv", "w").write("\n".join(lines) + "\n")
    open(prefix + ".summary.tsv", "w").write("\n".join(summary) + "\n")
    print("\n".join(summary))
