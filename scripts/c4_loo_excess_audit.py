#!/usr/bin/env python3
"""C4 production-v1 LOO: where the excess above the availability floor sits (evaluation only).

    scripts/c4_loo_excess_audit.py --campaign <campaign dir> --reeval <re-evaluation dir> --out <dir>

Reads the frozen campaign (jobs, completeness, each call's block genotypes and pooled candidate table) and the
re-evaluation written by tests/genotype_c4_loo_excess_reeval.sh (per-block called alleles and each truth
haplotype's closest, floor, allele). Writes:

  loo_excess_audit.tsv          one row per LOO call: distance, floor, excess, its stratum (0, 1-10, 11-100,
                                101-1000, >1000 bp), and the blocks carrying the excess
  loo_excess_blocks.tsv         one row per block with excess > 0 in a call above 10 bp: called pair, the
                                truth-exact pair (the allele equal to each truth haplotype, when one exists),
                                the floor pair, lengths, distance/floor/excess, genotype support, GQ, phase
                                support, and whether the floor pair occurs at that block in the union or
                                among the candidates the repair scored
Undefined truth blocks stay undefined. Nothing here feeds a call.
"""
import argparse
import csv
import pathlib
from collections import defaultdict


def kv(path):
    out = {}
    try:
        for line in open(path):
            parts = line.rstrip("\n").split("\t", 1)
            if len(parts) == 2 and parts[0] not in out:
                out[parts[0]] = parts[1]
    except FileNotFoundError:
        pass
    return out


def table(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def stratum(excess):
    if excess == 0:
        return "0"
    if excess <= 10:
        return "1-10"
    if excess <= 100:
        return "11-100"
    if excess <= 1000:
        return "101-1000"
    return ">1000"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign", required=True)
    parser.add_argument("--reeval", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    campaign, reeval, out = pathlib.Path(args.campaign), pathlib.Path(args.reeval), pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    jobs = {r["job"]: r for r in table(campaign / "jobs.tsv")}
    calls = [r for r in table(campaign / "completeness.tsv")
             if r["campaign"] == "loo" and r["outcome"] in ("strict_call", "genotype_equivalent_call")]
    run_rows, block_rows = [], []
    for entry in sorted(calls, key=lambda r: (int(r["coverage"]), r["subject"])):
        job = entry["job"]
        prefix = jobs[job]["out_prefix"]
        summary = kv(reeval / f"{job}.summary.tsv")
        blocks = table(reeval / f"{job}.blocks.tsv")
        excess = int(summary["excess_over_floor_total"])
        carrying = [b for b in blocks if b["excess"] not in ("NA", "") and int(b["excess"]) > 0]
        run_rows.append({"job": job, "donor": entry["subject"], "coverage": entry["coverage"],
                         "call_class": entry["outcome"], "phase_support": entry["phase_support"],
                         "distance": summary["unordered_distance_total"], "floor": summary["availability_floor_total"],
                         "excess": excess, "stratum": stratum(excess),
                         "blocks_with_excess": ",".join(f"{b['block_id']}:{b['excess']}" for b in carrying) or "none",
                         "undefined_truth_blocks": summary.get("undefined_truth_blocks", "NA")})
        if excess <= 10:
            continue
        genotype = {r["block_id"]: r for r in table(prefix + ".block_genotypes.tsv")}
        # Which unordered allele pairs occur, per block, in the union and among repair-scored candidates.
        in_union, in_repair = defaultdict(set), defaultdict(set)
        order = [b["block_id"] for b in blocks]
        for cand in table(prefix + ".mosaic_candidates.tsv"):
            h1 = cand["hap1_alleles"].split(",")
            h2 = cand["hap2_alleles"].split(",")
            target = in_repair if cand["source"] == "repair" else in_union
            for i, block_id in enumerate(order):
                target[block_id].add(tuple(sorted((int(h1[i]), int(h2[i])))))
        for b in carrying:
            floor_pair = tuple(sorted((int(b["truth1_closest_allele"]), int(b["truth2_closest_allele"]))))
            exact = []
            for side in ("1", "2"):
                exact.append(b[f"truth{side}_closest_allele"] if b[f"truth{side}_closest_distance"] == "0" else "none")
            g = genotype.get(b["block_id"], {})
            block_rows.append({
                "job": job, "donor": entry["subject"], "coverage": entry["coverage"], "block_id": b["block_id"],
                "kind": b["kind"], "called_pair": f"{b['called_allele1']}/{b['called_allele2']}",
                "called_bp": f"{b['called1_bp']}/{b['called2_bp']}",
                "truth_exact_pair": "/".join(exact), "truth_bp": f"{b['truth1_bp']}/{b['truth2_bp']}",
                "floor_pair": f"{floor_pair[0]}/{floor_pair[1]}",
                "floor_distances": f"{b['truth1_closest_distance']}/{b['truth2_closest_distance']}",
                "distance": b["unordered_distance"], "floor": b["availability_floor"], "excess": b["excess"],
                "genotype_support": g.get("genotype_support", "NA"), "conditional_GQ": g.get("conditional_GQ", "NA"),
                "phase_support": g.get("phase_support", "NA"),
                "floor_pair_in_union_at_block": floor_pair in in_union[b["block_id"]],
                "floor_pair_scored_by_repair_at_block": floor_pair in in_repair[b["block_id"]],
                "local_pair_rank": "NA (needs the truth-based diagnostic)",
            })
    for name, rows in (("loo_excess_audit.tsv", run_rows), ("loo_excess_blocks.tsv", block_rows)):
        with open(out / name, "w") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    counts = defaultdict(lambda: defaultdict(int))
    for r in run_rows:
        counts[r["coverage"]][r["stratum"]] += 1
    for cov in sorted(counts, key=int):
        print(cov, dict(counts[cov]))


if __name__ == "__main__":
    main()
