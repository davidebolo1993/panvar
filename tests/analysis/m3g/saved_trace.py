#!/usr/bin/env python3
"""M3-G saved-evidence trace of the ordinary refinement route at one block (TRUTH-GUIDED only through the named target
pair; reads saved caller outputs; computes nothing the caller did not record).

For each repair round r: the round's start incumbent (round 1: repair_anchor_*; round r: round r-1's result), and every
candidate the pool records as first scored in round r (source repair, scored_in round<r>) that differs from the start
incumbent only at BLOCK. Those are the round's scored single-block neighbours at BLOCK; their unordered pairs and
assignments are listed with their saved totals. Then: is TARGET (and CALLED) among them, in which assignment; and does the
saved pool contain TARGET at BLOCK anywhere (source, score). Stages the pool cannot show (generated-but-suppressed,
duplicate, rejected for prior) are reported only as the round's aggregate counters.

  saved_trace.py <run prefix> <block index> <block id> <target x,y> <called x,y> <out prefix>
"""
import csv
import sys

run, b, block_id, target, called, out = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]
target = tuple(sorted(map(int, target.split(","))))
called = tuple(sorted(map(int, called.split(","))))
st = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(run + ".mosaic_status.tsv") if "\t" in l}
rounds = list(csv.DictReader(open(run + ".repair_rounds.tsv"), delimiter="\t"))
cands = list(csv.DictReader(open(run + ".mosaic_candidates.tsv"), delimiter="\t"))
scores = {r["sequence_candidate_id"]: r for r in csv.DictReader(open(run + ".mosaic_scores.tsv"), delimiter="\t")}
P = lambda s: list(map(int, s.split(",")))
start = (P(st["repair_anchor_hap1"]), P(st["repair_anchor_hap2"]))
lines = ["round\tstart_ctx_pair\tmoved\trepaired_block\tblock_neighbours_scored\tdistinct_unordered_pairs\ttarget_scored\tcalled_scored\t"
         "best_block_neighbour_total_minus_start\tround_generated\tround_scored\tsuppressed_visited\tduplicate\tself_skipped\trejected_infinite_prior"]
detail = ["round\tcandidate_id\tpair\tassignment_hap1\tassignment_hap2\ttotal\ttotal_minus_start"]
for r in rounds:
    k = r["round"]
    h1, h2 = start
    start_total = float(r["start_total"])
    nb = []
    for c in cands:
        if c["source"] != "repair" or c["scored_in"] != f"round{k}":
            continue
        a, d = P(c["hap1_alleles"]), P(c["hap2_alleles"])
        others = [i for i in range(len(a)) if i != b]
        if all(a[i] == h1[i] and d[i] == h2[i] for i in others) or all(a[i] == h2[i] and d[i] == h1[i] for i in others):
            if (a[b], d[b]) != (h1[b], h2[b]):
                tot = float(scores[c["sequence_candidate_id"]]["score"])
                nb.append((tuple(sorted((a[b], d[b]))), a[b], d[b], tot, c["candidate_id"]))
    pairs = {x[0] for x in nb}
    best = max((x[3] for x in nb), default=None)
    lines.append("\t".join(map(str, [k, f"{h1[b]},{h2[b]}", r["moved"], r["repaired_block"], len(nb), len(pairs), target in pairs, called in pairs,
                                     f"{best - start_total:.6f}" if best is not None else "NA", r["neighbours_generated"], r["neighbours_scored"],
                                     r["suppressed_visited"], r["duplicate"], r["self_skipped"], r["rejected_infinite_prior"]])))
    for x in sorted(nb, key=lambda x: -x[3]):
        detail.append("\t".join(map(str, [k, x[4], f"{x[0][0]},{x[0][1]}", x[1], x[2], f"{x[3]:.6f}", f"{x[3] - start_total:.6f}"])))
    if r["moved"] == "true":
        start = (P(r["hap1_alleles"]), P(r["hap2_alleles"]))
anywhere = []
for c in cands:
    a, d = P(c["hap1_alleles"]), P(c["hap2_alleles"])
    if tuple(sorted((a[b], d[b]))) == target:
        anywhere.append(f"{c['candidate_id']}:{c['source']}:{c['scored_in']}:{scores[c['sequence_candidate_id']]['score']}")
src = {}
for c in cands:
    a, d = P(c["hap1_alleles"]), P(c["hap2_alleles"])
    for al in (a[b], d[b]):
        if al in target:
            src.setdefault((al, c["source"]), set()).add(c["candidate_id"])
open(out + ".rounds.tsv", "w").write("\n".join(lines) + "\n")
open(out + ".block_neighbours.tsv", "w").write("\n".join(detail) + "\n")
summary = [f"block\t{block_id} (index {b})", f"target_pair\t{target}", f"called_pair\t{called}",
           f"repair_final_incumbent\t{st['repair_final_incumbent_hap1']} | {st['repair_final_incumbent_hap2']}",
           f"pool_candidates_with_target_at_block\t{len(anywhere)}: " + "; ".join(anywhere),
           "pool_candidates_carrying_a_target_allele_at_block\t" + "; ".join(f"allele {al} source {s}: {len(v)}" for (al, s), v in sorted(src.items()))]
open(out + ".summary.tsv", "w").write("\n".join(summary) + "\n")
print("\n".join(summary))
print("\n".join(lines))
