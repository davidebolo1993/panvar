#!/usr/bin/env python3
"""M9-E bounded reachability audit of the three missed improving neighbours (E1, E3, E6 of M9-D). Read-only.

  m9e.py fixture <out>   context-aware inconsistency classification on synthetic rows (A vs reported winner; B carrier)
  m9e.py analyse <out>   terminal contexts and their final neighbourhood scans; the probe's sv:5 local-pair ranking
                         (reproduction gate: every recorded trunk move's local_pair_rank); rank, emission, cutoff gap and
                         ties of each frozen improving pair; marker availability from allele lengths (markers are syncmers
                         inside each block allele only) and emission tie classes; generation trace per assignment;
                         sequence-level scored-pool membership; the context-aware reclassification of M9-D's contrasts
"""
import csv
import hashlib
import os
import sys

mode, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
R = "/Users/davide.bolognini/panvar/tests/results"
C = f"{R}/m9c_acot_reliability_20261020/runs"
D = f"{R}/m9d_attribution_20261020/runs"
OUT = f"{R}/m9e_reachability_20261020/runs"
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l and not l.startswith("#")} if os.path.exists(p) else {}
labs = lambda s: [int(x) for x in s.split(",")]
key = lambda h1, h2: tuple(sorted((tuple(h1), tuple(h2))))
CASES = [("E1", "LOO.d1"), ("E3", "LOO.d3"), ("E6", "LOO.d6")]
K = 31


def write(name, rows):
    with open(f"{out}/{name}.tsv", "w") as f:
        if not rows:
            f.write("none\n")
            return
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


def classify(context, replacement_in_pool, replacement_total, context_total, winner_total, carrier_is_context):
    """Context-aware reading. An inconsistency is claimed only against the REPORTED WINNER: a hypothesis in the scored
    pool that scores above it. In context B the truth carrier is the context (in the pool, below the winner by
    construction of the winner); its beating a called-pair insertion is not by itself an inconsistency."""
    if context == "A":
        if replacement_total > winner_total:
            return "INCONSISTENCY (in pool, above the reported winner)" if replacement_in_pool else "SEARCH_GAP (improves the winner, never evaluated)"
        return "NO_IMPROVEMENT_OVER_WINNER"
    if carrier_is_context:
        if replacement_in_pool and replacement_total > winner_total:
            return "INCONSISTENCY (in pool, above the reported winner)"
        return "CONSISTENT (B carrier beats its called insertion; carrier not above the winner)" if context_total > replacement_total \
            else "CONSISTENT (called insertion beats the B carrier)"
    return "UNCLASSIFIED"


if mode == "fixture":
    cases = [("A replacement in pool above winner", ("A", True, -10.0, -12.0, -12.0, False), "INCONSISTENCY"),
             ("A replacement absent above winner", ("A", False, -10.0, -12.0, -12.0, False), "SEARCH_GAP"),
             ("A replacement below winner", ("A", False, -15.0, -12.0, -12.0, False), "NO_IMPROVEMENT_OVER_WINNER"),
             ("B carrier beats absent insertion, carrier below winner", ("B", False, -30.0, -20.0, -12.0, True), "CONSISTENT"),
             ("B called insertion absent, scores above winner", ("B", False, -10.0, -20.0, -12.0, True), "CONSISTENT"),
             ("B insertion in pool and above winner", ("B", True, -10.0, -20.0, -12.0, True), "INCONSISTENCY")]
    rows = [dict(case=n, got=classify(*a), expected=e, result="PASS" if classify(*a).startswith(e) else "FAIL") for n, a, e in cases]
    write("fixture", rows)
    print(open(f"{out}/fixture.tsv").read())
    sys.exit(0 if all(r["result"] == "PASS" for r in rows) else 1)

if mode == "analyse":
    blocks = []
    for r in T(f"{D}/models/LOO.block_alleles.tsv"):
        b, a = int(r["block_index"]), int(r["allele_index"])
        while len(blocks) <= b:
            blocks.append([])
        blocks[b].append(r["sequence"])
    term, gate, short, reach, trace, member = [], [], [], [], [], []
    m9d_contrasts = T(f"{D}/analyse/contrasts.tsv")
    m9d_rep = {(r["case"], r["context"]): r["label"] for r in T(f"{D}/analyse/reproduction.tsv")}
    for case, cell in CASES:
        d = f"{C}/call.{cell}"
        ids = [r["block_id"] for r in T(f"{d}/call.sequence_blocks.tsv")]
        b5, b6 = ids.index("sv:5"), ids.index("ctx:6")
        lg = T(f"{d}/call.locus_genotype.tsv")[0]
        win = key(labs(lg["hap1_block_path"]), labs(lg["hap2_block_path"]))
        rounds = T(f"{d}/call.repair_rounds.tsv")
        moved = [r for r in rounds if r["moved"] == "true"]
        trunk_end = key(labs(moved[-1]["hap1_alleles"]), labs(moved[-1]["hap2_alleles"])) if moved else None
        last = rounds[-1]
        branches = T(f"{d}/call.tie_branches.tsv")
        cand = T(f"{d}/call.mosaic_candidates.tsv")
        scored_in = {}
        for c in cand:
            scored_in[c["scored_in"]] = scored_in.get(c["scored_in"], 0) + 1
        row = dict(case=case, cell=cell, winner_equals_trunk_terminal=trunk_end == win, trunk_rounds=len(rounds), trunk_moves=len(moved),
                   trunk_final_round=last["round"], trunk_final_status=last["status"], trunk_final_moved=last["moved"],
                   trunk_final_generated=last["neighbours_generated"], trunk_final_scored=last["neighbours_scored"],
                   trunk_final_duplicate=last["duplicate"], trunk_final_suppressed_visited=last["suppressed_visited"],
                   trunk_final_self_skipped=last["self_skipped"], trunk_final_rejected_infinite_prior=last["rejected_infinite_prior"],
                   tie_branches=len(branches))
        for br in branches:
            t = key(labs(br["terminal_hap1"]), labs(br["terminal_hap2"]))
            row[f"branch{br['branch']}_terminal_equals_winner"] = t == win
            row[f"branch{br['branch']}_rounds"] = br["rounds"]
            row[f"branch{br['branch']}_stop"] = br["stop"]
            row[f"branch{br['branch']}_final_round_scored"] = scored_in.get(f"branch{br['branch']}_round{br['rounds']}", 0)
        term.append(row)
        # ---- ranking reproduction gate and the sv:5 ranks
        ranks = {}
        for r in T(f"{OUT}/ranks.{cell}/score.local_ranks.tsv"):
            ranks.setdefault(r["block_id"], {})[(int(r["allele_x"]), int(r["allele_y"]))] = r
        local_pairs = int(KV(f"{d}/call.mosaic_status.tsv").get("repair_local_pairs", "32"))
        for r in moved:
            p = (int(r["allele_x"]), int(r["allele_y"]))
            got = ranks.get(r["repaired_block"], {}).get(p, {}).get("rank", "absent")
            gate.append(dict(case=case, cell=cell, round=r["round"], block=r["repaired_block"], pair=f"{p[0]}/{p[1]}",
                             recorded_rank=r["local_pair_rank"], probe_rank=got, ok=got == r["local_pair_rank"]))
        sv5 = ranks.get("sv:5", {})
        ordered = sorted(sv5.values(), key=lambda r: int(r["rank"]))
        cutoff = ordered[local_pairs - 1]
        emis = [float(r["emission"]) for r in ordered]
        # ---- frozen improving alternatives from M9-D (context A truth replacements above the reported winner)
        ctx_label = m9d_rep[(case, "A")]
        improving = [r for r in m9d_contrasts if r["case"] == case and r["context"] == "A" and r["called"] == ctx_label and float(r["d_total"]) > 0]
        mem = {r["label"]: r for r in T(f"{D}/construct/membership.tsv")}
        cands = {}
        for line in open(f"{D}/construct/{case}.candidates.tsv"):
            f = line.rstrip("\n").split("\t")
            cands[f[0]] = (labs(f[1]), labs(f[2]))
        for imp in improving:
            lab = imp["truth_or_inserted"]
            h1, h2 = cands[lab]
            pair = tuple(sorted((h1[b5], h2[b5])))
            pr = sv5[pair]
            e = float(pr["emission"])
            reach.append(dict(case=case, improving=lab, d_total_vs_winner=imp["d_total"], sv5_pair=f"{pair[0]}/{pair[1]}", rank=pr["rank"], emission=e,
                              cutoff_rank=local_pairs, cutoff_emission=cutoff["emission"], gap_below_cutoff=float(cutoff["emission"]) - e,
                              in_local_list=pr["in_local_list"], pairs_tied_at_its_emission=sum(1 for x in emis if x == e),
                              pairs_tied_at_cutoff_emission=sum(1 for x in emis if x == float(cutoff["emission"])),
                              first_rank_at_cutoff_emission=next(int(r["rank"]) for r in ordered if float(r["emission"]) == float(cutoff["emission"])),
                              tie_rule="emission descending, then (allele_x, allele_y) ascending (genotype_repair.cpp CoordinateRepairSearch)"))
            # generation trace per assignment: the move sets sv:5 on both haplotypes from one local pair
            for asg, ok in (("this assignment", True),):
                in_list = pr["in_local_list"] == "true"
                stage = "excluded by the shortlist (never generated in any round or branch: the list is built once)" if not in_list \
                    else "in the shortlist: generated at the terminal scan unless self/visited/duplicate/prior-rejected (see membership)"
                trace.append(dict(case=case, improving=lab, assignment=f"hap1={h1[b5]} hap2={h2[b5]}", stage=stage))
            # sequence-level scored-pool membership (only candidates equal outside sv:5..ctx:6 can spell the same pair:
            # M9-D's inventory finds no within-block duplicates and cross-block equivalences only in sv:5..ctx:6)
            seg = lambda h: blocks[b5][h[b5]] + blocks[b6][h[b6]]
            outside = lambda h: tuple(h[:b5] + h[b6 + 1:])
            want = tuple(sorted((outside(h1), outside(h2))))
            checked = same = 0
            hits = []
            for c in cand:
                c1, c2 = labs(c["hap1_alleles"]), labs(c["hap2_alleles"])
                if tuple(sorted((outside(c1), outside(c2)))) != want:
                    continue
                checked += 1
                # pair the haplotypes by their outside labels, then compare the interval bytes (either pairing)
                m = (outside(c1) == outside(h1) and outside(c2) == outside(h2) and (seg(c1), seg(c2)) == (seg(h1), seg(h2))) or \
                    (outside(c1) == outside(h2) and outside(c2) == outside(h1) and (seg(c1), seg(c2)) == (seg(h2), seg(h1)))
                if m:
                    same += 1
                    hits.append(c["candidate_id"])
            member.append(dict(case=case, improving=lab, label_membership=mem[lab]["membership"], pool_candidates_equal_outside_interval=checked,
                               sequence_identical_in_pool=same, hits=",".join(hits) or "none"))
        # ---- short-allele / marker availability at sv:5
        called = tuple(sorted((labs(lg["hap1_block_path"])[b5], labs(lg["hap2_block_path"])[b5])))
        lens = [len(s) for s in blocks[b5]]
        distinct = sorted(set(emis), reverse=True)
        top_class = sum(1 for x in emis if x == distinct[0])
        no_marker_alleles = sum(1 for L in lens if L < K)
        listed = [(int(r["allele_x"]), int(r["allele_y"])) for r in ordered[:local_pairs]]
        for role, pair in (("called", called), ("improving", pair)):
            for a in pair:
                short.append(dict(case=case, role=role, allele=a, bp=lens[a], internal_31mers=max(0, lens[a] - K + 1),
                                  proposal_markers_possible=lens[a] >= K))
        short.append(dict(case=case, role="block sv:5 summary", allele="all", bp=f"{min(lens)}-{max(lens)}", internal_31mers=f"{no_marker_alleles}/{len(lens)} alleles shorter than 31 bp (no markers possible)",
                          proposal_markers_possible=f"{len(distinct)} distinct pair emissions over {len(emis)} pairs; top tie class {top_class}; "
                                                    f"local list = ranks 1-{local_pairs}: {';'.join(f'{x}/{y}' for x, y in listed[:12])}..."))
    # ---- context-aware reclassification of every M9-D contrast (E1-E6, C1-C2)
    recl = []
    for r in m9d_contrasts:
        case, cx = r["case"], r["context"]
        ctx_lab = m9d_rep[(case, cx)]
        sc = {x["label"]: x for x in T(f"{D}/score.{case}/score.scores_truth_free.tsv")}
        winner_total = float(sc[m9d_rep[(case, "A")]]["total"])
        a, b = r["truth_or_inserted"], r["called"]
        memb = {x["label"]: x for x in T(f"{D}/construct/membership.tsv")}
        if cx == "A":
            replacement = a if b == ctx_lab else b
            got = classify("A", not memb[replacement]["membership"].startswith("absent"), float(sc[replacement]["total"]), float(sc[ctx_lab]["total"]), winner_total, False)
        else:
            replacement = b if a == ctx_lab else a
            got = classify("B", not memb[replacement]["membership"].startswith("absent"), float(sc[replacement]["total"]), float(sc[ctx_lab]["total"]), winner_total, True)
        recl.append(dict(case=case, context=cx, contrast=f"{a} - {b}", m9d_reading=r["reading"], context_aware_reading=got))
    write("terminal_contexts", term)
    write("rank_gate", gate)
    write("reachability", reach)
    write("generation_trace", trace)
    write("sequence_membership", member)
    write("marker_availability", short)
    write("m9d_reclassified", recl)
    for n in ("terminal_contexts", "rank_gate", "reachability", "sequence_membership", "marker_availability"):
        print(f"## {n}\n" + open(f"{out}/{n}.tsv").read())
    print("reclassified:", {k: sum(1 for x in recl if x["context_aware_reading"].startswith(k)) for k in ("INCONSISTENCY", "SEARCH_GAP", "NO_IMPROVEMENT", "CONSISTENT")})
    sys.exit(0 if gate and all(g["ok"] for g in gate) else 1)
