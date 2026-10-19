#!/usr/bin/env python3
"""M4-B matched const/dist comparator (adapted from M3-F compare_m3f.py) (TRUTH-FREE; frozen before execution).

Per cell (run.<cell>.off vs run.<cell>.on, same binary, same inputs, profiles differing by repair_joint_search only):
status and call class, the selected unordered sequence pair (digest of the sorted called sequences, or of the winner's
unordered allele rows when no FASTA is written), block projections (unordered allele pairs per block), the winner group's
score and its fragment / exposure / prior terms, the joint decision / move / continuation / ties / skips, support v2
labels and phase support per block, the whole-caller wall time and peak RSS (ledger), the joint-pass seconds, retained
census bytes and footprint. A missing run or status is reported, never dropped.

Writes <prefix>.matched.tsv (cell, key, off, on, changed), <prefix>.blocks.tsv (per block, per arm) and
<prefix>.selection.sha256 (sha256 of both tables: the frozen predictions).

  compare_m3f.py <runs dir> <ledger.tsv> <out prefix> cell...
"""
import csv
import hashlib
import os
import sys

runs, ledger_path, lambda_path, prefix, cells = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:]
LAM = {r["cell"]: r for r in csv.DictReader(open(lambda_path), delimiter="\t")}
ledger = {r["unit"]: r for r in csv.DictReader(open(ledger_path), delimiter="\t")}
KEYS = ["status", "call_class", "no_call_reason", "search_conditional", "phase_support", "repair_rounds_run", "repair_moves", "repair_stop_reason",
        "repair_converged", "repair_final_incumbent_score", "winner_set_size", "winner_candidate_set_probability", "numeric_separation_gap",
        "fragment_term_gap_to_runner_up", "tie_branched", "tie_ties", "tie_branches_run", "tie_distinct_terminal_genotypes", "tie_incomplete_reason",
        "pool_identities", "block_support_labels_v2x", "phase_competitors_genotype_changing", "phase_competitors_genotype_preserving",
        "phase_competitors_mixed", "joint_search_status", "joint_search_decision", "joint_search_move", "joint_search_single_alternatives",
        "joint_search_catalogue_raw", "joint_search_catalogue_unique", "joint_search_winner_set_size", "joint_search_certified_gap",
        "joint_search_full_placements", "pool_joint_scored", "joint_continued", "joint_leg_rounds", "joint_leg_moves", "joint_leg_stop",
        "joint_leg_converged", "pool_joint_leg_scored", "seconds_total", "seconds_joint_search", "seconds_joint_leg",
        "resources_phys_footprint_peak_bytes", "copying_transition", "copying_length_bp", "locus_geometry_status", "locus_geometry_digest",
        "locus_geometry_frame_digest", "transition_vector_sha256", "transition_r_min", "transition_r_max", "transition_zero_distance_boundaries",
        "transition_zero_distance_informative_boundaries", "model_fingerprint", "repair_anchor_hap1", "repair_anchor_hap2", "pool_union_candidates", "pool_repair_scored"]


def table(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.exists(path) else []


def status(d):
    out = {}
    p = f"{d}/call.mosaic_status.tsv"
    if os.path.exists(p):
        for line in open(p):
            k, _, v = line.rstrip("\n").partition("\t")
            out.setdefault(k, v)
    return out


def winner(d):
    """The winner group's terms and its members' unordered allele rows (from the pooled tables)."""
    scores = [r for r in table(f"{d}/call.mosaic_scores.tsv") if r.get("in_winner_set") == "true"]
    if not scores:
        return {}, "none"
    ids = {r["sequence_candidate_id"] for r in scores}
    rows = set()
    for c in table(f"{d}/call.mosaic_candidates.tsv"):
        if c["sequence_candidate_id"] in ids:
            h1, h2 = c["hap1_alleles"].split(","), c["hap2_alleles"].split(",")
            rows.add(tuple(tuple(sorted((int(a), int(b)))) for a, b in zip(h1, h2)))
    best = max(scores, key=lambda r: float(r["score"]))
    terms = {k: best.get(k, "-") for k in ("score", "fragment_term", "exposure_term", "sequence_pair_prior", "numeric_error_bound")}
    terms["winner_groups"] = str(len(ids))
    return terms, hashlib.sha256(repr(sorted(rows)).encode()).hexdigest()[:16]


def seq_digest(d):
    p = f"{d}/call.mosaic_haplotypes.fa"
    if not os.path.exists(p):
        return "none"
    seqs = sorted("".join(l.strip() for l in chunk.split("\n")[1:]) for chunk in open(p).read().split(">")[1:])
    return hashlib.sha256("|".join(seqs).encode()).hexdigest()[:16]


def blocks(d):
    out = {}
    for r in table(f"{d}/call.block_genotypes.tsv"):
        out[r["block_id"]] = {k: r.get(k, "-") for k in ("allele1", "allele2", "genotype_support", "phase_support", "quality_filter", "genotype_support_v2x", "quality_filter_v2x",
                                                         "quality_reason_v2x", "support_gap_v2x", "genotype_filter")}
    return out


matched, block_rows = [], []
for cell in cells:
    off, on = (f"{runs}/run.{cell}.{arm}" for arm in ("const", "dist"))
    so, sn = status(off), status(on)
    for key in KEYS:
        a, b = so.get(key, "-"), sn.get(key, "-")
        matched.append((cell, key, a, b))
    (wo, ro), (wn, rn) = winner(off), winner(on)
    for k in ("score", "fragment_term", "exposure_term", "sequence_pair_prior", "numeric_error_bound", "winner_groups"):
        matched.append((cell, f"winner_{k}", wo.get(k, "-"), wn.get(k, "-")))
    matched.append((cell, "winner_unordered_rows_digest", ro, rn))
    matched.append((cell, "called_sequence_digest", seq_digest(off), seq_digest(on)))
    bo, bn = blocks(off), blocks(on)
    proj = lambda bs: hashlib.sha256(repr(sorted((k, tuple(sorted((v["allele1"], v["allele2"])))) for k, v in bs.items())).encode()).hexdigest()[:16]
    matched.append((cell, "block_projection_digest", proj(bo) if bo else "none", proj(bn) if bn else "none"))
    lam = LAM.get(cell, {})
    matched.append((cell, "frozen_geometry_digest_matches_dist_run", "-", str(sn.get("locus_geometry_digest") == lam.get("geometry_digest"))))
    matched.append((cell, "frozen_model_fingerprint_matches", str(so.get("model_fingerprint") == lam.get("model_fingerprint")), str(sn.get("model_fingerprint") == lam.get("model_fingerprint"))))
    matched.append((cell, "frozen_lambda_equals_recorded", "-", str(float(sn.get("copying_length_bp", "nan")) == float(lam.get("lambda_bp", "nan")) if lam else "-")))
    for arm in ("const", "dist"):
        u = ledger.get(f"run.{cell}.{arm}", {})
        matched.append((cell, f"ledger_{arm}_wall_seconds|peak_rss_gib|exit|timed_out",
                        "|".join(u.get(k, "-") for k in ("wall_seconds", "peak_rss_gib", "exit", "timed_out")), "-"))
    for b in sorted(set(bo) | set(bn), key=lambda x: (x.split(":")[0], int(x.rsplit(":", 1)[1]))):
        for arm, bs in (("const", bo), ("dist", bn)):
            v = bs.get(b, {})
            block_rows.append((cell, b, arm, *(v.get(k, "-") for k in ("allele1", "allele2", "genotype_support", "phase_support", "quality_filter", "genotype_support_v2x",
                                                                      "quality_filter_v2x", "quality_reason_v2x", "support_gap_v2x", "genotype_filter"))))
with open(f"{prefix}.matched.tsv", "w") as f:
    f.write("cell\tkey\tconst\tdist\tchanged\n")
    for c, k, a, b in matched:
        f.write(f"{c}\t{k}\t{a}\t{b}\t{'-' if b == '-' and k.startswith('ledger') else int(a != b)}\n")
with open(f"{prefix}.blocks.tsv", "w") as f:
    f.write("cell\tblock\tarm\tallele1\tallele2\tgenotype_support\tphase_support\tquality_filter\t"
            "genotype_support_v2x\tquality_filter_v2x\tquality_reason_v2x\tsupport_gap_v2x\tgenotype_filter\n")
    for r in block_rows:
        f.write("\t".join(str(x) for x in r) + "\n")
h = hashlib.sha256(open(f"{prefix}.matched.tsv", "rb").read() + open(f"{prefix}.blocks.tsv", "rb").read()).hexdigest()
open(f"{prefix}.selection.sha256", "w").write(h + "\n")
for cell in cells:
    get = {k: (a, b) for c, k, a, b in matched if c == cell}
    print(cell, "status", get["status"], "call", get["called_sequence_digest"], "rows", get["winner_unordered_rows_digest"],
          "joint", get["joint_search_status"], get["joint_search_decision"], "score", get["winner_score"])
print("selection digest", h)
