#!/usr/bin/env python3
"""M6-A summary (after the freeze; truth and structure enter here only).

Per cell, over EVERY candidate of the count ranking: the structural outcome of each ordered allele pair from units.py
(total content = units(h1) + units(h2) vs truth total; allocation = the unordered unit pair vs truth's), the best count
group's composition, the first group holding the truth's unit pair, whether the truth's unit pair exists in the catalogue
at all (missing-content check), and the caller's prediction. For the fragment-scored subset (select_eval.py): F, E,
prior, total and the sequence excess at the target block.
Frozen verdicts per cell:
  count-best structure  ALL_CORRECT (every best-group member has the truth unit pair) | TOTAL_ONLY (correct total, an
                        allocation differs or is mixed) | WRONG_TOTAL | MIXED_TOTAL | NOT_ASSESSED (LOW_INFORMATION)
  allocation resolved   false if the best group spans more than one unit allocation of the same total
  content available     whether any catalogue candidate carries the truth unit pair
  vs caller             IMPROVES (count-best ALL_CORRECT where the caller's structure is wrong) | PRESERVES (both correct)
                        | WORSENS (caller correct, count-best not ALL_CORRECT) | SAME_ERROR | DIFFERENT_ERROR

  summary.py <runs dir> <ledger> <out prefix> cell:block ...
"""
import csv
import sys

runs, ledger, prefix, specs = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l}
led = {r["unit"]: r for r in T(ledger)}
out = ["cell\tdepth\tcandidates\tgroups\tbest_group_size\tbest_structure\tallocation_resolved\ttruth_units\ttruth_pair_in_catalogue\tfirst_group_with_truth_pair\t"
       "caller_units\tcaller_structure\tcaller_group_rank\tcaller_loss_minus_best\tverdict_vs_caller\tcountfit_seconds\tcountfit_peak_footprint_gib"]
detail = ["cell\tlabel\tunits\tstructure\tcount_group\tloss_minus_best\tfragment\texposure\tprior\ttotal\ttotal_minus_caller\tblock_excess"]
for spec in specs:
    cell, block = spec.split(":")
    s = KV(f"{runs}/countfit.{cell}/cf.summary.tsv")
    cands = T(f"{runs}/countfit.{cell}/cf.candidates.tsv")
    units = {r["name"]: r["unit_copies"] for r in T(f"{runs}/units.{cell}/units.tsv")}
    tu = sorted(int(units[k]) for k in units if k.startswith("truth")) if all(units[k] != "LOW_INFORMATION" for k in units if k.startswith("truth")) else None
    U = lambda a: units.get(f"allele_{a}", "NA")
    def structure(a, b):
        ua, ub = U(a), U(b)
        if tu is None or "LOW_INFORMATION" in (ua, ub) or "NA" in (ua, ub):
            return "NOT_ASSESSED"
        p = sorted((int(ua), int(ub)))
        if p == tu:
            return "CORRECT"
        return "ALLOCATION_ERROR" if sum(p) == sum(tu) else "DOSAGE_ERROR"
    best = [r for r in cands if r["group_rank"] == "1"]
    bs = {structure(r["allele_h1"], r["allele_h2"]) for r in best}
    alloc = {tuple(sorted((U(r["allele_h1"]), U(r["allele_h2"])))) for r in best}
    best_struct = ("ALL_CORRECT" if bs == {"CORRECT"} else "NOT_ASSESSED" if "NOT_ASSESSED" in bs else
                   "TOTAL_ONLY" if bs <= {"CORRECT", "ALLOCATION_ERROR"} else "WRONG_TOTAL" if bs == {"DOSAGE_ERROR"} else "MIXED_TOTAL")
    resolved = len(alloc) == 1
    first_truth = next((r["group_rank"] for r in cands if structure(r["allele_h1"], r["allele_h2"]) == "CORRECT"), "absent")
    caller = next(r for r in cands if r["is_caller"] == "true")
    cs = structure(caller["allele_h1"], caller["allele_h2"])
    v = ("PRESERVES" if cs == "CORRECT" and best_struct == "ALL_CORRECT" else "WORSENS" if cs == "CORRECT" else
         "IMPROVES" if best_struct == "ALL_CORRECT" else "SAME_ERROR" if bs == {cs} else "DIFFERENT_ERROR")
    u = led.get(f"countfit.{cell}", {})
    out.append("\t".join(map(str, [cell, s["depth_haploid_kmer_mean"], len(cands), s["groups"], s["best_group_size"], best_struct, resolved, tu,
                                   first_truth != "absent", first_truth, f"{U(caller['allele_h1'])}+{U(caller['allele_h2'])}", cs, s["caller_group_rank"],
                                   s["caller_loss_minus_best"], v, s["seconds_total"], f"{int(s['phys_footprint_peak_bytes']) / 2**30:.2f}"])))
    sc = {r["label"]: r for r in T(f"{runs}/fragscore.{cell}/fs.scores.tsv")}
    blocks = {}
    for r in T(f"{runs}/fragscore.{cell}/fs.blocks.tsv"):
        if r["block_id"] == block:
            blocks[r["label"]] = r["block_excess"]
    ctot = float(next(r for k, r in sc.items() if "caller" in k)["total"])
    by_pair = {(r["allele_h1"], r["allele_h2"]): r for r in cands}
    for label, r in sc.items():
        tb = int(s["target_block_index"])
        a, b = r["hap1_alleles"].split(",")[tb], r["hap2_alleles"].split(",")[tb]
        cr = by_pair[(a, b)]
        detail.append("\t".join(map(str, [cell, label, f"{U(a)}+{U(b)}", structure(a, b), cr["group_rank"], cr["loss_minus_best"], r["fragment_term"],
                                          r["exposure_term"], r["prior"], r["total"], f"{float(r['total']) - ctot:.4f}", blocks.get(label, "NA")])))
open(prefix + ".tsv", "w").write("\n".join(out) + "\n\n" + "\n".join(detail) + "\n")
print("\n".join(out))
print()
print("\n".join(detail))
