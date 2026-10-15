#!/usr/bin/env python3
"""CYP2D6/D7 pilot reporting addendum 1: corrected denominators and callability beside conditional accuracy.

Reads only sealed package files (tables/, call_classification_a5.tsv, data/eval_a5/); writes <out dir>/*.tsv.
  denominators.tsv          per arm x coverage: block opportunities (cells x 19), block genotypes reported by calls,
                            truth-undefined among them, scored (= reported - truth-undefined), retained of reported and
                            of scored
  callability_accuracy.tsv  callability beside conditional and unconditional block outcomes and whole-sequence outcomes:
                            no-calls stay in the unconditional denominators
  recall_basis.tsv          the filter "recall": retained at-floor / at-floor among scored blocks of calls, with how many
                            of those at-floor blocks have a positive floor (at floor is not the same as correct)
  loo_spellability.tsv      per LOO cell: whether the truth diploid sequence is spellable from the LOO model's alleles at
                            every block (all availability floors 0, no truth-undefined block), i.e. whether an exact
                            whole-sequence call was attainable at all
Usage: addendum1.py <pilot package> <out dir>
"""
import csv
import os
import sys

PKG, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
BLOCKS = 19
ARMS = ["lzo", "mosaic", "loo"]


def table(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(name, rows):
    with open(os.path.join(OUT, name), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


cells = table(os.path.join(PKG, "tables", "cells.tsv"))
obs = [o for o in table(os.path.join(PKG, "tables", "block_observations.tsv")) if o["layer"] == "reported"]
und = [u for u in table(os.path.join(PKG, "tables", "truth_undefined_blocks.tsv")) if u["layer"] == "reported"]
wseq = {w["job"]: w for w in table(os.path.join(PKG, "tables", "whole_sequence.tsv"))}

den, ca, rb = [], [], []
for arm in ARMS:
    for cov in ("5", "30"):
        cs = [c for c in cells if c["arm"] == arm and c["coverage"] == cov]
        calls = [c for c in cs if c["exit"] == "0"]
        o = [x for x in obs if x["arm"] == arm and x["coverage"] == cov]
        u = [x for x in und if x["arm"] == arm and x["coverage"] == cov]
        reported = len(calls) * BLOCKS
        assert reported == len(o) + len(u), (arm, cov, reported, len(o), len(u))
        ret_scored = sum(x["genotype_filter"] == "PASS" for x in o)
        ret_und = sum(x["genotype_filter"] == "PASS" for x in u)
        at_floor = sum(float(x["excess"]) == 0 for x in o)
        den.append({"arm": arm, "coverage": cov, "cells": len(cs), "block_opportunities": len(cs) * BLOCKS,
                    "calls": len(calls), "block_genotypes_reported": reported, "truth_undefined": len(u),
                    "scored": len(o), "retained_of_reported": f"{ret_scored + ret_und}/{reported}",
                    "retained_of_scored": f"{ret_scored}/{len(o)}"})
        klass = [c["classification"] for c in cs]
        exact_seq = sum(wseq[c["job"]]["whole_sequence"].startswith("exact") and c["exit"] == "0" for c in cs)
        identified = sum(wseq[c["job"]]["whole_sequence"].startswith("exact (identified") for c in cs)
        ca.append({"arm": arm, "coverage": cov, "cells": len(cs),
                   "calls_strict": klass.count("CALL_strict_call"),
                   "calls_genotype_equivalent": klass.count("CALL_genotype_equivalent_call"),
                   "no_call_label_ambiguity": klass.count("SEQUENCE_IDENTICAL_BLOCK_LABEL_AMBIGUITY"),
                   "no_call_distinct_sequence_tie": klass.count("DISTINCT_SEQUENCE_TIE"),
                   "callability": f"{len(calls)}/{len(cs)}",
                   "blocks_at_floor_conditional_on_scored": f"{at_floor}/{len(o)}",
                   "blocks_at_floor_of_all_opportunities": f"{at_floor}/{len(cs) * BLOCKS}",
                   "strict_call_exact_sequence_of_all_cells": f"{exact_seq}/{len(cs)}",
                   "no_call_identified_sequence_exact_of_all_cells": f"{identified}/{len(cs)}"})
        ok = [x for x in o if float(x["excess"]) == 0]
        rb.append({"arm": arm, "coverage": cov, "scored_block_genotypes_of_calls": len(o),
                   "at_floor": len(ok), "at_floor_with_positive_floor": sum(float(x["floor"]) > 0 for x in ok),
                   "retained_at_floor": sum(x["genotype_filter"] == "PASS" for x in ok),
                   "recall_as_reported": f"{sum(x['genotype_filter'] == 'PASS' for x in ok)}/{len(ok)}",
                   "excluded_any": sum(x["genotype_filter"] != "PASS" for x in o),
                   "not_in_basis_no_call_cells": len(cs) - len(calls), "not_in_basis_truth_undefined": len(u)})
write("denominators.tsv", den)
write("callability_accuracy.tsv", ca)
write("recall_basis.tsv", rb)

loo = []
for c in [c for c in cells if c["arm"] == "loo"]:
    b = table(os.path.join(PKG, "data", "eval_a5", f"{c['job']}.reported.blocks.tsv"))
    undefined = [r["block_id"] for r in b if r["truth1_bp"] == "NA" or r["truth2_bp"] == "NA"]
    positive = [r["block_id"] for r in b if r["availability_floor"] not in ("NA", "0", "0.0")]
    loo.append({"job": c["job"], "classification": c["classification"],
                "truth_undefined_blocks": ",".join(undefined) or "none",
                "blocks_with_positive_floor": ",".join(positive) or "none",
                "truth_spellable_from_loo_model": "yes" if not undefined and not positive else "no",
                "whole_sequence": wseq[c["job"]]["whole_sequence"]})
write("loo_spellability.tsv", loo)
print("written", sorted(os.listdir(OUT)))
