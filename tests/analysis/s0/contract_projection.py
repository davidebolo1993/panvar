#!/usr/bin/env python3
"""Reporting-only projection of the sealed CYP2D6 pilot cells onto the sequence-first contract classes
(tests/genotype_sequence_label_reporting_design.md section 6). Reads sealed tables; writes nothing into the pilot.
Whole-sequence correctness columns are the sealed evaluator's, POST HOC: they must not set any confidence or class."""
import csv, collections, sys

PILOT = "/Users/davide.bolognini/panvar/tests/results/cyp2d6_shortread_pilot_v1/tables/whole_sequence.tsv"
CLASS = {  # current reporting classification -> (sequence_status, block reporting, prediction under the contract)
    "CALL_strict_call": ("RESOLVED", "all blocks", "retained (unchanged)"),
    "CALL_genotype_equivalent_call": ("GENOTYPE_RESOLVED_SEQUENCE_SET", "all blocks; phase sets", "representative pair, sequence confidence below block confidence"),
    "SEQUENCE_IDENTICAL_BLOCK_LABEL_AMBIGUITY": ("RESOLVED", "all but sv:3,sv:5; interval sv:3..sv:5 row", "NEW: sequence pair reported (was no-call)"),
    "DISTINCT_SEQUENCE_TIE": ("SEQUENCE_TIED", "blocks agreeing across the tie", "NEW: representative pair, low confidence (was no-call)"),
}


def main(out):
    rows = list(csv.DictReader(open(PILOT), delimiter="\t"))
    table = collections.Counter()
    with open(out, "w") as f:
        f.write("job\tcurrent\tsequence_status\tblocks_reported\tprediction\tpost_hoc_whole_sequence\n")
        for r in rows:
            c = r["classification"]
            key = c if c in CLASS else ("DISTINCT_SEQUENCE_TIE" if "TIE" in c.upper() else c)
            status, blocks, pred = CLASS.get(key, ("UNMAPPED", "?", "?"))
            f.write(f'{r["job"]}\t{c}\t{status}\t{blocks}\t{pred}\t{r["whole_sequence"]}\n')
            table[(c, status, r["whole_sequence"].split(" ")[0])] += 1
    for k, v in sorted(table.items()):
        print(*k, v, sep="\t")


if __name__ == "__main__":
    main(sys.argv[1])
