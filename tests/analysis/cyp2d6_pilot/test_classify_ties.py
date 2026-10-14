#!/usr/bin/env python3
"""Synthetic checks of classify_ties.py: the SEQUENCE_IDENTICAL_BLOCK_LABEL_AMBIGUITY label is given only to a tie whose
members spell one sequence with identical scores and differ only inside a degeneracy interval.

Model: three blocks A | B | C; block A alleles "AC", "ACGT"; block B "GTTT", "TT", "CCCC"; block C "G".
(0,0,0) spells AC+GTTT+G and (1,1,0) spells ACGT+TT+G: the same sequence, a cross-block label ambiguity over A..B.
(0,2,0) spells a distinct sequence. Usage: test_classify_ties.py <classify_ties.py>
"""
import os
import shutil
import subprocess
import sys
import tempfile

HEAD = "rank\tsequence_candidate_id\tcandidate_id\trepresentative\thap1_alleles\thap2_alleles\tsource\tscored_in\tscore\t" \
       "fragment_term\texposure_term\tsequence_pair_prior\tnumeric_error_bound\texact_path_prior\n"


def member(seq_id, cand, h1, h2, score="-10.5", frag="-8.0", expo="-2.0", prior="-0.5"):
    return f"1\t{seq_id}\t{cand}\tfalse\t{h1}\t{h2}\tpanel\tunion\t{score}\t{frag}\t{expo}\t{prior}\t0\t-1\n"


CASES = {
    "label_ambiguity": ([member(7, 1, "0,0,0", "0,0,0"), member(7, 2, "1,1,0", "0,0,0")],
                        "SEQUENCE_IDENTICAL_BLOCK_LABEL_AMBIGUITY"),
    "label_ambiguity_swapped_haplotypes": ([member(7, 1, "0,0,0", "0,2,0"), member(7, 2, "0,2,0", "1,1,0")],
                                           "SEQUENCE_IDENTICAL_BLOCK_LABEL_AMBIGUITY"),
    "distinct_sequences": ([member(7, 1, "0,0,0", "0,0,0"), member(8, 2, "0,2,0", "0,0,0")], "DISTINCT_SEQUENCE_TIE"),
    "same_id_distinct_spelling": ([member(7, 1, "0,0,0", "0,0,0"), member(7, 2, "0,2,0", "0,0,0")],
                                  "DISTINCT_SEQUENCE_TIE"),
    "score_differs": ([member(7, 1, "0,0,0", "0,0,0"), member(7, 2, "1,1,0", "0,0,0", frag="-8.1")],
                      "DISTINCT_SEQUENCE_TIE"),
    "difference_outside_interval": ([member(7, 1, "0,0,0", "0,0,0"), member(7, 2, "1,1,0", "0,0,0")],
                                    "DISTINCT_SEQUENCE_TIE"),
}


def build(root, name, members, interval):
    os.makedirs(os.path.join(root, "degeneracy", "models"), exist_ok=True)
    os.makedirs(os.path.join(root, "data", "runs_t"), exist_ok=True)
    with open(os.path.join(root, "degeneracy", "models", "original.block_alleles.tsv"), "w") as f:
        f.write("block_index\tblock_id\tkind\tallele_index\tlength\tsequence\n")
        for b, (bid, alleles) in enumerate((("A", ["AC", "ACGT"]), ("B", ["GTTT", "TT", "CCCC"]), ("C", ["G"]))):
            for a, s in enumerate(alleles):
                f.write(f"{b}\t{bid}\tSV\t{a}\t{len(s)}\t{s}\n")
    with open(os.path.join(root, "degeneracy", "degeneracy_summary.tsv"), "w") as f:
        f.write(f"model\tcoalescence_intervals\noriginal\t{interval}\n")
    with open(os.path.join(root, "jobs.tsv"), "a") as f:
        if f.tell() == 0:
            f.write("job\tarm\tsample\tcoverage\n")
        f.write(f"{name}\tlzo\tX\t30\n")
    prefix = os.path.join(root, "data", "runs_t", name)
    open(prefix + ".run.tsv", "w").write("job\texit\n" + f"{name}\t1\n")
    open(prefix + ".mosaic_status.tsv", "w").write(
        "key\tvalue\nstatus\tRESCORE_INCOMPLETE\ncall_class\tno_call\nno_call_reason\tdistinct_tie\nseconds_total\t1\n")
    open(prefix + ".equivalence_set.tsv", "w").write(HEAD + "".join(members))


def main():
    script = sys.argv[1]
    failures = 0
    for name, (members, expected) in CASES.items():
        root = tempfile.mkdtemp(prefix="classify_ties_")
        try:
            build(root, name, members, "C..C" if name == "difference_outside_interval" else "A..B")
            out = os.path.join(root, "out.tsv")
            subprocess.run([sys.executable, "-I", script, root, "runs_t", out, name], check=True, capture_output=True)
            got = open(out).read().splitlines()[1].split("\t")[-1]
        finally:
            shutil.rmtree(root)
        ok = got == expected
        failures += not ok
        print(("PASS " if ok else "FAIL ") + f"{name}: {got} (expected {expected})")
    print(f"{len(CASES)} cases, {failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
