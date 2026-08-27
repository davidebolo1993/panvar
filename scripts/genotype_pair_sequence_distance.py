#!/usr/bin/env python3
"""Sequence distance from a called haplotype pair to the sample's true haplotypes.

WHY THIS EXISTS. The non-mosaic ceiling in genotype_pair_ceiling.py maximises the number of blocks
carrying the EXACT truth allele label. A fragment likelihood does not optimise labels; it optimises
how well the candidate sequence explains the reads. Under leave-one-out the sample is off-panel, so
the pair that best explains the reads need not be the pair matching the most labels, and treating the
label ceiling as something the likelihood ought to reach silently assumes they are the same objective.

This measures the other one. If the caller's pair is CLOSER IN SEQUENCE while matching fewer labels,
the label metric is what needs fixing, not the model.

Usage: genotype_pair_sequence_distance.py <truth1.fa> <truth2.fa> <cand1.fa> <cand2.fa> [label]
Prints: label, TOTAL distance under the best assignment, the aligned-NM part, the unaligned part,
per-haplotype totals, per-haplotype aligned bp. The split is reported because a total that is mostly
unaligned bases means something different from one that is mostly mismatches.
"""
import subprocess
import sys
import os
import re

MM2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "external", "minimap2", "minimap2")


def _uncovered(length, intervals):
    """Bases of a sequence no alignment covers, from merged half-open intervals."""
    covered = 0
    last = -1
    for a, b in sorted(intervals):
        a = max(a, last)
        if b > a:
            covered += b - a
            last = b
    return max(0, length - covered)


def distance(ref, qry):
    """Edit distance between two haplotype FASTAs that PENALISES UNALIGNED SEQUENCE.

        D = NM over aligned blocks + unaligned truth bp + unaligned called bp

    Summing NM over minimap2's reported alignments alone is not safe: bases that never aligned are
    simply absent from the total, so a candidate missing 10-20 kb can score BETTER than one that
    represents the whole locus, purely because the missing part was never compared. Every base of
    both sequences has to be accounted for, either as aligned (and then charged its edits) or as
    unaligned (and then charged in full).

    Returns (total, aligned_bp, nm, unaligned_truth, unaligned_called).
    """
    out = subprocess.run([MM2, "-cx", "asm20", "--secondary=no", ref, qry],
                         capture_output=True, text=True)
    nm = blocks = 0
    qlen = tlen = 0
    qiv, tiv = [], []
    for line in out.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        m = re.search(r"NM:i:(\d+)", line)
        if m is None:
            continue
        nm += int(m.group(1))
        blocks += int(f[10])
        qlen = int(f[1])
        tlen = int(f[6])
        qiv.append((int(f[2]), int(f[3])))
        tiv.append((int(f[7]), int(f[8])))
    if qlen == 0 or tlen == 0:
        # Nothing aligned at all: charge both sequences in full rather than reporting 0 edits.
        qn = _fasta_len(qry)
        tn = _fasta_len(ref)
        return qn + tn, 0, 0, tn, qn
    uq = _uncovered(qlen, qiv)
    ut = _uncovered(tlen, tiv)
    return nm + uq + ut, blocks, nm, ut, uq


def _fasta_len(path):
    n = 0
    for line in open(path):
        if line[0] != ">":
            n += len(line.strip())
    return n


def main():
    t1, t2, c1, c2 = sys.argv[1:5]
    label = sys.argv[5] if len(sys.argv) > 5 else "pair"
    # Both assignments of the two candidates to the two truth haplotypes; the caller does not phase,
    # so scoring only one orientation would penalise a correct call for the order it was written in.
    a = [distance(t1, c1), distance(t2, c2)]
    b = [distance(t1, c2), distance(t2, c1)]
    # Select the homologue assignment on the FULL distance, not on aligned edits. Selecting on
    # aligned edits alone lets an assignment win by aligning less of the sequence.
    best = a if (a[0][0] + a[1][0]) <= (b[0][0] + b[1][0]) else b
    total = best[0][0] + best[1][0]
    nm = best[0][2] + best[1][2]
    unal = best[0][3] + best[0][4] + best[1][3] + best[1][4]
    print(f"{label}\t{total}\t{nm}\t{unal}\t{best[0][0]}\t{best[1][0]}"
          f"\t{best[0][1]}\t{best[1][1]}")


if __name__ == "__main__":
    main()
