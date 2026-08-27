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
Prints: label, total EDITS under the best assignment, total divergence, per-haplotype edits,
per-haplotype aligned bp. Edits are the primary figure: they are integers and they are what the
release goal is stated in.
"""
import subprocess
import sys
import os
import re

MM2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "external", "minimap2", "minimap2")


def divergence(ref, qry):
    """Aligned-base divergence between two haplotype FASTAs.

    From minimap2's NM:i: (an INTEGER edit count) over total aligned block length, not from de:f:.
    de:f: is printed to four decimals, which at the divergences seen here (2e-4 to 6e-4) quantises
    the answer to one or two steps and cannot separate two candidates that differ by a few dozen
    edits over 213 kb. Returns (divergence, aligned_bp, edits); aligned_bp matters because a low
    divergence over a small aligned fraction is not a good reconstruction and must not read as one.
    """
    out = subprocess.run([MM2, "-cx", "asm20", "--secondary=no", ref, qry],
                         capture_output=True, text=True)
    edits = den = 0
    for line in out.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 11:
            continue
        m = re.search(r"NM:i:(\d+)", line)
        if m is None:
            continue
        edits += int(m.group(1))
        den += int(f[10])
    return (edits / den if den else 1.0), den, edits


def main():
    t1, t2, c1, c2 = sys.argv[1:5]
    label = sys.argv[5] if len(sys.argv) > 5 else "pair"
    # Both assignments of the two candidates to the two truth haplotypes; the caller does not phase,
    # so scoring only one orientation would penalise a correct call for the order it was written in.
    a = [divergence(t1, c1), divergence(t2, c2)]
    b = [divergence(t1, c2), divergence(t2, c1)]
    best = a if (a[0][2] + a[1][2]) <= (b[0][2] + b[1][2]) else b
    total_edits = best[0][2] + best[1][2]
    total_div = best[0][0] + best[1][0]
    print(f"{label}\t{total_edits}\t{total_div:.6f}\t{best[0][2]}\t{best[1][2]}"
          f"\t{best[0][1]}\t{best[1][1]}")


if __name__ == "__main__":
    main()
