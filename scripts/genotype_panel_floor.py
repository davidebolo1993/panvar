#!/usr/bin/env python3
"""The best sequence reconstruction ANY pair of panel haplotypes could give, and the excess over it.

    E*_panel = min over pairs (a,b) of  min[ d(a,T1)+d(b,T2), d(a,T2)+d(b,T1) ]
             = min_h d(h,T1)  +  min_h d(h,T2)

The two terms are independent -- a and b range over the same panel and may coincide -- so this is one
pass over the panel, not a pair search: 2 alignments of the whole panel rather than O(n^2) pairs.

This is the number that separates PANEL LIMITATION from INFERENCE FAILURE, which is what the
exact-label ceiling could not do. E_excess = E_call - E*_panel is the part a better caller could still
recover; E*_panel itself is the part only a better panel could.

Truth is always the two held-out assembled haplotypes. It is deliberately NOT the panel pair sharing
the most markers: scoring against that would reward a caller for reproducing its own representation's
biases and would make discarded information invisible.

Usage: genotype_panel_floor.py <truth1.fa> <truth2.fa> <panel_all.fa> [label]
"""
import os
import re
import subprocess
import sys

MM2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "external", "minimap2", "minimap2")


def _uncovered(length, intervals):
    covered, last = 0, -1
    for a, b in sorted(intervals):
        a = max(a, last)
        if b > a:
            covered += b - a
            last = b
    return max(0, length - covered)


def per_query_distance(truth_fa, panel_fa):
    """distance(h, T) for every haplotype in panel_fa, in ONE minimap2 pass.

    Same accounting as the pair oracle: aligned NM plus every unaligned base of both sequences, so a
    haplotype cannot look good by aligning only part of itself.
    """
    out = subprocess.run([MM2, "-cx", "asm20", "--secondary=no", truth_fa, panel_fa],
                         capture_output=True, text=True)
    nm, qiv, tiv, qlen = {}, {}, {}, {}
    tlen = 0
    for line in out.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        m = re.search(r"NM:i:(\d+)", line)
        if m is None:
            continue
        q = f[0]
        nm[q] = nm.get(q, 0) + int(m.group(1))
        qlen[q] = int(f[1])
        tlen = int(f[6])
        qiv.setdefault(q, []).append((int(f[2]), int(f[3])))
        tiv.setdefault(q, []).append((int(f[7]), int(f[8])))
    return {q: nm[q] + _uncovered(qlen[q], qiv[q]) + _uncovered(tlen, tiv[q]) for q in nm}


def exact_distance(binary, fa_a, fa_b):
    """True global (Needleman-Wunsch) edit distance, via `panvar genotype-frag --exact-distance`."""
    out = subprocess.run([binary, "genotype-frag", "--exact-distance", fa_a, fa_b],
                         capture_output=True, text=True)
    try:
        return int(out.stdout.strip())
    except ValueError:
        return None


def write_one(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n{seq}\n")


def load_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name is not None:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    t1, t2, panel = sys.argv[1:4]
    label = sys.argv[4] if len(sys.argv) > 4 else "panel_floor"
    d1 = per_query_distance(t1, panel)
    d2 = per_query_distance(t2, panel)
    if not d1 or not d2:
        print(f"{label}\tNA\tNA\tNA\tNA\tNA")
        return
    b1 = min(d1, key=d1.get)
    b2 = min(d2, key=d2.get)

    # EXACT mode. The approximate floor above is alignment-derived: minimap2 decides what to align and
    # anything it declines is charged as whole unaligned bases. Safe as a bound, but it is not a
    # distance, and at an array it is the number a model would be judged against. So shortlist the
    # nearest few candidates approximately, compute a true global edit distance for those, and report
    # both -- the difference is the approximation's error, measured rather than assumed.
    if os.environ.get("EXACT_BINARY"):
        binary = os.environ["EXACT_BINARY"]
        topk = int(os.environ.get("EXACT_TOPK", "5"))
        tmp = os.environ.get("EXACT_TMP", "/tmp")
        panel_seqs = load_fasta(panel)
        best = []
        for idx, (truth_fa, d) in enumerate(((t1, d1), (t2, d2))):
            exact = {}
            for c in sorted(d, key=d.get)[:topk]:
                if c not in panel_seqs:
                    continue
                cf = os.path.join(tmp, f"exact_cand_{idx}.fa")
                write_one(cf, c, panel_seqs[c])
                e = exact_distance(binary, truth_fa, cf)
                if e is not None:
                    exact[c] = e
            if exact:
                bc = min(exact, key=exact.get)
                best.append((exact[bc], bc))
            else:
                best.append((None, None))
        if all(b[0] is not None for b in best):
            print(f"{label}\t{d1[b1] + d2[b2]}\t{best[0][0] + best[1][0]}"
                  f"\t{d1[b1]}\t{d2[b2]}\t{best[0][0]}\t{best[1][0]}"
                  f"\t{b1}\t{b2}\t{best[0][1]}\t{best[1][1]}")
            return

    print(f"{label}\t{d1[b1] + d2[b2]}\t{d1[b1]}\t{d2[b2]}\t{b1}\t{b2}")


if __name__ == "__main__":
    main()
