#!/usr/bin/env python3
"""The best any pair of COMPLETE panel haplotypes could do, given this panel and this truth.

This is the ceiling of a non-mosaic model. It is deliberately NOT the panel's ceiling: a mosaic
model may take a different haplotype at every block, and under that model every block whose truth
allele exists in the panel is reachable. The gap between the two is exactly what a mosaic layer is
worth, so comparing a whole-haplotype genotyper against THIS number rather than against the panel
says whether its evidence model is losing or whether the missing mosaic layer is.

Usage: genotype_pair_ceiling.py <alleles.tsv> <genotypes.tsv> <locus> <loo> <donor>
Writes the same per-block rows the A/B harness collects, with arm=ceiling.
"""
import sys
import collections


def main():
    al, gt, locus, loo, donor = sys.argv[1:6]
    hap = collections.defaultdict(dict)
    with open(al) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            h, b, _kind, a = line.rstrip("\n").split("\t")
            hap[h][int(b)] = int(a)

    rows = [line.rstrip("\n").split("\t") for line in open(gt)]
    col = {n: i for i, n in enumerate(rows[0])}
    truth, kind, nall = {}, {}, {}
    for r in rows[1:]:
        b = int(r[col["block_index"]])
        kind[b] = r[col["block_kind"]]
        nall[b] = r[col["n_alleles"]]
        t1, t2 = int(r[col["truth1"]]), int(r[col["truth2"]])
        if t1 >= 0 and t2 >= 0:
            truth[b] = (min(t1, t2), max(t1, t2))

    def agree(x, y, b, t):
        if b not in hap[x] or b not in hap[y]:
            return False
        return (min(hap[x][b], hap[y][b]), max(hap[x][b], hap[y][b])) == t

    names = list(hap)
    best, best_n = (names[0], names[0]) if names else (None, None), -1
    for i, x in enumerate(names):
        for y in names[i:]:
            n = sum(1 for b, t in truth.items() if agree(x, y, b, t))
            if n > best_n:
                best_n, best = n, (x, y)
    if best[0] is None:
        return
    x, y = best

    out = []
    for b in sorted(kind):
        if b not in truth:
            out.append(f"{locus}\t{loo}\t{donor}\tceiling\t{b}\t{kind[b]}\t{nall[b]}\t0\tNA")
        else:
            ok = 1 if agree(x, y, b, truth[b]) else 0
            out.append(f"{locus}\t{loo}\t{donor}\tceiling\t{b}\t{kind[b]}\t{nall[b]}\t1\t{ok}")
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
