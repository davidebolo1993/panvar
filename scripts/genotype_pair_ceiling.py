#!/usr/bin/env python3
"""The best any pair of COMPLETE panel haplotypes could do, given this panel and this truth.

This is the ceiling of a non-mosaic model. It is deliberately NOT the panel's ceiling: a mosaic
model may take a different haplotype at every block, and under that model every block whose truth
allele exists in the panel is reachable. The gap between the two is exactly what a mosaic layer is
worth, so comparing a whole-haplotype genotyper against THIS number rather than against the panel
says whether its evidence model is losing or whether the missing mosaic layer is.

LIMIT, stated because the number is easy to over-read: this maximises the count of blocks carrying
the EXACT truth allele label. It does not minimise sequence edit distance, and it ignores blocks whose
exact label is absent from the reduced panel even where a near-identical allele exists. So it is an
exact-label, single-pair ceiling and NOT a certified sequence-reconstruction ceiling. The release
goal is stated in sequence distance, so this number does not stand in for it.

Usage: genotype_pair_ceiling.py <alleles.tsv> <genotypes.tsv> <locus> <loo> <donor> [optimum.tsv]
Writes the same per-block rows the A/B harness collects, with arm=ceiling. With a sixth argument,
also writes the optimal pair, the number of tied optima and per-block reachability across all of them.
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
    if not names:
        return
    # EVERY optimum, not the first one found. Which blocks an optimal pair gets right differs between
    # tied optima, so reporting the profile of one of them silently picks a winner among equals -- and
    # the count of ties is itself the answer to "how identifiable is this at all".
    best_n = -1
    optima = []
    for i, x in enumerate(names):
        for y in names[i:]:
            n = sum(1 for b, t in truth.items() if agree(x, y, b, t))
            if n > best_n:
                best_n, optima = n, [(x, y)]
            elif n == best_n:
                optima.append((x, y))

    # The arm's own per-block profile comes from ONE optimum (the first), so it stays a real pair a
    # real caller could have chosen. Reachability across ALL optima is reported separately, because a
    # block reachable by some optimum but not this one is a different kind of miss.
    x, y = optima[0]
    reach = {b: any(agree(u, v, b, t) for u, v in optima) for b, t in truth.items()}

    out = []
    for b in sorted(kind):
        if b not in truth:
            out.append(f"{locus}\t{loo}\t{donor}\tceiling\t{b}\t{kind[b]}\t{nall[b]}\t0\tNA")
        else:
            ok = 1 if agree(x, y, b, truth[b]) else 0
            out.append(f"{locus}\t{loo}\t{donor}\tceiling\t{b}\t{kind[b]}\t{nall[b]}\t1\t{ok}")
    sys.stdout.write("\n".join(out) + "\n")

    if len(sys.argv) > 6:
        # Side channel for the gap decomposition: which pair to probe, and how identifiable the
        # optimum is. Written to a file rather than stdout so the row stream stays machine-readable.
        # EVERY optimum, not just the first. Probing one of several tied optima and reporting its
        # rank understates the likelihood whenever a different optimum ranks better -- the same shape
        # of selective comparison this project has had to retract before.
        with open(sys.argv[6], "w") as fh:
            fh.write("donor\thap1\thap2\tn_optima\tblocks_exact\tblocks_representable\tblocks_reachable_any_optimum\n")
            for u, v in optima:
                fh.write(f"{donor}\t{u}\t{v}\t{len(optima)}\t{best_n}\t{len(truth)}\t{sum(reach.values())}\n")


if __name__ == "__main__":
    main()
