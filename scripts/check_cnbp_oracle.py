#!/usr/bin/env python3
"""Assert FORMAT:CNBP against the allele VCF, which is an exact oracle for it.

CNBP claims to be the signed net bases a haplotype carries across a copy-number module relative to the
reference. The allele VCF spells every allele of every bubble explicitly, so for a haplotype whose GT
names allele k at that bubble, `len(allele k) - len(REF allele)` IS that quantity. The two are computed
by different code -- one sums node lengths over the widest source..sink span, the other spells the
walk -- so agreeing to the base is a real check on the span logic and on repeated-node multiplicity,
which are the two places this is easy to get subtly wrong.

It is an internal check, not a measure of biological truth: both sides read the same graph walks. Use
compare_svim_asm.py for the external comparison.

  check_cnbp_oracle.py <region.vcf> <alleles.vcf> [--verbose]

Exits non-zero if any comparison disagrees.
"""
import collections
import gzip
import sys


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_vcf(path):
    hdr, rows = None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
            elif not line.startswith("#"):
                rows.append(line.rstrip("\n").split("\t"))
    return hdr, rows


def info_of(field):
    d = {}
    for kv in field.split(";"):
        k, sep, v = kv.partition("=")
        d[k] = v if sep else True
    return d


def bubble_of(rec_id):
    return rec_id.split("_", 1)[0]          # bubble<N>_<TYPE>_<node>


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    verbose = "--verbose" in sys.argv
    if len(args) != 2:
        sys.exit(__doc__)
    rhdr, rrows = read_vcf(args[0])
    ahdr, arows = read_vcf(args[1])
    if rhdr is None or ahdr is None:
        sys.exit("check_cnbp_oracle: a VCF has no #CHROM header line")

    allele_rec = {bubble_of(r[2]): r for r in arows}

    n_dup = n_cmp = n_exact = n_unmatched = 0
    err = collections.Counter()
    examples = []
    for r in rrows:
        if info_of(r[7]).get("SVTYPE") != "DUP":
            continue
        n_dup += 1
        a = allele_rec.get(bubble_of(r[2]))
        if a is None:
            n_unmatched += 1
            continue
        lens = [len(a[3])] + [len(x) for x in a[4].split(",")]
        agt = a[8].split(":").index("GT")
        hap_allele = {ahdr[i]: a[i].split(":")[agt] for i in range(9, len(ahdr))}
        rkeys = r[8].split(":")
        if "CNBP" not in rkeys:
            continue
        ci = rkeys.index("CNBP")
        for i in range(9, len(rhdr)):
            v = r[i].split(":")
            if len(v) <= ci or v[ci] == ".":
                continue
            gt = hap_allele.get(rhdr[i])
            if gt is None or not gt.isdigit() or int(gt) >= len(lens):
                continue
            truth = lens[int(gt)] - lens[0]
            diff = int(v[ci]) - truth
            n_cmp += 1
            err[diff] += 1
            if diff == 0:
                n_exact += 1
            elif len(examples) < 8:
                examples.append((r[2], rhdr[i], int(v[ci]), truth, diff))

    pct = 100.0 * n_exact / n_cmp if n_cmp else 0.0
    print(f"DUP records {n_dup}  (allele record missing for {n_unmatched})")
    print(f"comparisons {n_cmp}  exact {n_exact} ({pct:.2f}%)")
    if verbose or n_exact != n_cmp:
        top = sorted(err.items(), key=lambda kv: -kv[1])[:8]
        print("CNBP minus allele delta:", ", ".join(f"{d}:{c}" for d, c in top))
        for e in examples:
            print(f"  {e[0]} {e[1]}: CNBP={e[2]} truth={e[3]} diff={e[4]}")
    if n_cmp == 0:
        sys.exit("check_cnbp_oracle: nothing to compare -- was --allele-vcf given to call?")
    sys.exit(0 if n_exact == n_cmp else 1)


main()
