#!/usr/bin/env python3
"""Assert FORMAT:CNBP against the allele VCF, which is an exact oracle for it.

CNBP claims to be the signed net bases a haplotype carries across a copy-number module relative to the
reference. The allele VCF spells every allele of every bubble explicitly, so for a haplotype whose GT
names allele k at that bubble, `len(allele k) - len(REF allele)` IS that quantity. The two are computed
by different code -- one sums node lengths over the widest source..sink span, the other spells the
walk -- so agreeing to the base is a real check on the span logic and on repeated-node multiplicity,
which are the two places this is easy to get subtly wrong.

It is an internal check, not a measure of biological truth: both sides read the same graph walks. Use
compare_svim_asm.py for the external comparison, and note its own ceiling -- it too starts from the
graph's haplotypes.

Everything it cannot compare is an error, not a smaller denominator. A silently shrinking denominator is
how a check reports success for work it never did: the run that compares nothing passes loudest.

  check_cnbp_oracle.py <region.vcf> <alleles.vcf> [--verbose] [--allow-skips]

Exits non-zero if any comparison disagrees, or if any expected comparison could not be made.
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
    if hdr is None:
        sys.exit(f"check_cnbp_oracle: {path} has no #CHROM header line")
    return hdr, rows


def info_of(field):
    d = {}
    for kv in field.split(";"):
        k, sep, v = kv.partition("=")
        d[k] = v if sep else True
    return d


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    verbose = "--verbose" in sys.argv
    allow_skips = "--allow-skips" in sys.argv
    if len(args) != 2:
        sys.exit(__doc__)
    rhdr, rrows = read_vcf(args[0])
    ahdr, arows = read_vcf(args[1])

    # Keyed on INFO/BUBBLE_ID, which both files carry and which is what actually ties a region record to
    # its allele record. The record ID prefix happens to encode the bubble today; it is a display string,
    # and keying on it would break silently the moment the naming changes.
    allele_rec = {}
    for r in arows:
        bid = info_of(r[7]).get("BUBBLE_ID")
        if bid is None:
            sys.exit(f"check_cnbp_oracle: allele record {r[2]} has no BUBBLE_ID")
        if bid in allele_rec:
            sys.exit(f"check_cnbp_oracle: two allele records claim bubble {bid} "
                     f"({allele_rec[bid][2]} and {r[2]})")
        allele_rec[bid] = r

    n_dup = n_expected = n_compared = n_exact = 0
    skipped = collections.Counter()
    err = collections.Counter()
    examples = []

    for r in rrows:
        if info_of(r[7]).get("SVTYPE") != "DUP":
            continue
        n_dup += 1
        bid = info_of(r[7]).get("BUBBLE_ID")
        rkeys = r[8].split(":")
        if "CNBP" not in rkeys:
            skipped["region record has no CNBP in FORMAT"] += 1
            continue
        ci = rkeys.index("CNBP")

        # every haplotype this record reports a CNBP for is a comparison we OWE
        want = [i for i in range(9, len(rhdr))
                if len(r[i].split(":")) > ci and r[i].split(":")[ci] != "."]
        n_expected += len(want)

        a = allele_rec.get(bid)
        if a is None:
            skipped[f"no allele record for bubble {bid}"] += len(want)
            continue
        akeys = a[8].split(":")
        if "GT" not in akeys:
            skipped[f"allele record for bubble {bid} has no GT"] += len(want)
            continue
        agt = akeys.index("GT")
        lens = [len(a[3])] + [len(x) for x in a[4].split(",")]
        hap_allele = {ahdr[i]: a[i].split(":")[agt] for i in range(9, len(ahdr))}

        for i in want:
            sample = rhdr[i]
            gt = hap_allele.get(sample)
            if gt is None:
                skipped["sample absent from the allele VCF"] += 1
                continue
            if not gt.isdigit():
                # A haplotype with a CNBP crosses the bubble, so the allele VCF must know which allele
                # it carries. '.' here means the two files disagree about who traverses what.
                skipped["sample has no allele call where the region VCF gave a CNBP"] += 1
                continue
            idx = int(gt)
            if idx >= len(lens):
                skipped["allele index past the end of the ALT list"] += 1
                continue
            try:
                got = int(r[i].split(":")[ci])
            except ValueError:
                skipped["CNBP is not an integer"] += 1
                continue
            truth = lens[idx] - lens[0]
            diff = got - truth
            n_compared += 1
            err[diff] += 1
            if diff == 0:
                n_exact += 1
            elif len(examples) < 8:
                examples.append((r[2], sample, got, truth, diff))

    n_skipped = n_expected - n_compared
    print(f"DUP records {n_dup}")
    print(f"comparisons expected {n_expected}  compared {n_compared}  skipped {n_skipped}")
    pct = 100.0 * n_exact / n_compared if n_compared else 0.0
    print(f"exact {n_exact} of {n_compared} ({pct:.2f}%)")
    for reason, count in skipped.most_common():
        print(f"  skipped: {count:6d}  {reason}")
    if verbose or n_exact != n_compared:
        top = sorted(err.items(), key=lambda kv: -kv[1])[:8]
        print("CNBP minus allele delta:", ", ".join(f"{d}:{c}" for d, c in top))
        for e in examples:
            print(f"  {e[0]} {e[1]}: CNBP={e[2]} truth={e[3]} diff={e[4]}")

    if n_expected == 0:
        sys.exit("FAIL: no comparison was expected at all -- was --allele-vcf given to call, and does "
                 "the region VCF contain DUP records?")
    bad = []
    if n_exact != n_compared:
        bad.append(f"{n_compared - n_exact} disagreement(s)")
    if n_skipped and not allow_skips:
        bad.append(f"{n_skipped} expected comparison(s) could not be made")
    if bad:
        sys.exit("FAIL: " + "; ".join(bad))
    print("OK")


main()
