#!/usr/bin/env python3
"""Compare panvar copy-number calls against per-haplotype ground truth, for any of the test genes.

Ground-truth formats differ by gene, selected with --mode:

  gene-count : an annotation BED (molecule, gene, start, end, strand); expected CN = number of rows
               whose gene is in --genes, per molecule.   (C4: --genes C4A,C4B ; CYP2D6: --genes CYP2D6)
  direct     : a 2-column table  <name> <tab> <count>  (GSTM1: gstm1.bed).
  repeats    : a header table with a name column and a count column
               (LPA: lpa.repeats.tsv, --name-col hapl --count-col copies).

Called CN is read from the panvar `call` region VCF: the per-sample CN of the largest DUP record (the
gene's copy-number module). Optionally also checks panphorte's copies.tsv (total + long/short for C4).

Reports exact and within-±1 concordance plus the called/expected CN distributions, so a gene either
"matches ground truth" or the mismatch is visible.

Examples:
  compare_copy_number.py --mode gene-count --truth c4.bed     --genes C4A,C4B --vcf c4/call/call.region.vcf
  compare_copy_number.py --mode gene-count --truth cyp2d6.bed --genes CYP2D6 --vcf cyp2d6/call/call.region.vcf
  compare_copy_number.py --mode direct     --truth gstm1.bed                  --vcf gstm1/call/call.region.vcf
  compare_copy_number.py --mode repeats    --truth lpa.repeats.tsv --name-col hapl --count-col copies \
                         --vcf lpa/call/call.region.vcf
"""
import argparse
import collections
import re
import sys


def strip_inv(name):
    """panphorte's reference flip may suffix reverse-complemented paths with _inv; drop it to match."""
    return re.sub(r"_inv$", "", name)


def load_truth(args):
    """Return {molecule -> expected copy number}."""
    exp = collections.Counter()
    if args.mode == "gene-count":
        genes = set(g.strip() for g in args.genes.split(",") if g.strip())
        with open(args.truth) as fh:
            fh.readline()  # header: molecule  gene  start  end  strand
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2:
                    continue
                if f[1] in genes:
                    exp[f[0]] += 1
        # molecules present in the BED but with zero target-gene rows still count as 0
        with open(args.truth) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f and f[0] not in exp:
                    exp[f[0]] = exp.get(f[0], 0)
    elif args.mode == "direct":
        with open(args.truth) as fh:
            first = fh.readline().rstrip("\n").split("\t")
            # header iff the 2nd column is not an integer
            if len(first) >= 2 and not first[1].lstrip("-").isdigit():
                pass  # header consumed
            elif len(first) >= 2:
                exp[first[0]] = int(first[1])
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2 and f[1].lstrip("-").isdigit():
                    exp[f[0]] = int(f[1])
    elif args.mode == "repeats":
        with open(args.truth) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ni, ci = hdr.index(args.name_col), hdr.index(args.count_col)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(ni, ci):
                    exp[f[ni]] = int(f[ci])
    return exp


def biggest_dup_cn(vcf_path):
    """Per-sample CN of the largest DUP record (the gene's copy-number module)."""
    samples, best = [], None
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if "SVTYPE=DUP" not in f[7]:
                continue
            svlen = abs(int(next(x for x in f[7].split(";") if x.startswith("SVLEN="))[6:]))
            if best is None or svlen > best[0]:
                best = (svlen, f)
    cn = {}
    if best is not None:
        for s, g in zip(samples, best[1][9:]):
            p = g.split(":")
            cn[strip_inv(s)] = int(p[1]) if len(p) > 1 and p[1] not in (".", "") else 0
    return cn


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", required=True, choices=["gene-count", "direct", "repeats"])
    ap.add_argument("--truth", required=True, help="ground-truth file (BED / 2-col / repeats table)")
    ap.add_argument("--vcf", required=True, help="panvar call region VCF")
    ap.add_argument("--genes", default="", help="comma list of target genes (gene-count mode)")
    ap.add_argument("--name-col", default="hapl", help="name column (repeats mode)")
    ap.add_argument("--count-col", default="copies", help="count column (repeats mode)")
    ap.add_argument("--label", default="", help="gene label for the report header")
    ap.add_argument("--dump-tsv", default="", help="append per-haplotype rows "
                    "(gene<TAB>sample<TAB>truth_cn<TAB>called_cn<TAB>is_reference) to this TSV, for plotting")
    ap.add_argument("--ref-substr", default="grch38,GRCh38,chm13,CHM13",
                    help="comma list; a sample whose name contains any is flagged is_reference=1")
    ap.add_argument("--offset", default="auto",
                    help="constant baseline offset between called and truth: an int, or 'auto' "
                         "(the most common called-minus-truth). Peak-multiplicity / coverage CN can sit "
                         "a fixed amount above the gene-count truth; after removing it the variation is "
                         "what should match. Use 0 for methods that report absolute CN (self-loop DUP).")
    args = ap.parse_args()

    exp = load_truth(args)
    cn = biggest_dup_cn(args.vcf)
    # compare on haplotypes present in BOTH the call and the ground truth
    shared = [s for s in cn if s in exp]
    label = args.label or args.mode
    print(f"== {label}: copy-number vs ground truth ==")
    print(f"  called haplotypes: {len(cn)} | ground-truth haplotypes: {len(exp)} | shared: {len(shared)}")
    if not shared:
        print("  NO SHARED HAPLOTYPES (name mismatch between VCF samples and ground truth)")
        sys.exit(1)
    n = len(shared)
    deltas = [cn[s] - exp[s] for s in shared]
    if args.offset == "auto":
        offset = collections.Counter(deltas).most_common(1)[0][0]
    else:
        offset = int(args.offset)
    exact = sum(1 for s in shared if cn[s] == exp[s])
    within1 = sum(1 for s in shared if abs(cn[s] - exp[s]) <= 1)
    exact_off = sum(1 for s in shared if cn[s] - offset == exp[s])
    within1_off = sum(1 for s in shared if abs(cn[s] - offset - exp[s]) <= 1)
    print(f"  exact CN == truth:  {exact}/{n} ({100*exact/n:.1f}%)")
    print(f"  CN within ±1:       {within1}/{n} ({100*within1/n:.1f}%)")
    print(f"  baseline offset (called-truth, mode): {offset:+d}")
    print(f"  exact after offset: {exact_off}/{n} ({100*exact_off/n:.1f}%)   "
          f"within ±1 after offset: {within1_off}/{n} ({100*within1_off/n:.1f}%)")
    print(f"  called CN dist: {dict(sorted(collections.Counter(cn[s] for s in shared).items()))}")
    print(f"  truth  CN dist: {dict(sorted(collections.Counter(exp[s] for s in shared).items()))}")
    if args.dump_tsv:
        refs = [r for r in args.ref_substr.split(",") if r]
        import os
        new = not os.path.exists(args.dump_tsv) or os.path.getsize(args.dump_tsv) == 0
        with open(args.dump_tsv, "a") as out:
            if new:
                # called_cn = absolute CN from the VCF (total collapsed-module count); baseline_offset =
                # the constant paralog baseline (called-minus-truth mode) so a plot can show the recovered
                # GENE copy number called_cn - baseline_offset against truth_cn.
                out.write("gene\tsample\ttruth_cn\tcalled_cn\tbaseline_offset\tis_reference\n")
            for s in shared:
                is_ref = 1 if any(r in s for r in refs) else 0
                out.write(f"{label}\t{s}\t{exp[s]}\t{cn[s]}\t{offset}\t{is_ref}\n")

    mism = [(s, cn[s], exp[s]) for s in shared if abs(cn[s] - offset - exp[s]) > 1][:10]
    if mism:
        print("  sample mismatches (|delta|>1), up to 10:")
        for s, c, e in mism:
            print(f"    {s}  called={c} truth={e}")


if __name__ == "__main__":
    main()
