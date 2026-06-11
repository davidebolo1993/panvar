#!/usr/bin/env python3
"""Compare panvar C4 copy-number calls against a per-molecule ground-truth BED.

Two ground-truth signals are checked:

  * total copy number   -- number of C4 (C4A/C4B) features per molecule.
  * long/short split     -- long modules carry the HERV-K insertion (gene length
                            > --long-min-bp); short modules do not.

panvar recovers total copy number from the DUP record of a panphorte-normalised
graph (run panphorte with a relaxed --min-similarity, e.g. 0.70, so the long and
short C4 modules collapse onto a single repeat unit). The long/short composition
is read straight from panphorte's copies.tsv (n_long / n_short columns), which are
derived from each copy's traversed length relative to the (long) reference unit.

Example:
  scripts/compare_c4_copy_number.py \
      --bed tests/real_data/c4.bed \
      --copies-tsv build/c4_70/panphorte.panphorte.copies.tsv \
      --vcf build/c4_70/call.region.vcf
"""
import argparse
import collections
import csv
import re


def strip_inv(name):
    """odgi flip may suffix reverse-complemented paths with _inv; drop it to match the BED."""
    return re.sub(r"_inv$", "", name)


def load_bed(path, long_min_bp):
    total = collections.Counter()
    long_ = collections.Counter()
    short = collections.Counter()
    with open(path) as fh:
        header = fh.readline()  # molecule  gene  start  end  strand
        for line in fh:
            mol, gene, start, end, strand = line.rstrip("\n").split("\t")
            if gene in ("C4A", "C4B"):
                total[mol] += 1
                if int(end) - int(start) > long_min_bp:
                    long_[mol] += 1
                else:
                    short[mol] += 1
    return total, long_, short


def biggest_dup_cn(vcf_path):
    """Per-sample CN of the largest DUP record (the C4 module locus)."""
    samples = []
    best = None
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
    ap.add_argument("--bed", required=True, help="ground-truth BED (molecule, gene, start, end, strand)")
    ap.add_argument("--copies-tsv", help="panphorte .copies.tsv (for total CN + long/short)")
    ap.add_argument("--vcf", help="panvar call region VCF (alternative total-CN source: biggest DUP)")
    ap.add_argument("--long-min-bp", type=int, default=17000, help="BED feature length above which a C4 module is 'long'")
    args = ap.parse_args()

    bed_total, bed_long, bed_short = load_bed(args.bed, args.long_min_bp)

    if args.copies_tsv:
        rows = list(csv.DictReader(open(args.copies_tsv), delimiter="\t"))
        has_ls = rows and "n_long" in rows[0]
        n = len(rows)
        ct_exact = ls_exact = 0
        for row in rows:
            mol = strip_inv(row["path_name"])
            copies = int(row["copies"])
            if copies == bed_total.get(mol, 0):
                ct_exact += 1
            if has_ls and int(row["n_long"]) == bed_long.get(mol, 0) and int(row["n_short"]) == bed_short.get(mol, 0):
                ls_exact += 1
        print(f"[copies.tsv] normalised molecules: {n}")
        print(f"  total CN == bed:      {ct_exact}/{n} ({100*ct_exact/n:.0f}%)")
        if has_ls:
            print(f"  long & short == bed:  {ls_exact}/{n} ({100*ls_exact/n:.0f}%)")
        else:
            print("  (copies.tsv has no n_long/n_short columns — re-run a panphorte build that emits them)")

    if args.vcf:
        cn = biggest_dup_cn(args.vcf)
        keys = list(cn.keys())
        if keys:
            exact = sum(1 for s in keys if cn[s] == bed_total.get(s, 0))
            within1 = sum(1 for s in keys if abs(cn[s] - bed_total.get(s, 0)) <= 1)
            print(f"[VCF DUP] samples with a CN field: {len(keys)}")
            print(f"  total CN == bed:      {exact}/{len(keys)} ({100*exact/len(keys):.0f}%)")
            print(f"  total CN within ±1:   {within1}/{len(keys)} ({100*within1/len(keys):.0f}%)")
            called = collections.Counter(cn[s] for s in keys)
            bed_d = collections.Counter(bed_total.get(s, 0) for s in keys)
            print(f"  called CN dist: {dict(sorted(called.items()))}")
            print(f"  bed   CN dist: {dict(sorted(bed_d.items()))}")
        else:
            print("[VCF DUP] no DUP record found")


if __name__ == "__main__":
    main()
