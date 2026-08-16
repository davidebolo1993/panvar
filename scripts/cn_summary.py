#!/usr/bin/env python3
"""Compact, diffable summary of what `call` reported for a results tree.

Written to answer one question across a pipeline change: did the copy-number calls move? Comparing two
runs by eye over six loci and thousands of records does not scale, and the interesting differences are
small. This prints one stable, sorted block per locus so `diff` does the work.

Deliberately excludes anything that is expected to move for uninteresting reasons -- node ids, bubble
ids and record ids are all reassigned by sorting and re-snarling, so a record is keyed by its reference
position and type rather than its name.

  cn_summary.py <results-dir> [locus ...]
"""
import collections
import os
import sys


def info_of(field):
    d = {}
    for kv in field.split(";"):
        k, sep, v = kv.partition("=")
        d[k] = v if sep else True
    return d


def read_vcf(path):
    hdr, rows = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
            elif not line.startswith("#"):
                rows.append(line.rstrip("\n").split("\t"))
    return hdr, rows


def summarize_region(path, out):
    hdr, rows = read_vcf(path)
    if hdr is None:
        out.append("  region VCF: unreadable")
        return
    by_type = collections.Counter(info_of(r[7]).get("SVTYPE", "?") for r in rows)
    out.append(f"  records {len(rows)}  " +
               " ".join(f"{k}={by_type[k]}" for k in sorted(by_type)))

    for r in sorted(rows, key=lambda r: (int(r[1]), info_of(r[7]).get("SVTYPE", ""))):
        inf = info_of(r[7])
        if inf.get("SVTYPE") != "DUP":
            continue
        keys = r[8].split(":")
        cn = collections.Counter()
        cnbp = {}
        if "CN" in keys:
            ci = keys.index("CN")
            bi = keys.index("CNBP") if "CNBP" in keys else None
            for i in range(9, len(hdr)):
                v = r[i].split(":")
                if len(v) <= ci or v[ci] == ".":
                    continue
                c = int(v[ci])
                cn[c] += 1
                if bi is not None and len(v) > bi and v[bi] != ".":
                    cnbp.setdefault(c, []).append(int(v[bi]))
        # POS is the anchor; keying on it rather than the record id keeps this stable across re-snarls
        out.append(f"  DUP at {r[1]}  REF_CN={inf.get('REF_CN','-')} RU_LEN={inf.get('RU_LEN','-')} "
                   f"SVLEN={inf.get('SVLEN','-')} AC={inf.get('AC','-')} AN={inf.get('AN','-')}")
        out.append("    CN " + " ".join(f"{c}:{n}" for c, n in sorted(cn.items())))
        if cnbp:
            med = {c: sorted(v)[len(v) // 2] for c, v in cnbp.items()}
            out.append("    median CNBP " + " ".join(f"{c}:{med[c]}" for c in sorted(med)))


def summarize_genes(path, out):
    if not os.path.exists(path):
        out.append("  dup_gene_cn.tsv: absent")
        return
    per = collections.Counter()
    miss = collections.Counter()
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        gi = head.index("genes") if "genes" in head else 3
        ci = head.index("cn") if "cn" in head else 4
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(gi, ci):
                continue
            if f[ci] in (".", ""):
                miss[f[gi]] += 1
            else:
                per[(f[gi], f[ci])] += 1
    out.append(f"  dup_gene_cn rows {sum(per.values()) + sum(miss.values())}")
    for gene in sorted({g for g, _ in per} | set(miss)):
        counts = " ".join(f"CN{c}:{n}" for (g, c), n in sorted(per.items()) if g == gene)
        nomiss = f"  missing:{miss[gene]}" if miss[gene] else ""
        out.append(f"    {gene:22s} {counts}{nomiss}")


def main():
    root = sys.argv[1]
    loci = sys.argv[2:] or sorted(
        d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)))
    out = []
    for locus in loci:
        out.append(f"### {locus}")
        call = os.path.join(root, locus, "call")
        region = os.path.join(call, "call.region.vcf")
        if os.path.exists(region):
            summarize_region(region, out)
        else:
            out.append("  region VCF: absent")
        summarize_genes(os.path.join(call, "call.dup_gene_cn.tsv"), out)
    print("\n".join(out))


main()
