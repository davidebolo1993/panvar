#!/usr/bin/env python3
"""Compare the variant BIMBAM written by `describe --only-variant` for the ORIGINAL and CURRENT VCFs.

Both VCFs go through the same describe binary, so any difference here is a difference in what call
wrote. A feature id is the VCF record ID, and IDs renumber (the _2/_3 suffixes follow emission
order), so rows are paired in two passes over what is still unpaired:

  1. same id, same alleles and dosage vector     identical
  2. same dosage vector under another id         RENAME (two ids that swap vectors are two renames)

Whatever is left is the ID-agnostic difference, and only it decides criticality: a leftover id on
both sides is CHANGED, the rest REMOVED or ADDED. A paired row whose alleles columns differ is a
CHANGED [alleles] item that is not CN-critical (its dosages are the same).

Rows are matched to samples through samples.txt.gz, so a reordered sample list is not a change.

  bimbamdiff.py <orig.describe_dir> <cur.describe_dir> [--label L] [--items out.tsv]

--items rows: kind (ADDED/REMOVED/RENAME/CHANGED), orig_id, cur_id, '-', fields, critical, detail,
pin (a hash of the complete alleles and dosage vectors involved, in sample-name order).
Exit 1 when the rows differ in any way other than renames.
"""
import argparse
import collections
import gzip
import hashlib
import sys


def load(d):
    rows = collections.OrderedDict()
    with gzip.open(f'{d}/haplotype/variant/bimbam_variant.bimbam.gz', 'rt') as fh:
        for line in fh:
            f = [x.strip() for x in line.rstrip('\n').split(',')]
            rows[f[0]] = (tuple(f[1:3]), tuple(f[3:]))
    with gzip.open(f'{d}/haplotype/variant/samples.txt.gz', 'rt') as fh:
        samples = [l.strip() for l in fh if l.strip()]
    return samples, rows


def ncar(x):
    return sum(1 for v in x if v not in ('0', 'NA'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('orig')
    ap.add_argument('cur')
    ap.add_argument('--label', default='')
    ap.add_argument('--items')
    a = ap.parse_args()
    sa, A = load(a.orig)
    sb, B = load(a.cur)
    items = []
    if sorted(sa) != sorted(sb):
        print(f'{a.label} bimbam: SAMPLE LISTS DIFFER ({len(sa)} vs {len(sb)})')
        items.append(('SAMPLES', '-', '-', '-', 'samples', 1, f'{len(sa)} vs {len(sb)}'))
    # sample-name order, so a pin does not depend on the order describe wrote
    common = sorted(set(sa) & set(sb))
    ia = [sa.index(s) for s in common]
    ib = [sb.index(s) for s in common]
    A = {k: (al, tuple(v[i] for i in ia)) for k, (al, v) in A.items()}
    B = {k: (al, tuple(v[i] for i in ib)) for k, (al, v) in B.items()}

    def pin(*rows):
        return hashlib.sha1('|'.join(f'{al}:{",".join(x)}' for al, x in rows).encode()).hexdigest()[:12]

    ua = [k for k in A if not (k in B and B[k] == A[k])]
    ub = {k for k in B if not (k in A and A[k] == B[k])}
    by_vec = collections.defaultdict(list)
    for k in B:
        if k in ub:
            by_vec[B[k][1]].append(k)
    rest = []
    for k in ua:
        cand = by_vec.get(A[k][1])
        if not cand:
            rest.append(k)
            continue
        j = k if k in cand else cand[0]   # the same id with other alleles first, then emission order
        cand.remove(j)
        ub.discard(j)
        if j != k:
            items.append(('RENAME', k, j, '-', 'ID', 0, 'same dosage vector', pin(A[k], B[j])))
        if A[k][0] != B[j][0]:
            items.append(('CHANGED', k, j, '-', 'alleles', 0, f'alleles {"/".join(A[k][0])}->{"/".join(B[j][0])}',
                          pin(A[k], B[j])))
    for k in rest:
        if k in ub:
            (ala, xa), (alb, xb) = A[k], B[k]
            ub.discard(k)
            f = ['alleles'] if ala != alb else []
            f.append('dosage')
            n = sum(1 for x, y in zip(xa, xb) if x != y)
            items.append(('CHANGED', k, k, '-', ','.join(f), 1, f'{n} haplotype dosages differ, '
                          f'nonzero {ncar(xa)}->{ncar(xb)}', pin(A[k], B[k])))
        else:
            items.append(('REMOVED', k, '-', '-', '-', 1, f'nonzero dosages={ncar(A[k][1])}', pin(A[k])))
    for k in B:
        if k in ub:
            items.append(('ADDED', '-', k, '-', '-', 1, f'nonzero dosages={ncar(B[k][1])}', pin(B[k])))
    va = collections.Counter(v for _, v in A.values())
    vb = collections.Counter(v for _, v in B.values())
    only_a, only_b = va - vb, vb - va
    nchg = sum(1 for x in items if x[0] == 'CHANGED')
    print(f'{a.label} bimbam: rows orig={len(A)} cur={len(B)} | by-id changed={nchg} '
          f'only_orig={sum(1 for k in A if k not in B)} only_cur={sum(1 for k in B if k not in A)} | '
          f'id-agnostic vectors_only_orig={sum(only_a.values())} vectors_only_cur={sum(only_b.values())}')
    for x in items:
        print(f'  {x[0]} {x[1]}' + (f' -> {x[2]}' if x[2] not in ('-', x[1]) else '') + f': {x[6]}')
    if a.items:
        with open(a.items, 'w') as fo:
            for x in items:
                fo.write('\t'.join(str(v) for v in x) + '\n')
    sys.exit(1 if any(x[0] != 'RENAME' for x in items) else 0)


if __name__ == '__main__':
    main()
