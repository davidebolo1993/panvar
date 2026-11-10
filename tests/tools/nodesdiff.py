#!/usr/bin/env python3
"""Compare the variant_nodes.tsv that `panvar call` writes next to the region VCF, ORIGINAL vs CURRENT.

One row per variant: variant_id, bubble_id, svtype, node_ids. IDs renumber, so rows are paired in
two passes over what is still unpaired, as bimbamdiff.py does:

  1. same id, same row        identical
  2. same bubble, svtype and node set under another id     RENAME (not a difference by itself)

What is left is a difference: an id on both sides is CHANGED (fields bubble_id, svtype, node_ids),
the rest REMOVED or ADDED. A node set is compared as a set; its order is reported separately
(fields node_order), since associate reads it as a set.

  nodesdiff.py <orig.variant_nodes.tsv> <cur.variant_nodes.tsv> [--label L] [--items out.tsv]

--items rows: kind, orig_id, cur_id, svtype, fields, critical, detail, pin (a hash of the complete
rows involved). Exit 1 when the tables differ in any way other than renames.
"""
import argparse
import collections
import hashlib
import sys


def h(s):
    return hashlib.sha1(str(s).encode()).hexdigest()[:12]


def load(path):
    with open(path) as fh:
        head = fh.readline().rstrip('\n')
        rows, dups = collections.OrderedDict(), []
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 4:
                continue
            if f[0] in rows:
                dups.append(f[0])
            rows[f[0]] = (f[1], f[2], f[3])
    return head, rows, dups


def key(r):
    return (r[0], r[1], frozenset(r[2].split(',')))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('orig')
    ap.add_argument('cur')
    ap.add_argument('--label', default='')
    ap.add_argument('--items')
    a = ap.parse_args()
    ha, A, da = load(a.orig)
    hb, B, db = load(a.cur)
    items = []   # (kind, orig_id, cur_id, svtype, fields, critical, detail, pin)
    if ha != hb:
        items.append(('COLUMNS', '-', '-', '-', 'header', 1, f'{ha} -> {hb}', h(f'{ha}|{hb}')))
    for side, d in (('orig', da), ('cur', db)):
        for k in d:
            items.append(('DUPID', k, '-', '-', 'id', 1, f'{side} lists {k} more than once', h(f'{side}|{k}')))

    ua = [k for k in A if not (k in B and B[k] == A[k])]
    ub = {k for k in B if not (k in A and A[k] == B[k])}
    by_key = collections.defaultdict(list)
    for k in B:
        if k in ub:
            by_key[key(B[k])].append(k)
    rest = []
    for k in ua:
        cand = by_key.get(key(A[k]))
        if not cand:
            rest.append(k)
            continue
        j = k if k in cand else cand[0]
        cand.remove(j)
        ub.discard(j)
        if j != k:
            items.append(('RENAME', k, j, A[k][1], 'ID', 0, 'same node set', h(f'{k}->{j}')))
        if A[k][2] != B[j][2]:
            items.append(('CHANGED', k, j, A[k][1], 'node_order', 0, 'same node set, different order',
                          h(f'{A[k]}|{B[j]}')))
    for k in rest:
        if k in ub:
            ub.discard(k)
            (ba, sa, na), (bb, sb, nb) = A[k], B[k]
            f, det = [], []
            if ba != bb:
                f.append('bubble_id')
                det.append(f'bubble {ba}->{bb}')
            if sa != sb:
                f.append('svtype')
                det.append(f'svtype {sa}->{sb}')
            xa, xb = set(na.split(',')), set(nb.split(','))
            if xa != xb:
                f.append('node_ids')
                det.append(f'nodes {len(xa)}->{len(xb)} (-{len(xa - xb)} +{len(xb - xa)})')
            items.append(('CHANGED', k, k, sa, ','.join(f), 0, '; '.join(det), h(f'{A[k]}|{B[k]}')))
        else:
            items.append(('REMOVED', k, '-', A[k][1], '-', 1, f'bubble {A[k][0]}, {len(A[k][2].split(","))} nodes',
                          h(A[k])))
    for k in B:
        if k in ub:
            items.append(('ADDED', '-', k, B[k][1], '-', 1, f'bubble {B[k][0]}, {len(B[k][2].split(","))} nodes',
                          h(B[k])))
    n = collections.Counter(x[0] for x in items)
    print(f'{a.label}: orig={len(A)} cur={len(B)} changed={n["CHANGED"]} renamed={n["RENAME"]} '
          f'removed={n["REMOVED"]} added={n["ADDED"]}')
    for x in items:
        print(f'  {x[0]} {x[3]} {x[1]}' + (f' -> {x[2]}' if x[2] not in ('-', x[1]) else '')
              + (f' [{x[4]}]' if x[0] == 'CHANGED' else '') + f': {x[6]}')
    if a.items:
        with open(a.items, 'w') as fo:
            for x in items:
                fo.write('\t'.join(str(v) for v in x) + '\n')
    sys.exit(1 if any(x[0] != 'RENAME' for x in items) else 0)


if __name__ == '__main__':
    main()
