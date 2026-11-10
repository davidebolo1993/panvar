#!/usr/bin/env python3
"""Count bubbles with more than one copy-number (SVTYPE=CNV, or DUP in older VCFs) record in a call region VCF.

A bubble's copy number is one quantity: two such records on one bubble (a module record next to a
UNIT one, say) count its extra copies twice in any dosage sum. Every bubble so flagged in the
CURRENT output is an item; the ORIGINAL output is counted for context only.

  cnrecords.py <cur.region.vcf> [--orig orig.region.vcf] [--label L] [--items out.tsv]

--items rows: MULTICN, '-', bubble<N>, DUP, records, 1, the record ids, pin (a hash of the sorted
ids and positions). Exit 1 if any bubble in the current VCF has more than one such record.
"""
import argparse
import collections
import hashlib
import sys


def dups_per_bubble(path):
    out = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 8:
                continue
            info = dict(kv.partition('=')[::2] for kv in f[7].split(';'))
            if info.get('SVTYPE') in ('CNV', 'DUP'):   # DUP: VCFs from before <CNV>
                out[info.get('BUBBLE_ID', '?')].append((f[2], f[1]))
    return {b: v for b, v in out.items() if len(v) > 1}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cur')
    ap.add_argument('--orig')
    ap.add_argument('--label', default='')
    ap.add_argument('--items')
    a = ap.parse_args()
    cur = dups_per_bubble(a.cur)
    orig = dups_per_bubble(a.orig) if a.orig else None
    items = []
    for b, recs in sorted(cur.items()):
        recs = sorted(recs)
        items.append(('MULTICN', '-', f'bubble{b}', 'DUP', 'records', 1,
                      f'{len(recs)} DUP records: ' + ', '.join(f'{i}@{p}' for i, p in recs),
                      hashlib.sha1(repr(recs).encode()).hexdigest()[:12]))
    print(f'{a.label}: bubbles_multi_cn cur={len(cur)}' + (f' orig={len(orig)}' if orig is not None else ''))
    for x in items:
        print(f'  {x[2]}: {x[6]}')
    if a.items:
        with open(a.items, 'w') as fo:
            for x in items:
                fo.write('\t'.join(map(str, x)) + '\n')
    sys.exit(1 if items else 0)


if __name__ == '__main__':
    main()
