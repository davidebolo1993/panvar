#!/usr/bin/env python3
"""Map a call VCF made on a renumber.py graph back to the original node and bubble ids.

Rewrites the record ID (bubble<b>_<TYPE>_<node>[_n], bubble<b>_<NAME>), BUBBLE_ID, START_NODE,
END_NODE, EVENT_NODES and EVENTID (bubble<b>_sub_<node>_<node>). Anything else is copied as is, so a
field that still differs after this is a real renumbering dependence.

  unmap_vcf.py <in.vcf> <renumbered.gfa.map> <out.vcf>
"""
import re
import sys

vcf, mp, out = sys.argv[1:4]
m, bm = {}, {}
with open(mp) as fh:
    for l in fh:
        o, n = l.rstrip('\n').split('\t')
        if o.startswith('B'):
            bm[n[1:]] = o[1:]
        else:
            m[n] = o


def node(x):
    return m.get(x, x)


def rid(v):
    mm = re.match(r'bubble(\d+)_(.*)$', v)
    if not mm:
        return v
    b, rest = bm.get(mm.group(1), mm.group(1)), mm.group(2)
    mm2 = re.match(r'(sub|[A-Z]+)_(\d+)(_\d+)?$', rest)
    if mm2 and mm2.group(1) != 'sub':
        rest = f'{mm2.group(1)}_{node(mm2.group(2))}{mm2.group(3) or ""}'
    else:
        mm3 = re.match(r'sub_(\d+)_(\d+)$', rest)
        if mm3:
            rest = f'sub_{node(mm3.group(1))}_{node(mm3.group(2))}'
    return f'bubble{b}_{rest}'


with open(vcf) as fh, open(out, 'w') as fo:
    for l in fh:
        if l.startswith('#'):
            fo.write(l)
            continue
        f = l.rstrip('\n').split('\t')
        f[2] = rid(f[2])
        kv = []
        for x in f[7].split(';'):
            k, eq, v = x.partition('=')
            if k in ('START_NODE', 'END_NODE', 'EVENT_NODES'):
                v = ','.join(node(y) for y in v.split(','))
            elif k == 'BUBBLE_ID':
                v = bm.get(v, v)
            elif k == 'EVENTID':
                v = rid(v)
            kv.append(k + ('=' + v if eq else ''))
        f[7] = ';'.join(kv)
        fo.write('\t'.join(f) + '\n')
