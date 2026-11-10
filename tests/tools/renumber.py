#!/usr/bin/env python3
"""Renumber every node id and bubble id, consistently in a GFA and its bubbles CSV.

Node ids go through a seeded random bijection onto a disjoint range (500000+), bubble ids through a
shuffle of 1..n, so nothing keeps its old number or its old relative order. The map is written to
<out.gfa>.map (old<TAB>new, bubble ids prefixed with B) for unmap_vcf.py.

  renumber.py <in.gfa> <in.bubbles.csv> <out.gfa> <out.bubbles.csv> <seed>
"""
import csv
import random
import re
import sys

csv.field_size_limit(10 ** 9)
gfa, bcsv, ogfa, obcsv, seed = sys.argv[1:6]
rng = random.Random(int(seed))
with open(gfa) as fh:
    ids = [l.split('\t', 2)[1] for l in fh if l.startswith('S\t')]
new = list(range(500000, 500000 + len(ids)))
rng.shuffle(new)
m = {o: str(n) for o, n in zip(ids, new)}
with open(gfa) as fh, open(ogfa, 'w') as out:
    for l in fh:
        f = l.rstrip('\n').split('\t')
        if f[0] == 'S':
            f[1] = m[f[1]]
        elif f[0] == 'L':
            f[1], f[3] = m[f[1]], m[f[3]]
        elif f[0] == 'P':
            f[2] = ','.join(m[s[:-1]] + s[-1] for s in f[2].split(','))
        elif f[0] == 'W':
            f[6] = re.sub(r'([<>])([^<>]+)', lambda x: x.group(1) + m[x.group(2)], f[6])
        out.write('\t'.join(f) + '\n')
with open(bcsv) as fh:
    rows = list(csv.DictReader(fh))
bids = [r['bubble_id'] for r in rows]
nb = list(range(1, len(bids) + 1))
rng.shuffle(nb)
bm = dict(zip(bids, map(str, nb)))
with open(obcsv, 'w', newline='') as out:
    w = csv.DictWriter(out, fieldnames=list(rows[0].keys()), quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
    w.writeheader()
    for r in rows:
        r['bubble_id'] = bm[r['bubble_id']]
        r['source'], r['sink'] = m[r['source']], m[r['sink']]
        r['inside_nodes'] = ';'.join(m[x] for x in r['inside_nodes'].split(';') if x)
        w.writerow(r)
with open(ogfa + '.map', 'w') as out:
    for o, n in m.items():
        out.write(f'{o}\t{n}\n')
    for o, n in bm.items():
        out.write(f'B{o}\tB{n}\n')
