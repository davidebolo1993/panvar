#!/usr/bin/env python3
"""Build a small GFA and its site table from a line spec, for call fixtures.

    spec_gfa.py <spec> <outdir>   writes <outdir>/g.gfa and <outdir>/g.bubbles.csv

Spec lines (# comments):
    N <name> <len> <seed>      a node of <len> pseudo-random bases
    N <name> =<SEQ>            a node with this sequence
    N <name> @<src>:<k>:<seed> a copy of node <src> with <k> bases changed
    P <name> <n1+,n2+,...>     a path; every adjacency it walks becomes an L line
    B <source> <sink> <a;b;c>  a site (source and sink forward) with these interior nodes
"""
import os
import random
import sys

spec, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
nodes, order, paths, sites = {}, [], [], []
for line in open(spec):
    t = line.split()
    if not t or t[0].startswith('#'):
        continue
    if t[0] == 'N':
        if t[2].startswith('='):
            nodes[t[1]] = t[2][1:]
        elif t[2].startswith('@'):
            src, k, seed = t[2][1:].split(':')
            s, r = list(nodes[src]), random.Random(int(seed))
            for _ in range(int(k)):
                i = r.randrange(len(s))
                s[i] = r.choice([c for c in 'ACGT' if c != s[i]])
            nodes[t[1]] = ''.join(s)
        else:
            r = random.Random(int(t[3]))
            nodes[t[1]] = ''.join(r.choice('ACGT') for _ in range(int(t[2])))
        order.append(t[1])
    elif t[0] == 'P':
        paths.append((t[1], t[2].split(',')))
    elif t[0] == 'B':
        sites.append((t[1], t[2], t[3]))
links = set()
for _, steps in paths:
    for a, b in zip(steps, steps[1:]):
        links.add((a[:-1], a[-1], b[:-1], b[-1]))
with open(os.path.join(out, 'g.gfa'), 'w') as f:
    f.write('H\tVN:Z:1.0\n')
    for n in order:
        f.write(f'S\t{n}\t{nodes[n]}\n')
    for a, ao, b, bo in sorted(links):
        f.write(f'L\t{a}\t{ao}\t{b}\t{bo}\t0M\n')
    for name, steps in paths:
        f.write(f'P\t{name}\t{",".join(steps)}\t*\n')
with open(os.path.join(out, 'g.bubbles.csv'), 'w') as f:
    f.write('bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,'
            'distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,'
            'max_inside_bp,inside_nodes\n')
    for i, (src, snk, inside) in enumerate(sites, 1):
        n = len(inside.split(';'))
        f.write(f'{i},{src},+,{snk},+,{n},{n + 2},{len(paths)},1,1,1,1,0,0,"{inside}"\n')
