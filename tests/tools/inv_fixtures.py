#!/usr/bin/env python3
"""Write the synthetic fixtures for tests/suites/call_invariants.sh.

Each fixture is a directory holding g.gfa and one or more bubbles CSVs written by hand, so the bubble
row is exactly the one under test (bubble would orient it to the reference, which is the point of
the reverse fixture). Node ids are numeric so tools/renumber.py and tools/unmap_vcf.py can map them.
Sequences are seeded pseudo-random, so no two nodes share content unless a fixture says so.

  rev    ref 1+2+4+5; an INS (node 3), a DEL (node 2), an INV (node 2-) and a substitution
         (node 2 -> node 6, one substitution). fwd.bubbles.csv has the row (1,4),
         rev.bubbles.csv the same bubble as (4,1)
  tie    ref 1+10+11+12+5 where 10, 11 and 12 spell the same 30 bp; hapA deletes 10+11 and hapB
         11+12, two 60 bp deletions with identical sequence 30 bp apart: they merge, and the
         representative is a size tie
  chain  three 100 bp insertions at one site on distinct nodes: B, A = B with its first 15 bases
         changed, C = B with its last 15 changed. A~B and B~C pass the 0.80 sequence identity,
         A~C (0.70) does not, so only single-linkage chaining puts all three in one record.
         g_ac.gfa drops B: the control that A and C really do not match directly
  sub    ref 1+2+4; hapSUB replaces node 2 by node 6, hapDEL deletes node 2. The two DELs merge;
         the INS is carried by hapSUB only. sub_delfirst.gfa / sub_subfirst.gfa differ only in
         P-line order
  loop   ref 1+7+2 with a 9 bp self-loop node 7; hapx10 walks node 7 ten times (+81 bp)

  inv_fixtures.py <out-dir>
"""
import os
import random
import sys

HDR = ('bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,'
       'path_support,distinct_alleles,ref_allele_support,alt_allele_support_max,'
       'alt_allele_support_min,min_inside_bp,max_inside_bp,inside_nodes\n')


def rs(rng, n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def mutate(rng, s, lo, hi):
    out = list(s)
    for i in range(lo, hi):
        out[i] = rng.choice([b for b in 'ACGT' if b != s[i]])
    return ''.join(out)


def write_gfa(path, seq, paths, extra_edges=()):
    edges = set(extra_edges)
    for p in paths.values():
        for a, b in zip(p, p[1:]):
            edges.add((a[:-1], a[-1], b[:-1], b[-1]))
    with open(path, 'w') as f:
        f.write('H\tVN:Z:1.0\n')
        for k in sorted(seq, key=int):
            f.write(f'S\t{k}\t{seq[k]}\n')
        for a, ao, b, bo in sorted(edges):
            f.write(f'L\t{a}\t{ao}\t{b}\t{bo}\t0M\n')
        for n, p in paths.items():
            f.write(f'P\t{n}\t{",".join(p)}\t*\n')


def write_csv(path, rows):
    with open(path, 'w') as f:
        f.write(HDR)
        for bid, src, snk, inside, npaths in rows:
            f.write(f'{bid},{src},+,{snk},+,{len(inside)},{len(inside) + 2},{npaths},{npaths},1,1,1,0,0,'
                    f'"{";".join(inside)}"\n')


def fx_rev(d, rng):
    seq = {'1': rs(rng, 40), '2': rs(rng, 60), '3': rs(rng, 70), '4': rs(rng, 40), '5': rs(rng, 25),
           '6': rs(rng, 55)}
    paths = {
        'ref#0#chr1:1001-1165': ['1+', '2+', '4+', '5+'],
        'hapINS#1#c1': ['1+', '2+', '3+', '4+', '5+'],
        'hapDEL#1#c1': ['1+', '4+', '5+'],
        'hapINV#1#c1': ['1+', '2-', '4+', '5+'],
        'hapSUB#1#c1': ['1+', '6+', '4+', '5+'],
    }
    write_gfa(f'{d}/g.gfa', seq, paths)
    write_csv(f'{d}/fwd.bubbles.csv', [(1, '1', '4', ['2', '3', '6'], 5)])
    write_csv(f'{d}/rev.bubbles.csv', [(1, '4', '1', ['2', '3', '6'], 5)])


def fx_tie(d, rng):
    unit = rs(rng, 30)
    seq = {'1': rs(rng, 40), '10': unit, '11': unit, '12': unit, '5': rs(rng, 40)}
    paths = {
        'ref#0#chr1:1-170': ['1+', '10+', '11+', '12+', '5+'],
        'hapA#1#c1': ['1+', '12+', '5+'],
        'hapB#1#c1': ['1+', '10+', '5+'],
        'hapR#1#c1': ['1+', '10+', '11+', '12+', '5+'],
    }
    write_gfa(f'{d}/g.gfa', seq, paths)
    write_csv(f'{d}/g.bubbles.csv', [(1, '1', '5', ['10', '11', '12'], 4)])


def fx_chain(d, rng):
    b = rs(rng, 100)
    seq = {'1': rs(rng, 50), '2': rs(rng, 50), '20': mutate(rng, b, 0, 15), '21': b,
           '22': mutate(rng, b, 85, 100)}
    paths = {
        'ref#0#chr1:1-100': ['1+', '2+'],
        'hapA#1#c1': ['1+', '20+', '2+'],
        'hapB#1#c1': ['1+', '21+', '2+'],
        'hapC#1#c1': ['1+', '22+', '2+'],
        'hapR#1#c1': ['1+', '2+'],
    }
    write_gfa(f'{d}/g.gfa', seq, paths)
    write_csv(f'{d}/g.bubbles.csv', [(1, '1', '2', ['20', '21', '22'], 5)])
    # the control: without the bridging B, A and C must stay two records
    del paths['hapB#1#c1']
    del seq['21']
    write_gfa(f'{d}/g_ac.gfa', seq, paths)
    write_csv(f'{d}/g_ac.bubbles.csv', [(1, '1', '2', ['20', '22'], 4)])


def fx_sub(d, rng):
    seq = {'1': rs(rng, 40), '2': rs(rng, 60), '4': rs(rng, 40), '6': rs(rng, 58)}
    ref = ('ref#0#chr1:1-140', ['1+', '2+', '4+'])
    dl = ('hapDEL#1#c1', ['1+', '4+'])
    sb = ('hapSUB#1#c1', ['1+', '6+', '4+'])
    rf = ('hapR#1#c1', ['1+', '2+', '4+'])
    write_gfa(f'{d}/sub_delfirst.gfa', seq, dict([ref, dl, sb, rf]))
    write_gfa(f'{d}/sub_subfirst.gfa', seq, dict([ref, sb, dl, rf]))
    write_csv(f'{d}/g.bubbles.csv', [(1, '1', '4', ['2', '6'], 4)])


def fx_loop(d, rng):
    seq = {'1': rs(rng, 60), '7': rs(rng, 9), '2': rs(rng, 60)}
    paths = {
        'ref#0#chr1:1-129': ['1+', '7+', '2+'],
        'hapx10#1#c1': ['1+'] + ['7+'] * 10 + ['2+'],
        'hapR#1#c1': ['1+', '7+', '2+'],
        'hapx1#2#c1': ['1+', '7+', '2+'],
    }
    write_gfa(f'{d}/g.gfa', seq, paths)
    write_csv(f'{d}/g.bubbles.csv', [(1, '1', '2', ['7'], 4)])


def main():
    out = sys.argv[1]
    for name, fn, seed in (('rev', fx_rev, 7), ('tie', fx_tie, 11), ('chain', fx_chain, 13),
                           ('sub', fx_sub, 17), ('loop', fx_loop, 19)):
        d = os.path.join(out, name)
        os.makedirs(d, exist_ok=True)
        fn(d, random.Random(seed))


if __name__ == '__main__':
    main()
