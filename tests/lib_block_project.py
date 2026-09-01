#!/usr/bin/env python3
"""Project named GFA paths onto the block chain, from the WALK.

    lib_block_project.py <gfa> <bubbles.csv> <name> [<name> ...]

Prints, per path and block:  name <TAB> block <TAB> kind <TAB> bp <TAB> md5

WHY THIS EXISTS. --force-haplotypes ADDS to the shortlist, it does not restrict it: forcing
HG00731/HG03248 at c4 still returned HG00096 as rank 1. Using it to project an arbitrary pair
therefore reports the CALLED pair's alleles under both labels, and every block compares equal --
measured, and it produced a table reading "11/11 exact" for a call sitting 32732 edits above the
floor. So the projection is computed here instead, from the graph, with no dependence on which pair
the caller preferred.

Blocks alternate flank/bubble: block 2i is the flank before bubble i+1, block 2i+1 is bubble i+1,
and there is a trailing flank. That is 2*nbubbles+1 blocks, matching what the caller emits (c4: 5
bubbles, 11 blocks, block 0 kind=flank).

A path that does not contain a bubble's source or sink cannot be segmented there; those blocks are
reported as ABSENT rather than guessed, which is the same rule the partial-frame work applies.
"""
import sys, csv, hashlib

def rc(s):
    return s[::-1].translate(str.maketrans('ACGTacgtN', 'TGCAtgcaN'))

def main():
    gfa, bub = sys.argv[1], sys.argv[2]
    want = set(sys.argv[3:])
    seq, paths = {}, {}
    for line in open(gfa):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'S':
            seq[f[1]] = f[2]
        elif f[0] == 'P' and f[1] in want:
            paths[f[1]] = f[2].split(',')
    order = []
    with open(bub) as fh:
        for row in csv.DictReader(fh):
            order.append((int(row['bubble_id']), row['source'], row['sink']))
    order.sort()
    for name in sys.argv[3:]:
        steps = paths.get(name)
        if steps is None:
            print('%s\tNA\tNOPATH\t\t' % name)
            continue
        ids = [s[:-1] for s in steps]
        def spell(lo, hi):           # steps[lo:hi]
            out = []
            for s in steps[lo:hi]:
                b = seq.get(s[:-1], '')
                out.append(b if s[-1] == '+' else rc(b))
            return ''.join(out).upper()
        # locate each bubble's source and sink, in walk order, monotonically
        marks, cur, ok = [], 0, True
        for _, src, snk in order:
            try:
                i = ids.index(src, cur)
                j = ids.index(snk, i + 1)
            except ValueError:
                marks.append(None); continue
            marks.append((i, j)); cur = j
        prev_end = 0
        for bi, m in enumerate(marks):
            if m is None:
                print('%s\t%d\tflank\tABSENT\t' % (name, 2 * bi))
                print('%s\t%d\tbubble\tABSENT\t' % (name, 2 * bi + 1))
                continue
            i, j = m
            fl = spell(prev_end, i + 1)          # flank up to and including the source
            al = spell(i + 1, j)                 # the allele: strictly inside
            for idx, kind, s in ((2 * bi, 'flank', fl), (2 * bi + 1, 'bubble', al)):
                print('%s\t%d\t%s\t%d\t%s' % (name, idx, kind, len(s),
                                              hashlib.md5(s.encode()).hexdigest()))
            prev_end = j
        tail = spell(prev_end, len(steps))
        print('%s\t%d\tflank\t%d\t%s' % (name, 2 * len(order), len(tail),
                                         hashlib.md5(tail.encode()).hexdigest()))

if __name__ == '__main__':
    main()
