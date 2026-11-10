#!/usr/bin/env python3
"""Count paths that traverse a bubble more than once.

call takes one walk per haplotype per bubble: bubble_steps keeps the WIDEST interval holding every
interior occurrence and module_span runs from the first source to the last sink. A haplotype that
traverses a bubble twice would therefore get one walk that includes the outside sequence between the
two traversals, and UNIT, MODULE and CNBP would count the interior twice. This is the guard
that says how often that happens on a graph; call itself flags it per bubble (MULTI_CROSSING).

A traversal is a maximal stretch of consecutive path steps that starts at one boundary node, ends at
the OTHER boundary node and has only bubble-interior nodes in between (none at all is allowed:
source,sink is a traversal). Either direction counts, and node orientation is ignored. A path with
two or more traversals of the same bubble is a multi-crosser and is reported. A path that visits a
boundary node more than twice without a second traversal (a touch, a hairpin back out, a detour
outside the bubble) is counted separately, for context only.

  crossings.py <graph.gfa> <bubbles.csv> [--label L] [--items out.tsv]

--items rows: CROSSING, bubble_id, '-', '-', 'paths', 1, detail, pin (a hash of every flagged path
with its traversal count). Exit 1 if any path is a multi-crosser.
"""
import argparse
import csv
import hashlib
import re
import sys

csv.field_size_limit(10 ** 9)


def load_paths(gfa):
    paths = {}
    with open(gfa) as fh:
        for line in fh:
            if line.startswith('P\t'):
                f = line.rstrip('\n').split('\t')
                paths[f[1]] = [s[:-1] for s in f[2].split(',') if s]
            elif line.startswith('W\t'):
                f = line.rstrip('\n').split('\t')
                paths['#'.join(f[1:4])] = re.findall(r'[<>]([^<>]+)', f[6])
    return paths


def traversals(steps, boundary, interior):
    """bubble_id -> (traversal count, boundary-node visits) along one path.

    The stretches walked from one bubble's boundary occurrences are disjoint (each ends at the first
    non-interior step), so this is one pass per path and bubble nesting level."""
    count, visits = {}, {}
    n = len(steps)
    for i, s in enumerate(steps):
        for bid, other in boundary.get(s, ()):
            visits[bid] = visits.get(bid, 0) + 1
            ins = interior[bid]
            j = i + 1
            while j < n and steps[j] in ins:
                j += 1
            if j < n and steps[j] == other:
                count[bid] = count.get(bid, 0) + 1
    return count, visits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gfa')
    ap.add_argument('bubbles')
    ap.add_argument('--label', default='')
    ap.add_argument('--items')
    a = ap.parse_args()
    paths = load_paths(a.gfa)
    # boundary node -> [(bubble_id, the other boundary)]; one scan per path instead of one per bubble
    boundary, interior, order = {}, {}, []
    with open(a.bubbles) as fh:
        for r in csv.DictReader(fh):
            bid, src, snk = r['bubble_id'], r['source'], r['sink']
            order.append(bid)
            if src == snk:
                continue
            boundary.setdefault(src, []).append((bid, snk))
            boundary.setdefault(snk, []).append((bid, src))
            interior[bid] = set(x for x in r.get('inside_nodes', '').split(';') if x)
    bad, touch = {}, {}
    total = 0
    for name, steps in paths.items():
        count, visits = traversals(steps, boundary, interior)
        total += sum(count.values())
        for bid, c in count.items():
            if c >= 2:
                bad.setdefault(bid, []).append((name, c))
        for bid, v in visits.items():
            if count.get(bid, 0) < 2 and v > 2:
                touch[bid] = touch.get(bid, 0) + 1
    items, multi_paths = [], set()
    for bid in order:
        if bid not in bad:
            continue
        b = sorted(bad[bid])
        multi_paths.update(n for n, _ in b)
        ex = ', '.join(f'{n} ({c} traversals)' for n, c in b[:3])
        items.append(('CROSSING', bid, '-', '-', 'paths', 1,
                      f'{len(b)} path(s) traverse more than once, e.g. {ex}',
                      hashlib.sha1(repr(b).encode()).hexdigest()[:12]))
    print(f'{a.label} crossings: bubbles={len(order)} paths={len(paths)} traversals={total} '
          f'bubbles_with_multi_crossing={len(items)} paths_multi_crossing={len(multi_paths)} '
          f'boundary_revisits_single_traversal={sum(touch.values())}')
    for x in items:
        print(f'  bubble {x[1]}: {x[6]}')
    if a.items:
        with open(a.items, 'w') as fo:
            for x in items:
                fo.write('\t'.join(map(str, x)) + '\n')
    sys.exit(1 if items else 0)


if __name__ == '__main__':
    main()
