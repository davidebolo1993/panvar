#!/usr/bin/env python3
"""rebuild_cases.py - write the small GFAs module_fold_stats.sh checks rebuild's repair and fold on.

Usage: rebuild_cases.py <out_dir>   (one <out_dir>/<case>.gfa per case; sequences are seeded, so the files
are the same on every run)

  unit2999, unit3000  a tandem unit just under and at the 3 kb module floor, in 1, 2 and 3 copies
  equal0, equal1000   two arrays, U and V, that change in opposite directions, so the haplotypes' lengths
                      match (equal1000: V's extra copies sit right after the reference's copy)
  diverged            the reference's two copies 92.5% identical, under the 95% that finds an array
  second_pass         hx starts 1.5 kb before a module and lacks one copy: the generator drops the deletion,
                      the fold does not take hx, and only a second pass writes the copy loss in
  inflate             a 1.5 kb deletion the generator drops, beside 100 single-base differences
  endins              a 3 kb insertion 1.5 kb before the haplotypes' end
  overrun_nospan      h1 lies outside a module and its walk runs past its end: a fold elsewhere must not pass it
  del_flat            two haplotypes delete 5 kb 1.5 kb before the end; nothing else varies
  adjacent            a 3 kb array U and a 20 kb array W 3 kb apart, changing in opposite directions: no node
                      between them can be cut, so U's cut runs into W's first copy
  adjacent_snp        adjacent, with each haplotype's own 3 SNPs in the spacer's last 150 bp, which W's cut takes in
  adjacent600         adjacent with a 600 bp spacer: W's anchor moves past U's copies, so W's cut would take them in
  fragment            a 4 kb array U and a 5 kb array V 10 kb apart; the reference holds U once and, 10 kb past V,
                      a 1.5 kb fragment of U at 97% identity
  fragment_half       fragment with a 2.1 kb fragment, over half the unit
  fragment_full       fragment with a full copy of U at 97% identity
  fragment_up         a 2.1 kb fragment of U 10 kb before U instead of past V
  fragment_long       fragment_half with h1 holding U 10 times, an array longer than the distance to the fragment
  fragment_pair       fragment with two full copies of U at 97% identity in place of the fragment
"""
import os
import random
import sys


def dna(n, seed):
    r = random.Random(seed)
    return ''.join(r.choice('ACGT') for _ in range(n))


def mutate(s, n, seed, lo=0, hi=None):
    # n single-base substitutions at random positions in [lo, hi)
    r = random.Random(seed)
    s = list(s)
    for i in r.sample(range(lo, len(s) if hi is None else hi), n):
        s[i] = r.choice([c for c in 'ACGT' if c != s[i]])
    return ''.join(s)


def write(path, haplotypes):
    with open(path, 'w') as out:
        out.write('H\tVN:Z:1.0\n')
        for i, (_, s) in enumerate(haplotypes, 1):
            out.write(f'S\t{i}\t{s}\n')
        for i, (name, _) in enumerate(haplotypes, 1):
            out.write(f'P\t{name}\t{i}+\t*\n')


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    for unit in (2999, 3000):
        a, u, b = dna(10000, 1), dna(unit, 2), dna(10000, 3)
        write(f'{out}/unit{unit}.gfa', [('ref', a + u * 2 + b), ('h1', a + u + b), ('h2', a + u * 3 + b)])
    for extra in (0, 1000):
        a, u, b, v, c, x = dna(30000, 1), dna(8000, 2), dna(30000, 3), dna(8000, 4), dna(30000, 5), dna(1000, 6)
        write(f'{out}/equal{extra}.gfa', [('ref', a + u * 3 + b + v + c),
                                          ('alt', a + u + b + v * 3 + (x if extra else '') + c)])
    a, u, b = dna(15000, 1), dna(4000, 2), dna(15000, 3)
    u2 = mutate(u, 300, 4)   # 300 of 4000 bases differ: 92.5% identity
    write(f'{out}/diverged.gfa', [('ref', a + u + u2 + b), ('h1', a + u + b), ('h2', a + u + u2 + u + b),
                                  ('h3', a + u + u2 + b)])
    # Drawn from one stream: which way the generator maps hx back depends on these exact sequences.
    rng = random.Random(11)
    def seq(n): return ''.join(rng.choice('ACGT') for _ in range(n))
    def mut(s, k):
        s = list(s)
        for i in rng.sample(range(len(s)), k): s[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[s[i]]
        return ''.join(s)
    f0, u, f2 = seq(15000), seq(8000), seq(10000)
    ins = seq(400)
    write(f'{out}/second_pass.gfa', [('ref', f0 + u + u + f2), ('h1', f0 + u + u + u + f2),
                                     ('h2', mut(f0, 20) + u + u + f2[:5000] + ins + f2[5000:]),
                                     ('h3', f0 + u + u + mut(f2, 20)), ('hx', f0[-1500:] + u + f2)])
    r = dna(50000, 71)
    alt = r[:1300] + r[2800:]
    alt = mutate(alt, 30, 72, 1300, 1600)    # a dense stretch beside the deletion
    alt = mutate(alt, 100, 73, 5000, 35000)  # and 100 single-base differences the rebuild collapses
    write(f'{out}/inflate.gfa', [('ref', r), ('alt', alt), ('branch', r[:25000] + dna(500, 76) + r[25000:])])
    a, b, ins = dna(40000, 81), dna(1500, 82), dna(3000, 83)
    write(f'{out}/endins.gfa', [('ref', a + b), ('h1', a + ins + b), ('h2', mutate(a, 20, 84) + ins + b),
                                ('h3', mutate(a, 40, 85) + b)])
    a, u, c, x, bl = dna(15000, 1), dna(8000, 2), dna(2000, 4), dna(100, 5), dna(40000, 13)
    write(f'{out}/overrun_nospan.gfa', [('ref', a + u + u + bl + c), ('h1', bl + c[:1400] + x),
                                        ('h2', a + u + bl + c), ('h3', a + u + u + u + bl + c)])
    a, d, b = dna(40000, 81), dna(5000, 82), dna(1500, 83)
    write(f'{out}/del_flat.gfa', [('ref', a + d + b), ('h1', a + b), ('h2', mutate(a, 20, 84) + b),
                                  ('h3', mutate(a, 40, 85) + d + b),
                                  ('br', a[:20000] + dna(500, 99) + a[20000:] + d + b)])
    a, u, s, w, c = dna(20000, 91), dna(3000, 92), dna(3000, 93), dna(20000, 94), dna(20000, 95)
    write(f'{out}/adjacent.gfa', [('ref', a + u * 2 + s + w * 2 + c), ('h1', a + u + s + w * 3 + c),
                                  ('h2', a + u * 3 + s + w + c), ('h3', a + u * 2 + s + w * 2 + c),
                                  ('h4', a + u * 3 + s + w * 3 + c)])
    def tail(seed): return mutate(s, 3, seed, 2850, 3000)
    write(f'{out}/adjacent_snp.gfa', [('ref', a + u * 2 + s + w * 2 + c), ('h1', a + u + tail(1) + w * 3 + c),
                                      ('h2', a + u * 3 + tail(2) + w + c), ('h3', a + u * 2 + tail(3) + w * 2 + c),
                                      ('h4', a + u * 3 + tail(4) + w * 3 + c)])
    s = dna(600, 93)
    write(f'{out}/adjacent600.gfa', [('ref', a + u * 2 + s + w * 2 + c), ('h1', a + u + s + w * 3 + c),
                                     ('h2', a + u * 3 + s + w + c), ('h3', a + u * 2 + s + w * 2 + c),
                                     ('h4', a + u * 3 + s + w * 3 + c)])
    a, u, b1, v, b2, c = dna(20000, 1), dna(4000, 2), dna(10000, 3), dna(5000, 4), dna(10000, 5), dna(20000, 6)
    for name, n in (('fragment', 1500), ('fragment_half', 2100), ('fragment_full', 4000)):
        frag = mutate(u[:n], n * 3 // 100, 7)
        write(f'{out}/{name}.gfa', [('ref', a + u + b1 + v * 2 + b2 + frag + c),
                                    ('h1', a + u * 3 + b1 + v * 2 + b2 + frag + c),
                                    ('h2', a + u + b1 + v * 3 + b2 + frag + c),
                                    ('h3', a + u * 2 + b1 + v + b2 + frag + c)])
    frag = mutate(u[:2100], 63, 7)
    write(f'{out}/fragment_long.gfa', [('ref', a + u + b1 + v * 2 + b2 + frag + c),
                                       ('h1', a + u * 10 + b1 + v * 2 + b2 + frag + c),
                                       ('h2', a + u + b1 + v * 3 + b2 + frag + c),
                                       ('h3', a + u * 2 + b1 + v + b2 + frag + c)])
    d1, d2 = mutate(u, 120, 8), mutate(u, 120, 9)
    write(f'{out}/fragment_pair.gfa', [('ref', a + u + b1 + v * 2 + b2 + d1 + d2 + c),
                                       ('h1', a + u * 3 + b1 + v * 2 + b2 + d1 + d2 + c),
                                       ('h2', a + u + b1 + v * 3 + b2 + d1 + d2 + c),
                                       ('h3', a + u * 2 + b1 + v + b2 + d1 + d2 + c)])
    frag, s = mutate(u[:2100], 63, 7), dna(10000, 8)
    write(f'{out}/fragment_up.gfa', [('ref', a + frag + s + u + b1 + v * 2 + b2 + c),
                                     ('h1', a + frag + s + u * 3 + b1 + v * 2 + b2 + c),
                                     ('h2', a + frag + s + u + b1 + v * 3 + b2 + c),
                                     ('h3', a + frag + s + u * 2 + b1 + v + b2 + c)])


if __name__ == '__main__':
    main()
