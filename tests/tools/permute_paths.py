#!/usr/bin/env python3
"""Rewrite a GFA with its P/W lines in another order; every other line keeps its place before them.

  permute_paths.py <in.gfa> <out.gfa> <rev|SEED>

'rev' reverses the path order, an integer seeds a shuffle. call must not depend on path order
beyond the sample-column order of its VCFs.
"""
import random
import sys

src, dst, mode = sys.argv[1:4]
other, paths = [], []
with open(src) as fh:
    for line in fh:
        (paths if line[:2] in ('P\t', 'W\t') else other).append(line)
if mode == 'rev':
    paths = paths[::-1]
else:
    random.Random(int(mode)).shuffle(paths)
with open(dst, 'w') as fo:
    fo.writelines(other + paths)
