#!/usr/bin/env python3
"""M3-B Amendment 1: split a frozen catalogue into two parts deterministically (catalogue order, first half / second half
of the non-incumbent rows). The incumbent row, with its expected total, heads both parts so each part must reproduce it.

  split_catalogue.py <catalogue.tsv> <out prefix>
"""
import sys

cat, prefix = sys.argv[1:3]
lines = [l for l in open(cat) if l.strip()]
inc, rest = lines[0], lines[1:]
assert inc.startswith("incumbent\t")
half = (len(rest) + 1) // 2
for i, chunk in enumerate((rest[:half], rest[half:]), start=1):
    open(f"{prefix}{i}.tsv", "w").write(inc + "".join(chunk))
    print(f"part{i}: {len(chunk)} candidates + incumbent")
