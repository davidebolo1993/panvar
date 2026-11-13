#!/usr/bin/env python3
"""Checks the order of the pairs whose mates were fetched at the end of a coordinate-sorted file: they close
the table, first those found through the index, in the order their kept reads came (the later record of the
pair), then those found by a second pass, in the order of their mates' records.

    fetch_order.py <table.tsv> <coordinate_sorted.sam> <found through the index> <found by a second pass>
Exits 1 when the order differs.
"""
import sys

table, sam, n_index, n_second = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
first, last = {}, {}
for i, line in enumerate(open(sam)):
    if line.startswith("@"):
        continue
    f = line.split("\t")
    if int(f[1]) & (256 | 512 | 2048):   # secondary, QC-failed, supplementary: skipped
        continue
    first.setdefault(f[0], i)
    last[f[0]] = i
rows = [line.split("\t")[0] for line in open(table)][1:]
tail = rows[len(rows) - n_index - n_second:]
kept = [last[n] for n in tail[:n_index]]
mates = [first[n] for n in tail[n_index:]]
if kept != sorted(kept) or mates != sorted(mates):
    print(f"fetched pairs out of order: through the index {kept}, by a second pass {mates}")
    sys.exit(1)
