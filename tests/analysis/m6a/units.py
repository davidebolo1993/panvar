#!/usr/bin/env python3
"""M6-A EVALUATION ONLY: repeat-unit copies per target-block allele and per truth haplotype, by M1-B's declared rule
(tests/analysis/m1b/structure.py, constants and anchor definition unchanged): unit source S0 = the shortest panel allele
>= 10,000 bp carried by a template; anchors = canonical 31-mers every 25 bp of S0, unique in S0 and present in >= 90% of
carried alleles; informative anchors vary in count (max >= 2); copies = median count over informative anchors; < 50
informative anchors -> LOW_INFORMATION. Never used for candidate membership or scoring.

  units.py <block-sequences prefix> <block id> <out tsv>
"""
import csv
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "m1b"))
import structure as S  # noqa: E402

prefix, block, out = sys.argv[1:4]
alleles = S.fasta(f"{prefix}.{block}.alleles.fa")
truth = S.fasta(f"{prefix}.{block}.truth.fa")
carried = {r["name"]: r["templates"] for r in csv.DictReader(open(f"{prefix}.block_sequences.tsv"), delimiter="\t") if r["block_id"] == block and r["source"] == "panel"}
used = [n for n in alleles if carried.get(n, "-") != "-"]
source = min((n for n in used if len(alleles[n]) >= S.MIN_SOURCE_BP), key=lambda n: (len(alleles[n]), n))
s0 = alleles[source]
pos_of, seen = {}, {}
for i in range(len(s0) - S.K + 1):
    k = S.canon(s0[i:i + S.K])
    seen[k] = seen.get(k, 0) + 1
    if i % S.STRIDE == 0:
        pos_of.setdefault(k, i)
anchors = {k: p for k, p in pos_of.items() if seen[k] == 1 and "N" not in k}
table = {n: S.counts(alleles[n], anchors) for n in alleles}
table.update({n: S.counts(s, anchors) for n, s in truth.items()})
conserved = [k for k in anchors if sum(table[n][k] >= 1 for n in used) >= S.CONSERVED * len(used)]
informative = [k for k in conserved if len({table[n][k] for n in used}) > 1 and max(table[n][k] for n in used) >= 2]
rows = ["name\tsource\tbp\tunit_copies\tsupport"]
for n in list(alleles) + list(truth):
    if len(informative) < S.MIN_INFORMATIVE:
        c, sup = "LOW_INFORMATION", 0.0
    else:
        vals = [table[n][k] for k in informative]
        c = statistics.median_low(vals)
        sup = sum(v == c for v in vals) / len(vals)
    rows.append(f"{n}\t{'truth' if n in truth else 'panel'}\t{len(alleles.get(n, truth.get(n, '')))}\t{c}\t{sup:.3f}")
open(out, "w").write("\n".join(rows) + f"\n# unit_source {source} ({len(s0)} bp); informative anchors {len(informative)}\n")
print(f"unit source {source} ({len(s0)} bp), informative anchors {len(informative)}; truth: " +
      ", ".join(f"{n} {r.split(chr(9))[3]}" for n, r in zip(list(truth), rows[-len(truth):])))
