#!/usr/bin/env python3
"""M6-A frozen, truth-free selection of the candidates the evaluation scores with the fragment model (score-set) after the
count-model freeze: the caller's prediction; every member of the count-best group, up to 20 in file order (a larger group
is reported as such: its remaining members are structurally evaluated from units but not fragment-scored); and the first
member of count groups 2-5. Haplotype paths = the frozen context with the target block's alleles substituted. The
caller's prediction carries its saved total as an expectation (it must reproduce).
Also writes the freeze digest over every cell's count outputs.

  select_eval.py <runs dir> <out dir> cell...
"""
import csv
import hashlib
import sys

runs, out, cells = sys.argv[1], sys.argv[2], sys.argv[3:]
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l}
h = hashlib.sha256()
for c in cells:
    d = f"{runs}/countfit.{c}"
    h.update(open(f"{d}/cf.candidates.tsv", "rb").read())
    h.update(open(f"{d}/cf.summary.tsv", "rb").read())
    s = KV(f"{d}/cf.summary.tsv")
    tb = int(s["target_block_index"])
    ctx = [s["context_hap1"].split(","), s["context_hap2"].split(",")]
    rows = list(csv.DictReader(open(f"{d}/cf.candidates.tsv"), delimiter="\t"))
    pick, seen = [], set()
    def add(r, label):
        key = (r["allele_h1"], r["allele_h2"])
        if key in seen:
            return
        seen.add(key)
        pick.append((label, r))
    caller = next(r for r in rows if r["is_caller"] == "true")
    best = [r for r in rows if r["group_rank"] == "1"]
    for i, r in enumerate(best[:20]):
        add(r, f"g1_m{i + 1}")
    for g in range(2, 6):
        first = next((r for r in rows if r["group_rank"] == str(g)), None)
        if first:
            add(first, f"g{g}_m1")
    add(caller, "caller")
    with open(f"{out}/eval_candidates.{c}.tsv", "w") as f:
        for label, r in pick:
            h1, h2 = list(ctx[0]), list(ctx[1])
            h1[tb], h2[tb] = r["allele_h1"], r["allele_h2"]
            is_caller = r["is_caller"] == "true"
            score = s["context_source"].rsplit("(score ", 1)[1].rstrip(")") if is_caller else "NA"
            f.write(f"{label}{'_caller' if is_caller and label != 'caller' else ''}\t{','.join(h1)}\t{','.join(h2)}\t{score}\n")
    print(f"{c}: best group size {len(best)} ({min(len(best), 20)} fragment-scored), {len(pick)} candidates selected")
open(f"{out}/count_selection.sha256", "w").write(h.hexdigest() + "\n")
print("count-model freeze digest", h.hexdigest())
