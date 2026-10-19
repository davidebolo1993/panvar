#!/usr/bin/env python3
"""M6-B: compare the exact-grouping count-fit reruns with the M6-A outputs (same dictionary, depth, loss and catalogue).
Per cell: every ordered candidate's loss (must be bit-identical as printed), group rank and group size; the summary keys;
the reported hash-collision buckets; and the frozen evaluation selection re-derived by M6-A's select_eval.py rule. Reports
every difference; presumes no collision.

  verify_groups.py <M6-A runs> <M6-B runs> <out dir> cell...
"""
import csv
import filecmp
import os
import subprocess
import sys

a_runs, b_runs, out, cells = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l}
lines = ["cell\tcandidates\tloss_differences\tgroup_rank_differences\tgroup_size_differences\tgroups_m6a\tgroups_exact\thash_collision_buckets\t"
         "caller_group_rank_m6a\tcaller_group_rank_exact\tbest_group_size_m6a\tbest_group_size_exact\tdepth_identical"]
bad = 0
for c in cells:
    A = {(r["allele_h1"], r["allele_h2"]): r for r in T(f"{a_runs}/countfit.{c}/cf.candidates.tsv")}
    B = {(r["allele_h1"], r["allele_h2"]): r for r in T(f"{b_runs}/countfit.{c}/cf.candidates.tsv")}
    assert set(A) == set(B), c
    dl = sum(A[k]["loss"] != B[k]["loss"] for k in A)
    dr = sum(A[k]["group_rank"] != B[k]["group_rank"] for k in A)
    ds = sum(A[k]["group_size"] != B[k]["group_size"] for k in A)
    sa, sb = KV(f"{a_runs}/countfit.{c}/cf.summary.tsv"), KV(f"{b_runs}/countfit.{c}/cf.summary.tsv")
    bad += dl + dr + ds + (sa["groups"] != sb["groups"]) + (sa["caller_group_rank"] != sb["caller_group_rank"]) + (sb["hash_collision_buckets"] != "0")
    lines.append("\t".join(map(str, [c, len(A), dl, dr, ds, sa["groups"], sb["groups"], sb["hash_collision_buckets"], sa["caller_group_rank"],
                                     sb["caller_group_rank"], sa["best_group_size"], sb["best_group_size"], sa["depth_haploid_kmer_mean"] == sb["depth_haploid_kmer_mean"]])))
# the frozen evaluation selection, re-derived from the exact-grouping outputs by M6-A's own rule
sel = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "m6a", "select_eval.py")
os.makedirs(f"{out}/sel", exist_ok=True)
subprocess.run([sys.executable, sel, b_runs, f"{out}/sel", *cells], check=True, capture_output=True)
same_sel = all(filecmp.cmp(f"{out}/sel/eval_candidates.{c}.tsv", f"{a_runs}/freeze/eval_candidates.{c}.tsv", shallow=False) for c in cells)
bad += 0 if same_sel else 1
lines += ["", f"evaluation selections identical to M6-A\t{same_sel}", f"VERDICT\t{'UNCHANGED' if bad == 0 else 'CHANGED'}"]
open(f"{out}/verify.tsv", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
sys.exit(0 if bad == 0 else 1)
