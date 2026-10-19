#!/usr/bin/env python3
"""M6-G prepare (panel sequences only; frozen rule, outcome-independent; no reads exist yet).

Library: the K = 2 period/long components with the largest founder weight among copies of M6-F's training-only
families (ties: more haplotypes, then smallest sha256), in M6-F's period/long reference coordinates (period_u13).
Windows: for L in {200, 600}, the first start s = 0, 50, 100, ... such that
  (a) neither component has an indel, deletion or uncertain correspondence within [s - 10, s + L + 10);
  (b) the two components differ at >= 1 position in [s, s + L);
  (c) the admitted residual sites a1 = s + floor(L/3), a2 = s + floor(2L/3) carry the same base in both components, and
      no component difference lies within 10 bp of either.
Admitted residual variation: at a1 and a2 only, the declared novel allele = the next base in ACGT (absent from both).
Weights: founder weights of the two components renormalised over the library (the restriction is explicit).
Residual divergence rate: eps_res = 1e-4 (M6-F's C1 selection in every period/long fold), distinct from the sequencing
error rate of the read model (0.002).

  prepare.py <out dir>
"""
import builtins
import hashlib
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "m6f"))
import m6f_core as core  # noqa: E402

M6F = "/Users/davide.bolognini/panvar/tests/results/m6f_copying_predictive_20261019/runs"
out = sys.argv[1]
ALLOW = {os.path.abspath(p) for p in [core.M1D_FA, f"{M6F}/prepare/copies.tsv", f"{M6F}/prepare/prepare.tsv", f"{M6F}/prepare/families.tsv", f"{M6F}/align/events.tsv"]}
_open = builtins.open


def guarded(path, mode="r", *a, **k):
    p = os.path.abspath(path)
    if "w" in mode or "a" in mode:
        assert p.startswith(os.path.abspath(out)), f"write outside the unit directory: {p}"
    else:
        assert p in ALLOW or p.startswith(os.path.abspath(out)), f"not an allowlisted input: {p}"
    return _open(path, mode, *a, **k)


builtins.open = guarded
EPS_RES, K, W = 1e-4, 2, 10
copies, refs, comps, families, evaluation = core.load(f"{M6F}/prepare", f"{M6F}/align/events.tsv")
S = "period/long"
seqs = core.fasta(core.M1D_FA)
pool = defaultdict(lambda: [0.0, set()])
for c in copies:
    if c["stratum"] == S and c["status"] == "TRAINING_ONLY" and float(c["weight"]) > 0:
        pool[c["unique_record"]][0] += float(c["weight"])
        pool[c["unique_record"]][1].add(c["template"])
lib = sorted(pool, key=lambda k: (-pool[k][0], -len(pool[k][1]), hashlib.sha256(seqs[k].encode()).hexdigest()))[:K]
wsum = sum(pool[k][0] for k in lib)
A = [comps[S][k] for k in lib]
L_ref = len(A[0].base)
BASES = "ACGT"
alt_of = {"A": "C", "C": "G", "G": "T", "T": "A"}
lines = [f"reference\t{refs[S]}", f"eps_res\t{EPS_RES}", f"components\t{','.join(lib)}",
         f"weights\t{','.join(f'{pool[k][0] / wsum:.17g}' for k in lib)}", f"founder_weights\t{','.join(f'{pool[k][0]:g}' for k in lib)}"]
ok = True
for L in (200, 600):
    found = None
    for s in range(0, L_ref - L, 50):
        lo, hi = max(0, s - W), min(L_ref, s + L + W)
        if any(a.bad[lo:hi].any() for a in A):
            continue
        diff = [x for x in range(s, s + L) if A[0].base[x] != A[1].base[x]]
        if not diff:
            continue
        a1, a2 = s + L // 3, s + (2 * L) // 3
        if A[0].base[a1] != A[1].base[a1] or A[0].base[a2] != A[1].base[a2]:
            continue
        if any(abs(d - a) <= W for d in diff for a in (a1, a2)):
            continue
        found = (s, diff, a1, a2)
        break
    if found is None:
        ok = False
        lines.append(f"window_{L}\tNONE (rule found no window)")
        continue
    s, diff, a1, a2 = found
    win = ["".join(BASES[b] for b in a.base[s:s + L]) for a in A]
    lines += [f"window_{L}_start\t{s}", f"window_{L}_component_differences\t{','.join(str(d - s) for d in diff)}",
              f"window_{L}_admitted_sites\t{a1 - s},{a2 - s}", f"window_{L}_novel_alleles\t{alt_of[win[0][a1 - s]]},{alt_of[win[0][a2 - s]]}"]
    for k, wseq in zip(lib, win):
        lines.append(f"window_{L}_seq_{k}\t{wseq}")
with open(f"{out}/library.tsv", "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(l if "_seq_" not in l else l[:60] + "..." for l in lines))
sys.exit(0 if ok else 1)
