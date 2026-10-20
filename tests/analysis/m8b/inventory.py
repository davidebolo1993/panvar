#!/usr/bin/env python3
"""M8-B corrected development/evaluation inventory and frozen per-locus splits (before any parameter scoring).

Development (previously used) donors = union of
  (a) authoritative sample-level records under tests/results (depth <= 5): directories named as a donor (per-sample
      reads/runs/LOO/LZO folders), LOO model files (<D>.loo.gfa, loo.<D>.gfa), mosaic plans (stage.log 'mosaic-truth A x B'),
      simulations (simulation.tsv 'donor', simulate.log '--donor'), campaign signatures ('donors' rows);
  (b) M8-A's DEVELOPMENT_DONORS.tsv (memory notes, tests/analysis, mosaic plans, M6-F evaluation families);
expanded to IGSR-pedigree families. Required checks: NA19705 and HG03520 must be development.
Evaluability rule (frozen): a locus is evaluated only with >= 4 uninspected families; ACOT evaluates the first 30
uninspected families by sha256("M8-B:" + family); otherwise the locus has no outer evaluation (INFEASIBLE) and only
training profiles are reported. Templates: the full models' active panels (M8-A support.full block_sequences).

  inventory.py <out dir>
"""
import csv
import hashlib
import os
import re
import sys
from collections import defaultdict

out = sys.argv[1]
R = "/Users/davide.bolognini/panvar/tests/results"
WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
PED = f"{R}/m6f_copying_predictive_20261019/inputs/20130606_g1k_3202_samples_ped_population.txt"
M8A = f"{R}/m8a_copying_calibration_20261020/support.full"
DONOR = re.compile(r"^(HG|NA)\d{5}$")
src = defaultdict(set)
for root, dirs, files in os.walk(R):
    depth = root[len(R):].count(os.sep)
    if depth >= 5:
        dirs[:] = []
    pkg = root[len(R) + 1:].split(os.sep)[0] if root != R else ""
    if pkg.startswith(("m8a_", "m8b_")):
        continue
    for x in dirs:
        if DONOR.match(x):
            src[x].add(f"per-sample directory ({pkg})")
    for f in files:
        m = re.match(r"^(?:loo\.)?((?:HG|NA)\d{5})(?:\.loo)?\.gfa$", f)
        if m and ("loo" in f):
            src[m.group(1)].add(f"LOO model ({pkg})")
        p = os.path.join(root, f)
        if f == "stage.log":
            for a, b in re.findall(r"mosaic-truth ([A-Z0-9]+) x ([A-Z0-9]+)", open(p, errors="ignore").read()):
                src[a].add(f"mosaic donor ({pkg})")
                src[b].add(f"mosaic donor ({pkg})")
        elif f == "simulation.tsv":
            for line in open(p, errors="ignore"):
                if line.startswith("donor\t") and DONOR.match(line.split("\t")[1].strip()):
                    src[line.split("\t")[1].strip()].add(f"simulation ({pkg})")
        elif f == "simulate.log":
            for d in re.findall(r"--donor ((?:HG|NA)\d{5})", open(p, errors="ignore").read()):
                src[d].add(f"simulation ({pkg})")
        elif f == "signature.tsv":
            for line in open(p, errors="ignore"):
                if line.startswith("donors\t"):
                    for d in line.split("\t")[1].split():
                        if DONOR.match(d):
                            src[d].add(f"campaign signature ({pkg})")
for line in open(f"{WT}/tests/analysis/m8a/DEVELOPMENT_DONORS.tsv"):
    if not line.startswith(("#", "donor\t")):
        src[line.split("\t")[0]].add("M8-A list")
ped = [l.split() for l in open(PED)][1:]
parent = {}


def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


for fam, s, fa, mo, *_ in ped:
    parent[find(s)] = find("FAM:" + fam)
    for p_ in (fa, mo):
        if p_ != "0":
            parent[find(s)] = find(p_)
famids = defaultdict(set)
for fam, s, *_ in ped:
    famids[find(s)].add(fam)
root_of = lambda d: find(d) if d in parent else "SELF:" + d
fam_name = lambda r_: "+".join(sorted(famids.get(r_, []))) or r_.replace("SELF:", "")
assert "NA19705" in src and "HG03520" in src, "required development donors missing"
os.makedirs(out, exist_ok=True)
with open(f"{out}/development_donors.tsv", "w") as f:
    f.write("donor\tsources\n")
    for d in sorted(src):
        f.write(f"{d}\t{'; '.join(sorted(src[d]))}\n")
summary = []
for locus in ("C4", "CYP2D6", "ACOT"):
    tab = [l.rstrip("\n").split("\t") for l in open(f"{M8A}/{locus}.block_sequences.tsv")][1:]
    templates = sorted({t for r in tab if r[1] == "panel" for t in r[5].split(",") if t and t != "-"})
    roots = defaultdict(list)
    for t in templates:
        roots[root_of(t.split("#")[0])].append(t)
    status = {}
    for r_, ts in roots.items():
        dns = {t.split("#")[0] for t in ts}
        dev_dns = {dn for dn in dns if dn in src}
        # whole pedigree family is development if ANY member (in the panel or not) was used
        fam_members = {s for fam_, s, *_ in ped if find(s) == r_} | dns
        used = {dn for dn in fam_members if dn in src}
        status[r_] = ("REFERENCE_GENOME" if any(dn.lower().startswith(("grch38", "chm13")) for dn in dns) else
                      "DEVELOPMENT" if used else "UNINSPECTED")
    un = sorted((r_ for r_ in roots if status[r_] == "UNINSPECTED"), key=lambda r_: hashlib.sha256(("M8-B:" + fam_name(r_)).encode()).hexdigest())
    if locus == "ACOT":
        evaluated = un[:30]
        verdict = "EVALUATED (30 by hash)" if len(un) >= 30 else ("EVALUATED" if len(un) >= 4 else "INFEASIBLE")
    else:
        evaluated = un if len(un) >= 4 else []
        verdict = "EVALUATED" if len(un) >= 4 else "INFEASIBLE (< 4 uninspected families)"
    with open(f"{out}/split.{locus}.tsv", "w") as f:
        for r_, ts in sorted(roots.items(), key=lambda x: fam_name(x[0])):
            for t in ts:
                f.write(f"{t}\t{fam_name(r_)}\t{'eval' if r_ in evaluated else 'train'}\t{status[r_]}\n")
    summary.append(dict(locus=locus, templates=len(templates), families=len(roots), development=sum(1 for v in status.values() if v == "DEVELOPMENT"),
                        reference_genome=sum(1 for v in status.values() if v == "REFERENCE_GENOME"), uninspected=len(un),
                        uninspected_families=",".join(fam_name(r_) for r_ in un[:12]) + ("..." if len(un) > 12 else ""), evaluated_families=len(evaluated),
                        evaluation=verdict))
with open(f"{out}/inventory_summary.tsv", "w") as f:
    keys = list(summary[0])
    w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(summary)
print(open(f"{out}/inventory_summary.tsv").read())
print("development donors:", len(src), "| NA19705:", sorted(src["NA19705"])[:4], "| HG03520:", sorted(src["HG03520"])[:4])
