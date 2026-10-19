#!/usr/bin/env python3
"""M6-F prepare (panel sequences + pedigree only; frozen with DESIGN.md).

  1. families: connected components of the IGSR 3,202-sample pedigree (FamilyID co-membership and parent links);
     donors absent from the pedigree are their own family
  2. status: TRAINING_ONLY for development families (DEVELOPMENT_DONORS.tsv, expanded to whole families), reference
     genomes (chm13v2, grch38) and the coordinate-reference donor HG00096; the others are ELIGIBLE
  3. evaluation: the first 8 eligible families by sha256("M6-F:" + family id) -- chosen here, before any scoring
  4. copies: one row per physical copy (component occurrence x carrier template) with role, length class and the
     family-aware weight (founder policy: 0 for a donor whose father and mother are both panel donors, else 1)
  5. coordinate references: per stratum (prefix; period long/short; tail long/short), the unique component with the
     largest positive weight among copies of NEVER-EVALUATED families (ties: more haplotypes, then smallest sha256)
  6. alignment pairs: every unique component of a stratum against that stratum's reference

  prepare.py <out dir>
"""
import builtins
import csv
import hashlib
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
M1D = "/Users/davide.bolognini/panvar/tests/results/m1d_unit_reachability_20261019/runs"
PED = "/Users/davide.bolognini/panvar/tests/results/m6f_copying_predictive_20261019/inputs/20130606_g1k_3202_samples_ped_population.txt"
DEV = os.path.join(HERE, "DEVELOPMENT_DONORS.tsv")
INV = f"{M1D}/inventory.a1/inventory"
BLOCK = f"{M1D}/seq.panel/seq.block_sequences.tsv"
out = sys.argv[1]
ALLOW = {os.path.abspath(p) for p in [PED, DEV, f"{INV}.tsv", f"{INV}.fa", BLOCK]}
_open = builtins.open


def guarded(path, mode="r", *a, **k):
    p = os.path.abspath(path)
    if "w" in mode or "a" in mode:
        assert p.startswith(os.path.abspath(out)), f"write outside the unit directory: {p}"
    else:
        assert p in ALLOW or p.startswith(os.path.abspath(out)), f"not an allowlisted input: {p}"
    return _open(path, mode, *a, **k)


builtins.open = guarded
LONG_BP, N_EVAL = 29000, 8
REFERENCE_GENOMES = {"chm13v2", "grch38"}
COORDINATE_DONORS = {"HG00096"}


def fasta(path):
    recs, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs[name] = []
        elif name:
            recs[name].append(line)
    return {k: "".join(v) for k, v in recs.items()}


comp = fasta(f"{INV}.fa")
rows = list(csv.DictReader(open(f"{INV}.tsv"), delimiter="\t"))
carriers = {r["name"]: [t for t in r["templates"].split(",") if t] for r in csv.DictReader(open(BLOCK), delimiter="\t")}
dev = {l.split("\t")[0] for l in open(DEV) if not l.startswith(("#", "donor\t"))}
ped = [l.split() for l in open(PED)][1:]
parent = {}


def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(a, b):
    parent[find(a)] = find(b)


parents_of = {}
for fam, s, fa, mo, *_ in ped:
    union(s, "FAM:" + fam)
    parents_of[s] = (fa, mo)
    for p in (fa, mo):
        if p != "0":
            union(s, p)
templates = sorted({t for r in rows for t in carriers[r["source_allele"]]})
donors = sorted({t.split("#")[0] for t in templates})
root = {d: (find(d) if d in parent else "SELF:" + d) for d in donors}
members = defaultdict(list)
for d in donors:
    members[root[d]].append(d)
fam_id = {}
for r_, ds in members.items():   # readable family id: the pedigree FamilyIDs of its members, else the donor
    ids = sorted({f for f, s, *_ in ped if s in ds})
    for d in ds:
        fam_id[d] = "+".join(ids) if ids else d
status, reason = {}, {}
for r_, ds in members.items():
    why = []
    if any(d in dev for d in ds):
        why.append("development:" + ",".join(sorted(d for d in ds if d in dev)))
    if any(d in REFERENCE_GENOMES for d in ds):
        why.append("reference_genome")
    if any(d in COORDINATE_DONORS for d in ds):
        why.append("coordinate_reference_donor")
    for d in ds:
        status[d] = "TRAINING_ONLY" if why else "ELIGIBLE"
        reason[d] = ";".join(why) if why else "-"
eligible = sorted({fam_id[d] for d in donors if status[d] == "ELIGIBLE"})
order = sorted(eligible, key=lambda f: hashlib.sha256(("M6-F:" + f).encode()).hexdigest())
evaluation = order[:N_EVAL]
for d in donors:
    if fam_id[d] in evaluation:
        status[d] = "EVALUATION"
in_panel = set(donors)
weight = {d: 0.0 if d in parents_of and all(p in in_panel for p in parents_of[d]) and "0" not in parents_of[d] else 1.0 for d in donors}

cls = lambda name, role: "-" if role == "prefix" else ("long" if len(comp[name]) >= LONG_BP else "short")
copies = []
for r in rows:
    pos = re.sub(r"^allele_\d+\.", "", r["component"])
    for t in carriers[r["source_allele"]]:
        d = t.split("#")[0]
        copies.append(dict(copy=f"{t}|{pos}", template=t, donor=d, family=fam_id[d], status=status[d], allele=r["source_allele"], position=pos,
                           role=r["type"], length_class=cls(r["unique_record"], r["type"]), stratum=f"{r['type']}/{cls(r['unique_record'], r['type'])}",
                           unique_record=r["unique_record"], weight=weight[d]))
with open(f"{out}/copies.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(copies[0]), delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(copies)
with open(f"{out}/families.tsv", "w") as f:
    f.write("donor\tfamily\tin_pedigree\tfather\tmother\tstatus\treason\tweight\teligible_rank\n")
    for d in donors:
        fa, mo = parents_of.get(d, ("-", "-"))
        rank = order.index(fam_id[d]) + 1 if fam_id[d] in order else "-"
        f.write(f"{d}\t{fam_id[d]}\t{int(d in parents_of)}\t{fa}\t{mo}\t{status[d]}\t{reason[d]}\t{weight[d]:g}\t{rank}\n")

strata = sorted({c["stratum"] for c in copies})
refs = {}
for s in strata:
    pool = defaultdict(lambda: [0.0, set()])
    for c in copies:
        if c["stratum"] == s and c["status"] == "TRAINING_ONLY" and c["weight"] > 0:
            pool[c["unique_record"]][0] += c["weight"]
            pool[c["unique_record"]][1].add(c["template"])
    refs[s] = min(pool, key=lambda k: (-pool[k][0], -len(pool[k][1]), hashlib.sha256(comp[k].encode()).hexdigest()))
fa = f"{INV}.fa"
with open(f"{out}/pairs.tsv", "w") as f:
    for s in strata:
        for k in sorted({c["unique_record"] for c in copies if c["stratum"] == s}, key=lambda k: (k.split("_u")[0], int(k.split("_u")[1]))):
            f.write(f"{s}:{k}\t{fa}\t{k}\t{fa}\t{refs[s]}\n")

checks = [
    ("pedigree confirms the three HGSVC trios (SH032, PR05, Y117)", all(any(f == x for f, *_ in ped) for x in ("SH032", "PR05", "Y117"))),
    ("trio children carry both parents in the pedigree", all(parents_of[c] == p for c, p in [("HG00514", ("HG00512", "HG00513")), ("HG00733", ("HG00731", "HG00732")), ("NA19240", ("NA19239", "NA19238"))])),
    ("at most 8 evaluation families, all eligible", 0 < len(evaluation) <= N_EVAL),
    ("every reference is carried by a training-only family", all(any(c["unique_record"] == refs[s] and c["status"] == "TRAINING_ONLY" for c in copies) for s in strata)),
    ("reference genomes and HG00096 are training-only", all(status[d] == "TRAINING_ONLY" for d in donors if d in REFERENCE_GENOMES | COORDINATE_DONORS)),
]
with open(f"{out}/prepare.tsv", "w") as f:
    for n, ok in checks:
        f.write(f"check\t{n}\t{'PASS' if ok else 'FAIL'}\n")
    f.write(f"families\t{len(members)} ({sum(1 for ds in members.values() if len(ds) > 1)} multi-donor: "
            + "; ".join(f"{fam_id[ds[0]]}={','.join(ds)}" for ds in members.values() if len(ds) > 1) + ")\n")
    f.write(f"not_in_pedigree\t{','.join(d for d in donors if d not in parents_of)}\n")
    f.write(f"eligible_families\t{len(eligible)}: {','.join(order)}\n")
    f.write(f"evaluation_families\t{','.join(evaluation)}\n")
    f.write(f"zero_weight_donors\t{','.join(d for d in donors if weight[d] == 0)}\n")
    for s in strata:
        cs = [c for c in copies if c["stratum"] == s]
        f.write(f"stratum\t{s}\treference {refs[s]} ({len(comp[refs[s]])} bp); copies {len(cs)}; unique {len({c['unique_record'] for c in cs})}; "
                f"evaluation copies {sum(1 for c in cs if c['status'] == 'EVALUATION')} in {len({c['family'] for c in cs if c['status'] == 'EVALUATION'})} families\n")
        f.write(f"reference\t{s}\t{refs[s]}\n")
        f.write(f"reference_multiplicity\t{s}\tevaluation copies identical to the reference (genuine multiplicity, kept): "
                f"{sum(1 for c in cs if c['status'] == 'EVALUATION' and c['unique_record'] == refs[s])}\n")
print(open(f"{out}/prepare.tsv").read())
sys.exit(0 if all(ok for _, ok in checks) else 1)
