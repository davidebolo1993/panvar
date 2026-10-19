#!/usr/bin/env python3
"""M6-E prepare (panel only): verify the M1-D inventory bindings, write the component/carrier table, choose the
reference components by the frozen panel-only rule and write the alignment pairs.

Reads ONLY the allowlisted inputs below (M1-D's panel inventory and the panel block sequences). Truth-selected
constructions, oracle patches, reads and evaluation outcomes are never opened (guarded).

  prepare.py <out dir>
"""
import builtins
import csv
import hashlib
import os
import re
import sys
from collections import Counter, defaultdict

M1D = "/Users/davide.bolognini/panvar/tests/results/m1d_unit_reachability_20261019/runs"
FROZEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "m1d", "INVENTORY_FROZEN.tsv")
INV = f"{M1D}/inventory.a1/inventory"
SEQ = f"{M1D}/seq.panel"
ALLOW = {os.path.abspath(p) for p in [FROZEN, f"{INV}.tsv", f"{INV}.fa", f"{INV}.anchors.tsv", f"{INV}.log.tsv",
                                       f"{SEQ}/seq.sv:7.alleles.fa", f"{SEQ}/seq.block_sequences.tsv", f"{SEQ}/run.log"]}
out = sys.argv[1]
_open = builtins.open


def guarded(path, mode="r", *a, **k):   # panel-only guard: inputs must be allowlisted; writes only under <out>
    p = os.path.abspath(path)
    if "w" in mode or "a" in mode:
        assert p.startswith(os.path.abspath(out)), f"write outside the unit directory: {p}"
    else:
        assert p in ALLOW or p.startswith(os.path.abspath(out)), f"not an allowlisted panel-only input: {p}"
    return _open(path, mode, *a, **k)


builtins.open = guarded
LONG_BP = 29000   # length classes: short (< 29 kb) vs long; the inventory has no component between 26.5 and 32.5 kb


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


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


checks = []
frozen = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(FROZEN) if l.startswith("sha256.")}
for ext in ["tsv", "fa", "anchors.tsv", "log.tsv"]:
    checks.append((f"inventory.{ext} digest equals INVENTORY_FROZEN", sha(f"{INV}.{ext}") == frozen[f"sha256.{ext}"]))
runlog = open(f"{SEQ}/run.log").read()
checks.append(("block sequences written from the NA19650 LOO model", "graphs/NA19650.loo.gfa" in runlog and "model fingerprint 14851130f453a071" in runlog))
rows = list(csv.DictReader(open(f"{INV}.tsv"), delimiter="\t"))
comp = fasta(f"{INV}.fa")
alleles = fasta(f"{SEQ}/seq.sv:7.alleles.fa")
block = list(csv.DictReader(open(f"{SEQ}/seq.block_sequences.tsv"), delimiter="\t"))
carriers = {r["name"]: [t for t in r["templates"].split(",") if t] for r in block}
templates = sorted({t for ts in carriers.values() for t in ts})
checks.append(("no NA19650 template in the active panel", not any(t.startswith("NA19650#") for t in templates)))
included = sorted({r["source_allele"] for r in rows}, key=lambda a: int(a.split("_")[1]))
checks.append(("allele_75 excluded from the inventory (Amendment 1) and reported", "allele_75" not in included and "allele_75" in alleles))
by_allele = defaultdict(list)
for r in rows:
    by_allele[r["source_allele"]].append(r)
recon_ok = True
for a in included:
    parts = sorted(by_allele[a], key=lambda r: int(r["begin"]))
    seq = "".join(comp[r["unique_record"]] for r in parts)
    recon_ok &= seq == alleles[a] and all(comp[r["unique_record"]] == alleles[a][int(r["begin"]):int(r["end"])] for r in parts)
    recon_ok &= all(hashlib.sha256(comp[r["unique_record"]].encode()).hexdigest() == r["sha256"] for r in parts)
checks.append(("116 included alleles reconstruct byte for byte from their components", recon_ok and len(included) == 116))
checks.append(("carriers in the inventory equal the panel block carriers", all(sorted(r["templates"].split(",")) == sorted(carriers[r["source_allele"]]) for r in rows)))

# component table: one row per unique component sequence, with carriers by physical copy, haplotype and donor
uc = defaultdict(lambda: {"copies": 0, "haps": set(), "donors": set(), "alleles": set(), "positions": Counter(), "type": None})
copy_rows = []
for r in rows:
    u = uc[r["unique_record"]]
    u["type"] = r["type"]
    pos = re.sub(r"^allele_\d+\.", "", r["component"])
    for t in carriers[r["source_allele"]]:
        u["copies"] += 1
        u["haps"].add(t)
        u["donors"].add(t.split("#")[0])
        copy_rows.append((t, t.split("#")[0], r["source_allele"], pos, r["type"], r["unique_record"], int(r["source_units"])))
    u["alleles"].add(r["source_allele"])
    u["positions"][pos] += 1
cls = lambda name: "long" if len(comp[name]) >= LONG_BP else "short"
with open(f"{out}/components.tsv", "w") as f:
    f.write("unique_record\ttype\tbp\tlength_class\tcopies\thaplotypes\tdonors\talleles\tpositions\n")
    for k in sorted(uc, key=lambda k: (uc[k]["type"], int(k.split("_u")[1]))):
        u = uc[k]
        f.write(f"{k}\t{u['type']}\t{len(comp[k])}\t{cls(k) if u['type'] != 'prefix' else '-'}\t{u['copies']}\t{len(u['haps'])}\t{len(u['donors'])}\t{len(u['alleles'])}\t"
                + ",".join(f"{p}:{n}" for p, n in sorted(u["positions"].items())) + "\n")
with open(f"{out}/copies.tsv", "w") as f:
    f.write("template\tdonor\tallele\tposition\ttype\tunique_record\tunits_in_haplotype\n")
    for c in copy_rows:
        f.write("\t".join(str(x) for x in c) + "\n")


def reference(kind, restrict=None):   # frozen rule: most physical copies, then haplotypes, then majority class, then sha
    cand = [k for k in uc if uc[k]["type"] == kind and (restrict is None or cls(k) == restrict)]
    major = Counter(cls(k) for k in cand for _ in range(uc[k]["copies"])).most_common(1)[0][0] if kind != "prefix" else None
    return min(cand, key=lambda k: (-uc[k]["copies"], -len(uc[k]["haps"]), 0 if kind == "prefix" or cls(k) == major else 1,
                                    hashlib.sha256(comp[k].encode()).hexdigest()))


refs = {"prefix": reference("prefix"), "period": reference("period"), "tail": reference("tail"), "period_short": reference("period", "short")}
fa = f"{INV}.fa"
with open(f"{out}/pairs.tsv", "w") as f:
    for kind in ["prefix", "period", "tail"]:
        for k in sorted((k for k in uc if uc[k]["type"] == kind), key=lambda k: int(k.split("_u")[1])):
            f.write(f"{kind}:{k}\t{fa}\t{k}\t{fa}\t{refs[kind]}\n")
    f.write(f"cross:tailref_vs_periodref\t{fa}\t{refs['tail']}\t{fa}\t{refs['period']}\n")
    for k in sorted((k for k in uc if uc[k]["type"] == "period"), key=lambda k: int(k.split("_u")[1])):
        f.write(f"period_short:{k}\t{fa}\t{k}\t{fa}\t{refs['period_short']}\n")
dh = Counter(t.split("#")[0] for t in templates)
multi = {a: ts for a, ts in carriers.items() if len(ts) > 1 and a in included}
with open(f"{out}/prepare.tsv", "w") as f:
    for name, ok in checks:
        f.write(f"check\t{name}\t{'PASS' if ok else 'FAIL'}\n")
    for k, v in refs.items():
        f.write(f"reference\t{k}\t{v} ({len(comp[v])} bp, {cls(v) if k != 'prefix' else '-'}, {uc[v]['copies']} copies, {len(uc[v]['haps'])} haplotypes)\n")
    f.write(f"panel\ttemplates\t{len(templates)} (included alleles: {sum(len(carriers[a]) for a in included)})\n")
    f.write(f"panel\tdonors\t{len(dh)}; haplotypes per donor {dict(Counter(dh.values()))}\n")
    f.write(f"panel\tcopies\t{len(copy_rows)} components ({Counter(c[4] for c in copy_rows)})\n")
    for a, ts in sorted(multi.items(), key=lambda x: int(x[0].split('_')[1])):
        donors = Counter(t.split("#")[0] for t in ts)
        f.write(f"multiplicity\t{a}\t{len(ts)} carriers; {'same donor (both haplotypes)' if max(donors.values()) > 1 else 'different donors'}: {','.join(ts)}\n")
    f.write(f"multiplicity\tunique components in >1 allele\t{sum(1 for u in uc.values() if len(u['alleles']) > 1)}\n")
    f.write(f"multiplicity\tHG00096 templates\t{','.join(t for t in templates if t.startswith('HG00096#'))}\n")
print(open(f"{out}/prepare.tsv").read())
sys.exit(0 if all(ok for _, ok in checks) else 1)
