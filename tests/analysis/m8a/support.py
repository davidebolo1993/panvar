#!/usr/bin/env python3
"""M8-A panel-only feasibility: family-aware leave-out inventory and held-out support under the copying model.

  support.py full <unit dir>      block-sequences on the full models (all blocks); per locus: templates, donors, families,
                                  development status; per (family, haplotype, block): label support, sequence support
  support.py loo <unit dir>       bounded frame-leakage check: for uninspected eligible donors that have an existing LOO
                                  model, project their true haplotypes onto the LOO block chain (block-sequences
                                  --truth-gfa) and test whether each truth block sequence is spelled by an allele carried
                                  by a training template (outside the donor's family)

Definitions (DESIGN.md): label support = the haplotype's own allele at the block is carried by some template outside its
family; sequence support = some template outside the family carries an allele with identical bytes (label-ambiguity
aware). A haplotype is FULLY SUPPORTED when every block has sequence support; otherwise the copying model (no mutation
model) assigns it probability 0 and it is reported, never dropped. Nearest-panel projections are never used.
"""
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict

mode, d = sys.argv[1], sys.argv[2]
R = "/Users/davide.bolognini/panvar/tests/results"
WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
PROBE = f"{WT}/build/genotype_probe"
PED = f"{R}/m6f_copying_predictive_20261019/inputs/20130606_g1k_3202_samples_ped_population.txt"
DEV = f"{WT}/tests/analysis/m8a/DEVELOPMENT_DONORS.tsv"
M4A = f"{R}/m4a_distance_prior_20261019/runs"
LOCI = {
    "C4": dict(gfa=f"{R}/genotype_c4_production_11fe681/graphs/c4.gfa", ref="HG00096#1#haplotype1-0000024:31848049-32074186", curation=None,
               geom=f"{M4A}/geometry.c4/geom.blocks.tsv", loo=lambda dn: f"{R}/genotype_c4_production_11fe681/graphs/{dn}.loo.gfa"),
    "CYP2D6": dict(gfa=f"{R}/cyp2d6_shortread_pilot_v1/data/graphs/cyp2d6.gfa", ref="grch38#1#chr22:42031864-42245566",
                   curation="/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv", geom=f"{M4A}/geometry.cyp_curated/geom.blocks.tsv",
                   loo=lambda dn: f"{R}/cyp2d6_shortread_pilot_v1/data/graphs/loo.{dn}.gfa"),
    "ACOT": dict(gfa=f"{R}/genotype_perf_bounded_20261015/inputs/acot/acot.gfa", ref="GRCh38#0#chr14:73416151-73668669", curation=None,
                 geom=f"{M4A}/geometry.acot/geom.blocks.tsv", loo=lambda dn: f"{R}/genotype_perf_bounded_20261015/inputs/acot/{dn}.loo.gfa"),
}
REFERENCE_GENOMES = ("grch38", "chm13", "GRCh38", "CHM13")
T = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


def run(cmd, log):
    with open(log, "a") as f:
        f.write("$ " + " ".join(cmd) + "\n")
        f.flush()
        return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT).returncode


def fasta(path):
    recs, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs[name] = []
        elif name:
            recs[name].append(line)
    return {k: "".join(v).upper() for k, v in recs.items()}


def families():
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
        for p in (fa, mo):
            if p != "0":
                parent[find(s)] = find(p)
    ids = defaultdict(set)
    for fam, s, *_ in ped:
        ids[find(s)].add(fam)
    return (lambda dn: find(dn) if dn in parent else "SELF:" + dn), (lambda root: "+".join(sorted(ids.get(root, []))) or root.replace("SELF:", ""))


def donor_of(template):
    return template.split("#")[0]


def write(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)


root_of, fam_name = families()
dev = {l.split("\t")[0] for l in open(DEV) if not l.startswith(("#", "donor\t"))}
log = f"{d}/commands.log"

if mode == "full":
    inv, sup, fam_rows = [], [], []
    for locus, L in LOCI.items():
        blocks = [r["block"] for r in T(L["geom"])]
        cmd = [PROBE, "block-sequences", "--out", f"{d}/{locus}", "--reference", L["ref"], "--model-gfa", L["gfa"], "--blocks", ",".join(blocks)]
        if L["curation"]:
            cmd += ["--panel-curation", L["curation"]]
        assert run(cmd, log) == 0, f"block-sequences failed for {locus}"
        tab = [l.rstrip("\n").split("\t") for l in open(f"{d}/{locus}.block_sequences.tsv")][1:]
        seqs = {b: fasta(f"{d}/{locus}.{b}.alleles.fa") for b in blocks}
        path = defaultdict(dict)             # template -> block -> allele name
        carriers = defaultdict(set)          # (block, allele) -> templates
        for row in tab:
            b, src, idx, name = row[0], row[1], row[2], row[3]
            if src != "panel":
                continue
            for t in [x for x in row[5].split(",") if x and x != "-"]:
                if b in path[t]:
                    path[t][b] = "AMBIGUOUS"
                else:
                    path[t][b] = name
                carriers[(b, name)].add(t)
        templates = sorted(path)
        incomplete = [t for t in templates if len(path[t]) != len(blocks) or "AMBIGUOUS" in path[t].values()]
        dup_seq = sum(len(v) - len(set(v.values())) for v in seqs.values())
        donors = sorted({donor_of(t) for t in templates})
        roots = defaultdict(list)
        for dn in donors:
            roots[root_of(dn)].append(dn)
        status = {}
        for r_, dns in roots.items():
            st = ("REFERENCE_GENOME" if any(dn.startswith(REFERENCE_GENOMES) for dn in dns) else
                  "INSPECTED_DEVELOPMENT" if any(dn in dev for dn in dns) else "UNINSPECTED")
            for dn in dns:
                status[dn] = st
            fam_rows.append(dict(locus=locus, family=fam_name(r_), donors=",".join(dns), templates=sum(1 for t in templates if donor_of(t) in dns),
                                 status=st, related_donors_in_panel=len(dns)))
        identical_paths = len(templates) - len({tuple(path[t].get(b) for b in blocks) for t in templates})
        inv.append(dict(locus=locus, blocks=len(blocks), templates=len(templates), donors=len(donors), families=len(roots),
                        multi_donor_families=sum(1 for v in roots.values() if len(v) > 1),
                        uninspected_families=sum(1 for r_, v in roots.items() if status[v[0]] == "UNINSPECTED"),
                        inspected_families=sum(1 for r_, v in roots.items() if status[v[0]] == "INSPECTED_DEVELOPMENT"),
                        reference_genome_families=sum(1 for r_, v in roots.items() if status[v[0]] == "REFERENCE_GENOME"),
                        templates_incomplete_or_ambiguous=len(incomplete), duplicate_allele_sequences_within_blocks=dup_seq,
                        templates_with_an_identical_allele_path_to_another=identical_paths))
        for r_, dns in roots.items():
            fam_t = {t for t in templates if donor_of(t) in dns}
            train = [t for t in templates if t not in fam_t]
            for t in sorted(fam_t):
                if t in incomplete:
                    sup.append(dict(locus=locus, family=fam_name(r_), status=status[dns[0]], template=t, assessable="no (incomplete/ambiguous block projection)"))
                    continue
                lab_miss, seq_miss, lab_only_miss, miss_blocks = 0, 0, 0, []
                for b in blocks:
                    a = path[t][b]
                    lab = any(x in train for x in carriers[(b, a)])
                    s = seqs[b].get(a)
                    seq = lab or any(seqs[b].get(path[x][b]) == s for x in train)
                    lab_miss += not lab
                    seq_miss += not seq
                    lab_only_miss += (not lab) and seq
                    if not seq:
                        miss_blocks.append(b)
                sup.append(dict(locus=locus, family=fam_name(r_), status=status[dns[0]], template=t, assessable="yes",
                                blocks_label_unsupported=lab_miss, blocks_sequence_unsupported=seq_miss, label_only_unsupported=lab_only_miss,
                                fully_supported="yes" if seq_miss == 0 else "no", unsupported_blocks=",".join(miss_blocks),
                                identical_full_path_in_training="yes" if any(all(path[x][b] == path[t][b] for b in blocks) for x in train) else "no"))
    write(f"{d}/inventory.tsv", inv)
    write(f"{d}/families.tsv", fam_rows)
    write(f"{d}/support.tsv", sup)
    # summary by locus and status
    summ = []
    for locus in LOCI:
        for st in ("UNINSPECTED", "INSPECTED_DEVELOPMENT", "REFERENCE_GENOME"):
            rows = [r for r in sup if r["locus"] == locus and r["status"] == st]
            ass = [r for r in rows if r["assessable"] == "yes"]
            if not rows:
                continue
            summ.append(dict(locus=locus, status=st, haplotypes=len(rows), assessable=len(ass),
                             fully_supported=sum(1 for r in ass if r["fully_supported"] == "yes"),
                             identical_full_path_in_training=sum(1 for r in ass if r["identical_full_path_in_training"] == "yes"),
                             unsupported_blocks_total=sum(int(r["blocks_sequence_unsupported"]) for r in ass),
                             label_only_unsupported_blocks=sum(int(r["label_only_unsupported"]) for r in ass),
                             families=len({r["family"] for r in rows})))
    write(f"{d}/support_summary.tsv", summ)
    for p in ("inventory.tsv", "support_summary.tsv"):
        print(open(f"{d}/{p}").read())

if mode == "loo":
    full = T(f"{d.rsplit('/', 1)[0]}/support.full/support.tsv")
    fam = T(f"{d.rsplit('/', 1)[0]}/support.full/families.tsv")
    rows = []
    for locus, L in LOCI.items():
        un = [r for r in fam if r["locus"] == locus and r["status"] == "UNINSPECTED"]
        cands = [(r["family"], dn) for r in un for dn in r["donors"].split(",") if os.path.exists(L["loo"](dn))][:8]   # frozen cap: first 8 by family order
        for family, dn in cands:
            out = f"{d}/{locus}.{dn}"
            blocks = [r["block"] for r in T(L["geom"])]
            cmd = [PROBE, "block-sequences", "--out", out, "--reference", L["ref"], "--model-gfa", L["loo"](dn), "--blocks", ",".join(blocks),
                   "--truth-gfa", L["gfa"], "--donor", dn]
            if L["curation"]:
                cmd += ["--panel-curation", L["curation"]]
            code = run(cmd, log)
            if code != 0:
                rows.append(dict(locus=locus, family=family, donor=dn, status="LOO frame unavailable (block-sequences exit %d; see commands.log)" % code))
                continue
            tab = [l.rstrip("\n").split("\t") for l in open(f"{out}.block_sequences.tsv")][1:]
            loo_blocks = sorted({r[0] for r in tab}, key=lambda b: [x[0] for x in tab].index(b))
            fam_donors = next(r["donors"].split(",") for r in fam if r["locus"] == locus and r["family"] == family)
            train_alleles = defaultdict(set)
            for r in tab:
                if r[1] == "panel":
                    ts = [x for x in r[5].split(",") if x and x != "-" and donor_of(x) not in fam_donors]
                    if ts:
                        train_alleles[r[0]].add(r[3])
            for k in (1, 2):
                miss, noframe = [], []
                for b in loo_blocks:
                    al = fasta(f"{out}.{b}.alleles.fa")
                    trf = f"{out}.{b}.truth.fa"
                    tr = fasta(trf).get(f"truth{k}") if os.path.exists(trf) else None
                    if tr is None:
                        noframe.append(b)
                        continue
                    if not any(al.get(a) == tr for a in train_alleles[b]):
                        miss.append(b)
                full_row = [r for r in full if r["locus"] == locus and r["template"].startswith(dn + "#") and r["template"].split("#")[1] == str(k)]
                rows.append(dict(locus=locus, family=family, donor=dn, haplotype=k, loo_blocks=len(loo_blocks), truth_blocks_without_frame=len(noframe),
                                 loo_frame_unsupported_blocks=len(miss), loo_frame_fully_supported="yes" if not miss and not noframe else "no",
                                 unsupported=",".join(miss), full_frame_fully_supported=full_row[0].get("fully_supported", "NA") if full_row else "NA"))
    write(f"{d}/loo_support.tsv", rows)
    print(open(f"{d}/loo_support.tsv").read())
