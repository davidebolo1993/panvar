#!/usr/bin/env python3
"""M7-B canonical pools (the frozen M7-A manifest; nothing added or removed).

  pools.py <unit dir>
Per cell: spell D/X candidates (genotype_probe spell, model frame), copy O records (M1-D/M1-E FASTA, model frame), check
orientation and identity, canonicalise haplotypes by sequence and diplotypes by unordered pair, write mapping.tsv, pool
FASTAs, --priors files and Panvar candidate files.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m7b_common import *  # noqa: E402,F401,F403

d = sys.argv[1]
log = f"{d}/commands.log"
man = table(f"{M7A}/CANDIDATES.tsv")
assert sha(f"{M7A}/CANDIDATES.tsv").startswith("a8cd9516"), "M7-A manifest changed"
checks, mapping = [], []


def kmers(s, k=31, step=1):
    return {s[i:i + k] for i in range(0, len(s) - k + 1, step)}


for cell, c in CELLS.items():
    rows = [r for r in man if r["cell"] == cell]
    L = LOCI[c["locus"]]
    ref = spell_gfa_path(L["gfa"], L["ref_path"])
    refk = kmers(ref)
    # D / X: spell
    cand = f"{d}/{cell}.spell_candidates.tsv"
    with open(cand, "w") as f:
        for i, r in enumerate(rows):
            if r["representation"] == "panel_allele_paths":
                f.write(f"cand{i}\t{r['hap1_alleles']}\t{r['hap2_alleles']}\n")
    cmd = [SCORING_PROBE, "spell", "--out", f"{d}/{cell}.spell", "--reference", L["ref_path"], "--model-gfa", c["model"], "--candidates", cand, "--block", c["block"]]
    if c["curation"]:
        cmd += ["--panel-curation", c["curation"]]
    code = run(cmd, log)
    assert code == 0, f"spell failed for {cell}"
    sp = fasta(f"{d}/{cell}.spell.fa")
    spelled_tab = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(f"{d}/{cell}.spell.spelled.tsv")}
    seqs = {}   # row index -> (seq1, seq2, block interval info)
    for i, r in enumerate(rows):
        if r["representation"] == "panel_allele_paths":
            seqs[i] = (sp[f"cand{i}.1"], sp[f"cand{i}.2"], spelled_tab.get(f"cand{i}.1"), spelled_tab.get(f"cand{i}.2"))
        else:
            fa = fasta(r["fasta"])
            r1, r2 = r["records"].split("|")
            seqs[i] = (fa[r1], fa[r2], None, None)
    # identity check: M1-D's reported.1/.2 equal the spelled reported haplotypes (C4 LOO)
    if cell == "c4_loo_na19650_10":
        m1d = fasta(f"{R}/m1d_unit_reachability_20261019/runs/construct.c4_10.a2/candidates.fa")
        rep = next(i for i, r in enumerate(rows) if "reported" in r["labels"].split("|"))
        ok = {m1d["reported.1"], m1d["reported.2"]} == {seqs[rep][0], seqs[rep][1]}
        checks.append(dict(cell=cell, check="M1-D reported.1/.2 == spelled reported haplotypes", ok=ok))
    # orientation and N checks (sampled 31-mers, every 97th position)
    worst = 1.0
    for i, (s1, s2, _, _) in seqs.items():
        for s in (s1, s2):
            assert "N" not in s, f"{cell} row {i} contains N"
            sample = [s[j:j + 31] for j in range(0, len(s) - 31, 97)]
            fwd = sum(1 for k in sample if k in refk) / len(sample)
            rev = sum(1 for k in sample if revcomp(k) in refk) / len(sample)
            worst = min(worst, fwd - rev)
    checks.append(dict(cell=cell, check="every haplotype: forward 31-mer containment in the locus reference minus reverse-complement containment (worst)",
                       value=f"{worst:.3f}", ok=worst > 0.5))
    # canonicalise
    hap_id, hap_seq = {}, {}
    for i in sorted(seqs, key=lambda i: ("DXO".index(rows[i]["tier"]), i)):
        for s in seqs[i][:2]:
            h = hashlib.sha256(s.encode()).hexdigest()
            if h not in hap_id:
                hap_id[h] = f"{cell}.h{len(hap_id):02d}"
                hap_seq[hap_id[h]] = s
    dip_of = {}
    for i, r in enumerate(rows):
        a, b = (hap_id[hashlib.sha256(s.encode()).hexdigest()] for s in seqs[i][:2])
        dip = ",".join(sorted((a, b)))
        dip_of[i] = dip
        mapping.append(dict(cell=cell, row=i, tier=r["tier"], labels=r["labels"], sources=r["sources"], representation=r["representation"],
                            hap1_alleles=r["hap1_alleles"], hap2_alleles=r["hap2_alleles"], target_alleles=r["target_alleles"], unit_copies=r["unit_copies"],
                            unit_class=r["unit_class"], unit_total=r["unit_total"], also_truth_selected=r.get("also_truth_selected", ""),
                            caller_scored=r.get("caller_scored", ""), saved_expected_total=r["saved_expected_total"],
                            hap1_sha256=hashlib.sha256(seqs[i][0].encode()).hexdigest()[:16], hap2_sha256=hashlib.sha256(seqs[i][1].encode()).hexdigest()[:16],
                            canonical_hap1=a, canonical_hap2=b, canonical_diplotype=dip,
                            block_interval_1="-" if seqs[i][2] is None else f"{seqs[i][2][3]}-{seqs[i][2][4]}",
                            block_interval_2="-" if seqs[i][3] is None else f"{seqs[i][3][3]}-{seqs[i][3][4]}"))
    pools = {"D": ["D"]}
    if any(r["tier"] == "X" for r in rows):
        pools["DX"] = ["D", "X"]
    if any(r["tier"] == "O" for r in rows):
        pools["DO"] = ["D", "O"]
    for name, tiers in pools.items():
        dips = sorted({dip_of[i] for i, r in enumerate(rows) if r["tier"] in tiers})
        haps = sorted({h for dip in dips for h in dip.split(",")})
        with open(f"{d}/{cell}.pool.{name}.fa", "w") as f:
            for h in haps:
                f.write(f">{h}\n{hap_seq[h]}\n")
        with open(f"{d}/{cell}.priors.{name}.tsv", "w") as f:
            for dip in dips:
                f.write(f"{c['locus']}\t{dip}\t0\n")
        checks.append(dict(cell=cell, check=f"pool {name}", value=f"{len(dips)} canonical diplotypes, {len(haps)} canonical haplotypes "
                           f"(from {sum(1 for r in rows if r['tier'] in tiers)} frozen candidates)", ok=True))
    # Panvar candidate files: D and X by allele paths (label = row id), O explicit by copying the M1-D/M1-E lines
    for name, tier in (("D", "D"), ("X", "X")):
        sel = [(i, r) for i, r in enumerate(rows) if r["tier"] == tier]
        if sel:
            with open(f"{d}/{cell}.panvar.{name}.tsv", "w") as f:
                for i, r in sel:
                    f.write(f"row{i}\t{r['hap1_alleles']}\t{r['hap2_alleles']}\t{r['saved_expected_total']}\n")
    osel = [(i, r) for i, r in enumerate(rows) if r["tier"] == "O" and r["representation"] == "explicit_sequences"]
    if osel:
        with open(f"{d}/{cell}.panvar.O.tsv", "w") as f:
            for i, r in osel:
                src = r["fasta"].replace("/candidates.fa", "/candidates.explicit.tsv")
                label = r["labels"].split("_", 1)[1]
                line = next(l for l in open(src) if l.split("\t")[0] == label)
                parts = line.rstrip("\n").split("\t")
                parts[0] = f"row{i}"
                parts[-1] = "NA"
                f.write("\t".join(parts) + "\n")
    # block sequences for the evaluation-only truth distance (D, X)
    with open(f"{d}/{cell}.blocks.fa", "w") as f:
        for i, (s1, s2, t1, t2) in seqs.items():
            if t1 is None:
                continue
            for k, (s, t) in enumerate(((s1, t1), (s2, t2)), 1):
                b, e = int(t[3]), int(t[4])
                f.write(f">row{i}.{k}\n{s[b:e]}\n")
write_tsv(f"{d}/mapping.tsv", mapping)
write_tsv(f"{d}/pools_checks.tsv", checks)
print(open(f"{d}/pools_checks.tsv").read())
sys.exit(0 if all(r["ok"] for r in checks) else 1)
