#!/usr/bin/env python3
"""M7-B one Locityper fixed-pool run.

  loc.py <cell> <pool D|DX|DO> <seed> <unit dir>
target (database = the pool's canonical haplotypes), verify the database, genotype with --priors restricting the genotype
set, then verify: options == priors, every genotype has a stage-2 (annealing) row and 20 attempts, read-pair counts.
"""
import gzip
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m7b_common import *  # noqa: E402,F401,F403

cell, pool, seed, d = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
c = CELLS[cell]
for r in table(f"{OLD}/runs/provenance/lineage.tsv"):
    if r["lineage"] == c["lineage"]:
        assert r["verdict"] == "PASS", f"lineage {c['lineage']} unresolved: comparison stopped"
log = f"{d}/commands.log"
P = f"{OLD}/runs/pools"
assert open(f"{OUT}/runs/preproc.{cell}/calibration_ok").read().strip() == "PASS", "calibration failed validation: comparison stopped"
refdir = f"{OUT}/runs/reference.{c['locus']}"
loc_len = int({l.split("\t")[0]: l.split("\t")[1] for l in open(f"{refdir}/locus.tsv")}["length"])
contig = LOCI[c["locus"]]["contig"]
pool_fa = f"{P}/{cell}.pool.{pool}.fa"
priors = f"{P}/{cell}.priors.{pool}.tsv"
n_gt = sum(1 for _ in open(priors))
checks = []
code = run([LOCITYPER, "target", "-d", f"{d}/db", "-r", f"{refdir}/ref.fa", "-j", f"{refdir}/ref.jf", "-l", c["locus"], f"{contig}:1-{loc_len}",
            pool_fa, "-@", 4, "--jellyfish", JELLYFISH], log)
assert code == 0, "locityper target failed"
db = fasta(f"{d}/db/loci/{c['locus']}/haplotypes.fa.gz")
want = fasta(pool_fa)
checks.append(dict(check="database haplotypes == pool (names and sequences); nothing discarded",
                   value=f"{len(db)} in database, {len(want)} in pool", ok=dict(db) == dict(want) and not os.path.exists(f"{d}/db/loci/{c['locus']}/discarded_haplotypes.txt")))
code = run([LOCITYPER, "genotype", "-i", f"{c['reads']}/reads_1.fq", f"{c['reads']}/reads_2.fq", "-p", f"{OUT}/runs/preproc.{cell}/preproc", "-d", f"{d}/db",
            "-o", f"{d}/gt", "-@", 4, "-s", seed, "--debug", 1, "--basis", "none", "--priors", priors, "--filt-diff", "1e9", "--prob-thresh", "-1e9",
            "--dont-skip", "-S", "greedy:i=5k,a=1", "-S", f"anneal:i={n_gt},a=20", "--strobealign", STROBEALIGN, "--samtools", SAMTOOLS], log)
assert code == 0, "locityper genotype failed"
gdir = f"{d}/gt/loci/{c['locus']}"
res = json.load(gzip.open(f"{gdir}/res.json.gz", "rt"))
opts = res.get("options", [])
want_gts = {l.split("\t")[1] for l in open(priors)}
got = {o["genotype"] for o in opts}
canon = lambda g: ",".join(sorted(g.split(",")))
checks.append(dict(check="res.json options == priors set", value=f"{len(got)} options, {len(want_gts)} declared",
                   ok={canon(g) for g in got} == want_gts and all(o.get("lik_mean") is not None for o in opts)))
sol = subprocess.run([BROTLI, "-dc", f"{gdir}/sol.csv.br"], capture_output=True, text=True).stdout.splitlines()
stage2 = {canon(l.split("\t")[1]) for l in sol if l.split("\t")[0] == "2"}
checks.append(dict(check="every canonical genotype has a stage-2 (annealing) row in sol.csv.br", value=f"{len(stage2)}/{len(want_gts)}", ok=stage2 == want_gts))
attempts, acct, hdr = {}, {}, None
for fn in sorted(os.listdir(gdir)):
    if fn.startswith("sol_ext") and fn.endswith(".br"):
        for l in subprocess.run([BROTLI, "-dc", f"{gdir}/{fn}"], capture_output=True, text=True).stdout.splitlines():
            f = l.split("\t")
            if f and f[0] == "stage":
                hdr = f
                continue
            if f and f[0] == "2" and len(f) >= 3:
                g = canon(f[1])
                attempts[g] = attempts.get(g, 0) + 1
                if hdr and len(f) == len(hdr):   # per-attempt read accounting at the scoring stage
                    for k in ("total_reads", "unmapped", "out_of_bounds"):
                        if k in hdr:
                            acct.setdefault((g, k), []).append(float(f[hdr.index(k)]))
checks.append(dict(check="20 annealing attempts per genotype (sol_ext)", value=f"{sorted(set(attempts.values())) if attempts else 'no sol_ext file'}",
                   ok=bool(attempts) and set(attempts.values()) == {20} and set(attempts) == want_gts))
text = open(log).read()
reads = [l.strip() for l in text.splitlines() if re.search(r"(recruit|read pairs|reads|unmapped|Mapp|used)", l, re.I) and not l.startswith("$")][:20]
mean = lambda v: sum(v) / len(v) if v else ""
rows = [dict(genotype=canon(o["genotype"]), lik_mean_log10=o["lik_mean"], lik_sd_log10=o.get("lik_sd"), prob=o.get("prob"), log10_prob=o.get("log10_prob"),
             scored_reads_mean=mean(acct.get((canon(o["genotype"]), "total_reads"), [])), unmapped_mean=mean(acct.get((canon(o["genotype"]), "unmapped"), [])),
             out_of_bounds_mean=mean(acct.get((canon(o["genotype"]), "out_of_bounds"), [])))
        for o in opts]
write_tsv(f"{d}/options.tsv", rows)
with open(f"{d}/run_summary.tsv", "w") as f:
    f.write(f"cell\t{cell}\npool\t{pool}\nseed\t{seed}\ngenotypes_declared\t{len(want_gts)}\ntotal_reads\t{res.get('total_reads')}\n"
            f"top_genotype\t{res.get('genotype')}\nquality\t{res.get('quality')}\nunexplained_reads\t{res.get('unexpl_reads')}\n"
            f"warnings\t{res.get('warnings')}\n")
    for l in reads:
        f.write(f"log\t{l}\n")
# target-window GC support (windows with weight >= 0.001 against the cell's background per-GC window counts)
sup = {int(r["gc"]): int(r["nwindows"]) for r in table(f"{OUT}/runs/preproc.{cell}/gc_support.tsv")}
win = [l.rstrip("\n").split("\t") for l in gzip.open(f"{gdir}/windows.bed.gz", "rt") if not l.startswith("#")]
used = [(int(w[3]), float(w[7])) for w in win if float(w[7]) >= 0.001]
weak = [g for g, _ in used if sup.get(g, 0) < 100]
with open(f"{d}/target_gc_support.tsv", "w") as f:
    f.write(f"target_windows_used\t{len(used)}\nwindows_in_gc_bins_with_<100_background_windows\t{len(weak)}\nfraction_weak\t{len(weak) / len(used) if used else float('nan'):.4f}\n"
            f"target_gc_range\t{min(g for g, _ in used) if used else 'NA'}-{max(g for g, _ in used) if used else 'NA'}\n"
            f"weak_gc_values\t{','.join(str(g) for g in sorted(set(weak)))}\n")
write_tsv(f"{d}/checks.tsv", checks)
print(open(f"{d}/checks.tsv").read())
print(open(f"{d}/run_summary.tsv").read()[:2500])
sys.exit(0 if all(r["ok"] for r in checks) else 1)
