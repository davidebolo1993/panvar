#!/usr/bin/env python3
"""Synthetic LPA-like example for the SAMPLE-LEVEL k-mer GWAS demo.

Models the real pipeline: a pangenome HAPLOTYPE PANEL (the graph paths) plus a DIPLOID COHORT
of individuals, each assigned two haplotypes (a cosigt-style samples.tsv). GWAS then tests the
SAMPLE phenotype against the SAMPLE genotype = aggregate over its two haplotypes.

Three loci (between single-copy anchors L0..L3):
  KIV-2 VNTR (CAUSAL, multiplicity-only): L0 - K (Ck K)x(CN-1) - L1, folded via a K-Ck-K cycle
      (no self-loop). Per-haplotype CN in [1..6], ALWAYS >= 1 -> the unit k-mers are present in
      EVERY haplotype, so presence/absence sees no contrast; only the COUNT (= copy number)
      carries signal. Diploid sample dosage = CN_A + CN_B (the real KIV-2 biology vs Lp(a)).
  Decoy A (NEG control, presence/absence, unlinked): L1 - Sa - [DA?] - L2.
  Decoy B (POS control, presence/absence, phenotype-linked): L2 - Sb - [DB?] - L3 -> recovered
      by BOTH methods (shows P/A is fine for indel/SNP-like variants, just blind to pure CNV).

Sample phenotype:
  lpa_continuous = BASE - SLOPE*(CN_A+CN_B) + DB_EFFECT*has_DB_sample + N(0, SIGMA)   (inverse)
  case_binary    = 1 if lpa_continuous > median (high Lp(a) = case)

Writes (committed):
  synthetic/lpa.gfa, synthetic/lpa.snarls.jsonl,
  synthetic/samples.tsv      (cosigt-style: sample <tab> hapA <tab> hapB),
  synthetic/phenotypes.tsv   (per SAMPLE),
  synthetic/truth.tsv        (causal nodes).
"""
import os
import random
import statistics

SEED = 17
N_HAPLOTYPES = 40          # the pangenome haplotype panel (graph paths)
N_SAMPLES = 150            # diploid individuals genotyped against the panel
CN_CHOICES = [1, 2, 3, 4, 5, 6]
BASE, SLOPE, DB_EFFECT, SIGMA = 160.0, 11.0, 18.0, 8.0   # Lp(a)-like arbitrary units

random.seed(SEED)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "synthetic")
os.makedirs(OUT, exist_ok=True)

names, seqs, order = {}, {}, []
def mk(name, length):
    if name in names:
        return names[name]
    nid = str(len(order) + 1)
    names[name] = nid
    seqs[nid] = "".join(random.choice("ACGT") for _ in range(length))
    order.append(nid)
    return nid

L = [mk(f"L{i}", 40) for i in range(4)]      # L0..L3 flank anchors
K  = mk("K", 80); Ck = mk("Ck", 15)          # KIV-2 unit + spacer (causal)
Sa = mk("Sa", 50); DA = mk("DA", 60)         # decoy A (neg)
Sb = mk("Sb", 50); DB = mk("DB", 60)         # decoy B (pos)

def step(name, rev=False):
    return (names[name], rev)

def build_hap(cn, has_da, has_db):
    p = [step("L0"), step("K")]
    for _ in range(cn - 1):
        p += [step("Ck"), step("K")]
    p += [step("L1"), step("Sa")]
    if has_da:
        p.append(step("DA"))
    p += [step("L2"), step("Sb")]
    if has_db:
        p.append(step("DB"))
    p.append(step("L3"))
    return p

# ---- haplotype panel (graph paths) ----
haps = {}            # hap_path_name -> (cn, has_da, has_db)
hap_paths = {}       # hap_path_name -> walk
hap_names = []
for h in range(N_HAPLOTYPES):
    if h == 0:
        cn, da, db = 3, True, True   # hap00 = call --reference-path: carries every locus
    else:
        # CN cycles 1..6; DA is assigned ORTHOGONAL to CN (alternates within each CN stratum)
        # so it is a genuine null (uncorrelated with CN and hence with the phenotype). DB is
        # the phenotype-linked decoy (it enters the phenotype directly, so correlation with CN
        # is irrelevant); ~half the haplotypes carry it.
        cn = (h % 6) + 1
        da = ((h // 6) % 2 == 0)
        db = ((h // 3) % 2 == 0)
    name = f"hap{h:02d}#0#chrL:1-3000"
    haps[name] = (cn, da, db)
    hap_paths[name] = build_hap(cn, da, db)
    hap_names.append(name)

# ---- diploid cohort: each sample = 2 haplotypes from the panel ----
samples = []         # (sample_id, hapA, hapB)
prows = []           # (sample_id, kiv2_dosage, has_db_sample, lpa)
for i in range(N_SAMPLES):
    a, b = random.choice(hap_names), random.choice(hap_names)
    cn_a, _, db_a = haps[a]
    cn_b, _, db_b = haps[b]
    dosage = cn_a + cn_b                       # diploid KIV-2 copy number
    has_db = 1 if (db_a or db_b) else 0
    lpa = BASE - SLOPE * dosage + DB_EFFECT * has_db + random.gauss(0, SIGMA)
    sid = f"sample{i:03d}"
    samples.append((sid, a, b))
    prows.append([sid, dosage, has_db, lpa])

median_lpa = statistics.median(r[3] for r in prows)
for r in prows:
    r.append(1 if r[3] > median_lpa else 0)

# ---- write GFA (haplotype panel) ----
edges = set()
for p in hap_paths.values():
    for (na, ra), (nb, rb) in zip(p, p[1:]):
        edges.add((na, "-" if ra else "+", nb, "-" if rb else "+"))
with open(os.path.join(OUT, "lpa.gfa"), "w") as f:
    f.write("H\tVN:Z:1.0\n")
    for nid in order:
        f.write(f"S\t{nid}\t{seqs[nid]}\n")
    for a, ao, b, bo in sorted(edges):
        f.write(f"L\t{a}\t{ao}\t{b}\t{bo}\t0M\n")
    for name, p in hap_paths.items():
        walk = ",".join(f"{nid}{'-' if rev else '+'}" for nid, rev in p)
        f.write(f"P\t{name}\t{walk}\t*\n")

with open(os.path.join(OUT, "lpa.snarls.jsonl"), "w") as f:
    for a, b in [(L[0], L[1]), (L[1], L[2]), (L[2], L[3])]:
        f.write('{"start": {"node_id": "%s"}, "end": {"node_id": "%s"}, "type": 1, '
                '"start_end_reachable": true, "directed_acyclic_net_graph": true}\n' % (a, b))

with open(os.path.join(OUT, "samples.tsv"), "w") as f:
    f.write("sample\thaplotype_1\thaplotype_2\n")
    for sid, a, b in samples:
        f.write(f"{sid}\t{a}\t{b}\n")

with open(os.path.join(OUT, "phenotypes.tsv"), "w") as f:
    f.write("sample\tkiv2_dosage\thas_db\tlpa_continuous\tcase_binary\n")
    for sid, dosage, has_db, lpa, case in prows:
        f.write(f"{sid}\t{dosage}\t{has_db}\t{lpa:.4f}\t{case}\n")

with open(os.path.join(OUT, "truth.tsv"), "w") as f:
    f.write("locus\tnodes\tmodel\tlinked\texpect\n")
    f.write(f"KIV2_VNTR\t{names['K']},{names['Ck']}\tmultiplicity\tyes\tcount-significant, P/A-not\n")
    f.write(f"decoy_A\t{names['DA']}\tpresence_absence\tno\tnon-significant (both)\n")
    f.write(f"decoy_B\t{names['DB']}\tpresence_absence\tyes\tsignificant (both)\n")

print(f"wrote {OUT}/lpa.gfa  nodes={len(order)} edges={len(edges)} haplotypes={len(hap_paths)}")
print(f"  samples={len(samples)} (diploid) | KIV-2 unit node={names['K']} decoyA={names['DA']} decoyB={names['DB']}")
print(f"  KIV-2 sample dosage range: {min(r[1] for r in prows)}..{max(r[1] for r in prows)}")
print(f"  cases (high Lp(a)): {sum(r[4] for r in prows)}/{len(prows)} (median Lp(a)={median_lpa:.1f})")
