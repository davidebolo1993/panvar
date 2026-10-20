"""M7-B shared configuration (frozen with DESIGN.md)."""
import csv
import hashlib
import os
import shlex
import subprocess
from collections import OrderedDict

R = "/Users/davide.bolognini/panvar/tests/results"
OUT = f"{R}/m7b_fixed_pool_20261020"
WT = "/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references"
M7A = f"{WT}/tests/analysis/m7a"
SCORING_PROBE = f"{WT}/build/genotype_probe"
C4_SIM_CANDIDATES = [f"{WT}/build/genotype_probe", "/Users/davide.bolognini/panvar-worktrees/mosaic-11fe681/build/genotype_probe"]
CYP_SIM = f"{R}/cyp2d6_shortread_pilot_v1/bin/genotype_probe"
ENV = "/Users/davide.bolognini/miniconda3/envs/locityper/bin"
LOCITYPER, STROBEALIGN, SAMTOOLS = f"{ENV}/locityper", f"{ENV}/strobealign", f"{ENV}/samtools"
JELLYFISH = "/Users/davide.bolognini/miniconda3/pkgs/kmer-jellyfish-2.3.1-py312pl5321h1f8763b_6/bin/jellyfish"
BROTLI = "/opt/homebrew/bin/brotli"
CURATION = "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"
SEED1, SEED2, PREPROC_SEED = 20261020, 20261021, 20261020
BG_PARTS = [(f"{R}/genotype_pilot_acot_2x_d24e471/panel.sorted.gfa", "GRCh38#0#chr14:73416151-73668669"),
            (f"{R}/genotype_pilot_gstm1_2x_d24e471/panel.sorted.gfa", "GRCh38#0#chr1:109561616-109802331")]
SPACER = 1000
C4REF = "HG00096#1#haplotype1-0000024:31848049-32074186"
CYPREF = "grch38#1#chr22:42031864-42245566"
LOCI = {"C4": dict(gfa=f"{R}/genotype_c4_production_11fe681/graphs/c4.gfa", ref_path=C4REF, contig="c4ref"),
        "CYP2D6": dict(gfa=f"{R}/cyp2d6_shortread_pilot_v1/data/graphs/cyp2d6.gfa", ref_path=CYPREF, contig="cypref")}
CELLS = OrderedDict([
    ("c4_lzo_hg00171_5", dict(locus="C4", lineage="C4", donor="HG00171", depth=5, seed=20260916, block="sv:7", t=7,
        model=f"{R}/genotype_c4_production_11fe681/graphs/c4.gfa", reads=f"{R}/genotype_c4_production_11fe681/reads/HG00171/cov5",
        run=f"{R}/m4b_distance_prior_20261019/runs/run.c4_lzo_hg00171_5.const/call", curation=None,
        truth_block=f"{R}/m6a_count_model_20261019/runs/blockseq.c4_lzo_hg00171_5/seq.sv:7.truth.fa", control=True)),
    ("cyp_lzo_hg02769_5", dict(locus="CYP2D6", lineage="CYP", donor="HG02769", depth=5, seed=20260921, block="sv:9", t=9,
        model=f"{R}/cyp2d6_shortread_pilot_v1/data/graphs/cyp2d6.gfa", reads=f"{R}/cyp2d6_shortread_pilot_v1/data/reads/HG02769/cov5",
        run=f"{R}/m4b_distance_prior_20261019/runs/run.cyp_lzo_hg02769_5.const/call", curation=CURATION,
        truth_block=f"{R}/m6a_count_model_20261019/runs/blockseq.cyp_lzo_hg02769_5/seq.sv:9.truth.fa", control=True)),
    ("c4_loo_na19650_10", dict(locus="C4", lineage="C4", donor="NA19650", depth=10, seed=20260916, block="sv:7", t=7,
        model=f"{R}/genotype_c4_production_11fe681/graphs/NA19650.loo.gfa", reads=f"{R}/genotype_c4_production_11fe681/reads/NA19650/cov10",
        run=f"{R}/m1b_structural_20261019/runs/call.c4_10/call", curation=None,
        truth_block=f"{R}/m6a_count_model_20261019/runs/blockseq.c4_loo_na19650_10/seq.sv:7.truth.fa", control=False)),
    ("cyp_loo_na19650_5", dict(locus="CYP2D6", lineage="CYP", donor="NA19650", depth=5, seed=20260921, block="sv:9", t=9,
        model=f"{R}/cyp2d6_shortread_pilot_v1/data/graphs/loo.NA19650.gfa", reads=f"{R}/cyp2d6_shortread_pilot_v1/data/reads/NA19650/cov5",
        run=f"{R}/m1b_structural_20261019/runs/call.cyp_5/call", curation=CURATION,
        truth_block=f"{R}/m6a_count_model_20261019/runs/blockseq.cyp_loo_na19650_5/seq.sv:9.truth.fa", control=False)),
])
LINEAGE_DONORS = {"C4": ["NA19650", "HG00171"], "CYP": ["NA19650", "HG02769"]}


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def table(path):
    return list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))


def fasta(path):
    recs, name = OrderedDict(), None
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = lambda p: gzip.open(p, "rt")
    for line in opener(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs[name] = []
        elif name is not None:
            recs[name].append(line)
    return OrderedDict((k, "".join(v).upper()) for k, v in recs.items())


COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def spell_gfa_path(gfa, path_name):
    """Spell a GFA P-line path, reverse-complementing minus steps."""
    segs, steps = {}, None
    for line in open(gfa):
        if line.startswith("S\t"):
            f = line.rstrip("\n").split("\t")
            segs[f[1]] = f[2]
        elif line.startswith("P\t"):
            f = line.rstrip("\n").split("\t")
            if f[1] == path_name:
                steps = f[2].split(",")
    if steps is None:
        raise KeyError(f"no path {path_name} in {gfa}")
    out = []
    for st in steps:
        s = segs[st[:-1]]
        out.append(s if st[-1] == "+" else revcomp(s))
    return "".join(out).upper()


def region_length(path_name):
    a, b = path_name.rsplit(":", 1)[1].split("-")
    return int(b) - int(a)


def run(cmd, log, cwd=None):
    with open(log, "a") as f:
        f.write("$ " + " ".join(shlex.quote(str(c)) for c in cmd) + "\n")
        f.flush()
        p = subprocess.run([str(c) for c in cmd], stdout=f, stderr=subprocess.STDOUT, cwd=cwd)
    return p.returncode


def write_tsv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(rows)
