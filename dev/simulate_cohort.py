#!/usr/bin/env python3
"""Simulate a cohort over a graph's haplotypes, with a phenotype driven by a copy number.

Each sample belongs to a 1000 Genomes population and carries two distinct haplotypes drawn from that
population's paths, so allele frequencies differ between populations as in the panel. Covariates: a random
sex, and genotype PC1-PC5 drawn around the population's mean with its spread (population_pcs.tsv).

Phenotypes (--model):
  dosage    a gene's expression: copies over both haplotypes (from a pangene table, --gene) + ancestry
            offset + sex effect + noise
  isoform   a protein level each allele makes, lower with more repeat copies, as Lp(a) with KIV-2 (log
            scale): log(exp(-rate * copies_1) + exp(-rate * copies_2)) + ancestry offset + sex effect +
            noise, so the allele with fewer copies dominates

The ancestry offset confounds, as copies also differ between populations. Its continental part (random, or
--continent-offset) enters as a linear function of the PCs, so the PCs remove it; its population part
(--population-sd) enters as is, so some confounding remains. The noise makes the genetic part explain
--genetic-r2 over the whole cohort; within populations it explains less, and the log reports both.

Writes <out>.samples.tsv (sample, hap1, hap2), <out>.phenotype.tsv (sample, phenotype, sex, PC1-PC5) and
<out>.truth.tsv (sample, population, copies, each part of the phenotype).

  simulate_cohort.py --gfa graph.gfa --copies gstm.bed --gene GSTM1 \\
      --populations samples.tsv --population-pcs population_pcs.tsv -o out [--n 1000] [--seed 1]
  simulate_cohort.py --gfa lpa.gfa --copies lpa.repeats.tsv --model isoform --rate 0.1 \\
      --continent-offset AFR=0.9 ... -o out
"""
import math
import argparse
import collections
import gzip
import random
import re
import statistics


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def w_path_name(sample, hap, seqid, start, end):
    """A W line's path name as panvar gives it: sample#hap#seqid:start-end, the start made 1-based; no
    interval without a numeric start, and an end that is not a number written '*'."""
    numeric = lambda t: re.fullmatch(r"[0-9]{1,18}", t) is not None
    if not numeric(start):
        return "%s#%s#%s" % (sample, hap, seqid)
    return "%s#%s#%s:%d-%s" % (sample, hap, seqid, int(start) + 1, end if numeric(end) else "*")


def read_tsv(path):
    with _open(path) as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        return [dict(zip(hdr, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def solve(A, b):
    """A x = b by Gaussian elimination with partial pivoting."""
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[piv] = M[piv], M[c]
        for r in range(n):
            if r != c and M[c][c] != 0.0:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def main():
    ap = argparse.ArgumentParser(description="A cohort over a graph's haplotypes, its phenotype driven by a copy number.",
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=32, width=100))
    ap.add_argument("--gfa", required=True, metavar="GFA", help="graph with sample#hap#contig paths")
    ap.add_argument("--copies", required=True, metavar="TSV", help="pangene BED (with --gene), or haplotype and copies")
    ap.add_argument("--gene", help="gene counted in a pangene BED")
    ap.add_argument("--model", choices=("dosage", "isoform"), default="dosage", help="phenotype model (default dosage)")
    ap.add_argument("--rate", type=float, default=0.1, metavar="F", help="isoform: log drop per copy (default 0.1)")
    ap.add_argument("--populations", required=True, metavar="TSV", help="sample, population, superpopulation, sex")
    ap.add_argument("--population-pcs", required=True, metavar="TSV", help="PC1-PC5 mean and sd per population")
    ap.add_argument("-o", "--out-prefix", required=True, metavar="PREFIX", help="output prefix")
    ap.add_argument("--n", type=int, default=1000, help="samples (default 1000)")
    ap.add_argument("--genetic-r2", type=float, default=0.8, metavar="F", help="variance the genetic part explains (default 0.8)")
    ap.add_argument("--continent-sd", type=float, default=0.3, metavar="F", help="SD of continental offsets (default 0.3)")
    ap.add_argument("--continent-offset", action="append", default=[], metavar="C=X", help="fixed continental offset; repeatable")
    ap.add_argument("--population-sd", type=float, default=0.1, metavar="F", help="SD of population offsets (default 0.1)")
    ap.add_argument("--sex-effect", type=float, default=0.1, metavar="F", help="sex effect (default 0.1)")
    ap.add_argument("--seed", type=int, default=1, metavar="N", help="random seed (default 1)")
    a = ap.parse_args()
    if not 0.0 < a.genetic_r2 <= 1.0:
        ap.error("--genetic-r2 must be in (0, 1]")
    rng = random.Random(a.seed)

    paths = []
    with _open(a.gfa) as fh:
        for line in fh:
            if line.startswith("P\t"):
                paths.append(line.split("\t", 2)[1])
            elif line.startswith("W\t"):
                paths.append(w_path_name(*line.split("\t", 6)[1:6]))
    copies = collections.Counter()
    molecules = set()
    for row in read_tsv(a.copies):
        if "molecule" in row:
            if not a.gene:
                ap.error("a pangene table needs --gene")
            molecules.add(row["molecule"])
            if row["gene"] == a.gene:
                copies[row["molecule"]] += 1
        else:
            molecules.add(row["hapl"])
            copies[row["hapl"]] = int(row["copies"])
    pop = {r["sample"]: r for r in read_tsv(a.populations)}
    pcs = {r["population"]: r for r in read_tsv(a.population_pcs)}

    # Haplotypes usable for a population: a path of a labelled sample that the copy table covers.
    pool = collections.defaultdict(list)
    for p in paths:
        sample = p.split("#")[0]
        if "#" in p and sample in pop and p in molecules and pop[sample]["population"] in pcs:
            pool[pop[sample]["population"]].append(p)
    pops = sorted(q for q in pool if len(pool[q]) >= 2)
    weights = [len(pool[q]) for q in pops]

    samples = []
    for i in range(a.n):
        q = rng.choices(pops, weights)[0]
        h1, h2 = rng.sample(pool[q], 2)
        pc = [rng.gauss(float(pcs[q]["PC%d_mean" % k]), float(pcs[q]["PC%d_sd" % k])) for k in range(1, 6)]
        samples.append({"sample": "sim%05d" % (i + 1), "population": q,
                        "superpopulation": pcs[q]["superpopulation"], "hap1": h1, "hap2": h2,
                        "copies": copies[h1] + copies[h2], "sex": rng.randint(0, 1), "pcs": pc})
        if a.model == "dosage":
            samples[-1]["genetic"] = copies[h1] + copies[h2]
        else:
            samples[-1]["genetic"] = math.log(math.exp(-a.rate * copies[h1]) + math.exp(-a.rate * copies[h2]))

    sd_g = statistics.pstdev([s["genetic"] for s in samples]) or 1.0
    continents = sorted({s["superpopulation"] for s in samples})
    if a.continent_offset:
        fixed_offsets = {k: float(v) for k, v in (x.split("=", 1) for x in a.continent_offset)}
        offset = {c: fixed_offsets.get(c, 0.0) for c in continents}
    else:
        offset = {c: rng.gauss(0.0, a.continent_sd * sd_g) for c in continents}
    # Least squares of the continents' offsets on (1, PC1-PC5); each sample takes the fitted value.
    X = [[1.0] + s["pcs"] for s in samples]
    yv = [offset[s["superpopulation"]] for s in samples]
    coef = solve([[sum(r[i] * r[j] for r in X) for j in range(6)] for i in range(6)],
                 [sum(r[i] * v for r, v in zip(X, yv)) for i in range(6)])
    population_offset = {q: rng.gauss(0.0, a.population_sd * sd_g) for q in pops}
    for s, r in zip(samples, X):
        s["offset"] = sum(c * x for c, x in zip(coef, r))
        s["population_offset"] = population_offset[s["population"]]
        s["sex_part"] = a.sex_effect * sd_g * s["sex"]
    fixed = [s["genetic"] + s["offset"] + s["population_offset"] + s["sex_part"] for s in samples]
    var_g = sd_g ** 2
    # Noise so that the genetic part explains --genetic-r2 of genetic part + noise.
    var_noise = var_g * (1.0 - a.genetic_r2) / a.genetic_r2
    for s, f in zip(samples, fixed):
        s["noise"] = rng.gauss(0.0, var_noise ** 0.5)
        s["phenotype"] = f + s["noise"]

    out = a.out_prefix
    with open(out + ".samples.tsv", "w") as fh:
        fh.write("sample\thap1\thap2\n")
        for s in samples:
            fh.write("%s\t%s\t%s\n" % (s["sample"], s["hap1"], s["hap2"]))
    with open(out + ".phenotype.tsv", "w") as fh:
        fh.write("sample\tphenotype\tsex\t" + "\t".join("PC%d" % k for k in range(1, 6)) + "\n")
        for s in samples:
            fh.write("%s\t%.5f\t%d\t%s\n" % (s["sample"], s["phenotype"], s["sex"],
                                              "\t".join("%.6g" % v for v in s["pcs"])))
    with open(out + ".truth.tsv", "w") as fh:
        fh.write("sample\tpopulation\tsuperpopulation\tcopies\tgenetic\tcontinental_offset\tpopulation_offset"
                 "\tsex_part\tnoise\n")
        for s in samples:
            fh.write("%s\t%s\t%s\t%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n" % (
                s["sample"], s["population"], s["superpopulation"], s["copies"], s["genetic"], s["offset"],
                s["population_offset"], s["sex_part"], s["noise"]))
    ph = [s["phenotype"] for s in samples]
    r_g = statistics.correlation([s["genetic"] for s in samples], ph)
    r_c = statistics.correlation([s["copies"] for s in samples], ph)
    # Share of genetic part + noise the genetic part explains, overall and within populations.
    g = [s["genetic"] for s in samples]
    noise = statistics.pvariance([s["noise"] for s in samples])
    pop_mean = {q: statistics.fmean(s["genetic"] for s in samples if s["population"] == q)
                for q in {s["population"] for s in samples}}
    g_within = statistics.pvariance([s["genetic"] - pop_mean[s["population"]] for s in samples])
    share = statistics.pvariance(g) / (statistics.pvariance(g) + noise)
    share_within = g_within / (g_within + noise)
    print("samples: %d from %d populations; copies (sum of both haplotypes) %d-%d; r2 with the phenotype: "
          "genetic part %.2f, copies %.2f; genetic share of genetic part + noise: %.2f overall, %.2f within "
          "populations" % (
              len(samples), len({s["population"] for s in samples}), min(s["copies"] for s in samples),
              max(s["copies"] for s in samples), r_g * r_g, r_c * r_c, share, share_within))


if __name__ == "__main__":
    main()
