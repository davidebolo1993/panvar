# LPA GWAS example inputs (simulated)

Ready-to-use phenotype and cohort tables for the [GWAS example](../../../docs/gwas/example.md), so the real
`panvar` pipeline can be exercised end-to-end without first generating a cohort.

**These are simulated**, not a real study. A structured cohort of 10,000 diploid individuals is drawn from the
LPA test graph's haplotypes (`tests/real_data/lpa.gfa.gz`); the Lp(a) phenotype is generated with an inverse
KIV-2 copy-number effect plus a subpopulation confounder, age/sex, and noise (see `make_lpa_phenotype.py`,
`--seed 42`). The point is to demonstrate the method and its corrections, not to report a finding.

| file | contents |
|------|----------|
| `samples.tsv` | cosigt `sample <tab> haplotype_1 <tab> haplotype_2` (the diploid cohort over the graph's haplotypes) |
| `pheno.quant.tsv` | `sample, phenotype (log10 Lp(a)), Age, Sex, PC1..PC10` (with ~5% NA) |
| `pheno.binary.tsv` | same columns, phenotype = high-risk case/control (0/1) |
| `pheno.quant.nopc.tsv` | quantitative phenotype without the PC columns (the naive, structure-uncorrected analysis) |

Regenerate with: `python3 tests/gwas/make_lpa_phenotype.py <copies.tsv> <out> --n 10000 --seed 42 --sim-markers 0`
(`copies.tsv` comes from running `panphorte` on the LPA graph).
