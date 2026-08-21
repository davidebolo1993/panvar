# Association with panvar — a worked example

The association half of the pipeline, run on the output of the [walkthrough](walkthrough.md). Variant calling and feature description produce the genotype matrices; this page tests them.

## The cohort

`tests/gwas/make_lpa_phenotype.py` reads each real haplotype's copy number at the tandem array, splits the haplotypes into subpopulations with different array frequencies, pairs them into diploid samples, and simulates a quantitative phenotype driven by the diploid copy number plus a subpopulation offset. The offset is what makes structure correction necessary rather than decorative.

Two cohorts exist, and the numbers on this page come from the second:

- `tests/gwas/lpa/` is committed and runnable directly: a 6,000-sample cosigt table with quantitative and binary phenotypes, each with and without ancestry principal components.
- `scripts/regen_results.sh lpa` generates a smaller cohort under `results/real_data/lpa/gwas/real/` and runs everything below on it. That run produced every figure here: 300 samples, of which 285 have complete phenotype and covariate data and are used.

The genotype matrices live under `results/real_data/lpa/gwas/desc/`, one folder per substrate, at both haplotype and diploid-sample level.

## Test the variant unit

```bash
panvar associate \
  --genotypes results/real_data/lpa/gwas/desc/sample/variant/bimbam_variant.bimbam.gz \
  --samples   results/real_data/lpa/gwas/desc/sample/variant/samples.txt.gz \
  --feature-annot results/real_data/lpa/gwas/desc/sample/variant/feature_annot.variant.tsv.gz \
  --phenotype results/real_data/lpa/gwas/real/pheno.quant.tsv \
  -o results/real_data/lpa/gwas/associate/assoc_variant_quant
```

The unit is auto-detected as `variant` from the sidecar, so the test count is the number of structural-variant calls rather than the number of markers. Here 12 calls are tested after the frequency filter drops 5.

The array's duplication, `bubble6_DUP_5100`, comes out at `p = 2.9e-09`. Because the variant substrate carries the duplication's copy number as its dosage rather than a presence flag, that test is on the copy number itself.

![Region scan, variant substrate](img/assoc_variant_quant.manhattan.png)

## Test the graph and k-mer substrates

```bash
panvar associate \
  --genotypes results/real_data/lpa/gwas/desc/sample/graph/bimbam_graph.bimbam.gz \
  --samples   results/real_data/lpa/gwas/desc/sample/graph/samples.txt.gz \
  --feature-annot results/real_data/lpa/gwas/desc/sample/graph/feature_annot.graph.tsv.gz \
  --phenotype results/real_data/lpa/gwas/real/pheno.quant.tsv \
  -o results/real_data/lpa/gwas/associate/assoc_graph_quant
```

These substrates keep every node, edge and k-mer test, so they fine-map within the locus at the cost of many correlated tests.

![Region scan, graph substrate](img/assoc_graph_quant.manhattan.png)

The graph scan's top feature is node `5100` — the repeat-unit node itself, and the same node the variant-level record is anchored on, at the same p-value. The two substrates are describing one thing at different granularities.

![Region scan, k-mer substrate](img/assoc_kmers_quant.manhattan.png)

The k-mer substrate tells the same story more finely, its strongest marker reaching `p = 7.7e-10`, with the repeat-unit k-mers all peaking together.

Across the three, on the same cohort:

| substrate | tested | `Meff` | method | significant at `p_bonf_meff` | at `q_bh < 0.05` |
|---|---|---|---|---|---|
| variant | 12 | 10 | eigenvalue | 1 | 5 |
| graph | 2,115 | 9 | distinct bubbles | 184 | 977 |
| k-mers | 1,990 | 9 | distinct bubbles | 942 | 1,178 |

The variant unit is the honest denominator: one test per call. The other two are for locating the signal within the locus, not for counting how many independent things were found.

## Reading the corrections

- `p_bonf` is the raw Bonferroni correction over every feature tested. It carries a formal guarantee and is conservative when the features are correlated, which within one locus they always are.
- `p_bonf_meff` scales by `Meff` instead. In the variant unit `Meff` is the phenotype-blind eigenvalue estimate of the genotype correlation matrix, which never looks at the phenotype, so the threshold it implies cannot be circular. Here it reads 10 against 9 from LD clumping, and `meff_method` records which drove the threshold. In the feature units `Meff` is the number of distinct bubbles, which is a biological grouping rather than a statistical one.
- `q_bh` is the Benjamini-Hochberg false-discovery rate, and is the primary control.
- `p_conditional` and `cond_role` answer a different question: not whether a marker is associated, but whether it is associated independently of the others. The array's duplication is the sole independent signal here; every other feature collapses under conditioning on it.

## The pipeline, stage by stage

`scripts/plot_associate_pipeline.R` draws one facet per stage, recolouring the same markers by what survives each, so the funnel from every test down to the independent signal is visible in one figure.

1. `TEST` — every marker is tested and shown at its raw p-value.
2. `FILTER MAF` — markers below `--min-maf` are greyed out.
3. `CLUMP` — variant unit only: variants are grouped by genotype correlation into leads and shadows.
4. `CORRECT` — both thresholds are drawn, and each marker is coloured by the strictest it passes.
5. `CONDITION` — conditional p-values after the forward-stepwise stage: shadows collapse and only independent signals stay above the line.

![Association pipeline, variant tier](img/assoc_variant_quant.pipeline.png)

![Association pipeline, graph tier](img/assoc_graph_quant.pipeline.png)

![Association pipeline, k-mer tier](img/assoc_kmers_quant.pipeline.png)

The feature tiers have no `CLUMP` facet, since clumping applies to the variant unit only.

## Reading `lambda_gc` at one locus

`lambda_gc` reads the median chi-square on the assumption that most tests are null. That holds genome-wide and fails at a single locus, where a real signal and everything in linkage disequilibrium with it can be most of the tests. It reads 5.85 on the variant scan above and 10.13 on the graph scan — those measure the effect, not inflation. The run summary says which situation it is in and only calls it an inflation estimate when there are enough tests and few enough of them significant.

## Structure correction, on a panel where it can be judged

To show what ancestry covariates buy, the driver also builds a synthetic panel of mostly-null markers alongside the real copy-number dosage, so a genome-wide-like setting exists where `lambda_gc` does mean what it usually means.

```bash
# no ancestry covariates
panvar associate --genotypes <panel.bimbam.gz> --samples <...> --feature-annot <...> \
  --phenotype <pheno.quant.nopc.tsv> --min-maf 0.02 -o sim_naive
# the same phenotype, with the principal components as covariates
panvar associate --genotypes <panel.bimbam.gz> --samples <...> --feature-annot <...> \
  --phenotype <pheno.quant.tsv>      --min-maf 0.02 -o sim_pc
# a linear mixed model instead, with an external relationship matrix
panvar associate --genotypes <panel.bimbam.gz> --samples <...> --feature-annot <...> \
  --phenotype <pheno.quant.nopc.tsv> --model lmm --kinship <grm.tsv> --min-maf 0.02 -o sim_lmm
```

| analysis | `lambda_gc` | significant at `p_bonf_meff` |
|---|---|---|
| no ancestry covariates | 2.22 | 122 |
| with principal components | 1.02 | 1 |
| linear mixed model | 2.22 | 122 |

The principal components do the work: inflation collapses to nominal and the 122 apparent hits reduce to the one real signal.

The mixed model does not, on this panel. It runs with the supplied relationship matrix, but the fitted variance ratio comes out around 1.3e5, meaning essentially all variance is residual and the random effect is negligible — so it behaves like the uncorrected model and reproduces its numbers. That is a property of this synthetic panel, whose structure is carried by a handful of markers rather than spread across a genome, not evidence about mixed models generally. Treat the mixed-model path as experimental: its only external check is a correlation against an established implementation, which cannot detect a systematic difference in effect size, standard error or p-value.

## See also

- [associate](modules/associate.md) — units, options, output columns and limitations.
- [algorithms/associate.md](algorithms/associate.md) — which test is used where, and why.
- [walkthrough](walkthrough.md) — the calling and feature-description steps that produce these matrices.
