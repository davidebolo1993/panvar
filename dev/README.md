# Dev scripts

Scripts that check panvar against outside truth: truth sets, other callers, simulations. They are not part of the pipeline. Run them from the repository root; every argument is optional, with defaults in the repository (each script's header lists them), and an empty argument (`""`) keeps its default.

## Copy number against the truth sets

`dev/call_truth.sh` compares the example's copy numbers with the truth sets in `tests/real_data`, haplotype by haplotype: each locus (its `CNV` record) and each gene (its rows in `dup_gene_cn.tsv`) against pangene's gene counts, and the tandem repeat against a table of copy counts. `dev/compare_copy_number.py` writes one table per comparison (a gene with rows at several sites needs `--site <bubble_id>`), and `dev/plot_cn_correlation.R` plots them, one panel per locus or gene, in `tests/examples/call/copy_number`.

```bash
bash dev/call_truth.sh        # after scripts/call_examples.sh
# one comparison alone, and its plot
python3 dev/compare_copy_number.py --truth tests/real_data/paralogs/c4.bed --vcf tests/examples/call/c4/c4.region.vcf \
  --genes C4A,C4B --label c4 -o c4.cn.tsv
Rscript dev/plot_cn_correlation.R --table c4.cn.tsv --out c4.cn
```

## `call` against svim-asm

`dev/benchmark_callers.sh` calls svim-asm on the haplotypes of every graph of `tests/examples/benchmark` (`dev/assembly_variants.py --sv-tables`: minimap2 alignment of each haplotype to the reference, svim-asm per haplotype, truvari collapse across haplotypes), scores `call` and svim-asm against `present.tsv`, and plots precision, recall and F1 per graph (`callers.png`) and by event size (`callers.by_size.png`, with `dev/plot_sv_callers.R`). It needs svim-asm and truvari.

`dev/compare_sv_callers.py` scores each caller at two levels: each haplotype's own calls (`call`'s `carrier_events.tsv`; svim-asm per haplotype) and the merged output (`call`'s VCF, copy numbers as `CN` against `REF_CN`; svim-asm after truvari collapse). The matching is the script's own, following truvari bench's rules and defaults without its sequence comparison (truvari bench cannot read `present.tsv`, per-haplotype carriers or copy numbers). A call matches a truth event of `present.tsv` when it is on the same haplotype, within `--refdist` (500 bp), of a size within `--pctsize` (0.7), and in the same direction. Matching is one to one first; then, since the alignment splits a tandem gain into one event per copy, an event left over may claim several adjacent ones of the other side that together reach its size. Only calls inside the stretch benchmark scored count. In `call`'s own events, a carrier's events inside its copy stretch are part of it, as in benchmark; in the VCF, a copy-number record is a count, so the records inside it are calls of their own.

```bash
bash dev/benchmark_callers.sh        # after scripts/benchmark_examples.sh; svim-asm and truvari on PATH
bash dev/benchmark_callers.sh "" "" "" <svim-asm env>/bin <truvari env>/bin   # or their folders
# one graph alone: svim-asm's calls, both callers scored, the plot
python3 dev/assembly_variants.py --gfa tests/examples/refine/c4/c4.normalized.sorted.gfa --reference GRCh38 \
  --tools <svim-asm env>/bin --merge-tools <truvari env>/bin --out-dir svim/c4 --sv-tables svim/c4/c4
python3 dev/compare_sv_callers.py --label c4 --benchmark tests/examples/benchmark/c4/c4 --calls tests/examples/call/c4/c4 \
  --svim svim/c4/c4 -o c4.callers.tsv
Rscript dev/plot_sv_callers.R --out callers.png c4.callers.tsv
```

## Inputs of the `associate` example

`dev/associate_inputs.sh` makes the inputs of `scripts/associate_examples.sh` in `tests/synthetic_data/associate/<locus>/`: the cohort's samples, phenotypes and covariates, and each matrix at the sample and the haplotype level. It needs svim-asm and truvari. Per locus:

1. `dev/simulate_cohort.py` draws 1,000 samples, each with two haplotypes of one 1000 Genomes population, and covariates sex and PC1-PC5 around the population's values ([MAGE](https://doi.org/10.1038/s41586-024-07882-7) principal components, `tests/real_data/populations`). The phenotype is a copy number's effect: an expression following a gene's copies, or a protein's level in the blood set mostly by the allele with fewer repeat copies. To it are added an ancestry offset, a sex effect and noise. The offset has a continental part that is linear in the PCs, which they remove, and a smaller population part they capture only partly, so some confounding remains with them.
2. `describe --samples` writes the sample-level variant matrix of `call`'s records.
3. `dev/assembly_variants.py` aligns each haplotype to the reference with minimap2 and calls SNVs with paftools and SVs with svim-asm. The SVs are merged across haplotypes with bcftools merge and truvari collapse, at truvari's default thresholds but with no upper size limit (by default truvari leaves SVs over 50 kb unmerged). A haplotype reads `NA` where its alignment deletes the site.

The phenotype files are three: sex and the PCs (`cohort.phenotype.tsv`); sex only (`no_pcs`); and the PCs with the simulated genetic part and the population offsets the PCs miss (`null`). `cohort.truth.tsv` holds what was simulated per sample.

```bash
bash dev/associate_inputs.sh         # after scripts/refine_examples.sh and call_examples.sh; svim-asm and truvari on PATH
# the cohort alone
python3 dev/simulate_cohort.py --gfa tests/examples/refine/gstm/gstm.normalized.sorted.gfa \
  --copies tests/real_data/paralogs/gstm.bed --gene GSTM1 --populations tests/real_data/populations/samples.tsv \
  --population-pcs tests/real_data/populations/population_pcs.tsv -o cohort
```

## `associate` against the tools it follows

`dev/gemma/validate_gemma.sh` compares the linear and mixed models with GEMMA on the same BIMBAM files, in `tests/examples/gemma`. A mixed model needs relatives, so at each locus it pairs the example's samples into parents and gives each couple three children, one haplotype from each parent, with a phenotype following the copy-number record plus a background shared within families. Both models agree with GEMMA's to within 1e-4 of a standard error, and the variance share the kinship explains to 0.001. Rows with missing dosages are not compared: GEMMA fills a missing dosage with the row's mean, as panvar's mixed model does but its linear model does not. GEMMA reads a dosage above 2 as invalid and skips the row without a message, so it gets `describe --scale-dosage` matrices.

```bash
bash dev/gemma/validate_gemma.sh     # gemma on PATH; reads tests/synthetic_data/associate
python3 dev/gemma/make_families.py tests/synthetic_data/associate/gstm families/gstm   # the families alone
```

`dev/validate_associate/validate.sh` checks each method against the tool it follows, on data it simulates (`simulate.R`: 600 unrelated samples with sex and two PCs, 90 records from common to rare, among them copy numbers and records tagging one, with a measured and a case/control phenotype; and 150 families of four with their pedigree kinship). `reference.R` recomputes every result with the reference and prints the largest difference over the rows compared; the run fails if one is outside its tolerance. The table is written to `tests/examples/validate_associate/comparison.tsv`. It needs R with the SKAT, susieR and SPAtest packages, and GEMMA.

```bash
bash dev/validate_associate/validate.sh               # simulate.R, the runs, then reference.R
bash dev/validate_associate/validate.sh --calibrate   # also the null calibration below (much longer)
bash dev/validate_associate/calibrate.sh build/panvar tests/examples/validate_associate   # the calibration alone, after validate.sh
```

| method | reference | largest difference |
|---|---|---|
| linear model: `beta`, `se`, `p` | R `lm` | 1e-5 of a standard error; 1e-6 in −log10 p |
| dose classes: `p_genotypic`, `p_nonadditive` | R `anova` of the additive and class models | 1e-6 in −log10 p |
| `--rint` | R `qnorm((rank(y) - 0.5) / n)`, then `lm` | 2e-6 |
| `p_bonf`, `q_bh` | R `p.adjust` | 5e-6, relative |
| Meff | Li and Ji's formula on the correlation matrix's eigenvalues, in R | identical |
| logistic score test, normal tail | R `glm`, `anova(test = "Rao")` | 2e-6 in −log10 p |
| saddlepoint tail | SPAtest's ScoreTest_SPA | 5e-5 |
| counted tails | the carriers' outcomes enumerated in R (the hypergeometric without covariates) | 2e-6 |
| `log_or` | R `glm` | 3e-6 of a standard error |
| burden, SKAT, SKAT-O, ACAT-O | SKAT package, `SKAT()`; ACAT-O over its six p-values | burden and SKAT 1e-6; SKAT-O and ACAT-O 0.006 in −log10 p |
| fine-mapping | susieR `susie(L = 10, refine = TRUE)` | the same credible sets; their PIPs 1e-6; other PIPs up to 0.05 |
| conditional selection | forward selection with backward elimination written out with R `lm` | the same signals; conditional p 2e-6 |
| mixed model | GEMMA `-lmm 1` | 7e-6 of a standard error; 2e-6 in −log10 p |
| `--lmm-fast` | the EMMAX procedure written out in R: the null REML variance, then generalized least squares | 2e-6 |
| `lmm_h2` | GEMMA's pve in the null model | 0.001 |

- The SKAT package computes its tails by Davies' method at an absolute accuracy of 1e-6, a few percent of a p near 1e-5: the reference asks the same routine for 1e-12. SKAT-O's remaining 1% is the package's numerical integration.
- SPAtest reads genotypes 0 to 2, so copy numbers are rescaled to that range, which leaves a saddlepoint p unchanged.
- Case/control groups are compared where `associate` uses the plain asymptotic test: more than `--exact-carriers` carriers and no record past |z| = 2.
- Outside the credible sets, a weak effect can differ: `associate` searches each effect's prior variance over its whole range and keeps the effect when any variance raises the evidence, where susieR's optimiser can drop it, so a little probability spreads over rows in no set.
- GCTA-COJO works from summary statistics; on the samples' own data its procedure is ordinary regression with the selected records as covariates, which is what the R version does.

With `--calibrate`, `calibrate.sh` then runs null simulations for what no tool computes, and writes `tests/examples/validate_associate/calibration.txt`:

| on noise phenotypes | false hits | target |
|---|---|---|
| permutation threshold (200 runs, 200 shuffles each) | 5.5% of runs | 5% |
| logistic, normal or saddlepoint tail, 30% / 3% cases | 4.7% / 5.1% below 0.05; 0.9% / 0.9% below 0.01 | 5%, 1% |
| logistic, counted tail, 30% / 3% cases | 2.1% / 1.8% below 0.05; 0.3% / 0.2% below 0.01 (conservative, as counts of few outcomes are) | 5%, 1% |
| `--group`, 2% cases: burden / SKAT / SKAT-O | 4.7% / 9.5% / 7.8% below 0.05; 0.8% / 4.1% / 3.4% below 0.01 | 5%, 1% |
| `delta_r2` interval, 300 noise records against one | 0 of 100 exclude 0 | at most 5% |

## Read recruitment against simulations

`dev/recruit_simulate.sh` makes the read sets of the recruitment benchmark in `tests/examples/recruit_sim`, from the paralog and tandem-repeat graphs of `tests/real_data`. Per locus, a sample and its relatives are held out of the graph (`dev/recruit_sources.py held-out`), and the locus is indexed without them against two backgrounds: GRCh38 and CHM13 (the main arm), and GRCh38 alone (so CHM13 is a genome the index never saw). Reads come from the held-out sample's two haplotypes (the positives), from GRCh38 2 Mb around the locus, from the decoy windows, from 15 kb across each locus edge, and from 1x of GRCh38 and of CHM13 with the loci masked; each set in its own file, so the file is the truth. Short pairs come from wgsim and ART's HiSeq 2500 profile, long reads from pbsim3's ONT model and `dev/simulate_hifi.py`; ART pairs and ONT reads carry their true span in their names (`dev/sim_truth.py`). With `WG=1` it adds a whole-GRCh38 arm and whole-CHM13 negatives.

`dev/recruit_benchmark.sh` recruits every set in every arm and `dev/recruit_score.py` scores them: certain and recruited share of the locus fragments, fragments from elsewhere taken, and boundary reads by whether they lie inside, across or outside the edge, with seconds per set. `dev/recruit_calibrate.py` pools the negatives into the false-certain rate per technology at each private-site cut, with its 95% bound. `dev/private_elsewhere.cpp` finds a locus's private k-mers on chromosomes a background lacks, and `dev/relabel_index.cpp` relabels them shared, to emulate a whole-genome background from a partial one.

```bash
bash dev/recruit_simulate.sh "" "" GRCh38.fa.gz CHM13.fa.gz        # wgsim, art_illumina and pbsim3 on PATH (PBSIM_MODEL)
bash dev/recruit_benchmark.sh
python3 dev/recruit_calibrate.py tests/examples/recruit_sim/benchmark
```

`dev/recruit_seeds.sh` recruits every set with each screen (`recruit --screen`: closed syncmers, minimizers at the same density, every k-mer) and tables certain and recruited fragments and seconds per screen. `dev/recruit_baselines.sh` takes the same sets with BBDuk and mirabait (every 31-mer of the locus haplotypes as bait), and `dev/locityper_compare.sh` with Locityper's recruitment on the same background (its jellyfish counts of GRCh38); both table each tool beside recruit's certain and recruited fragments.

```bash
bash dev/recruit_seeds.sh tests/examples/recruit_seeds
bash dev/locityper_compare.sh tests/examples/recruit_locityper grch38.jf targets/      # locityper and jellyfish on PATH
bash dev/recruit_baselines.sh tests/examples/recruit_baselines targets/ "" tests/examples/recruit_locityper   # bbduk.sh (and mirabait)
```

The targets folder holds each locus's haplotypes as `<locus>.fa`.

## Recruitment on real reads

`dev/recruit_assembly_truth.sh` scores a real sample against its own assembly: it recruits the sample's reads (CRAM, BAM or FASTQ), finds each locus on the sample's two haplotypes by aligning the locus's GRCh38 interval to them, aligns every read to each haplotype on its own (aligned to both at once, a read the two share has MAPQ 0), and `dev/recruit_assembly_score.py` classes each fragment by where its mates align: in the locus (MAPQ 20 or more, or below), elsewhere, or unaligned, against recruit's fate. The indexes must be built from graphs without the sample and its relatives (`dev/recruit_sources.py held-out`). `dev/recruit_real.sh` is the lighter check: a sample's GRCh38 alignments around each locus and in its decoy windows, recruited and compared with where the aligner put them (`dev/recruit_real_score.py`). `dev/recruit_formats.sh` aligns a FASTQ pair and writes it in every layout (BAM sorted by coordinate with and without index, by name, CRAM, unaligned BAM): each must give the FASTQ's table.

```bash
THREADS=16 CRAM_REFERENCE=GRCh38_full_analysis_set_plus_decoy_hla.fa bash dev/recruit_assembly_truth.sh truth \
  sample.final.cram sample.paternal.fa.gz sample.maternal.fa.gz GRCh38_full_analysis_set_plus_decoy_hla.fa index/*.plx
bash dev/recruit_formats.sh reads_1.fq.gz reads_2.fq.gz GRCh38.fa.gz formats index/*.plx       # minimap2 and samtools on PATH
```

## Copy number from certain reads

`dev/recruit_gene_cn.py` reads each gene's copy number in a sample from its certain reads alone (`recruit --reads-out --certain-only`): the k-mers whose count in the graph's annotated haplotypes equals each haplotype's copies of the gene (in at least 90% of them) read one copy each, so their median depth over that of the backbone k-mers (once in every haplotype), times two, is the diploid copy number. The truth is the sample's own rows in the gene tables of `tests/real_data` (or the repeat's copy counts). `dev/recruit_hard_cases.sh` runs it on the simulated samples and on a sample without a gene that most haplotypes carry, held out of its graph, which must read 0.

```bash
bash dev/recruit_hard_cases.sh tests/examples/recruit_cn
python3 dev/recruit_gene_cn.py --gfa held_out.gfa --sample <name> --genes tests/real_data/paralogs/<locus>.bed --reads certain_reads/*.fq.gz
```

