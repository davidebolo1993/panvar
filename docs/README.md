# panvar Documentation

This directory contains module-focused documentation for the `panvar` CLI.

## Modules

1. `docs/modules/bubble.md`
   Module 1 (`panvar bubble`): bubble-site extraction/refinement from any GFA via internal sort + cactus snarl finding
2. `docs/modules/panphorte.md`
   Module 2 (`panvar panphorte`): normalize tandem-repeat bubbles into a compact, copy-number-explicit GFA, internally re-sorted + re-snarled for `call`
3. `docs/modules/call.md`
   Module 3 (`panvar call`): graph-native structural variant calling (DEL/INS/INV/DUP) into a multi-sample VCF
4. `docs/modules/describe.md`
   Module 4 (`panvar describe`): per-bubble k-mer feature tables / sample-level GWAS inputs
5. `docs/modules/inspect.md`
   Utility (`panvar inspect`): clustering, path FASTA and node/edge traversal matrices for one or all called bubbles


## Guides

- `docs/gwas_example.md`
  Worked **pangenome association** (LPA KIV-2 copy number → Lp(a)): concepts, sample-level testing via cosigt aggregation, multiplicity vs presence/absence on both the k-mer and node/edge substrates, Manhattan/QQ, and traceback to variants.
- `docs/algorithm_example.md`
  Tiny hand-traced datasets that walk through the internals end to end: bubble normalization, the panphorte approximate collapse, the `call` event typing / merge / copy-number arithmetic, and the `describe` syncmer sampling and discriminative filter. Read this to see *how* a result is computed, not just *what* the command prints.
- `docs/references.md`
  Background literature for the concepts used here (closed syncmers, cactus/ultrabubble snarls, Lp(a)/KIV-2 copy-number biology, copy-number vs presence/absence association).

## Reproducing the results

The example commands write under `results/` (gitignored). Regenerate everything — every gene region
(`lpa`, `c4`, `gstm1`, `cyp2d6`) plus the synthetic smoke, with plots and copy-number validation — with:

```bash
conda activate base            # for Rscript (plots); optional
bash scripts/regen_results.sh  # or: scripts/regen_results.sh lpa c4 gstm1 cyp2d6 synthetic
```

Real inputs are the gzipped graphs in `tests/real_data/*.gfa.gz`; copy-number ground truth is in
`tests/real_data/{c4,cyp2d6,gstm1}.bed` and `lpa.repeats.tsv`.

### Copy number: one method per locus topology

How a locus is represented in the pangenome decides how `call` reads its copy number — and which graph
it reads. The driver picks the right recipe per gene:

| Region | Topology | Call substrate | CN method | Concordance vs ground truth |
|--------|----------|----------------|-----------|------------------------------|
| **LPA** (KIV-2) | tandem repeat | `panphorte` graph | `--cn-from-multiplicity` (REP self-loop) | **465/465 = 100%** |
| **C4** (RCCX) | PGGB-collapsed paralog | `bubble` graph | `--cn-from-coverage` (full-walk) | **131/131 = 100%** |
| **GSTM1** | deletion/CNV (segdup) | `bubble` graph | `--cn-from-multiplicity` (peak) | **159/159 = 100%** (constant paralog baseline) |
| **CYP2D6** | PGGB-collapsed paralog | `bubble` graph | `--cn-from-coverage` (full-walk) | concordant vs D6+D7; residual = unannotated CYP2D8P/hybrid |

The principle: PGGB collapses **identical** copies onto shared nodes (copy number = node multiplicity);
`panphorte` collapses a **variable tandem** into one REP node. So tandem loci are called on the
`panphorte` graph, while PGGB-collapsed paralog clusters are called on the unfolded `bubble` graph (where
the multiplicity is intact). The regen driver writes a per-gene concordance dotplot
(`results/cn_correlation.png`, called CN from the VCF vs the BED, faceted by gene, reference haplotype
highlighted). See [modules/call.md](modules/call.md) for the full method description.

![Copy-number concordance: panvar calls vs ground truth](../results/cn_correlation.png)

Each point is a haplotype's recovered **gene** copy number (the absolute VCF count minus the constant
folded-paralog baseline, annotated per facet) vs the BED truth. LPA and C4 have baseline 0 and sit exactly
on `y = x`; GSTM1 (+2: its GSTM2–5 segdup paralogs) and CYP2D6 (+1: CYP2D8P) report the total collapsed
module, so a constant baseline is removed. CYP2D6's residual scatter is the variable, unannotated
CYP2D8P/2D7-hybrid copies (see call.md).

### Gene annotation (`--gtf`)

An optional reference-coordinate GTF (Ensembl/GENCODE) projects gene names onto the graph across the
pipeline. It needs a **PanSN** reference path (`sample#hap#contig:start-end`, so chromosome + absolute
start are known); convert older underscore-named graphs first with `scripts/pansn_rename.py`. lncRNAs are
skipped. When supplied:

- `bubble --gtf` and `panphorte --gtf` each write `<prefix>.bandage_genes.csv` (`Name,Colour,Gene`) to
  highlight genes per bubble in Bandage — both, because panphorte's collapse renumbers nodes.
- `call --gtf` adds `INFO=GENES` to every record, writes `<prefix>.node_genes.tsv`, and a per-gene DUP
  copy-number table `<prefix>.dup_gene_cn.tsv` that **separates** paralogs the graph folds together
  (reliable rows, e.g. CYP2D6 vs CYP2D7/2D8P) and honestly **collapses** near-identical ones it can't
  resolve (unreliable rows with the module total, e.g. C4A;C4B). See [modules/call.md](modules/call.md).
- the GWAS traceback (`scripts/gwas_demo.py --node-genes <call.node_genes.tsv>`) names the gene behind an
  association. `describe` itself is k-mer based and does **not** consume the GTF.

## Doc Structure Convention

Each module page follows the same structure:

- what it does
- required inputs
- key options
- outputs
- algorithm overview
- runnable example (reads a `tests/real_data/*.gfa.gz` graph, writes under `results/real_data/<region>/`)
