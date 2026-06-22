# panvar Documentation

Documentation for the `panvar` CLI, organized into three subfolders plus a shared bibliography.

| folder / file | what's in it |
|---------------|--------------|
| **[modules/](modules/)** | one short **usage** page per module: what it does, inputs, key options, outputs, a runnable example. Start here to *run* a module. |
| **[algorithms/](algorithms/)** | the matching **mechanism + worked traces** (and deep option semantics) for each module. Read these to understand *how* a result is computed, or when a usage page links a term here. |
| **[gwas/](gwas/)** | the association workflow: [gwas/primer.md](gwas/primer.md) (GWAS concepts from scratch) and [gwas/example.md](gwas/example.md) (a runnable LPA walk-through + GEMMA validation). |
| **[references.md](references.md)** | tools (GitHub) and papers behind each module. |

## Pipeline

`bubble → panphorte → call → describe → associate`, plus the `inspect` utility.

1. [modules/bubble.md](modules/bubble.md) — bubble-site extraction from a GFA (internal sort + snarl finder).
2. [modules/panphorte.md](modules/panphorte.md) — normalize tandem-repeat bubbles into a copy-number-explicit GFA.
3. [modules/call.md](modules/call.md) — graph-native SV calling (DEL/INS/INV/DUP) into a multi-sample VCF.
4. [modules/describe.md](modules/describe.md) — per-bubble k-mer / node-edge features + BIMBAM dosage genotypes.
5. [modules/associate.md](modules/associate.md) — GWAS on the describe genotypes (GLM or LMM, MAF filter, region-wide multiple testing, λ, gene column).
- [modules/inspect.md](modules/inspect.md) — clustering, path FASTA, and node/edge matrices for a called bubble.

## Reproducing the results

The per-gene driver scripts regenerate everything under `results/` (gitignored), data only (plotting
commands are included but commented):

```bash
PYTHON=~/miniconda3/bin/python scripts/genes/lpa.sh      # bubble→inspect→panphorte→inspect→call→describe→associate
scripts/genes/c4.sh ; scripts/genes/gstm1.sh ; scripts/genes/cyp2d6.sh
```

Real inputs are the gzipped graphs in `tests/real_data/*.gfa.gz`; copy-number ground truth is in
`tests/real_data/{c4,cyp2d6,gstm1}.bed` and `lpa.repeats.tsv`. R plots use the `scripts/plot_*.R` helpers
(need `Rscript` + `ggplot2`).

## Copy number: one method per locus topology

How a locus is represented in the pangenome decides how `call` reads its copy number — and which graph it
reads. This table is the canonical reference (the module pages link here):

| Region | Topology | Call substrate | CN method | Concordance vs truth |
|--------|----------|----------------|-----------|----------------------|
| **LPA** (KIV-2) | tandem repeat | `panphorte` graph | `--cn-from-multiplicity` (REP self-loop) | **465/465 = 100%** |
| **C4** (RCCX) | PGGB-collapsed paralog | `bubble` graph | `--cn-from-coverage` (full-walk) | **131/131 = 100%** |
| **GSTM1** | deletion/CNV (segdup) | `bubble` graph | `--cn-from-multiplicity` (peak) | **159/159 = 100%** |
| **CYP2D6** | PGGB-collapsed paralog | `bubble` graph | `--cn-from-coverage` (full-walk) | concordant vs D6+D7; residual = unannotated CYP2D8P/hybrid |

The principle: PGGB collapses **identical** copies onto shared nodes (copy number = node multiplicity);
`panphorte` collapses a **variable tandem** into one REP node. So tandem loci are called on the `panphorte`
graph, PGGB-folded paralog clusters on the unfolded `bubble` graph. Mechanics + traces:
[algorithms/call.md](algorithms/call.md#copy-number--three-ways).

## Gene annotation (`--gtf`)

An optional reference-coordinate GTF projects gene names onto the graph (needs a **PanSN** reference path;
lncRNAs skipped). `bubble`/`panphorte` write Bandage gene CSVs; `call --gtf` adds `INFO=GENES`,
`node_genes.tsv`, and the per-gene DUP table `dup_gene_cn.tsv` (separating resolvable paralogs, collapsing
near-identical ones). `node_genes.tsv` is the node→gene bridge for the GWAS traceback —
`associate --node-genes` joins it to emit a `gene` column. Details:
[modules/call.md](modules/call.md#gene-annotation---gtf), trace:
[algorithms/call.md](algorithms/call.md#gene-annotation-trace---gtf).
