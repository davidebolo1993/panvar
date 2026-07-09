# `panvar` documentation

Documentation for the `panvar` CLI, organized into three subfolders plus a shared bibliography.

| folder / file | what's in it |
|---------------|--------------|
| **[modules/](modules/)** | one short usage page per module: what it does, inputs, key options, outputs. |
| **[algorithms/](algorithms/)** | the matching worked traces for each module.|
| **[walkthrough.md](walkthrough.md)** | the whole pipeline on one real locus (**LPA**), step by step with commands and plots. |
| **[gwas/](gwas/)** | the association (GWAS) worked example on the LPA output: a simulated Lp(a) phenotype, multiple-testing corrections, structure control. |
| **[references.md](references.md)** | tools and papers behind each module. |

## Pipeline

`bubble → panphorte → call → describe → associate`, plus the `inspect` utility.

1. [modules/bubble.md](modules/bubble.md) — bubble-site extraction from a GFA.
2. [modules/panphorte.md](modules/panphorte.md) — normalize tandem-repeat bubbles into a copy-number-explicit GFA.
3. [modules/call.md](modules/call.md) — graph-native SV calling (DEL/INS/INV/DUP) into a multi-sample VCF.
4. [modules/describe.md](modules/describe.md) — per-bubble variant, k-mer and node-edge features.
5. [modules/associate.md](modules/associate.md) — GWAS on the describe genotypes.
- [modules/inspect.md](modules/inspect.md) — clustering, path FASTA, and node/edge matrices for a called bubble.
- [modules/benchmark.md](modules/benchmark.md) — round-trip QV / reconstruction identity of the caller's own output.
