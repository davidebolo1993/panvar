# `panvar` documentation

Documentation for the `panvar` CLI, organized into three subfolders plus a shared bibliography.

| folder / file | what's in it |
|---------------|--------------|
| **[modules/](modules/)** | one short usage page per module: what it does, inputs, key options, outputs, a runnable example. Start here to *run* a module. |
| **[algorithms/](algorithms/)** | the matching worked traces for each module. Read these to understand how a result is computed |
| **[gwas/](gwas/)** | the association workflow with a runnable LPA walk-through |
| **[references.md](references.md)** | tools (GitHub) and papers behind each module. |

## Pipeline

`bubble → panphorte → call → describe → associate`, plus the `inspect` utility.

1. [modules/bubble.md](modules/bubble.md) — bubble-site extraction from a GFA.
2. [modules/panphorte.md](modules/panphorte.md) — normalize tandem-repeat bubbles into a copy-number-explicit GFA.
3. [modules/call.md](modules/call.md) — graph-native SV calling (DEL/INS/INV/DUP) into a multi-sample VCF.
4. [modules/describe.md](modules/describe.md) — per-bubble variant, k-mer and node-edge features.
5. [modules/associate.md](modules/associate.md) — GWAS on the describe genotypes.
- [modules/inspect.md](modules/inspect.md) — clustering, path FASTA, and node/edge matrices for a called bubble.
