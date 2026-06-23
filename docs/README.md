# panvar Documentation

Documentation for the `panvar` CLI, organized into three subfolders plus a shared bibliography.

| folder / file | what's in it |
|---------------|--------------|
| **[modules/](modules/)** | one short **usage** page per module: what it does, inputs, key options, outputs, a runnable example. Start here to *run* a module. |
| **[algorithms/](algorithms/)** | the matching **mechanism + worked traces** (and deep option semantics) for each module. Read these to understand *how* a result is computed, or when a usage page links a term here. |
| **[gwas/](gwas/)** | the association workflow: [gwas/example.md](gwas/example.md) — a runnable LPA walk-through that explains each concept and correction in context, with GEMMA validation. |
| **[references.md](references.md)** | tools (GitHub) and papers behind each module. |

## Pipeline

`bubble → panphorte → call → describe → associate`, plus the `inspect` utility.

1. [modules/bubble.md](modules/bubble.md) — bubble-site extraction from a GFA (internal sort + snarl finder).
2. [modules/panphorte.md](modules/panphorte.md) — normalize tandem-repeat bubbles into a copy-number-explicit GFA.
3. [modules/call.md](modules/call.md) — graph-native SV calling (DEL/INS/INV/DUP) into a multi-sample VCF.
4. [modules/describe.md](modules/describe.md) — per-bubble k-mer / node-edge features + BIMBAM dosage genotypes.
5. [modules/associate.md](modules/associate.md) — GWAS on the describe genotypes (GLM or LMM, MAF filter, region-wide multiple testing, λ, gene column).
- [modules/inspect.md](modules/inspect.md) — clustering, path FASTA, and node/edge matrices for a called bubble.
