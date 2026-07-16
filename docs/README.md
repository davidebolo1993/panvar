# `panvar` documentation

Documentation for the `panvar` CLI, organized into three subfolders plus a shared bibliography.

| folder / file | what's in it |
|---------------|--------------|
| **[modules/](modules/)** | one short usage page per module: what it does, inputs, key options, outputs. |
| **[algorithms/](algorithms/)** | the matching worked traces for each module.|
| **[walkthrough.md](walkthrough.md)** | the whole pipeline on one real locus (LPA), step by step with commands and plots. |
| **[gwas/example.md](gwas/example.md)** | the association (GWAS) worked example on the LPA output: a simulated Lp(a) phenotype, multiple-testing corrections, structure control. |
| **[references.md](references.md)** | tools and papers behind each module. |

## Pipeline

`bubble → panphorte → [refine] → call → describe → associate`, plus the `inspect` and `benchmark` utilities.

1. [modules/bubble.md](modules/bubble.md) — bubble-site extraction from a GFA.
2. [modules/panphorte.md](modules/panphorte.md) — normalize tandem-repeat bubbles into a copy-number-explicit GFA.
3. [modules/refine.md](modules/refine.md) — *(opt-in)* POA-realign bubble interiors to remove graph-builder alignment artifacts (DUP-safe).
4. [modules/call.md](modules/call.md) — graph-native SV calling (DEL/INS/INV/DUP) into a multi-sample VCF.
5. [modules/describe.md](modules/describe.md) — per-bubble variant, k-mer and node-edge features.
6. [modules/associate.md](modules/associate.md) — GWAS on the describe genotypes.
- [modules/inspect.md](modules/inspect.md) — clustering, path FASTA, and node/edge matrices for a called bubble.
- [modules/benchmark.md](modules/benchmark.md) — round-trip reconstruction identity of the caller's own output.
