# panvar Documentation

This directory contains module-focused documentation for the `panvar` CLI.

## Modules

1. `docs/modules/bubble.md`
   Module 1 (`panvar bubble`): site extraction/refinement from precomputed `vg snarls`
2. `docs/modules/inspect.md`
   Utility (`panvar inspect`): path FASTA and node/edge traversal matrices for one or all called bubbles
2b. `docs/modules/panphorte.md`
   Utility (`panvar panphorte`): normalize tandem-repeat bubbles into a compact GFA
3. `docs/modules/call.md`
   Module 2 (`panvar call`): graph-native structural variant calling (DEL/INS/INV/DUP) into a multi-sample VCF
4. `docs/modules/describe.md`
   Module 3 (`panvar describe`): per-bubble k-mer feature tables + pyseer-ready / sample-level GWAS inputs
5. `docs/modules/example.md`
   Worked end-to-end example

## Guides

- `docs/gwas_example.md`
  Worked **pangenome k-mer GWAS** (LPA KIV-2 copy number → Lp(a)): concepts, sample-level testing via
  cosigt aggregation, multiplicity vs presence/absence, Manhattan/QQ, and traceback to variants.
- `docs/presentation.md`
  Slide-ready tube-map schematics (panphorte normalization scenarios; called SV types in the graph).

## Doc Structure Convention

Each module page follows the same structure:

- what it does
- required inputs
- outputs
- algorithm overview
- key options
- runnable example
