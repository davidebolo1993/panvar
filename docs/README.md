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
  Worked **pangenome k-mer GWAS** (LPA KIV-2 copy number → Lp(a)): concepts, sample-level testing via cosigt aggregation, multiplicity vs presence/absence, Manhattan/QQ, and traceback to variants.
- `docs/algorithm_example.md` XXX


## Doc Structure Convention

Each module page follows the same structure:

- what it does
- required inputs
- outputs
- algorithm overview
- key options
- runnable example
