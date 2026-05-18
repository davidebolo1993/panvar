# panvar Documentation

This directory contains module-focused documentation for the `panvar` CLI.

## Modules

1. `docs/modules/bubble.md`  
   Module 1 (`panvar bubble`): site extraction/refinement from precomputed `vg snarls`
2. `docs/modules/allele.md`  
   Module 2 (`panvar allele`): allele extraction and clustering per bubble
3. `docs/modules/call.md`  
   Module 3 (`panvar call`): SV calling from clustered alleles
4. `docs/modules/describe.md`  
  Module 4 (`panvar describe`): per-bubble haplotype feature tables for downstream association workflows
5. `docs/modules/example.md`  
  Worked sequence-level examples of clustering/calling/merge/debug behavior

## Doc Structure Convention

Each module page follows the same structure:

- what it does
- required inputs
- outputs
- algorithm overview
- key options
- runnable example
