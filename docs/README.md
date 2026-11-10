# `panvar` documentation

Documentation for the `panvar` CLI.

## Pipeline

`[rebuild] → bubble → panphorte → refine → call → describe → associate`, plus the `inspect` and `benchmark` utilities.

### Main modules

1. [rebuild.md](rebuild.md) — re-induce a locus graph too fragmented to decompose
2. [bubble.md](bubble.md) — bubble-site extraction from a GFA
3. [panphorte.md](panphorte.md) — normalize tandem-repeat bubbles
4. [refine.md](refine.md) — POA-realign bubble interiors
5. [call.md](call.md) — graph-native SV calling (DEL/INS/INV/SUB, CNV)
6. [describe.md](describe.md) — per-bubble variant, k-mer and node/edge features.
7. [associate.md](associate.md) — association testing on the `describe` features

### Utilities

- [inspect.md](inspect.md) — clustering, path FASTA, and node/edge matrices for a called bubble.
- [benchmark.md](benchmark.md) — haplotypes' large variation call's output accounts for
