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

### Read recruitment

`background → recruit`: count the genome once, then sort each sample's reads to the locus graphs without aligning.

1. [background.md](background.md) — every k-mer of a set of genomes counted once, with anchors on a reference
2. [recruit.md](recruit.md) — locus indexes against the background, and each read called certain, ambiguous or unrecruited per locus

### Utilities

- [inspect.md](inspect.md) — clustering, path FASTA, and node/edge matrices for a called bubble.
- [benchmark.md](benchmark.md) — haplotypes' large variation call's output accounts for
