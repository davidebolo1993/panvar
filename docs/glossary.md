# Glossary

Terms that mean something specific in `panvar`, and are easy to read as meaning something adjacent.
The authoritative definitions are the VCF header lines and the algorithm pages; this is the short form
and a pointer.

## Sites and support

**bubble / site** — a snarl in the graph: a pair of boundary nodes and everything between them, where
haplotypes diverge. `bubble` emits top-level snarls only, and the sites it emits are pairwise disjoint.

**interior span** — the sequence a path carries strictly between the boundaries. This, not divergence
from the reference, is what `--min-variant-bp` / `--max-variant-bp` measure, which is why they are also
spelled `--min-interior-bp` / `--max-interior-bp`. A 1 bp substitution inside a 1 kb allele has a 1 kb
interior span.

**`path_support`** — how many panel paths **traverse** the site. On a fully typed panel nearly every
haplotype crosses nearly every bubble, so this is close to the panel size and says little about any
particular allele. Use `alt_allele_support_max` / `--min-alt-support` for the question a support filter
is usually being asked.

**`distinct_alleles`** — the number of distinct walks across the site, which is what varies between
loci: C4 has a bubble with 127 distinct alleles across 131 paths.

## Copy number

The fields differ by **route**, and reading one route's field under another route's meaning is the most
common mistake here. `CN_METHOD` names the route.

**`REF_CN`** — the reference's copy number of the unit. `REF_CN_SOURCE` says whether that was an exact
self-loop traversal count (`REP_TRAVERSAL`) or a heuristic (`MAX_NODE_MULTIPLICITY`).

**`FORMAT:CN`** — this haplotype's integer copy number.

**`RU_LEN`** — repeat-unit length, one copy, in bp. Emitted **only** for `CN_METHOD=REP`, where the unit
is a literal repeat and `(CN − REF_CN) × RU_LEN` really is the haplotype's size change.

**`CN_UNIT_BP`** — the calibration constant a collapsed module's dosage was divided by
(`CN_METHOD=MODULE_BP` only). It is the **shared** per-copy content, not a whole copy, so
`(CN − REF_CN) × CN_UNIT_BP` understates a carrier — at GSTM1 by about 2.4×, and across the reference
loci by anywhere from ×1.01 to ×36.75. It is deliberately not called `RU_LEN`.

**`FORMAT:CNBP`** — the linear bp this haplotype gains (+) or loses (−) across the module against the
reference, from the spelled walk. **This is the per-haplotype event size on every route**, and on a
collapsed module it is the only one: `SVLEN` is absent there because no single record-level size
describes carriers spanning −32,030 to +18,504 bp. Checked externally against svim-asm to a median of
2–6 bp.

**`CN_MODULE_REF_BP`** — the reference bp across the whole bubble interior, the interval `CNBP` is
summed over. `END − POS == CN_MODULE_REF_BP` on a module record, which is the cheapest check that it is
placed correctly.

**`CNR_MARGIN`** — 0.5 minus the distance from a whole number of units. Near 0.5 the integer `CN` is
unambiguous; near 0 the rounding was a coin flip.

## Reconstruction levels

`benchmark` reports four, and only the last one reconstructs the emitted VCF. The first three implant
the haplotype's own true block and are therefore **ceilings**.

| level | reconstructs | bounds |
|---|---|---|
| `graph` | reference walk plus the haplotype's steps at every block sharing a node with **any** call at that bubble | can the graph hold this haplotype, and did the caller flag the divergent blocks |
| `called` | the same, restricted to blocks a **specific** record is attributed to and that reach `--min-sv-bp` | what the retained records would reach if each reproduced its block exactly |
| `carrier` | `called`, plus the requirement that this haplotype's `GT` names a record overlapping the block | the same ceiling after carrier assignment |
| `genotype` | the reference plus **only** the edits this haplotype's `GT` names | what a consumer of the VCF actually gets |

`graph ≥ called ≥ carrier` are nested. **`genotype` is not bounded by them** — it applies every record
the haplotype carries rather than only attributed eligible blocks — so the last two loss terms are
signed.

**Every QV figure this project reported before the module review pass is the `graph` column.**

**`called`** means a record shares at least one node with the block. It is attribution, not coverage:
it does not establish that the record spans the block or represents it correctly. `missed` is therefore
a firm negative and `called` an upper bound on discovery.

**baseline** — the plain reference with no edits applied. It is the denominator of `gap_closed`, and the
metric's own correctness check: a haplotype genotyped as carrying nothing must reconstruct
byte-identically to it.

## The two VCFs

**region VCF** (`call.region.vcf`) — the interpreted, human-readable output. Records are merged for
readability, which means a merged record hands every carrier the **representative's** sequence, and a
`DUP` is reconstructed by tiling an inferred span. Reconstructs 39–97% of the reference-to-truth
distance depending on locus.

**allele VCF** (`call --allele-vcf`) — one record per bubble with every distinct allele spelled out and
each haplotype's `GT` indexing its own. Reconstructs its haplotypes with **0 bp residual on all six
reference loci**. That demonstrates the representation is lossless; it is a serialization ceiling, not a
statement about call sensitivity.

## Thresholds

**`--min-sv-bp`** — the size floor for an event to count as callable. It must match between `call` and
`benchmark`, or variation that was correctly not called is charged as a miss.

**`--cluster-similarity`** — the walk-identity threshold `inspect` clusters at, derived from a shingle
Jaccard as `2J/(1+J)`.

**`--min-similarity`** — the identity at which `panphorte` treats a block as a copy of the unit. `1.0` is
exact and sequence-preserving; below that the collapse of divergent copies is lossy.
