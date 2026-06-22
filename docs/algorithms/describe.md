# describe — algorithm & worked example

Mechanism and a hand-traced example for **Module 4**. For usage/flags see [modules/describe.md](../modules/describe.md);
citations: [references.md](../references.md#describe).

## Terms

- **k-mer** — a length-`k` DNA substring; **canonicalised** to `min(forward, reverse-complement)` of its
  2-bit code (`A=0,C=1,G=2,T=3`), so a k-mer and its reverse complement count as one feature.
- **closed syncmer** — a sampling rule that keeps a k-mer only when its smallest internal `s`-mer sits at
  either end; it samples sequence space evenly and degrades gracefully under substitutions (one base change
  disturbs only the overlapping syncmers).
- **node/edge dosage** — how many times a path traverses a graph node / an oriented step-to-step edge
  (`from±>to±`); a graph-local, coordinate-keyed feature complementary to k-mers.

## K-mer encoding & syncmer selection

Each k-mer is encoded as a 2-bit integer and canonicalised; the matrix uses compact names `K1, K2, …`, the
feature map resolves them back to the DNA string + node provenance.

**Canonicalisation** (k=4): `GACT = 2 0 1 3 = 135`; revcomp `AGTC = 0 2 3 1 = 45`; canonical = `min = AGTC`.

**Closed-syncmer selection** (k=4 → `s = max(1, min(11, (k+2)/3)) = 2`). Keep when the smallest internal
2-mer is at offset 0 or `k−s=2`; drop if it's in the middle. On `GACACG`:

```text
k-mer   2-mers (off 0,1,2)   min at        keep?
GACA    GA AC CA             AC @ 1         DROP (middle)
ACAC    AC CA AC             AC @ 0         KEEP (end)
CACG    CA AC CG             AC @ 1         DROP (middle)
```

So `GACACG` contributes one syncmer: `ACAC`. `--feature-mode all` instead keeps every canonical k-mer
(`GACA, ACAC, CACG`). The default `k=31` fits assembled haplotypes (no error-tolerance reason to shorten;
each marker is tested locally so genome-wide uniqueness is irrelevant).

## Counting + the discriminative filter

Each retained feature is counted **per haplotype** (counts = multiplicity, so copy-number expansions stay
faithful), then kept by a two-part rule (`--min-paths N`, default 1):

1. **Copy-number features always kept** — if the count *varies among carriers* (`min_nonzero ≠ max`) it
   tracks copy number, informative even if present everywhere.
2. **Otherwise a symmetric minor-presence (MAF) cut** — drop when `min(present, absent) ≤ N` (removes
   no-contrast features: all-equal, singletons, all-but-one).

Worked example over 5 haplotypes (`--min-paths 1`):

```text
feature   counts             present/absent   min_nz≠max?   verdict
K1=ACAC   1,1,3,1,0          4/1              3≠1 yes       KEEP (copy-number: 3× in hap3)
K2        1,1,1,1,1          5/0              1=1 no        DROP (no contrast)
K3        2,0,0,0,0          1/4              2=2 no        DROP (singleton)
K4        1,1,1,1,0          4/1              1=1 no        DROP (all-but-one)
```

`--min-paths 0` disables rule 2 (keeps every non-constant feature). The same filter applies to node/edge
dosage features.

## Node provenance

Because k-mers are spelled from the bubble's node sequences, each kept k-mer records the set of graph nodes
its occurrences touch (in the feature map). So a significant marker maps back to a node/edge in the bubble —
and to the variant `call` makes there. The pooled cohort survivors are written as BIMBAM mean-genotype
dosage (`bimbam_{kmers,graph}.bimbam.gz`).

## Node + edge dosage (complementary substrate)

Built from the same canonical walks as the k-mers:

- **node dosage** — a *descriptive* traversal count, **not** a CN call: count > 1 can be a tandem (`…A,A…`)
  *or* a node revisited elsewhere (`…A,B,A…`); for short nodes it carries no CN meaning. Treat as
  presence/abundance.
- **edge dosage** — count of each oriented `from±>to±` transition; adjacency-aware, so it is the better
  tandem signal (a tandem block repeats an edge; scattered reuse does not) and an inversion appears as
  distinct (orientation-flipped) edge features.

The discriminative filter drops constant features (e.g. source/sink nodes). Node dosages match `inspect`'s
`node_counts.tsv` totals by construction. Under `--variant-nodes` both substrates are confined to the
variant nodes (+ `--variant-flank-bp`): a k-mer is masked unless it lies on a variant node, a node feature
is counted only if it is a variant node, an edge only if **both** endpoints are.
