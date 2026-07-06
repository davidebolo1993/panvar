# Module `describe` - algorithm

Mechanism and a hand-traced example for the **describe** module. For usage/flags see [modules/describe.md](../modules/describe.md); References in [references.md](../references.md#describe).

## Terms

- **k-mer** — a length-`k` DNA substring; **canonicalised** to `min(forward, reverse-complement)` of its 2-bit code (`A=0,C=1,G=2,T=3`), so a k-mer and its reverse complement count as one feature.
- **closed syncmer** — a sampling rule that keeps a k-mer only when its smallest internal `s`-mer sits at either end; it samples sequence space evenly and degrades gracefully under substitutions (one base change disturbs only the overlapping syncmers).
- **node/edge dosage** — how many times a path traverses a graph node / an oriented step-to-step edge (`from±>to±`); a graph-local, coordinate-keyed feature complementary to k-mers.
- **variant dosage** — one genotype per SV call read straight from `call`'s VCF: copy number for a `DUP`, allele presence otherwise. The variant, not the k-mer, is the multiple-testing unit.

## Graph features

### K-mers 

#### Enconding and selection

Each k-mer is packed into a 2-bit integer (`A=0, C=1, G=2, T=3`) and canonicalised — replaced by whichever of the k-mer and its reverse complement has the smaller code — so a sequence and its reverse strand count as one feature. At `k=4`, for instance, `GACT` encodes to 135 and its reverse complement `AGTC` to 45, so the canonical form is `AGTC`. The BIMBAM matrix names features compactly (`K1, K2, …`), and the feature map keeps the resolution back to the DNA string and its node provenance.

By default `describe` does not keep every k-mer but only the closed syncmers: a k-mer is retained when its smallest internal `s`-mer sits at one of the two ends, and dropped when that minimiser falls in the middle (here `s = max(1, min(11, (k+2)/3))`, so `s = 2` at `k=4`). This samples the sequence evenly and is robust to substitutions, since a single base change only disturbs the few syncmers that overlap it. Scanning the three 4-mers of `GACACG`:

```text
k-mer   2-mers (off 0,1,2)   min at        keep?
GACA    GA AC CA             AC @ 1         DROP (middle)
ACAC    AC CA AC             AC @ 0         KEEP (end)
CACG    CA AC CG             AC @ 1         DROP (middle)
```

so `GACACG` contributes the single syncmer `ACAC`; passing `--feature-mode all` instead keeps every canonical k-mer (`GACA, ACAC, CACG`). The default `k=31` suits assembled haplotypes — there is no read-error reason to shorten it, and because each marker is tested locally, genome-wide k-mer uniqueness does not matter.

#### Counting and filtering

Each retained feature is then counted **per haplotype**, and because the count is a true multiplicity rather than a presence flag, a copy-number expansion is recorded faithfully (three copies of a k-mer show as 3, not 1). A feature is kept by a two-part rule (`--min-paths N`, default 1):

1. A feature whose count varies among its carriers (`min_nonzero ≠ max`) is always kept: that variation is a copy-number signal, and it is informative even when the feature is present in every haplotype.
2. Otherwise the feature must clear a symmetric minor-presence (MAF) cut — it is dropped when `min(present, absent) ≤ N`, which removes the features with no contrast to test: those present everywhere, singletons, and all-but-one.

Worked example over 5 haplotypes (`--min-paths 1`):

```text
feature   counts             present/absent   min_nz≠max?   verdict
K1=ACAC   1,1,3,1,0          4/1              3≠1 yes       KEEP (copy-number: 3× in hap3)
K2        1,1,1,1,1          5/0              1=1 no        DROP (no contrast)
K3        2,0,0,0,0          1/4              2=2 no        DROP (singleton)
K4        1,1,1,1,0          4/1              1=1 no        DROP (all-but-one)
```

`--min-paths 0` disables rule 2 (keeps every non-constant feature). The same filter applies to node/edge dosage features.

#### Provenance

Because k-mers are spelled from the bubble's node sequences, each kept k-mer records the set of graph nodes its occurrences touch (in the feature map). So a significant marker maps back to a node/edge in the bubble — and to the variant `call` makes there. The pooled cohort survivors are written as BIMBAM mean-genotype dosage (`bimbam_{kmers,graph}.bimbam.gz`).

## Node and edge dosage

Built from the same canonical walks as the k-mers:

- **node dosage** — a descriptive traversal count: count > 1 can be a tandem (`…A,A…`) or a node revisited elsewhere (`…A,B,A…`); for short nodes it carries no CN meaning. Treat as presence/abundance.
- **edge dosage** — count of each oriented `from±>to±` transition; adjacency-aware, so it is the better tandem signal (a tandem block repeats an edge; scattered reuse does not) and an inversion appears as distinct (orientation-flipped) edge features.

The discriminative filter drops constant features (e.g. source/sink nodes). Node dosages match `inspect`'s `node_counts.tsv` totals by construction. Under `--variant-nodes` both substrates are confined to the variant nodes (+ `--variant-flank-bp`): a k-mer is masked unless it lies on a variant node, a node feature is counted only if it is a variant node, and an edge only if **both** endpoints are.

## Variant features

The k-mer and graph substrates are spelled from the graph itself. The third substrate (`--variant-vcf`) is different: it reads `call`'s region VCF back and emits one genotype row per SV call, so the feature that gets tested is the variant `call` already typed rather than a fragment of it. Its samples (VCF columns) are the haplotype paths, and each record contributes one BIMBAM row whose dosage is read off `FORMAT`:

- a `DUP` carries a per-sample `CN`, so its dosage is that copy number — the honest signal for a copy-number locus;
- a multiallelic record (`NALLELES > 1`) expands into one row per ALT allele (`<id>_a1`, `<id>_a2`, …), each an indicator of whether the haplotype carries that allele;
- a `DEL`/`INS`/`INV` uses the `0/1` genotype directly, a presence dosage.

A haplotype whose genotype is `.` did not traverse the bubble and becomes `NA`, exactly as in the other substrates. Reading one VCF record over four haplotypes (`h4` misses the two graph-collapsed bubbles):

```text
VCF record (FORMAT)          h1      h2      h3      h4      →  BIMBAM row(s)
DUP   V1   GT:CN           1:3     0:2     1:4     .:.       →  V1      3,  2,  4,  NA
DEL   V2   GT              1       0       1       0         →  V2      1,  0,  1,  0
multi V3   GT (2 alleles)  1       2       0       .         →  V3_a1   1,  0,  0,  NA
                                                                V3_a2   0,  1,  0,  NA
```

The `DUP` keeps a reference-like sample's actual copy number (`h2` shows `2`, not `0`), so a loss and a gain are both testable dosages rather than a folded presence flag.

Unlike the k-mer/graph substrates, the variant layer runs no discriminative filter — it emits every call and leaves frequency pruning to `associate --min-maf`. The point of the layer is the multiple-testing denominator: the many k-mers and nodes inside one event are correlated, so testing them all inflates the count of independent tests; one dosage per variant is the statistically honest GWAS unit. The sidecar `feature_annot.variant.tsv.gz` carries each row's `svtype`, `gene`, `AF`, and `AN` for the traceback, and `--samples` sums a sample's haplotype dosages into the diploid `bimbam_variant.samples.bimbam.gz` (CN = CN_A + CN_B). `--scale-dosage` applies here too, min-max mapping each variant's dosage to `0..2`.
