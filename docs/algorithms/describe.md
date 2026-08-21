# Module `describe` - algorithm

Mechanism for the `describe` module. For usage/flags see [modules/describe.md](../modules/describe.md); References in [references.md](../references.md#describe).

`describe` turns a called locus into association-ready genotype features across three substrates: k-mers and node/edge dosages spelled from the graph, and one dosage per structural-variant call read from `call`'s VCF. Copy-number variation remains a multiplicity in every substrate where it is available rather than being reduced to presence or absence.

## How it works

### 1. Take each haplotype's canonical walk

The graph-derived substrates start from each haplotype's walk across a bubble, canonicalized to run from one boundary to the other so opposite-strand crossings are comparable. A haplotype taking the direct route between the boundaries, carrying nothing between them, is a traverser like any other and is counted as one. The variant substrate later reads the corresponding genotypes from the VCF.

### 2. Spell k-mer features from the walk

Each k-mer is packed into a 2-bit integer (`A=0, C=1, G=2, T=3`) and canonicalised — replaced by whichever of the k-mer and its reverse complement has the smaller code — so a sequence and its reverse strand count as one feature. At `k=4`, for instance, each base is a 2-bit digit (`A=00, C=01, G=10, T=11`) concatenated left-to-right into one integer, first base most-significant: `GACT` = `10 00 01 11` in binary = `2·4³ + 0·4² + 1·4 + 3` = 135, and its reverse complement `AGTC` = `00 10 11 01` = `0·4³ + 2·4² + 3·4 + 1` = 45, so the canonical form is `AGTC` (the smaller code). Per-bubble files use compact names (`K1, K2, …`) and map them back to sequence and node provenance; pooled BIMBAM rows use the canonical DNA sequence itself so the same marker has one identifier across bubbles.

By default `describe` does not keep every k-mer but only the closed syncmers: a k-mer is retained when its smallest internal s-mer sits at one of the two ends, and dropped when that minimiser falls in the middle (here `s = max(1, min(11, (k+2)/3))`, so `s = 2` at `k=4`). This samples the sequence evenly and is robust to substitutions, since a single base change only disturbs the few syncmers that overlap it. Scanning the three 4-mers of `GACACG`:

```text
k-mer   2-mers (off 0,1,2)   min at        keep?
GACA    GA AC CA             AC @ 1         DROP (middle)
ACAC    AC CA AC             AC @ 0         KEEP (end)
CACG    CA AC CG             AC @ 1         DROP (middle)
```

so `GACACG` contributes the single syncmer `ACAC`; passing `--feature-mode all` instead keeps every canonical k-mer (`GACA, ACAC, CACG`). The default `k=31` suits assembled haplotypes.

### 3. Count and filter to the discriminative features

Each retained feature is then counted per haplotype, and because the count is a true multiplicity rather than a presence flag, a copy-number expansion is recorded faithfully (three copies of a k-mer show as 3, not 1). A feature is kept by a two-part rule (`--min-paths N`, default 1):

1. A feature whose count varies among its carriers is always kept: that variation is a copy-number signal, and it is informative even when the feature is present in every haplotype.
2. Otherwise the feature must clear a symmetric minor-presence cut — it is dropped when `min(n_present, n_absent) ≤ N`, where `N` is `--min-paths`. This removes features with no usable contrast: those present everywhere and, at the default `N=1`, singletons and all-but-one features.

Worked example over 5 haplotypes (`--min-paths 1`):

```text
feature   counts             present/absent   min_nz≠max?   verdict
K1=ACAC   1,1,3,1,0          4/1              3≠1 yes       KEEP (copy-number: 3× in hap3)
K2        1,1,1,1,1          5/0              1=1 no        DROP (no contrast)
K3        2,0,0,0,0          1/4              2=2 no        DROP (singleton)
K4        1,1,1,1,0          4/1              1=1 no        DROP (all-but-one)
```

`--min-paths 0` disables rule 2 (keeps every non-constant feature). The same filter applies to node/edge dosage features.

### 4. Record where each feature came from

Because k-mers are spelled from the bubble's node sequences, each kept k-mer records the set of graph nodes its occurrences touch (in the feature map). So a significant marker maps back to a node/edge in the bubble — and to the variant `call` makes there. The pooled cohort survivors are written as BIMBAM mean-genotype dosage (`bimbam_{kmers,graph}.bimbam.gz`).

### 5. Take node and edge dosage from the same walks

Built from the same canonical walks as the k-mers:

- node dosage: a descriptive traversal count: count > 1 can be a tandem (`…A,A…`) or a node revisited elsewhere (`…A,B,A…`); for short nodes it carries no `CN` meaning. Treat as presence/abundance.
- edge dosage: count of each oriented `from±>to±` transition; adjacency-aware, so it is the a better tandem signal compared to node dosage. An inversion appears as distinct (orientation-flipped) edge features.

The discriminative filter drops constant features (e.g. source/sink nodes). Node dosages match `inspect`'s `node_counts.tsv` totals by construction. Under `--variant-nodes` both substrates are confined to the variant nodes, widened by `--variant-flank-bp`: a k-mer is masked unless it lies on a variant node or within the flank, a node feature is counted only if it is a variant node or the flank reaches it, and an edge only if both endpoints are.

The flank has two granularities, deliberately. For k-mers it admits exactly N bases at the node end facing the variant node, the rest being masked to `N` so the scanner resets; for node and edge dosage it admits the whole node, because a node dosage is a property of the whole node and there is no partial-node count to give. So the same flag selects more nodes than bases. Keeping a neighbour entire for the k-mers as well is what made `--variant-flank-bp 30` admit a 100 kb node and every k-mer in it.

A traversal does not require an interior node. A path taking the direct source→sink edge — a pure deletion, or the short side of an insertion — is a traverser of the bubble and receives a dosage of `0` at the nodes it skips, not `NA`. Requiring an interior node meant such a path was not observed at all, and the node that discriminates the deletion was then discarded as non-discriminative because only its carriers had been seen: the site produced no features whatever.

### 6. Read the variant substrate from the call output

The k-mer and graph substrates are spelled from the graph itself. The third substrate (`--variant-vcf`) is different: it reads `call`'s region VCF back and emits one genotype row per SV call, so the feature that gets tested is the variant `call` already typed rather than a fragment of it. Its samples (VCF columns) are the haplotype paths, and each record contributes one BIMBAM row whose dosage is read off `FORMAT`:

- a `DUP` carries a per-sample `CN`, so its dosage is that copy number rather than a presence indicator;
- a multiallelic record (`NALLELES > 1`) expands into one row per ALT allele (`<id>_a1`, `<id>_a2`, …), each an indicator of whether the haplotype carries that allele;
- a `DEL`/`INS`/`INV` uses the `0/1` genotype directly, a presence dosage.

A haplotype whose genotype is `.` did not traverse the bubble and becomes `NA`, exactly as in the other substrates. `NA` means unobserved, and is kept distinct from `0` throughout: a pooled feature is finite only when every bubble contributing it is observable on that path, and a diploid sample only when every assigned haplotype traverses — a partial sum reported as a whole one would bias dosage downward. A `DUP` must carry a usable `CN`; falling back to `GT` presence would substitute a 0/1 indicator for a copy-number dosage. Reading one VCF record over four haplotypes (`h4` misses the two graph-collapsed bubbles):

```text
VCF record (FORMAT)          h1      h2      h3      h4      →  BIMBAM row(s)
DUP   V1   GT:CN           1:3     0:2     1:4     .:.       →  V1      3,  2,  4,  NA
DEL   V2   GT              1       0       1       0         →  V2      1,  0,  1,  0
multi V3   GT (2 alleles)  1       2       0       .         →  V3_a1   1,  0,  0,  NA
                                                                V3_a2   0,  1,  0,  NA
```

The `DUP` keeps a reference-like sample's actual copy number (`h2` shows `2`, not `0`), so a loss and a gain are both testable dosages rather than a folded presence flag.

Unlike the k-mer/graph substrates, the variant layer runs no discriminative filter — it emits every call and leaves frequency pruning to `associate`. It provides a coarser testing unit than the many correlated k-mers and nodes inside an event, although nearby or paired variant records can still be correlated. The sidecar `feature_annot.variant.tsv.gz` (in the `variant/` folder) carries each row's `svtype`, `gene`, `AF`, and `AN` for the traceback, and `--samples` sums a sample's haplotype dosages into the diploid `sample/variant/bimbam_variant.bimbam.gz`.

## Worked trace

The steps below follow the six above, one for one. One bubble, boundaries `a` and `b`, with an interior node `v` that three haplotypes carry and three delete outright:

```text
h1-h3: a v b        carry the interior; h1 is the reference
h4-h6: a   b        take the direct route, deleting it
```

1. Take each haplotype's canonical walk. All six cross the bubble, including the three that carry nothing between the boundaries. Counting only the interior-node carriers would leave the feature that distinguishes them observed on one side only.

2. Spell k-mer features from the walk. Each walk is spelled and scanned, so the carriers yield k-mers spanning `v` and its junctions with the flanks, while the deleting haplotypes yield k-mers spanning the junction the deletion creates, which the carriers do not have.

3. Count and filter to the discriminative features. A k-mer carried by all six is constant and is dropped; the ones spanning `v` and the ones spanning the deletion junction each split the panel three against three and are kept.

4. Record where each feature came from. Each kept k-mer records the nodes its occurrences touch and the bubble it belongs to, so a significant marker can be traced back to a site and, through the call output, to a variant.

5. Take node and edge dosage from the same walks. Node `v` reads 1 for the carriers and 0 for the others; the edge from `a` straight to `b` reads 1 for the deleting haplotypes and 0 for the carriers. A haplotype not reaching the bubble at all would read as missing rather than zero.

6. Read the variant substrate from the call output. With a VCF supplied, the deletion at this bubble also appears as one dosage row taken from each haplotype's genotype, which is the same split expressed at the level of the call rather than of the graph.
