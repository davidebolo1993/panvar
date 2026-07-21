# Module `rebuild` - algorithm

Mechanism for the `rebuild` module. For usage/flags see [modules/rebuild.md](../modules/rebuild.md); References in [references.md](../references.md#rebuild).

## Terms

- **degree** — for a node, the number of distinct link neighbours counting both ends.
- **hub** — a node whose degree is at or above `--hub-degree`. A cluster of hubs is the signature of a locus resolved by transitive closure rather than by colinear alignment.
- **richness** — how much distinct sequence a haplotype carries, used to order haplotypes for progressive construction.
- **seed/backbone** — the first haplotype added. It becomes the coordinate backbone every later haplotype is aligned against, so the choice is not neutral.
- **augmentation** — inserting the parts of a haplotype that did not align as new nodes and edges.


A graph builder that computes the transitive closure of all pairwise alignments will merge every near-identical copy of a repeat unit into one node. In a low-complexity locus that produces a small set of nodes with very high degree, and a decomposition with no clean nesting: snarl finding returns a single large site rather than a set of variant sites, so nothing downstream can call it.

Progressive construction does not have this failure mode. Each haplotype is chained colinearly through the graph built so far, so a repeat copy either follows an existing path or branches locally; there is no global closure step to collapse distinct copies onto one node. The locus is therefore reconstructible from the same haplotype sequences into a graph that does decompose. The cost is that progressive generation works above a size threshold, so variation below it is not represented.

## How it works

### 1. Gate 

One pass over the input GFA. Degree is computed per node as the number of distinct link neighbours, and the graph is called pathological when

```text
#nodes with degree >= --hub-degree   >=   --min-hubs
```

Two further signals are computed and logged but do not decide: the maximum degree, and node density (nodes per kb of the longest haplotype span).
A graph that fails the gate is copied to the output unchanged. This is what makes the module safe to put in front of a pipeline unconditionally: loci that were already decomposable are not touched, so their downstream output is unchanged.

### 2. Order

Progressive construction is order-sensitive. The seed becomes the coordinate backbone, and early haplotypes shape what later ones can align to, so haplotypes are added most-complete-first.

Completeness is measured by k-mer richness. For each haplotype the spelled sequence is scanned with a rolling 2-bit encoding, giving the number of distinct k-mers and the number of k-mer positions. Haplotypes are then ordered lexicographically, descending:

```text
(distinct k-mers, total k-mers)
```

Diversity decides first: a haplotype carrying many copies of one unit is not rewarded for the copies, only for the distinct sequence it contributes. Abundance discriminates only between haplotypes that carry the same amount of distinct sequence, where the one with more k-mer mass — more repeat copies of what it already has — is preferred. That secondary key is not a formality: on real cohorts roughly half to four-fifths of haplotypes share a distinct-k-mer count with at least one other, so it decides a large share of the order.

The two keys are deliberately not blended into a weighted sum. A sum lets abundance outrank diversity, so a haplotype with less distinct sequence but more repeat copies could seed the graph, which is the opposite of what the ordering is for.

Richness and length are correlated but not equivalent, which is the reason for measuring richness rather than sorting by length: richness discounts a haplotype that is long because a repeat expanded in it.

### 3. Generate 

Haplotypes are written out in richness order and handed to the generator, which reads the first as the seed graph and then, for each subsequent haplotype in turn, aligns it against the graph built so far and augments the graph with the segments that did not align. Only differences of at least `--min-var` are inserted; smaller ones are left at the backbone allele. This is the step that is delegated rather than reimplemented, and passing the order is how this module steers it.

### 4. Recover walks

The generated graph carries no path information, so each haplotype is mapped back to it and its walk is read out directly, in memory. Where a haplotype produces more than one chain the longest is taken, rather than concatenating chains, because stitching two chains together can produce a walk that is not a coherent traversal.

Each step of a walk is an oriented vertex, so a haplotype that traverses a segment on the reverse strand is recorded as such. Preserving that orientation through to the `L` and `P` lines is what allows an inversion to survive into the emitted graph as an inversion bubble; flattening every link to the forward strand would silently discard that entire class of variant.

A haplotype that fails to map produces no `P` line, so the number of recovered walks is reported and any shortfall warned about. Nothing here treats the reference as special: it is spelled, ordered and mapped like every other haplotype, and it seeds the graph only if it happens to be the richest. The reference is supplied later, to `bubble`, which sorts the graph along it.

### 5. Emit

The output is a plain GFA: `S` lines, `L` lines with both orientations, and one `P` line per recovered haplotype. Segments are renumbered to consecutive integers starting at 1, with the generator's original name kept as an `SN:Z:` tag, because tooling that stores node ids as integers will refuse a non-numeric segment name outright.

## Worked trace

Four haplotypes over one locus. `A`, `B`, `C` are shared segments; `U` is a repeat unit; `V` is a segment unique to one haplotype.

```text
h1 :  A  U U U U  B      C          (longest, but three of its four U are repeat copies)
h2 :  A  U U      B  V   C
h3 :  A  U U U    B      C
h4 :  A  U        B      C
```

1. Gate. Suppose the input graph has 14 nodes of degree ≥ 50 against a threshold of 10: pathological, so the rebuild proceeds. Had it been under the threshold, the input would have been copied out and the run would end here.

2. Order. Distinct k-mer content is what separates the haplotypes, not length. `h1` is the longest but its extra length is repeat copies, which add no distinct k-mers; `h2` is shorter but carries `V`, which nothing else has. So `h2` has the most distinct sequence and seeds. `h1` and `h3` then tie on diversity: both span the same k-mers, including the ones crossing a `U`-`U` junction, and differ only in how many copies they carry. The tie falls to total k-mers, putting `h1` ahead. `h4` has a single `U`, so it never spans a `U`-`U` junction and has strictly fewer distinct k-mers than either. The order is `h2, h1, h3, h4`.

3. Generate. `h2` becomes the seed graph, `A U U B V C`. `h1` aligns along it and contributes two further copies of `U`, augmented as a branch. `h3` aligns entirely onto what exists — one path through three copies — and adds nothing. `h4` likewise. Only `h1` augmented, and only where it exceeded the threshold; a haplotype differing from the backbone by a single base would have added nothing at all.

4. Recover. Each haplotype is mapped back and its walk read out, giving four `P` lines over the shared segment set.

5. Emit. Segments are renumbered `1..n` and links written with orientation.

The result is that the four haplotypes share `A`, `B`, `C` and the copies of `U`, differing by copy count and by the presence of `V` — a graph in which the copy-number site and the insertion site are separate, nested, decomposable bubbles. Each of those is a bubble `bubble` can find and `call` can type, and the copy-number one is what `panphorte` folds into a repeat-unit self-loop.
