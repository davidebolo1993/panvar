# Module `rebuild` - algorithm

Mechanism for the `rebuild` module. For usage/flags see [modules/rebuild.md](../modules/rebuild.md); References in [references.md](../references.md#rebuild).

pggb induces its graph with seqwish, whose transitive closure over all-pairs alignments collapses sequences that recur across the locus — repeat copies, or the same short segment reused elsewhere — onto shared nodes. Where that recurrence is dense the shared nodes reach very high degree and this complicates downstream variant calling, especially of duplications. Progressive graph construction does not have this failure mode. Each haplotype is threaded colinearly through the graph built so far, so recurring sequence either follows an existing path or branches locally. The locus is therefore reconstructible from the same haplotype sequences into a graph that does decompose. The cost is that progressive generation works above a size threshold, so variation below it is not represented.

## How it works

### 1. Gate

One pass over the input GFA decides whether the graph needs rebuilding. For each node it counts the node's degree. An `L` line joins one end of a segment to one end of another, so a node has links off both of its ends; its degree is the number of distinct other nodes those links reach, pooled across both ends (a neighbour reached by two links, or off both ends, still counts once). A node whose degree is at least `--hub-degree` is a hub — the mark of seqwish having merged sequence that recurs across the locus onto one shared node. The graph is called pathological when

```text
#nodes with degree >= --hub-degree   >=   --min-hubs
```

Two further signals are computed and logged (but do not decide the rebuilding): the maximum degree, and node density (nodes per kb of the longest haplotype span). A graph that fails the gate is copied to the output unchanged, so a locus that was already decomposable is never touched.

### 2. Order

Progressive construction is order-sensitive. The first haplotype added — the seed — becomes the coordinate backbone every later haplotype is aligned against, and early haplotypes shape what later ones can align to, so haplotypes are added most-complete-first. Completeness is measured by k-mer richness. For each haplotype the spelled sequence is scanned with a rolling 2-bit encoding, giving the number of distinct k-mers and the number of k-mer positions. Haplotypes are then ordered lexicographically, descending:

```text
(distinct k-mers, total k-mers)
```

Diversity decides first: a haplotype carrying many copies of one unit is not rewarded for the copies, only for the distinct sequence it contributes. Abundance discriminates only between haplotypes that carry the same amount of distinct sequence, where the one with more k-mer mass — more repeat copies of what it already has — is preferred. The two keys are deliberately not blended into a weighted sum. A sum lets abundance outrank diversity, so a haplotype with less distinct sequence but more repeat copies could seed the graph, which is the opposite of what the ordering is for.


### 3. Generate 

Haplotypes are written out in richness order and handed to the generator, which reads the first as the seed graph and then, for each subsequent haplotype in turn, aligns it against the graph built so far and augments the graph — inserting the segments that did not align as new nodes and edges. Only differences of at least `--min-var` are inserted; smaller ones are left at the backbone allele. This is the step that is delegated to the minigraph library.

### 4. Recover walks

The generated minigraph graph carries no path information, so each haplotype is mapped back to it and its walk is read out directly. Mapping returns one or more chains — each a run of graph vertices covering a contiguous stretch of the haplotype, in order. A haplotype that maps cleanly gives a single chain spanning it end to end; it splits into several only when a stretch of it — divergent, rearranged, or missing as an assembly gap — leaves no single graph path to thread it onto. Separate chains need not be joined by a link in the graph, so concatenating them into one `P` line could assert an edge that does not exist — a walk that cannot actually be traversed. The single longest chain (the one covering the most of the haplotype) is therefore taken as its walk. Each step of a walk is an oriented vertex, so a haplotype that traverses a segment on the reverse strand is recorded as such. Preserving that orientation through to the `L` and `P` lines is what allows an inversion to survive into the emitted graph as an inversion bubble.


### 5. Emit

The output is a plain GFA: `S` lines, `L` lines with both orientations, and one `P` line per recovered haplotype. Segments are renumbered to consecutive integers starting at 1, with the generator's original name kept as an `SN:Z:` tag, because tooling that stores node ids as integers will refuse a non-numeric segment name outright.

## Worked trace

Five haplotypes over one locus. `A`, `B`, `C` are shared segments; `U` is a repeat unit; `V` is a segment unique to one haplotype. `h5` comes from a fragmented assembly whose middle (where `B` sits) is an assembly gap, so its sequence jumps from the repeat array straight to `C`.

```text
h1 :  A  U U U U  B       C          (longest, but three of its four U are repeat copies)
h2 :  A  U U      B  V    C
h3 :  A  U U U    B       C
h4 :  A  U        B       C
h5 :  A  U U    [ gap ]   C          (assembly gap where B would be)
```

1. Gate. Suppose the input graph has 14 nodes of degree ≥ 50 against a threshold of 10: pathological, so the rebuild proceeds. Had it been under the threshold, the input would have been copied out and the run would end here.

2. Order. Distinct k-mer content is what separates the haplotypes, not length. `h1` is the longest but its extra length is repeat copies, which add no distinct k-mers; `h2` is shorter but carries `V`, which nothing else has. So `h2` has the most distinct sequence and seeds. `h1` and `h3` then tie on diversity: both span the same k-mers, including the ones crossing a `U`-`U` junction, and differ only in how many copies they carry. The tie falls to total k-mers, putting `h1` ahead. `h4` and `h5` carry the least and follow. The order is `h2, h1, h3, h4, h5`.

3. Generate. `h2` becomes the seed graph, `A U U B V C`. `h1` aligns along it and contributes two further copies of `U`, augmented as a branch. `h3`, `h4` and `h5` align onto what already exists and add nothing — `h5`'s gap has no sequence to insert, and its `A U U` and `C` both match nodes already there. Only `h1` augmented, and only where it exceeded the threshold; a haplotype differing from the backbone by a single base would have added nothing at all.

4. Recover. `h1`–`h4` each map back as a single chain and read out cleanly, giving four `P` lines over the shared segment set. `h5` does not: the graph reaches `C` only through `B` (paths `…U → B → C` or `…U → B → V → C`), but `h5` has no `B`, so there is no `U → C` edge to carry it across the gap. Mapping therefore returns two chains — one over the `A U U` prefix, one over the lone `C` — with nothing joining them. Concatenating into `A U U C` would assert a `U → C` link the graph does not contain, so the longer chain wins: `A U U` covers more of `h5` than `C` does, so `h5`'s `P` line is `A U U` and the `C` tail is dropped. Its recovered coverage falls below 1, which the run log reports. (This fragmented case is uncommon; a well-assembled haplotype maps as one chain.)

5. Emit. Segments are renumbered `1..n` and links written with orientation; five `P` lines are written, `h5`'s being partial.

The result is that the haplotypes share `A`, `B`, `C` and the copies of `U`, differing by copy count and by the presence of `V` — a graph in which the copy-number site and the insertion site are separate, nested, decomposable bubbles. Each of those is a bubble `bubble` can find and `call` can type, and the copy-number one is what `panphorte` folds into a repeat-unit self-loop.
