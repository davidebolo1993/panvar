# Module `rebuild` - algorithm

Mechanism for the `rebuild` module. For usage/flags see [modules/rebuild.md](../modules/rebuild.md); references in [references.md](../references.md#rebuild).

A graph induced by transitive closure over all-pairs alignments can collapse sequence that recurs across a locus — repeat copies, or the same short segment reused elsewhere — onto shared nodes. Where that recurrence is dense the shared nodes reach very high degree, which makes downstream variant calling hard, particularly for duplications. Progressive graph construction instead threads each haplotype colinearly through the graph built so far, so recurring sequence tends to follow an existing path or branch locally. This can produce a graph that decomposes more cleanly, at the cost of representing only variation above a size threshold.

## How it works

### 1. Gate

One pass over the input decides whether the graph needs rebuilding. A link joins one end of a segment to one end of another, so a node has links off both ends and degree is counted per end: a node's degree is the larger of its two, not the two pooled. Pooling conflates two different shapes, since twenty-five neighbours on each side is a clean two-sided branch while fifty on one side is the tangle the gate exists to find, and pooled they read the same. Self-loops are counted separately, because that is how a folded tandem array appears and it is not pathology.

A node whose degree is at least `--hub-degree` is a hub. The graph is called pathological when

```text
#nodes with degree >= --hub-degree   >=   --min-hubs
```

Maximum degree and node density are logged alongside but do not decide anything. A graph that passes the gate is copied to the output unchanged, so a locus that was already decomposable is never touched.

### 2. Order

Progressive construction is order-sensitive. The first haplotype added becomes the coordinate backbone every later one is aligned against, so haplotypes are added richness-first: the number of distinct k-mers decides first and the total number breaks ties.

Diversity deciding first is the point. A haplotype carrying many copies of one unit is not rewarded for the copies, only for the distinct sequence it contributes; abundance separates only haplotypes that carry the same distinct content. Combining the two into a single weighted score would let abundance outrank diversity and seed the graph with a haplotype carrying less distinct sequence.

### 3. Generate

Haplotypes are written out in that order and handed to the generator, which reads the first as the seed graph and then, for each one in turn, aligns it against the graph built so far and inserts the segments that did not align as new nodes and edges. Only differences of at least `--min-var` are inserted; smaller ones stay at the backbone allele.

Two length thresholds apply in sequence, and the first is the one that bites. An alignment whose block length falls under `--min-align-len` is discarded whole, before any of its differences are examined, so `--min-var` never sees them. The generator's own default for that threshold assumes whole-chromosome input, which a locus graph cannot reach — every alignment is then rejected and the result collapses to the seed sequence alone. Setting it to `0` scales the threshold from the haplotype lengths instead, capped so no haplotype is excluded by a bar it could never clear.

### 4. Recover walks

The generated graph carries no path information, so each haplotype is mapped back to it and its walk read out. Mapping returns one or more chains, each a run of graph vertices covering a contiguous stretch of the haplotype in order. A haplotype that maps cleanly gives one chain spanning it end to end; it splits only where a stretch — divergent, rearranged, or missing as an assembly gap — leaves no single graph path to thread it onto.

Separate chains need not be joined by a link, so concatenating them would assert an edge that does not exist and produce a walk nobody can traverse. The chain with the most matching bases is therefore selected; identity and mapping quality break ties deterministically, and the rest is dropped. Each step is an oriented vertex, so a haplotype traversing a segment on the reverse strand is recorded as doing so, which is what lets an inversion survive into the emitted graph as an inversion bubble.

### 5. Check and accept

A rebuilt graph is only useful if the haplotypes come back. Each recovered walk is re-spelled from the rebuilt graph and compared against the original haplotype, giving four separate numbers rather than one: how much of the haplotype lies between its first and last aligned base, how much is covered by bases that actually match, the identity within the aligned region, and the identity of the re-spelled walk. The first hides internal gaps and the last does not, so a threshold set on either alone would mean something different.

The rebuild is accepted only if every haplotype is recovered, every consecutive step is joined by a link that exists, and each walk clears `--min-matched-cover` and `--min-recovered-identity`. A named `--reference-path` must clear them too, since everything downstream is expressed relative to it. If any check fails the rebuilt graph is discarded and the input is written unchanged, with the reason reported. `--allow-loss` accepts anyway and records what was violated.

Acceptance is about fidelity only. Whether the rebuild actually reduced tangling is measured separately, by re-running the gate's own degree count on the result, and reported rather than enforced.

### 6. Emit

The output is a plain GFA: `S` lines, `L` lines with both orientations, and one `P` line per recovered haplotype. Segments are renumbered to consecutive integers, keeping the generator's original name as a tag, because tooling that stores node ids as integers refuses a non-numeric segment name. The audit is written after the graph, so it always describes a file that exists.

## Worked trace

The steps below follow the six above, one for one. Five haplotypes over one locus: `A`, `B`, `C` are shared segments, `U` is a repeat unit, `V` is unique to one haplotype. `h5` comes from a fragmented assembly whose middle is a gap, so its sequence jumps from the repeat array straight to `C`.

```text
h1 :  A  U U U U  B       C          longest, but three of its four U are repeat copies
h2 :  A  U U      B  V    C
h3 :  A  U U U    B       C
h4 :  A  U        B       C
h5 :  A  U U    [ gap ]   C          assembly gap where B would be
```

1. Gate. Suppose the input has 14 nodes of degree at least 50, against a threshold of 10: pathological, so the rebuild proceeds. Under the threshold, the input would be copied out and the run would end here.

2. Order. Distinct content separates the haplotypes, not length. `h1` is longest but its extra length is repeat copies, which add no distinct k-mers; `h2` is shorter but carries `V`, which nothing else has, so `h2` seeds. `h1` and `h3` then tie on diversity, spanning the same k-mers and differing only in copy count, and the tie falls to total k-mers, putting `h1` ahead. The order is `h2, h1, h3, h4, h5`.

3. Generate. `h2` becomes the seed graph `A U U B V C`. `h1` aligns along it and contributes two further copies of `U` as a branch. `h3`, `h4` and `h5` align onto what exists and add nothing: `h5`'s gap has no sequence to insert, and its `A U U` and `C` both match nodes already present. A haplotype differing from the backbone by a single base would also have added nothing.

4. Recover walks. `h1` to `h4` each map as a single chain and read out cleanly. `h5` does not: the graph reaches `C` only through `B`, but `h5` has no `B`, so no edge carries it across the gap. Mapping returns two chains, one over `A U U` and one over the lone `C`, with nothing joining them. Concatenating would assert a link the graph does not contain, so the longer chain wins and `h5`'s walk is `A U U`.

5. Check and accept. `h1` to `h4` re-spell their haplotypes within the bounds. `h5` does not: its walk covers only the prefix, so its matched cover falls well below `--min-matched-cover` and it is recorded as `low_cover`. One haplotype failing is enough, so the rebuild is refused and the input graph is written unchanged. Passing `--allow-loss` would accept it instead, with `h5` recorded as the haplotype that was not recovered.

6. Emit. On an accepted run, segments are renumbered from 1 and links written with orientation, one `P` line per haplotype, followed by the audit giving each haplotype's cover, identity and status.

Had every haplotype been recovered, the result would share `A`, `B`, `C` and the copies of `U`, differing by copy count and by the presence of `V` — a graph in which the copy-number site and the insertion site are separate, nested, decomposable bubbles.
