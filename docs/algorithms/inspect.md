# Module `inspect` - algorithm

Mechanism for the `inspect` utility's `--cluster` and `--genotype-table` modes. For usage/flags see [modules/inspect.md](../modules/inspect.md); references in [references.md](../references.md#inspect).

`--cluster` groups haplotypes crossing a bubble by how similarly they traverse it and selects one representative per connected group. Identical walks are collapsed before comparison; additional walks can cluster when their similarity reaches the chosen threshold. A walk is a sequence of oriented node steps, compared through shared windowed content.

## How it works

### 1. Collapse identical walks

Haplotypes spelling the same step sequence fold to one distinct walk before anything is compared, so the work scales with the number of distinct alleles rather than with the panel.

### 2. Shingle each walk

Each distinct walk is broken into shingles, windows of consecutive steps. The multiset of its shingles is the walk's fingerprint. Multiplicity is kept rather than discarded: the second occurrence of a shingle is a different element from the first, which is what lets two haplotypes carrying the same repeat unit at different copy numbers be told apart. As plain sets they would be identical.

### 3. Score similarity

Similarity is derived from the shingle Jaccard as `2J / (1 + J)`.

Where a walk is short enough that its sketch holds every shingle it has, `J` is computed directly as intersection over union, since the sketch is the multiset and there is nothing to estimate. Where a walk is long enough to be truncated, `J` is estimated from the smallest values of the two sketches' union. Taking the union first matters: each sketch holds the smallest hashes of its own shingles, so a shingle both walks carry can sit inside one sketch and below the other's cutoff, and comparing the stored sketches directly is biased by exactly that, increasingly so as the two walks differ in length.

A walk too short to shingle at all falls back to a length-weighted Jaccard over its oriented node tokens.

### 4. Cluster

Every pair at or above `--cluster-similarity` is connected, and the clusters are the connected components of that graph. Membership is therefore transitive: two walks below the threshold can end up together through a chain of intermediates.

### 5. Pick a representative

Each cluster is represented by the walk minimizing its maximum, then mean, distance to the others, with the better-supported walk breaking ties. Among the haplotypes realizing that walk, the lexicographically smallest name is reported, so the output does not depend on the order of records in the input file.

## Worked trace

The steps below follow the five above, one for one. Four haplotypes cross one bubble; `a` and `b` are the boundaries and `(u, v)` is the repeat unit:

```text
h1 : a u v u v b              unit (u,v) twice
h2 : a u v u v b              identical to h1
h3 : a u v u v u v b          unit (u,v) three times
h4 : a x y b                  a different allele
```

1. Collapse identical walks. `h1` and `h2` spell the same sequence and fold into one distinct walk carrying two haplotypes. `h3` and `h4` are distinct, so three walks go forward.

2. Shingle each walk, here in windows of three steps:

```text
h1 : (a,u,v) (u,v,u) (v,u,v) (u,v,b)
h3 : (a,u,v) (u,v,u) (v,u,v) (u,v,u) (v,u,v) (u,v,b)
h4 : (a,x,y) (x,y,b)
```

   As plain sets `h1` and `h3` hold exactly the same shingles, so a set comparison would call them identical and merge the two-copy and three-copy alleles. Keeping counts separates them: the two middle shingles occur once in `h1` and twice in `h3`.

3. Score similarity. Both walks are short enough to be held completely, so the Jaccard is computed rather than estimated: the intersection over minimum counts is 4 and the union over maximum counts is 6, giving `J = 0.667` and an identity of `0.80`. `h4` shares no shingle with either, so its identity to both is 0.

4. Cluster. At `--cluster-similarity 0.90` nothing links: `0.80` is below the threshold, so the two-copy and three-copy alleles stay apart rather than chaining into one band, and `h4` stays alone. At a threshold of `0.79` the first two would join.

5. Pick a representative. Each cluster here has one walk, so each represents itself; `h1`'s cluster reports the lexicographically smaller of the two haplotype names carrying it.

Resulting `clusters.tsv`:

```text
cluster_id  n_paths  representative_path  members
0           2        h1                   h1;h2
1           1        h3                   h3
2           1        h4                   h4
```

## Genotype table

`--genotype-table` turns each bubble into a small genotype matrix over its distinct alleles.

1. Spell each crossing path. The walk is the same `source → sink` walk the per-path FASTA spells, so a haplotype crossing on the opposite strand is reverse-complemented into the bubble's orientation and a pure deletion is the source and sink alone.
2. Collapse identical sequences. Haplotypes spelling the same bases share one allele, whatever walk produced them.
3. Number the alleles by descending haplotype count, ties by the first path carrying them in GFA order, so `allele1` is the most common and the numbering never depends on hashing.
4. Orient. With `--reference-path`, the reference's step indices through the bubble are compared: running backwards means the reference crosses `sink → source`, and every allele is reverse-complemented once, so all alleles share the reference strand. Dedup is unaffected, since a reverse complement is a bijection.
5. Fill rows. A haplotype row marks the one allele it spells. A sample row, from a cosigt table, adds its haplotypes' indicators, so the row sums to the sample's ploidy. A haplotype that does not cross the bubble has no allele to mark, and its row, or any sample row using it, is `NA` throughout: a partial sum would read as a real low count.

Trace, on the fixture in `tests/inspect_stats.sh`: `full` and `full2` walk `1+ 2+ 3+`, `rev` walks `3- 2- 1-`, `del` walks `1+ 3+`, and `other` never reaches the bubble.

```text
alleles.fa.gz                    genotypes.tsv.gz (haplotypes)    with --samples
>allele1 full,full2,rev          haplotype allele1 allele2        sample  allele1 allele2
AAAAAAAAAACCCCCCCCCCGGGGGGGGGT   full      1       0              het     1       1      (full, del)
>allele2 del                     full2     1       0              same    2       0      (full, full)
AAAAAAAAAAGGGGGGGGGT             del       0       1              shared  2       0      (full, rev)
                                 rev       1       0              miss    NA      NA     (full, other)
                                 other     NA      NA
```

With `--reference-path rev` the same alleles are written reverse-complemented (`allele1` becomes `ACCCCCCCCCGGGGGGGGGGTTTTTTTTTT`), because `rev` crosses `sink → source`.
