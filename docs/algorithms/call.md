# Module `call` - algorithm

Mechanism for the `call` module. For usage/flags see [modules/call.md](../modules/call.md); references in [references.md](../references.md#call).

`call` types structural variants by comparing walks rather than sequences. Within a bubble, each haplotype's route from one boundary to the other is a sequence of nodes, and so is the reference's; the difference between the two routes is the variant. Working in node space means an event is defined by which nodes a haplotype does or does not visit, which is exact and cheap, and it lets a tandem array be typed by how many times a node is traversed instead of as a pile of insertions.

## How it works

### 1. Group alleles

Haplotypes taking identical routes across the bubble carry the same allele. They are grouped and the allele is diffed once, with genotypes expanded back over the group afterwards, so cost scales with distinct alleles rather than panel size.

### 2. Fold copy-number loops

Consecutive repeats of one node — the self-loop a folded tandem array leaves behind — collapse to a single token before diffing. Extra copies then surface as copy number rather than as a spurious insertion of the repeated sequence.

### 3. Diff at anchors into typed events

Both routes are split at the nodes they share, and only the gap between two consecutive shared nodes is compared. Each gap becomes one typed event: a deletion where only the reference has nodes there, an insertion where only the haplotype does, an inversion where the haplotype's run is the reference's run reverse-complemented, or a substitution where both are present at once, emitted as a co-located deletion and insertion sharing one `EVENTID`.

### 4. Coalesce within a haplotype

One real event can be split across gaps by the graph's node boundaries. Consecutive same-type events within `--merge-distance-bp` are therefore joined into one, measured either in reference coordinates or in the haplotype's own sequence, since a large intervening deletion can push two parts of one insertion far apart in reference space while leaving them adjacent in the haplotype.

### 5. Merge across haplotypes

The same variant found in many haplotypes should be one record, not many. Events are clustered by single-linkage: two link when they are the same type, lie within a position window, and either share enough nodes or align well enough in sequence.

The window is `--merge-distance-bp` widened by the smaller event's size, so one large indel whose breakpoint wanders between haplotypes still reads as one site. The node test is a length-weighted Jaccard against `--merge-jaccard`; where that fails, sequence identity against `--merge-seq-identity` catches an allele threaded through different nodes, length-gated by `--merge-size-ratio` so two events of very different sizes do not merge on a shared substring. The largest member represents the record, and the merge evidence is recorded in `NMERGED`, `SVLEN_RANGE`, `MERGE_JACCARD`, `MERGE_SEQID` and `MERGE_DIAMETER`.

Single-linkage is transitive, so a chain of pairwise-similar events can join two members that are not similar to each other. `MERGE_DIAMETER` reports how far apart the most distant members are, which is what makes that visible.

Copy-number records do not take this path; they merge on a shared repeat-unit node.

### 6. Force-call and filter

Once a record's representative is fixed, every haplotype that was not already a carrier is re-tested against it and added if it matches. This recovers carriers whose own event fell below the size threshold, which would otherwise be genotyped as reference. A record survives only if its representative reaches `--min-sv-bp` and its carriers reach `--min-haplotypes` and `--min-alt-af`.

### 7. Place the record

`POS` is the base before the event rather than its first affected base, which is the symbolic convention: `REF` then carries one real reference base and `ALT` is the symbol. The event occupies `POS+1` onward, and `END` is the last reference base it spans.

| record | `POS` | `END` | `END − POS` |
|---|---|---|---|
| `DEL` / `INV` | base before the event | last affected base | `SVLEN` in absolute value |
| `INS` | base before the insertion | `POS` | 0, since an insertion spans no reference |
| `DUP`, `CN_SCOPE=REPEAT_UNIT` | last base of the upstream flank | last base of the first copy | `RU_LEN` |
| `DUP`, `CN_SCOPE=COLLAPSED_MODULE` | last base of the near boundary | base before the far boundary | `CN_MODULE_REF_BP` |

The two duplication rows differ because the routes count different things. A repeat unit is one node, so its span is that node's length; a collapsed module's span is the bubble interior, which is what `FORMAT:CNBP` sums over, and both boundaries are excluded from that sum.

The identity in the last row holds while every boundary is visited once, and `CN_SPAN_AMBIGUOUS` marks the records where it does not. `POS` and `END` come from the widest span, first source occurrence to last sink; `CN_MODULE_REF_BP` and `CNBP` sum node length by traversal count over the module's nodes. Where a boundary recurs, the span encloses the boundary visits lying between the outermost ones and the node sum does not, so the two differ by exactly that sequence. Treat the identity as a placement check on unflagged records, and a flagged one as a deliberate span choice rather than a uniquely resolved interval.

A deletion or inversion starting at the region's first base has no preceding base to anchor on and keeps its original anchor.

## Worked trace

The steps below follow the seven above, one for one. Reference interior `R = A B C D E F` in single-letter node tokens. Six haplotypes cross the bubble:

```text
ref  : A B C D E F
hX   : A B   D E F          C deleted
hW   : A B   D E F          C deleted, same route as hX
hZ   : A     D E F          B and C deleted
hI   : A B P Q C D E F      P, Q inserted after B
hS   : A B C X E F          D replaced by X
hV   : A B d c E F          C D inverted; lowercase is the reverse-complement run
```

1. Group alleles. `hX` and `hW` take the same route, so they are one allele and are diffed once. The other four are distinct, giving five alleles for six haplotypes.

2. Fold copy-number loops. No node repeats consecutively here, so nothing folds and every route is diffed as written.

3. Diff at anchors into typed events. The tokens appearing once in both routes are the anchors, and only the gap between consecutive anchors is compared. `hX` and `hW` give a deletion of `C`; `hZ` gives one deletion of `B C`, since both fall in a single gap; `hI` gives an insertion of `P Q`; `hS`'s gap holds reference-only `D` and haplotype-only `X` at once, so it is a substitution, emitted as a paired deletion and insertion sharing an `EVENTID`; `hV`'s gap is the reference run `C D` traversed reverse-complemented, so it is an inversion.

4. Coalesce within a haplotype. Each haplotype produced its events in a single gap, so there is nothing adjacent to join.

5. Merge across haplotypes. `hX` and `hW` have identical node sets and become one record with two carriers. `hZ` shares only `C` with them, which is below `--merge-jaccard`, so it stays separate. The substitution pair and the inversion match nothing else and remain singletons.

6. Force-call and filter. Every non-carrier is re-tested against each surviving representative; none matches here. Records below `--min-sv-bp` or `--min-haplotypes` would be dropped at this point.

7. Place the record. Each record is anchored on the base before its first affected base, and `END` set from the table above.

```text
SVTYPE=DEL  EVENT_NODES=C     AC=2   carriers hX, hW
SVTYPE=DEL  EVENT_NODES=B,C   AC=1   carrier  hZ
SVTYPE=INS  EVENT_NODES=P,Q   AC=1   carrier  hI
SVTYPE=DEL  EVENT_NODES=D     AC=1   carrier  hS   EVENTID=e1
SVTYPE=INS  EVENT_NODES=X     AC=1   carrier  hS   EVENTID=e1
SVTYPE=INV  EVENT_NODES=C,D   AC=1   carrier  hV
```

## Copy number

Copy number is read off the walk rather than inferred from coverage, by one of three routes. `CN_METHOD` records which was used, and `CN_SCOPE` says whether a copy means a copy of a repeat unit or of a whole collapsed module.

`REP` applies where the array has been folded to a self-loop. The copy number is the number of traversals, which is exact, and `RU_LEN` is a real per-copy length, so `(CN − REF_CN) × RU_LEN` gives the haplotype's size change.

`MODULE_BP` applies where several paralogous copies were collapsed onto shared nodes rather than folded into a unit. There is no single repeat unit to count, so the bases a haplotype carries across the module are divided by a unit calibrated from the reference. `CN_UNIT_BP` is that calibration constant and describes only the content the copies share, so it understates a carrier's real gain or loss; `FORMAT:CNBP` is the per-haplotype size. `--cn-unit-spacing` takes the step from the spacing of the panel's own values instead of from the reference, and `--max-cn-model-residual` refuses the call when the values do not fall near a consistent ladder.

`PEAK` applies where neither of the above does: copy number is taken from the highest multiplicity any interior node reaches. One short node visited several times can set it, so absolute copy number on this route is a heuristic and the records say so with `CN_CONFIDENCE=HEURISTIC`.

Without `--cn`, only the self-loop route runs, since it needs no inference. A folded extra copy then surfaces as an ordinary insertion, flagged `INS_SUBTYPE=DUP` under `--classify-ins`, and a lost copy as a deletion, neither carrying a per-sample copy number.

## Gene annotation

With `--gtf` and a PanSN reference path, gene intervals are projected onto reference nodes, so each record can name the genes it overlaps in `INFO/GENES` and `<prefix>.node_genes.tsv` maps every node to its genes.

For a duplication, the module's total copy number can sometimes be split per gene. Where the paralogs are divergent enough, each is given the k-mers private to it within the locus and a haplotype's dosage for that gene is read from how many of them it carries. Where they are too similar for any k-mer to be private, no split is attempted and the module's total stands. The result is `<prefix>.dup_gene_cn.tsv`, whose `reliable` column marks which rows rest on enough private markers to be trusted.
