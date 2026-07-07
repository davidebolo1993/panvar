# Module `call` - algorithm

Mechanism for the `call` module. For usage/flags see [modules/call.md](../modules/call.md); References in [references.md](../references.md#call).

## Terms

- **anchor** — a node token that occurs exactly once in *both* the reference and the haplotype walk, so it must correspond to the same position in each. The anchors the two walks share are chained in the order they appear (a monotonic subsequence), pinning the walks together at those points. 
- **gap block** — a maximal stretch between two anchors where reference and haplotype differ; each becomes one typed event (`DEL` / `INS` / `INV` / `substitution`).
- **full walk** — the widest `source→sink` span with every repeated traversal counted (a node visited twice contributes its length twice); contrast the minimal-span walk that visits each distinct node once.
- **peak multiplicity** — the visit count of a haplotype's single busiest node in the bubble.

## Algorithm

For each bubble:

1. Group identical `source→sink` walks into distinct alleles; call once per allele (genotypes expand by membership);
2.  Fold a node's consecutive self-repeats (a `REP` self-loop) into one alignment token, so extra copies surface as copy number, not a spurious `INS`/ `DEL`;
3. Split both walks at their shared anchors and align only the segments between them. Each gap block becomes one typed event: a `DEL` where only the reference has nodes, an `INS` where only the haplotype does, an `INV` where the haplotype's run is the reverse-complement node-walk of the reference's, or a substitution — a co-located `DEL` and `INS` that share one `EVENTID`.
4. Merge consecutive same-type events whose gap is ≤ `--merge-distance-bp`, measured in either reference space or the haplotype's own sequence space.
5. Merge across haplotypes into one record per site, by transitive single-linkage (connected components), seeding from the events down to `--rescue-min-bp`. Two events link when they are the same type, sit within a position window, and either their node sets overlap or their sequences align. The window is `--merge-distance-bp` widened by the smaller event's size, so one large indel whose breakpoint wanders by kilobases across haplotypes still reads as one site; the node/sequence gate below keeps that fuzziness from over-merging genuinely distinct events. The overlap test is a length-weighted node-set Jaccard ≥ `--merge-jaccard`; if that misses, the fallback is sequence identity ≥ `--merge-seq-identity` (length-gated by `--merge-size-ratio`). The largest member represents the record and `MERGE_*` records the evidence. Copy-number (`DUP`) records do not use this path — they merge separately on a shared `REP` self-loop node.
6. Re-test every non-carrier against its own diff, adding it as a carrier when it matches (sub-threshold events get `GT=1` instead of `0`); one pass (representative is fixed). Keep records whose representative reaches `--min-sv-bp` and carrier count reaches `--min-haplotypes`.

#### Merge keys — Jaccard vs sequence identity

The two link tests in step 5 fail on opposite kinds of event, which is why both exist:

- same content, different nodes — node Jaccard low, sequence identity high: one allele threaded through different graph nodes. The sequence gate rescues it — its main reason to exist.
- same nodes, poorly-aligning sequence — Jaccard high, sequence identity low: a shared backbone with a large internal indel or low-complexity content. The Jaccard gate rescues it (length-weighted, order- and orientation-blind).
- different sizes — only the sequence path is size-gated (`--merge-size-ratio`); Jaccard has no size gate.

So lower `--merge-jaccard` to merge events that share a backbone, and lower `--merge-seq-identity` / `--merge-size-ratio` to merge similar content threaded through different nodes.

## Worked trace (non-`DUP` events)

Reference interior `R = A B C D E F` — single-letter node tokens (a `REP` self-loop, if present, is dropped from this walk and typed separately as a `DUP`, so extra copies are never mistyped as an insertion). Six haplotypes cross the bubble:

```text
ref  : A B C D E F
hX   : A B   D E F          C deleted
hW   : A B   D E F          C deleted (same walk as hX)
hZ   : A     D E F          B and C deleted
hI   : A B P Q C D E F      P, Q inserted after B
hS   : A B C X E F          D replaced by X   (substitution)
hV   : A B d c E F          C D inverted      (lowercase = reverse-complement run, order reversed)
```

1. Diff each haplotype against the reference. The tokens that appear once in both walks are the anchors (`A B … E F`, in order); only the gap between two consecutive anchors is aligned, and each gap becomes one typed event:
   - `hX`/`hW`: the gap between anchors `B` and `D` holds only the reference's `C`: `DEL{C}`.
   - `hZ`: the gap between `A` and `D` holds the reference's `B C`: one `DEL{B,C}` (both fall in a single gap block).
   - `hI`: the gap after `B` holds the haplotype-only `P Q`: `INS{P,Q}`.
   - `hS`: the gap between `C` and `E` has reference-only `D` and haplotype-only `X` at once: it is a substitution, emitted as a co-located `DEL{D}` and `INS{X}` sharing one `EVENTID`.
   - `hV`: the gap between `B` and `E` is the reference run `C D` traversed as its reverse-complement node-walk: `INV{C,D}`.
2. Coalesce within a haplotype. Same-type events that sit in different gap blocks but within `--merge-distance-bp` are joined into one, measured either in reference space or in the haplotype's own sequence space (the latter catches insertions that a large intervening deletion would otherwise push apart). None of these trigger it.
3. Cluster across haplotypes by transitive single-linkage. `hX:DEL{C}` and `hW:DEL{C}` have identical node sets (Jaccard 1.0), so they become one record with `AC = 2`; `hZ:DEL{B,C}` shares only `C` (Jaccard `1/2 < --merge-jaccard`), so it stays its own record; `hS`'s substitution pair and `hV`'s `INV{C,D}` match nothing else, so each stays a singleton. Events also link when their sequences match (not just their nodes), which rescues one allele threaded through different nodes — see [Merge keys](#merge-keys--jaccard-vs-sequence-identity).
4. Force-call, then filter. Every non-carrier is re-tested against its own diff and added if it matches — so a sub-threshold copy of an event gets `GT = 1` rather than `0` — and a record is kept only if its representative reaches `--min-sv-bp` and its carriers reach `--min-haplotypes`.

Resulting records:

```text
SVTYPE=DEL  EVENT_NODES=C     AC=2   carriers hX, hW
SVTYPE=DEL  EVENT_NODES=B,C   AC=1   carrier  hZ
SVTYPE=INS  EVENT_NODES=P,Q   AC=1   carrier  hI
SVTYPE=DEL  EVENT_NODES=D     AC=1   carrier  hS   EVENTID=e1
SVTYPE=INS  EVENT_NODES=X     AC=1   carrier  hS   EVENTID=e1   (paired with the DEL above)
SVTYPE=INV  EVENT_NODES=C,D   AC=1   carrier  hV
```

Without `--cn`: the `DEL` / `INS` / `INV` records above are emitted either way — copy number is the only thing that changes. The always-on self-loop route still calls a genuine tandem a `DUP` (a `REP` self-loop needs no flag), but the coverage and peak-multiplicity routes are off: a folded extra copy then surfaces as an ordinary `INS` (flagged `INS_SUBTYPE=DUP` under `--classify-ins`) and a lost copy as a `DEL`, both without a per-sample `CN`. Add `--cn` to turn those into `DUP` records carrying `CN` / `REF_CN`.

## Worked trace (`DUP` events)

`call --cn` reads copy number straight off the walk and emits it as a `DUP` record carrying a per-sample `CN`. How it reads the number depends on how the locus is folded in the graph. 

#### A tandem repeat (self-loop `REP`)

`panphorte` has folded the tandem array into one `REP` node with a self-loop, so a walk's copy number is simply how many times it loops that node (`L`/`R` are flanks):

```text
reference : L REP REP REP R          loops REP ×3
hapX      : L REP REP R              loops REP ×2
hapY      : L REP REP REP REP R      loops REP ×4
```

1. The unit is already a single `REP` self-loop node (≥ `--min-sv-bp`), so its integer loop count is the copy number, exactly. The reference loops it three times, so `REF_CN = 3`.
2. Count each haplotype's loops: `hapX` traverses `REP` twice (`CN = 2`, one copy lost) and `hapY` four times (`CN = 4`, one gained); a haplotype at three loops is reference-like.
3. Emit one `DUP` at the bubble with `REF_CN = 3` and a per-sample `CN` from the loop count. `GT = 1` marks only the carriers whose `CN` differs from `REF_CN`, so `AC`/`AF` stay meaningful.

```text
SVTYPE=DUP  REF_CN=3  RU_LEN=<unit bp>     hapX CN=2   hapY CN=4   (reference-like samples CN=3)
```

This is the only reading that fires without `--cn` too: the self-loop route is always on, so a `REP` tandem is a `DUP` regardless of the flag.

#### A collapsed paralog cluster (coverage)

Here `pggb` collapsed near-identical paralog copies onto shared nodes, and the reference itself traverses them more than once. No single node can show a loss, so we measure the total sequence over the folded nodes and calibrate one copy from the reference:

```text
folded nodes (reference visits each ≥2×):  reference total = 10,000 bp,  ref folds 2×
hapA : 15,000 bp over the folded nodes
hapB :  5,000 bp
hapC : 10,000 bp
```

1. Calibrate the one-copy unit from the reference: `one-copy bp = reference folded-node bp / ref_fold = 10,000 / 2 = 5,000`, and `REF_CN = ref_fold = 2`.
2. Divide each haplotype's folded-node bp by that unit: `hapA` gives 3 (a gain), `hapB` gives 1 (a loss), `hapC` gives 2 (reference-like). One measure recovers both losses and gains.
3. Emit one `DUP` with `REF_CN = 2` and per-sample `CN` — the absolute total of all the collapsed paralogs together. Measuring over the folded nodes only is what keeps an edit to unique content in the bubble (interstitial sequence, a single-visit paralog) from being misread as copy number.

```text
SVTYPE=DUP  REF_CN=2      hapA CN=3   hapB CN=1   hapC CN=2
```

#### A single-copy-reference duplication (peak multiplicity)

The reference here visits each inside node once — it does not fold — so there is no reference fold for coverage to calibrate against. This reading covers a haplotype that folds an extra copy back onto those nodes (a duplication where the reference is the single/short allele):

```text
reference : each inside node ×1            ref_peak = 1
hapD      : busiest inside node ×2         one extra copy folded back
hapE      : each inside node ×1            no extra copy
```

1. `REF_CN = ref_peak = 1`. There is no unit to divide by, so we compare peak multiplicities directly.
2. Take each haplotype's busiest inside node: `hapD` peaks at 2 (`> 1`) and is interpreted as a gain; `hapE` stays at 1. Using the global peak, not a per-node excess, keeps scattered cluster background from inflating the call.
3. Emit a `DUP` for `hapD` with `REF_CN = 1`, `CN = 2`, and `SVLEN = Σ node_len × excess` (the duplicated bp). This reading is gains only — a loss leaves no extra fold to count.

```text
SVTYPE=DUP  REF_CN=1      hapD CN=2   (hapE is not a carrier)
```


## Gene annotation trace (`--gtf`)

With `--gtf`, `call` records which genes each variant touches and splits a `CN` call that spans several genes into per-gene `CN` — from the GTF (Gene Transfer Format) alone, by k-mer dosage. The reference path's PanSN (Pangenome Sequence Naming) name (`sample#hap#chrom:start-end`) supplies the chromosome and start coordinate the projection needs.

First the genes are projected onto the reference nodes: walking the reference and accumulating node lengths gives each node a coordinate span, intersected with the GTF gene intervals to tag every node with the gene (or genes) it falls in. Where the reference revisits the same nodes at two coordinates — a folded paralog — that node is tagged with both genes, which is precisely why graph multiplicity alone cannot tell the paralogs apart. When a `DUP` overlaps two or more genes, each gene's discriminative sketch is its merged coding sequence (CDS; the k-mers private to it versus its paralogs), and per-haplotype copy number is the dosage of those private k-mers over the module sub-walk.


## Worked trace (assign genes in `DUP` events)

#### Divergent paralogs (private-k-mer dosage)

A folded cluster carries two coding paralogs `A` and `B` plus a pseudogene `P`. Coding sequence is where paralogs differ, so each gene's sketch is its merged CDS — except `P`, a gene with no CDS, which falls back to its gene span:

```text
gene   marker   private k-mers (unique vs the other paralogs)
A      CDS      1200
B      CDS       800
P      span     4800     (no CDS -> gene span)
```

1. Build each gene's canonical k-mer set from its marker and keep the k-mers private to it (absent from every paralog in the cluster). These sets discriminate from the reference alone.
2. Scan a haplotype's module sub-walk and count each gene's private k-mers. A gene present in `c` copies carries each private k-mer about `c` times, so `dosage = hits / private-set size ≈ c`, rounded to the copy number.
3. Emit one row per gene with its evidence. Keeping a CDS-less gene `P` in the set matters even when its own copy number is not separately validated: its private k-mers absorb its locus, so one of its copies is not mis-counted as a coding paralog.

```text
one haplotype:  A  hits=990/1200  dosage=0.83  -> CN=1
                B  hits=560/800   dosage=0.70  -> CN=1
                P  hits=5500/4800 dosage=1.15  -> CN=1
```


#### Near-identical pair (per-site consensus)

Pooled dosage suffices for divergent paralogs but blurs on near-identical pairs: two paralogs `A` and `B` that differ at only a handful of coding sites, where gene conversion makes individual copies mosaics — so two converted `A` copies count like one `A` + one `B`. When the pairwise k-mer Jaccard flags a pair as near-identical, `call` splits it by per-site consensus. Take a haplotype whose module total (from coverage) is 2 copies:

```text
per-site A-allele vs B-allele k-mer fraction, one haplotype (module total = 2):
  site 1  A=1.0 B=0.0
  site 2  A=1.0 B=0.0
  site 3  A=0.5 B=0.5   <- a copy gene-converted at this site
  ...     (more sites at A=1.0 B=0.0)
```

1. Align the two coding sequences once (reference-side) to enumerate the divergent columns; at each column collect the k-mers carrying the A-allele and the B-allele.
2. Per haplotype, count each site's A-allele vs B-allele k-mers. Because every copy carries some allele at every site, `A + B` sums to the module total, so each site independently reports the split as the fraction `A/(A+B)`.
3. Take the median of that fraction across all sites: the converted site (`0.5`) is an outlier out-voted by the clean sites (`1.0`), so the median is `1.0`. Split the reliable module total by it → `A = round(2 × 1.0) = 2`, `B = 0`. A pooled count would have read the mosaic as one copy each.

```text
SVTYPE=DUP ...  gene A CN=2   gene B CN=0     (per-site recovers the mosaic; pooled dosage would miss)
```

Every `dup_gene_cn.tsv` row carries `hits`, `priv_kmers`, and `dosage`, so the split is auditable and the copy number carries through the whole spectrum — a deletion reads `0`, a single copy `1`, a duplication `2`. A gene with no CDS is sketched from its gene span instead of its coding sequence; a gene with no private k-mers at all (indistinguishable from a paralog) is reported as the module total with `reliable=0`.

