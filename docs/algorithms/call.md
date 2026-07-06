# Module `call` - algorithm

Mechanism and a hand-traced example for the `call` module. For usage/flags see [modules/call.md](../modules/call.md); References in [references.md](../references.md#call).

## Terms

- **anchor** — a node token unique in *both* the reference walk and a haplotype walk; a monotonic subsequence of anchors aligns the two walks, and only the segments between anchors are aligned.
- **gap block** — a maximal stretch between two anchors where reference and haplotype differ; each becomes one typed event (DEL/INS/INV/substitution).
- **full walk** — the widest `source→sink` span with every repeated traversal counted (a node visited twice contributes its length twice); contrast the **minimal-span** walk that visits each distinct node once.
- **peak multiplicity** — the visit count of a haplotype's single busiest node in the bubble.

## Algorithm

For each bubble:

1. Group identical `source→sink` walks into distinct alleles; call once per allele (genotypes expand by membership).
2.  Fold a node's consecutive self-repeats (a `REP` self-loop) into one alignment token, so extra copies surface as copy number, not a spurious INS/DEL.
3. Split both walks at their shared anchors and align only the segments between them. Each gap block becomes one typed event: a `DEL` where only the reference has nodes, an `INS` where only the haplotype does, an `INV` where the haplotype's run is the reverse-complement node-walk of the reference's, or a substitution — a co-located `DEL` and `INS` that share one `EVENTID`.
4. Merge consecutive same-type events whose gap is ≤ `--merge-distance-bp`, measured in either reference space or the haplotype's own sequence space.
5. Merge across haplotypes. Transitive single-linkage (connected components), seeding from events down to `--rescue-min-bp`; two events link on same type + a position window + (Jaccard ≥ `--merge-jaccard` or sequence identity ≥ `--merge-seq-identity`, gated by `--merge-size-ratio`). Largest member represents; `MERGE_*` records the evidence. Copy-number records merge separately on shared `REP` identity.
6. Re-test every non-carrier against its own diff, adding it as a carrier when it matches (sub-threshold events get `GT=1` instead of `0`); one pass (representative is fixed). Keep records whose representative reaches `--min-sv-bp` and carrier count reaches `--min-haplotypes`.

## Worked trace

Reference interior `R = A B C D E F` — single-letter node tokens (a `REP` self-loop, if present, is dropped from this walk and typed separately as a `DUP`, so extra copies are never mistyped as an insertion). Four haplotypes cross the bubble:

```text
ref  : A B C D E F
hX   : A B   D E F          C missing
hW   : A B   D E F          C missing (same walk as hX)
hZ   : A     D E F          B and C missing
hI   : A B P Q C D E F      P, Q inserted after B
```

1. Diff each haplotype against the reference. Anchoring on the tokens the two walks share and aligning only the segments between anchors: `hX` and `hW` each read a `DEL{C}`, `hZ` reads a single `DEL{B,C}` (both missing nodes fall in one gap block), and `hI` reads an `INS{P,Q}` anchored after `B`. A haplotype run that reverse-complemented the reference's would instead type as an `INV`, and a co-located reference-only + haplotype-only pair as a substitution — a `DEL` and `INS` sharing one `EVENTID`.
2. Coalesce within a haplotype. Same-type events that sit in *different* gap blocks but within `--merge-distance-bp` are joined into one, measured either in reference space or in the haplotype's own sequence space (the latter catches insertions that a large intervening deletion would otherwise push apart). None of these four trigger it.
3. Cluster across haplotypes by transitive single-linkage. `hX:DEL{C}` and `hW:DEL{C}` have identical node sets (Jaccard 1.0), so they become one record with `AC = 2`; `hZ:DEL{B,C}` shares only `C` (Jaccard `1/2 < --merge-jaccard`), so it stays its own record. Events also link when their **sequences** match, which rescues one allele threaded through different nodes.
4. Force-call, then filter. Every non-carrier is re-tested against its own diff and added if it matches — so a sub-threshold copy of an event gets `GT = 1` rather than `0` — and a record is kept only if its representative reaches `--min-sv-bp` and its carriers reach `--min-haplotypes`.

Resulting records:

```text
SVTYPE=DEL  EVENT_NODES=C     AC=2   carriers hX, hW
SVTYPE=DEL  EVENT_NODES=B,C   AC=1   carrier  hZ
SVTYPE=INS  EVENT_NODES=P,Q   AC=1   carrier  hI
```

## Copy number

`call --cn` reads copy number straight off the walk and emits it as a `DUP` record carrying a per-sample `CN`. How it reads the number depends on how the locus is folded in the graph. 

### Worked trace — a tandem repeat (self-loop `REP`)

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

### Worked trace — a collapsed paralog cluster (coverage)

Here `pggb` collapsed near-identical paralog copies onto shared nodes, and the reference itself traverses them ≥ 2×. No single node can show a loss, so we measure the total sequence over the folded nodes and calibrate one copy from the reference:

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

### Worked trace — a single-copy-reference duplication (peak multiplicity)

The reference here visits each inside node once — it does not fold — so there is no reference fold for coverage to calibrate against. This reading covers a haplotype that folds an extra copy back onto those nodes (a duplication where the reference is the single/short allele):

```text
reference : each inside node ×1            ref_peak = 1
hapD      : busiest inside node ×2         one extra copy folded back
hapE      : each inside node ×1            no extra copy
```

1. `REF_CN = ref_peak = 1`. There is no unit to divide by, so we compare peak multiplicities directly.
2. Take each haplotype's busiest inside node: `hapD` peaks at 2 (`> 1`) and is interpreted as a gain; `hapE` stays at 1. Using the global peak, not a per-node excess, keeps scattered cluster background from inflating the call.
3. Emit a `DUP` for `hapD` with `REF_CN = 1`, `CN = 2`, and `SVLEN = Σ node_len × excess` (the duplicated bp). This reading is **gains only** — a loss leaves no extra fold to count.

```text
SVTYPE=DUP  REF_CN=1      hapD CN=2   (hapE is not a carrier)
```

### Worked trace — copy number alongside sequence events

A `DUP` never hides the bubble's genuinely-novel `DEL`/`INS`/`INV` calls. Take the coverage bubble above and add a haplotype `hapF` carrying a 300 bp novel insertion inside the module:

1. The coverage `DUP` is emitted exactly as before (the total module `CN` per sample).
2. `hapF`'s insertion has different carriers and a different size from the duplicated content, so it is kept as its own `INS` record next to the `DUP`.
3. Only the copy itself is not re-emitted, to avoid counting it twice: a copy-number-loss `DEL` (≥ half a unit) is already the DUP's reduced `CN`, and the duplicated-content `INS` (the extra copy's sequence, carried by the gain samples at ~one-copy size) is already its raised `CN` — so both fold into the `DUP` instead of appearing separately.

```text
SVTYPE=DUP  REF_CN=2               (module CN per sample)
SVTYPE=INS  SVLEN=300  hapF        (kept — a novel insertion, not the duplicated copy)
```

**Without `--cn`.** Only the always-on self-loop `REP` reading fires: a tandem is still a `DUP`, but a folded extra copy (the coverage and peak cases) instead surfaces as a dup-like `INS` (`INS_SUBTYPE=DUP` under `--classify-ins`) with no per-sample copy number, and a loss surfaces as a `DEL`. Run `call --cn` to get `CN`/`REF_CN`; the sequence-resolved `DEL`/`INS`/`INV` are present either way.

## Merge keys — Jaccard vs sequence identity

Two events merge across haplotypes on the position window and either node-set Jaccard (`--merge-jaccard`) or sequence identity (`--merge-seq-identity`; tried only if Jaccard misses). They fail in orthogonal ways:

- **same content, different nodes**: Jaccard low, sequence high (one allele threaded through different graph nodes — a microsatellite tangle). The sequence gate rescues it (its main reason to exist).
- **same nodes, poorly-aligning sequence**: Jaccard high, sequence low (shared backbone with a big internal indel / low-complexity content). Jaccard rescues it (length-weighted, order/orientation-blind).
- **different sizes**: only the sequence path is gated by `--merge-size-ratio`; Jaccard has no size gate.

So lower `--merge-jaccard` to merge events sharing a backbone; lower `--merge-seq-identity`/`--merge-size-ratio` to merge similar content threading different nodes.


## Gene annotation trace (`--gtf`)

With `--gtf`, `call` records which genes each variant touches and, for a copy-number call spanning several genes, splits the module total into per-gene copy numbers — **from the GTF alone, by k-mer dosage, with no per-haplotype alignment**. The reference path's PanSN name (`sample#hap#chrom:start-end`) supplies the chromosome and start coordinate the projection needs.

First the genes are projected onto the reference nodes: walking the reference and accumulating node lengths gives each node a coordinate span, which is intersected with the GTF gene intervals to tag every node with the gene (or genes) it falls in. Where the reference revisits the same nodes at two coordinates — a folded paralog — that node is tagged with both genes, which is precisely why graph multiplicity on its own cannot tell the paralogs apart.

When a `DUP` overlaps two or more genes, `call` builds each gene's discriminative sketch from the GTF: the merged **coding sequence** (the span for a gene without CDS, e.g. a pseudogene), reduced to the canonical k-mers **private** to that gene versus its paralogs. Coding sequence is where paralogs differ, so this sketch discriminates from the reference alone. Each haplotype's module sub-walk is then scanned and each private k-mer counted; a gene present in `c` copies carries each of its private k-mers about `c` times, so the per-copy dosage (`hits ÷ private-set size`) rounds to the copy number. A gene with no private k-mers (indistinguishable from a paralog) reports the module total with `reliable=0`.

Pooled dosage suffices for divergent paralogs (CYP2D6 vs 2D7), but blurs on **near-identical** pairs like C4A/C4B: they differ at only a handful of coding sites, and gene conversion makes individual copies mosaics, so two converted C4A copies count like one C4A + one C4B. For a pair the k-mer Jaccard flags as near-identical, `call` instead splits by **per-site consensus**: it aligns the two coding sequences once (reference-side) to enumerate the divergent columns, and at each such site counts the haplotype's A-allele vs B-allele k-mers. Because every copy carries some allele at every site, the two counts sum to the module total, and the **median of the per-site A-fraction across all sites** rejects the odd converted site as an outlier — recovering the per-copy majority vote that a pooled count loses. The module total (the reliable coverage count) is then split by that median.

Each `dup_gene_cn.tsv` row carries the evidence — `hits`, `priv_kmers`, `dosage` — so a call reads directly (e.g. `C4A hits=47/509, dosage=0.09 → 0` and `C4B hits=210/492, dosage=0.43 → 1` on an A0B1 haplotype), and the copy number carries through the whole spectrum (a deletion → `0`, single copy → `1`, duplication → `2`).

