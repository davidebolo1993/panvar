# Module `benchmark` - algorithm

Mechanism for `benchmark`. For commands and outputs see [modules/benchmark.md](../modules/benchmark.md); references are listed in [references.md](../references.md).

## How it works

### 1. Build a truth-event ledger

For every reference-traversed bubble, the reference walk is compared with each haplotype walk. Both are sequences of oriented graph steps. A longest-common-subsequence alignment supplies shared anchors; each non-identical interval between anchors becomes one truth event.

For each event:

- `ref_bp` and `hap_bp` are the sequence lengths on the two walks;
- `size_bp = max(ref_bp, hap_bp)`, matching the size rule used by `call` for a replacement;
- a call is attributed only when one of its nodes lies inside the event;
- the event becomes `called`, `missed` or `below_threshold` according to `--min-sv-bp`.

Sizing events from graph walks avoids an ambiguity of sequence alignments: the same edit can be split into several short mismatch runs by different equally good alignments.

### 2. Build reconstruction ceilings

The reference is used as a template. At selected truth events, its steps are replaced with the haplotype's true steps:

| reconstruction | truth event is pasted when |
|----------------|----------------------------|
| `graph` | it overlaps any node touched by a call in that bubble |
| `called` | it is attributed to a record and reaches the size threshold |
| `carrier` | it meets `called`, and this haplotype carries an overlapping record |

These levels deliberately use the true haplotype sequence. They separate discovery, attribution and genotyping from VCF representation, but cannot prove that REF/ALT is correct.

An additional in-scope floor pastes every event at or above the threshold. Its remaining error is entirely sub-threshold variation and is used to partition the final residual.

### 3. Reconstruct from the VCF

With `--vcf`, sample names are joined exactly to graph paths and records carried by each haplotype are applied to the reference in descending coordinate order.

| record | operation |
|--------|-----------|
| explicit ALT | replace REF with the selected allele |
| `DEL` | remove `|SVLEN|` bases after the anchor |
| `INS` | insert `INSSEQ` after the anchor |
| `INV` | reverse-complement `POS+1..END` |
| `DUP` | tile the inferred reference span using `CN` or the `CNBP` length change |

Explicit alleles, deletions, insertions and inversions are sequence-resolved. Duplication projection is marked `heuristic`. Records that cannot be applied are counted by reason rather than silently ignored.

The allele VCF follows a different contract: one explicit record represents all alleles of one bubble. It can test serialization accurately, but it has no one-to-one join to the merged call records in `variant_nodes.tsv`; per-call carrier attribution is therefore not reported in this mode.

### 4. Score and explain the residual

Every reconstruction is globally aligned to the true haplotype sequence. The primary quantity is edit distance (`delta`); identity and QV are derived presentations of the same comparison.

All cross-level totals use the same set of haplotype-bubble observations. This prevents an apparently better level from arising only because hard observations were omitted.

The VCF residual is partitioned into:

1. `out_of_scope` — sub-threshold truth;
2. `discovery_or_attribution` — eligible truth not covered by a record;
3. `carrier_missed` — a record exists, but the haplotype is not called as a carrier;
4. `representation` — the carried record does not reproduce the truth;
5. `false_positive_damage` — the VCF edits a haplotype without eligible truth there.

The five terms sum to the VCF residual. The last two may be signed because applying a real VCF edit can occasionally improve or worsen a region relative to a truth-pasting ceiling.

The unchanged reference supplies the baseline for:

```text
variation_recovered = (baseline_delta - level_delta) / baseline_delta
gap_closed           = (baseline_delta - region_vcf_delta)
                       / (baseline_delta - graph_delta)
```

`gap_closed` is undefined when the baseline and graph ceiling have the same error.

## Worked trace

The steps below follow the four above, one for one. Threshold `--min-sv-bp 50`. One haplotype `H`
against the reference, five bubbles, and the records `call` emitted:

```text
bubble  reference walk   H's walk      divergent stretches
1       1+,2+,3+,4+      1+,3+,4'+     node 2 deleted, 60 bp   -- and node 4 replaced, 20 bp
2       5+,6+,7+         5+,7+         node 6 deleted, 80 bp
3       8+,9+,10+        8+,10+        node 9 deleted, 70 bp
4       11+,12+,13+      11+,13+       node 12 deleted, 100 bp
5       14+,15+          14+,15+       none

ID   BUBBLE  SVTYPE  SVLEN  nodes      GT(H)
v1   1       DEL     -60    2,3        1        node 3 is shared with the 20 bp stretch
v3   3       DEL     -70    9          0        emitted, but H is not genotyped as a carrier
v4   4       DEL     -88    12         1        encodes 88 bp of a 100 bp deletion
v5   5       INS     +55    14         1        no true difference here
                                                bubble 2 has no record at all
```

1. Build a truth-event ledger. Lining each walk up against the reference gives six divergent
   stretches. Sized from the walks: 60, 20, 80, 70, 100, and nothing at bubble 5. The 20 bp stretch is
   `below_threshold`. Bubble 2's 80 bp stretch is `missed`, since no record's node lies inside it. The
   60, 70 and 100 bp stretches are `called`. `v5` creates no truth event, because there is no
   difference for it to describe.

2. Build reconstruction ceilings. `graph` pastes any stretch touching a node some call touches: `v1`
   holds node 3, which the 20 bp stretch also uses, so that stretch is pasted too even though nothing
   was called on it. `called` requires attribution and the size threshold, so it pastes the 60, 70 and
   100 bp stretches and leaves the 20 bp one. `carrier` drops the 70 bp stretch as well, because
   `v3` has `GT=0` for this haplotype.

3. Reconstruct from the VCF. Only the edits `H`'s genotypes name are applied: `v1` reproduces its
   deletion exactly, `v4` removes 88 of 100 bp and leaves 12, and `v5` inserts 55 bp where nothing is
   missing. `v3` is not carried, so bubble 3 stays reference. Bubbles with no record stay reference.

4. Score and explain the residual. Each reconstruction is aligned to `H`'s true sequence. The
   stretches here are disjoint indels, so their edit distances add:

```text
stretch          baseline  graph  called  carrier  region_vcf  loss term
60 bp deletion         60      0       0        0           0  -
20 bp replacement      20      0      20       20          20  out_of_scope
80 bp deletion         80     80      80       80          80  discovery_or_attribution
70 bp deletion         70      0       0       70          70  carrier_missed
100 bp deletion       100      0       0        0          12  representation
55 bp insertion         0      0       0        0          55  false_positive_damage
total delta           330     80     100      170         237
```

   The five terms sum to `20 + 80 + 70 + 12 + 55 = 237`, exactly the region-VCF residual, and the
   headline figures follow:

```text
variation_recovered(graph)      = (330 - 80)  / 330 = 75.8%
variation_recovered(called)     = (330 - 100) / 330 = 69.7%
variation_recovered(carrier)    = (330 - 170) / 330 = 48.5%
variation_recovered(region_vcf) = (330 - 237) / 330 = 28.2%
gap_closed                      = (330 - 237) / (330 - 80) = 37.2%
```

A high `graph` or `called` ceiling therefore proves nothing on its own: carrier assignment, record
representation and false-positive edits only enter at the later levels.
