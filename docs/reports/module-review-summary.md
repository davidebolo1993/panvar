# Module review: corrected measurements

Dated 2026-08-21, at the release review.

This is the record of numbers this project reported and later had to correct. It is history, not
unfinished work -- `docs/review_followups.md` holds only the latter, and its own introduction says
completed work should not stay there. It is kept because several of these figures were quoted for
months before they were checked, and anyone reading an older note, plot or commit message needs to be
able to find out which reading was wrong and why.

Each entry states the old number, the corrected one, and the reason the first was wrong.

## Benchmark

### Corrected measurements — what the old numbers meant

- **Every QV figure recorded for this project before the ledger is the `graph` column**, which is an
  optimistic upper bound: a block is substituted when it shares a node with ANY call at the bubble, and
  no genotype is read at all. It is kept, unchanged and relabelled. The `called` reconstruction beside
  it is the strict form.
- **`over_threshold_bp` read 0 at all six reference loci and could not do otherwise.** The residual was
  split by contiguous alignment-run length, and edlib's co-optimal edit path distributes a gap over
  chance matches: a clean 60 bp deletion of non-repetitive sequence returns as fourteen runs of 1–10
  bases. Any claim resting on that split, and on the carrier truth flag built from it, is withdrawn.
- **"Recall is 100% everywhere; FN = 0 at all six loci" is withdrawn.** It followed from the truth flag
  almost never firing. On the walk-derived ledger, ACOT has 145 genotyping false negatives and CYP2D6
  has 3. One ACOT case is a 56,889 bp divergent block at bubble 7 whose haplotype is genotyped `0` at
  every one of the nine records there.
- **Carrier precision was understated, substantially.** ANKRD36C 3.9% → 69.2%, CYP2D6 55.7% → 100%,
  LPA 77.4% → 94.7%, ACOT 76.1% → 98.6%.
- **No above-threshold truth event is left untouched by some record** at any of the six loci: 0 `missed`.
  That is a discovery *ceiling*, not a proof of correct discovery — `called` means a record shares at
  least one node with the block, not that it covers or represents it. The firm half of the result is
  the negative: had a record been absent entirely, the event would have shown as `missed`.
- **Where the loss lives is representation.** The genotype residual is partitioned exactly, over one
  common observation set, into five consecutive steps. Across the six loci: `discovery_or_attribution`
  is **0 bp everywhere** (the retained calls restore exactly as much as a perfect in-scope caller
  would), `carrier_missed` is 0 at four loci and at most 15.6 kb, `false_positive_damage` is small and
  sometimes negative, and `representation` runs 43 kb to 8.57 Mb. Out-of-scope sub-threshold variation
  is 13–78 kb and is reported separately, since `call` was never asked to emit it.

  Three accounting faults were found in the first version of that partition and are fixed: comparative
  totals were taken over different populations (a partial VCF join put `carrier` ABOVE `called` and
  produced a loss of minus 20 bp); false-positive carrier damage was charged to representation; and
  below-threshold variation was charged to discovery. Each has a fixture.


## Describe

### Corrected measurements

- **A pure deletion produced no features at all.** All three step-building sites required an interior
  bubble node, so a path taking the direct source→sink edge was not counted as a traverser; the node
  that discriminates the deletion was then dropped as non-discriminative *because only its carriers had
  been observed*. Real effect, like-for-like: cyp2d6 +245 path observations (+38%), lpa +496 (+13%),
  c4 +26, and c4's 1340 "NA" cells were dropped traversals rather than missing data. Any dosage matrix
  built before this is missing its deletion carriers.
- **Both association gates re-run and pass.** The LPA KIV-2 signal is unmoved (bubble 7, beta −0.0193,
  p_bonf 6.1e-06) on a 13% larger substrate, and a matched null (same flags, only the code differing)
  gives type-I 0.05052→0.05082 at 0.05 and 0.00081→0.00081 at 0.001, lambda 0.9094→0.9116, with the
  worst per-feature KS p and the maxT threshold bit-identical.
- **The ~23% per-feature non-uniformity in that null is pre-existing**, present identically before and
  after, and belongs to `associate`'s recorded limits rather than to this pass.

## Bubble

- **Retained sites are pairwise disjoint, and the enclosing snarl was the wrong thing to keep.** The
  recorded policy was to retain the enclosing snarl and drop the smaller ones. Measured at ANKRD36C,
  it is the opposite: one enclosing site closes 3.98% of the reconstruction gap, the ten smaller ones
  close 89.54%, at the same record count. `smaller` is the behaviour. The allele VCF is lossless under
  both, so this affects only the interpreted output, and only ANKRD36C -- the other five loci emit no
  overlapping sites and are byte-identical.

## Inspect

- **The sketch Jaccard was biased low, and the bias grew with the length ratio.** `|shared| / |union|`
  over two independently truncated bottom-k sketches read 0.0920 against a true 0.1667 at a 3x size
  gap -- 45% low, and a 3x size gap is exactly two haplotypes carrying the same repeat unit at
  different copy numbers. Replaced by the union-bottom-k estimator: 0.1684. 23 of 26 cluster files
  changed, in both directions.
- **Complete sketches were estimated rather than compared.** A walk of 512 shingles or fewer has an
  untruncated sketch, which IS its shingle multiset, and the estimator still subsampled it: 0.5 where
  0.667 was in memory. Corrected 2026-08-21. It fires on no reviewed locus -- every walk at those
  bubbles runs to tens of thousands of steps -- so the fixture is the documented worked example.

## Bubble, again: the 64-endpoint cap

- **A repeated boundary was silently clipped at 64 occurrences.** The interval search enumerated the
  endpoints after each start and stopped at 64, so a path revisiting a boundary more often reported a
  truncated interior: a 70-repeat path read 64. Corrected 2026-08-21, exactly, by observing that the
  interior count is monotone under containment -- so the widest interval attains the maximum and the
  shortest one attaining it is the tightest window around the interior steps. Byte-identical on all
  six reference loci, because a top-level snarl boundary is a cut node each path crosses once; it was
  reachable through `--snarls-in` and through any externally supplied bubble CSV.
