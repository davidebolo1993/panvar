# Which estimator should report array copy number? Pre-registration

Date: 2026-08-23. Frozen before the extended cohort is run. No production source changes.

## Why this exists

On simulated chimeras at LPA KIV-2, the continuous `mass_bp` estimator reached 0.42 repeat units of
length error while every pair-selection method exceeded 1.5, which suggested replacing `called_bp`
with `mass_bp` as the copy-number answer at arrays.

On 31 REAL held-out individuals that reverses. Per array block, mean absolute length error:

| block | n | `called_bp` | `mass_bp` | `mass_bp` closer |
|---|---:|---:|---:|---:|
| LPA 6 | 16 | 6 bp | 940 | 0/16 |
| LPA 13 (KIV-2) | 16 | 4,859 | 4,899 | 6/16 |
| LPA 18 | 16 | 2 | 1,025 | 0/16 |
| LPA 20 | 16 | 1 | 2,265 | 0/16 |
| CYP2D6 2 | 15 | 4 | 4,904 | 0/15 |
| CYP2D6 9 | 15 | 817 | 7,282 | 0/15 |
| CYP2D6 10 | 15 | 8 | 5,138 | 0/15 |

Pair selection is essentially exact at every array EXCEPT KIV-2, where the two estimators tie at
0.88 units. The simulated result did not generalise to real individuals, so the question is reopened
with power rather than settled by either sample.

## Frozen method

For every panel individual having both haplotypes in the locus graph and a 1000G high-coverage CRAM:

1. `samtools view` the locus +/- 20 kb from the remote CRAM (region only, never the whole file).
2. Collate to paired FASTQ.
3. `panvar genotype` with DEFAULT options, that individual's two assembly haplotypes passed to both
   `--exclude-haplotypes` and `--truth-haplotypes`.
4. Record, per block: `true_bp`, `called_bp`, `mass_bp`, `identity`, `best_identity`, `block_class`.

No flag sweeps, no tuning, no per-sample choices, no HMM changes. Samples already scored are not
re-run. A sample whose region pull fails after three attempts is skipped and reported as skipped.

## Declared outcome measures

Per array block, over all individuals:

- **A** mean `|called_bp - true_bp|` against mean `|mass_bp - true_bp|`
- **B** the fraction of observations where `mass_bp` is strictly closer
- **C** the identity headroom, `best_identity - identity`

## Declared decision rule

`mass_bp` replaces `called_bp` as the reported copy number at a block class only if, at that class,
**B > 0.60 AND its mean error is lower**. Otherwise `called_bp` remains the copy-number answer and
`mass_bp` stays a diagnostic.

## Declared prediction

On the n=16/n=15 data the rule is failed everywhere: B is 0.00 at six of seven array blocks and 0.375
at KIV-2. **The prediction is that the extended cohort will NOT pass the rule at any block**, and that
KIV-2 will remain near 0.9 units for both estimators.

If that prediction holds, the conclusion is that array LENGTH at KIV-2 is a joint limitation of both
estimators and is not fixable by choosing between them — which retires a line of work rather than
opening one, and redirects effort to the identity headroom (measure C), which on the current cohort is
0.0272 at LPA and 0.0153 at CYP2D6 excluding the one unrepresentable individual.

A result that contradicts the prediction is the interesting one and would be reported as such.

## Scope limits, stated now

Two loci, one ancestry-skewed cohort, 30x short reads, region-mapped only (no unmapped-read rescue),
and a path-holdout rather than a graph rebuild — the last measured at 0.003% node and 0.17% edge
leakage for one individual, so believed small but verified only once.
