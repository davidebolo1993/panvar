# Pre-registration: where do the shared-marker counts come from?

Written before any attribution was computed. Follows `CYP2D6_ORACLE_RESULT.md`, which established
that unfiltered markers restore identifiability under ideal counts (24/24) and fail on observed reads
(0/24), with 48% of the shared-marker mass unexplained by blocks 3 and 5.

The question is **not** "is there extra mass" -- that is measured. It is **where it comes from** and
**whether it is what breaks the decision.**

## Check 1: full-haplotype multiplicity. The decisive one.

Scan the two complete truth haplotype sequences directly, per donor, per marker code:

```
total_truth_multiplicity   = occurrences of the marker in both truth haplotypes, whole locus
target_multiplicity        = what blocks 3 and 5 account for
outside_target_multiplicity = total - target
```

Then re-run A3 and A4 on observed counts using **`outside_target_multiplicity` measured from the
sequence**, not reconstructed from block annotations.

This is decisive because it removes the block decomposition from the accounting entirely:

| outcome | reading |
|---|---|
| observed recovery **succeeds** | the block model was missing real marker sources; build a locus-wide shared-marker factor |
| observed recovery **still fails** | multiplicities were already right; the fault is the count likelihood, depth model or sampling model |

## Check 2: which markers actually drive the wrong decision

Pooled excess mass proves nothing about causation -- it could sit entirely on markers that do not
separate the candidates. For the certified pair against the observed optimum, per marker:

observed count; predicted count under the certified pair; predicted under the competitor; the
resulting delta log-likelihood; source category; filter fate.

Ranked by |delta LL|. **We need to know whether the wrong call is caused by five markers, twenty, or
diffuse inflation.** Those imply different fixes and the pooled number cannot distinguish them.

## Full count ledger

Per donor, per marker, partition the observed count into: block 3; block 5; other bubble blocks;
backbone and flanks; additional occurrences within the same allele; canonical-sequence collisions;
sequencing-error matches; unexplained. Because reads are simulated, each read's source haplotype and
position are known, so this should close exactly apart from sequencing error.

## Already excluded by construction, not to be re-tested

- **One observed count shared across marker slots.** Slot-to-code is strictly 1:1 over all 36 markers.
- **Repeated occurrences inside one allele.** Per-allele multiplicity reaches 10 and is represented.

## Controls, declared now

- Score every unique marker **sequence** once; never once per slot.
- Repeat over several read seeds, so systematic misspecification is separable from sampling variance.
- Keep the 3-pair equivalence class intact; never force one representative.
- The three donors with no representable truth at block 5 get a separate certified-floor analysis and
  are not counted as failures.

## What would close this route

If check 1 fails and check 2 shows diffuse inflation rather than a handful of decisive markers, then
the restored markers cannot be used on observed reads under any count model this caller can build, and
CYP2D6 marker restoration is closed in favour of fragment evidence.

## No production change

Not before check 1 passes and check 2 identifies the driving markers.
