# Fragment-level genotyping: a working prototype, and what it measures

Built to run beside `panvar genotype`, not instead of it. Nothing in production moved: `genotype-frag`
is a separate subcommand behind the same `PANVAR_ENABLE_EXPERIMENTAL_GENOTYPE` gate.

This follows `FRAGMENT_EVIDENCE_PREREGISTRATION.md` (arm F4) and the architecture review that proposed
replacing the evidence model outright. It does not settle that proposal. It answers a narrower
question that had to come first, and the answer changed the plan twice.

## What was built

`src/genotype_fragments.cpp`, two modes over the same fragment likelihood:

| mode | candidate | unit of observation |
|---|---|---|
| block-local (default) | one block's alleles, in a flanked context | physical fragment |
| `--haplotype-mode` | a pair of COMPLETE panel haplotypes | physical fragment |

Both score `log L = SUM over fragments log( 0.5 P(frag|a) + 0.5 P(frag|b) )`, with `P(frag|candidate)`
from both mates' infix alignment edit distance and the implied insert length. Each fragment enters
once. Mates are joined by read name, so interleaved and split R1/R2 input both work and neither has
to be declared.

The starting premise held: `count_reads` discards the read name and the mate at the first step, so
production genuinely has no access to fragment linkage, and every mechanism built to compensate --
confinement, over-expected filtering, marker clumps, the `rho` discount, the adjacency channel -- is
reconstructing something the representation threw away.

## Three defects found by building it, each measured

**1. Block-local scoring is structurally blocked at cyp2d6, and not by its scoring.** Under
leave-ZERO-out -- reads simulated from haplotypes that ARE in the panel -- the truth pair at block 5
ranked 14 of 105. The reason is in the sequence, not the model: panel allele 7 (556 bp) CONTAINS the
truth label's allele 0 (261 bp), and both begin at the same position in the sample's own haplotype.
The remaining 295 bp of the sample's real sequence lives in the NEIGHBOURING block. Any block-local
context has to glue some flank there; the sample's flank is unknown, which is the genotyping problem;
and a fragment spanning the junction then correctly prefers whichever candidate supplies more real
sequence before the guess begins.

**The evidence that decides the block sits across the block boundary.** No reweighting reaches that,
which is consistent with section 4 of `GENOTYPE_STATUS_2026-08-25.md`, where nothing did.

**2. A Gaussian insert-size term has unbounded influence.** NA18939's own haplotype 1 carries 13.6 kb
the panel's typical haplotype does not. Its two mates anchored to different copies, the implied insert
came out at 13.6 kb, and `z^2/2` alone cost 35,000 nats -- so the sample's own haplotype scored 15x
worse than an unrelated one. Fixed two ways: mates are now placed JOINTLY over their candidate
anchors, and the insert term is mixed against a uniform so a discordant pair costs a bounded amount.
This is the same unbounded-influence defect the status doc records for the negative-binomial emission,
in a new channel.

**3. A shortlist-wide anchor cap made the answer depend on the shortlist size.** Identical reads put
NA18939's truth at rank 2 with `--max-haplotypes 48` and rank 1 with 96, purely because the larger
shortlist pushed more syncmers past a shared occurrence cap. The cap is now per haplotype.

## Result: 10 donors, simulated 30x, paired per block

### Leave-zero-out -- the disqualifying test

| stratum | prod | haplotype mode |
|---|---:|---:|
| all blocks | **190/190 (100%)** | **190/190 (100%)** |
| bubble | 90/90 | 90/90 |
| backbone | 80/80 | 80/80 |

The prototype passes. Before the three fixes above it did not: it lost 8 blocks, all on one donor.

### Leave-one-out

| locus | stratum | prod | haplotype mode | **non-mosaic ceiling** |
|---|---|---:|---:|---:|
| cyp2d6 | all | **84/123 (68%)** | 78/123 (63%) | **104/123 (85%)** |
| cyp2d6 | bubble | 41/68 (60%) | 40/68 (59%) | 59/68 (87%) |
| gstm1 | all | 86/149 (58%) | 87/149 (58%) | **117/149 (79%)** |
| gstm1 | **bubble** | 56/87 (64%) | **62/87 (71%)** | 75/87 (86%) |

Per donor, cyp2d6: the prototype fixes 6 blocks and breaks 12, negative on 5 donors and positive on 3.
gstm1: fixes 17, breaks 16, positive on 5 donors and negative on 3.

**The prototype does not beat production, and on cyp2d6 it loses.** Its only clear gain is gstm1's
bubble blocks, 64% -> 71%, which is where the alleles shorter than k live -- the case a marker model
cannot see at all and a fragment spanning the allele and both its boundaries can.

## The measurement that decides what to do next

The **non-mosaic ceiling** column is the best any pair of COMPLETE panel haplotypes could reach
against this truth, computed exhaustively over all pairs (`scripts/genotype_pair_ceiling.py`). It is
not the panel's ceiling: a mosaic model may take a different haplotype at every block, and under that
model every representable block is reachable.

Two things follow, and the second reverses what I expected:

1. **The ceiling is far above production** -- 85% against 68% at cyp2d6, 79% against 58% at gstm1. A
   whole-haplotype-pair model that reached its own ceiling would beat the current caller
   substantially, **with no mosaic layer at all**.

2. **The prototype is far below that ceiling** -- 63% against 85%. So the missing mosaic layer is NOT
   what is costing it. Twenty-two points are available inside the architecture it already has, and the
   fragment likelihood is not collecting them.

I expected the opposite: that a non-mosaic model would be capped near production and that the gap was
the mosaic layer. That is refuted. **Building the factor graph now would be building the wrong thing.**

## What this says about the proposed rewrite

The architectural claim -- score physical fragments against complete haplotypes, not marker counts
against blocks -- is supported. Finding 1 is independent evidence for its central point, arrived at
from the sequence rather than from the argument.

The proposed ORDER is not supported. The review's own POC 1 gate ("substantially higher call rate")
would have been read as a pass on gstm1 bubble blocks and a fail on cyp2d6, with no way to tell which
mattered. The ceiling column is what makes the result actionable, and it says the next work is inside
the whole-haplotype likelihood, not past it.

Specifically **not** yet done, and now known not to be the priority: the block factor graph, mosaic
transitions, the HMM.

## Next, in order

1. **Close the gap to the non-mosaic ceiling.** 22 points at cyp2d6, 21 at gstm1. Per-block, the
   ceiling knows which pair is right and the prototype does not: that difference is directly
   inspectable with `--top-pairs` and `.hap_scores.tsv`, which report each haplotype's sequence
   likelihood, placed fragments, covered windows and coarse-shortlist score separately.
2. **The depth channel is real and is not yet a safe default.** `--coverage-weight` scores each
   haplotype's own length in windows, which is the only thing that can see copy number -- reads from a
   duplicated segment align perfectly to a single copy. At weight 1 it reuses the fragments the
   sequence term already used, and it broke leave-zero-out on NA18939. Default 0, same conclusion this
   project reached for `--edge-weight`, and for the same reason.
3. **Only then** ask whether mosaics are worth building, by re-measuring the gap to the ceiling and
   the gap from the ceiling to full representability.

## Caveats this rests on

- 10 donors, one seed, two loci, simulated reads (`wgsim`: uniform error, no GC or fragment bias).
  The pre-registration asks for multiple seeds and 28 donors; this is neither.
- Per-block observations within a donor are correlated; the per-donor net counts are reported for
  that reason and should be read before the aggregate.
- The prototype has no chain, no linkage and no HMM, so production is being compared against a model
  with a whole layer missing. That is deliberate, and the ceiling column is what makes the comparison
  interpretable anyway.
- `--haplotype-mode` costs about 9 s per cyp2d6 donor at 30x against about 4 s for production.

## Reproducing

```
cmake -S . -B build_genotype_review -DPANVAR_ENABLE_EXPERIMENTAL_GENOTYPE=ON
cmake --build build_genotype_review --target panvar -j8

# leave-zero-out first, always: a failure there is an implementation defect, not a panel limit
OUT=/tmp/ab LOO=0 ARMS="prod ceiling hap" tests/genotype_frag_ab.sh cyp2d6 10
OUT=/tmp/ab LOO=1 ARMS="prod ceiling hap" tests/genotype_frag_ab.sh cyp2d6 10
```

Both arms see identical reads, an identical panel and an identical exclusion, so a difference between
them is the evidence model and nothing else. Per-block rows land in `$OUT/<locus>/blocks.tsv`.
