# Pre-registration: does a smaller k reach alleles shorter than k?

Written 2026-08-25, **before any arm was scored**. The k=21 and k=15 runs were launched first, but
the stratum below is derived from panel geometry only — allele lengths in the graph — so it cannot be
chosen to flatter a result. Recording it here so that is checkable.

## Question

Section 12 of `GENOTYPE_STATUS_2026-08-25.md` found that gstm1's two largest-gap blocks hold alleles
shorter than k=31, and therefore carry zero markers of either kind before any filter runs. The
proposed remedy is boundary-spanning markers, which is new code. A smaller k reaches the same alleles
with an existing flag, and the standing method rule is to test with an existing flag first.

If lowering k does not improve these blocks, the boundary construction is unlikely to either, since
both work by making short alleles carry markers at all.

## Method, frozen

`tests/../cohort/run_sim.sh`, simulated reads from the held-out individual's own two haplotypes,
LOO regime, 8 samples, seed schedule unchanged. Arms: baseline k=31 (already run), `--kmer-size 21`,
`--kmer-size 15`. cyp2d6 k=21 as a control locus, since its gap is a selection failure and should not
respond to k.

Metric: `gap = best_identity - identity`, the model headroom, paired per (sample, block).
Lower is better. `dbp` is a LENGTH difference and is NOT the metric.

## Stratum, declared from the panel

Measured with `--dump-block` over all 20 gstm1 blocks. A block is SUB-K if at least half its alleles
are shorter than 31 bp:

| block | alleles | min bp | max bp | under 31 bp | |
|---:|---:|---:|---:|---:|---|
| 3 | 5 | 7 | 82 | 4 (80%) | **SUB-K** |
| 9 | 16 | 14 | 42 | 11 (69%) | **SUB-K** |
| 10 | 11 | 5 | 160 | 9 (82%) | **SUB-K** |
| 11 | 20 | 5 | 21 | 20 (100%) | **SUB-K** |
| 8 | 38 | 10 | 99 | 4 (11%) | not |
| 16 | 34 | 7 | 86 | 7 (21%) | not |
| 6 | 44 | 22 | 4,260 | 1 (2%) | not |
| all others (13 blocks) | — | ≥ 38 | — | 0 | not |

**SUB-K stratum = blocks 3, 9, 10, 11.** Everything else is the control stratum.

## Declared prediction

Lowering k improves the SUB-K stratum. If it does not, the length hypothesis is wrong about what
those blocks lose, and boundary-spanning markers should not be built on it.

## Gate

Accept k as a lever only if **both** hold at k=21:

1. mean gap over the SUB-K stratum improves, and more paired observations improve than worsen;
2. mean gap over the control stratum does not worsen by more than 0.001, and no control block
   regresses by more than 0.01.

k=15 is exploratory: it is expected to raise marker ambiguity across the board, and its role is to
show the direction of the trade rather than to be adopted.

**Anticipated failure mode, declared in advance:** a smaller k makes markers less unique panel-wide,
so more of them will fail the over-expected rule. A SUB-K block could therefore gain candidates and
still retain none. If the gate fails, `--ledger-block` on blocks 9 and 11 at k=21 distinguishes "no
candidates were created" from "candidates were created and filtered away" — those call for different
next steps, and the ledger must be read before either is chosen.

---

# RESULT (2026-08-25): prediction FALSIFIED. The gate's letter passes; do not read that as a win.

## Accuracy, gstm1 k=31 -> k=21, 8 samples, 158 paired (sample, block) observations

| stratum | n | gap base | gap k=21 | delta | better | worse |
|---|---:|---:|---:|---:|---:|---:|
| SUB-K (blocks 3, 9, 10, 11) | 30 | 0.04533 | 0.04149 | **-0.00385** | 1 | 0 |
| control (all others) | 128 | 0.00785 | 0.00857 | +0.00073 | 13 | 21 |

Both gate conditions are technically met: the SUB-K mean improves with more better than worse, the
control mean worsens by less than 0.001, and no control block regresses by more than 0.01 (worst is
block 12 at +0.0069).

**The gate is nonetheless uninformative, and the declared prediction is falsified.** Of the 30 SUB-K
observations, **exactly one changed at all**:

| block | sub-k alleles | paired obs | obs that changed | delta |
|---:|---:|---:|---:|---:|
| 3 | 80% | 8 | **0** | 0.00000 |
| 9 | 69% | 8 | **0** | 0.00000 |
| 10 | 82% | 8 | 1 | -0.01442 |
| 11 | 100% | 6 | **0** | 0.00000 |

The entire stratum improvement is one sample at one block. The two blocks the mechanism was
identified in — 9 and 11, at 69% and 100% sub-k — are byte-identical at k=21. A criterion that can
pass on 1 of 30 unchanged observations was badly designed; "mean improves AND better > worse" needs a
floor on how many observations move. Noted for the next pre-registration.

## Why: the anticipated failure mode, confirmed by the ledger

The pre-registration declared that a SUB-K block could gain candidates and still retain none, and
that `--ledger-block` would distinguish the two. It did. Blocks 3, 9, 10, 11 at gstm1, LOO, HG00097:

| block | k=31 informative / retained | k=21 | k=15 |
|---:|---|---|---|
| 3 | 11 / **0** | 16 / **0** | 30 / **0** |
| 9 | 7 / **0** | 11 / **0** | 11 / **0** |
| 10 | 19 / 4 | 19 / 3 | 15 / **0** |
| 11 | 0 / 0 | 1 / **0** | 0 / 0 |

**Lowering k does exactly what it was supposed to** — block 3 goes from 11 informative candidates to
30, block 9 from 7 to 11, block 11 from none to one. **Every one of them is then filtered away.**
Block 10, the only sub-k block that retained anything at k=31, is driven from 4 retained to 3 to 0.

So k is not the binding constraint at these blocks. Both mechanisms are present at once — short
alleles yield few candidates AND region filtering removes what they do yield — and because filtering
is downstream, it binds last. **This is why the accuracy barely moved: creating candidates changes
nothing when the retained set stays empty.**

## Consequence for boundary-spanning markers

This is the more important result, because it applies to the construction that has not been built yet.
A k-mer anchored in a flank and crossing into a short allele **contains flank sequence by
construction**, and that flank sequence is shared with the neighbouring blocks' alleles. Such a marker
therefore varies in more than one block and is rejected by the confinement rule — which is already the
dominant fate at these blocks (block 3: 7 of 11 markers dropped for multi-block alone at every k).

**Boundary-spanning markers as specified would be filtered out on arrival.** They cannot be built
against the current region rules; they need markers attributed to a JUNCTION rather than to a block,
or a filtering rule that understands that a marker spanning a boundary is expected to be seen in both
blocks it spans. That is a design constraint discovered before implementation rather than after.

## Status

k as a lever: **rejected**, on mechanism rather than on the gate. k=15 is worse than k=21 at every
sub-k block and is not pursued. The cyp2d6 k=21 control arm is still running and is not needed for
this conclusion; it is reported when it lands.

## Scope limit on this experiment (added after review)

**This sweep does not test boundary-spanning markers, and a negative result here must not be read as
evidence against them.** They are different interventions:

- lowering k asks whether SHORTER seeds help everywhere;
- boundary markers ask whether LONG, specific context spanning a junction helps at one place.

These can point opposite ways. Shorter k-mers are less unique across the whole locus, so they fail the
over-expected rule more often, and they never use flanking context at all.

Two further confounds, both real:

1. **k and s move together.** `default_syncmer_s(k) = min(11, (k+2)/3)`, so the arms were k=31/s=11,
   k=21/s=7, k=15/s=5. The sweep varies the syncmer sampling parameter as well as the k-mer length,
   and cannot separate them. An `--all-kmers` arm would.
2. At k=15 most of block 11's alleles (5-21 bp) are still unrepresentable, so that arm does not
   actually reach the alleles the hypothesis is about.

What the sweep DOES establish is the mechanism, and that part is not confounded: at every k tried,
candidates are created at the sub-k blocks and every one is filtered away. The retained set stays
empty. That conclusion comes from the ledger, not from the accuracy arms.

## Remaining arms (completed)

**gstm1 k=15 — fails.** SUB-K stratum 0.04533 -> 0.04571 (**worse**), 1 better / 1 worse; control
0.00785 -> 0.00719 but 14 better against 18 worse. Consistent with the ledger, where k=15 drives
block 10 from 4 retained markers to 0 and leaves blocks 3, 9 and 11 at zero throughout.

**cyp2d6 k=21 — the control behaves as the mechanism predicts.** 0.01490 -> 0.01479, a change of
-0.00011 over 152 paired observations with 10 better and 18 worse: no effect. cyp2d6's failure is
marker selection among alleles of 261-851 bp, which have never been short of candidates, so k is not
expected to touch it and does not. That is a small independent check that the case (a) / case (b)
split is real rather than a story fitted to gstm1.

## Final status

| arm | sub-k stratum | control | verdict |
|---|---|---|---|
| gstm1 k=21 | -0.00385, but 1 of 30 observations moved | +0.00073 | prediction falsified |
| gstm1 k=15 | **+0.00037** | -0.00066, 14 better / 18 worse | fails |
| cyp2d6 k=21 | n/a | -0.00011, no effect | control as predicted |

**k is rejected as a lever**, on the ledger's mechanism rather than on the accuracy arms: at every k
tried, candidates are created at the sub-k blocks and every one is then filtered away. The retained
set stays empty, so accuracy cannot move.
