# Mosaic genotyper: current state

The first document to read. `genotype_mosaic.md` is the contract; `genotype_mosaic_checkpoint.md` is the
archive of how we got here and is consulted only when a specific round's reasoning is needed. This file
carries only what is true now.

## Model and output contract

The locus is a chain of blocks alternating SV bubbles (>= 50 bp) with context blocks. A panel template is
one graph path spanning the chain in a consistent orientation; a template that does not reach both
terminal anchors has a **partial frame** and fails model construction closed. The only sanctioned way a
template leaves the panel is an explicit `--panel-curation` manifest, which is recorded in the status
with the reason each template would otherwise have carried and a digest of the decision.

Candidates are unordered pairs of whole haplotypes: every panel diplotype, plus up to K off-panel mosaics
assembled by the beam from alleles already in the graph. A mosaic need not exist as any graph path.

Run statuses and what each emits:

| status | meaning | emits |
|---|---|---|
| `COMPLETE` | one diplotype certified above every rival beyond floating-point error | genotypes, block genotypes, scores, candidates, phased haplotype FASTA |
| `BLOCK_CONSENSUS` | scoring finished; no whole-locus winner separable, but tied candidates agree on some blocks | genotypes and block genotypes, with disagreeing blocks as explicit `NOCALL` rows; **no** phased FASTA |
| `MODEL_INCOMPLETE` | the structural model could not be built | status only |
| `CANDIDATE_INCOMPLETE` | proposal stopped at a declared budget | status only |
| `SPELLING_INCOMPLETE` | a candidate could not be spelled | status only |
| `RESCORE_INCOMPLETE` | a resource cap was reached, or a tie with no block agreement | status only |

Block consensus is **refused**, leaving `RESCORE_INCOMPLETE` with no calls, after model failure, resource
refusal, an incomplete candidate set, truncated scoring, or missing provenance. Phase is always
`UNRESOLVED` under consensus, and no unique phased haplotype is ever emitted there. Confidence is
conditional on the scored candidate set, never global, unless the search was exhaustive and untruncated.

**A K-truncated prefix is never grounds for consensus.** Candidates below K were constructed and ordered
by proposal mass but never exactly scored, so their order proves nothing about whether one of them wins
or disagrees blockwise. When the primary ranking ties, the search was exhaustive, and the beam retained
candidates outside K, the complete retained beam is scored exactly and the output is derived from that
ranking alone -- candidate and score tables included. It needs no second scoring pass and no new
parameter: the shared pass already covers the full beam whenever the search is exhaustive.

That is a widened pass, not a free one. The pass runs before the primary ranking exists, so it cannot
know whether the ranking will tie: on an **exhaustive** search it scores the whole retained beam even
when the primary winner turns out to be unique and no audit is ever attempted. The extra work is bounded
by the beam size. Zero added work was *observed* on the five C4 plans, and only because every one of them
is non-exhaustive, so the widened pass never applied there -- that measurement does not generalise to an
exhaustive run.

Eight facts are reported separately, and `candidate_set_conditional` is not redefined -- it still
describes the prefix: `primary_prefix_k`, `full_beam_candidates`, `proposal_search_exhaustive`,
`primary_prefix_truncated`, `full_beam_audit_attempted`, `candidates_exactly_audited`, `audit_completed`,
`final_candidate_scope` (the scope of the set that actually determined the call).

A refused run reports the work it performed: `observed_piece_hits`, `observed_physical_origins`, the
placement workers actually used, and the phase timings up to the refusal.

## Frozen defaults

State limit 1,024 (production) and 4,096 (audit only). `max_piece_hits` 2e9, `max_physical_origins` 5e8,
`max_state_classes` 2e8, `max_mass_cells` 2e8, `piece_hit_publish_batch` 4,096. Mismatch band 6, minimum
piece 12 bp, k-mer 31, syncmer 15, min variant 50 bp, switch probability 0.01. Simulation: seed 20260916
salted per donor or plan, 100 bp reads, insert 350 +/- 40, 4 sigmas, substitution error 0.002, outlier
weight 0.001. C4 reference `HG00096#1#haplotype1-0000024:31848049-32074186`.

## Trusted scientific results

- **CYP2D6, curated**: excluding the one partial template `NA18989#1#haplotype1-0000031:36617797-36825011`
  makes the model buildable; 19 of 19 blocks exact against HG00096 truth, every block read-supported.
  Phase `UNRESOLVED` -- exact *unordered genotype* accuracy, not phased accuracy.
- **C4 synthetic mosaic**: the beam recovers a genuinely off-panel mosaic. `bridgeable` placed the truth
  at full-beam rank 6 and called 11 of 11 blocks exactly, for haplotypes in no panel row.
- **Where the truth is lost when it is lost**: `one_switch` is a beam-width case (absent at 1,024,
  rank 6 and never lost at 4,096). `every_block` (10 switches) is crowded out at both widths, with the
  loss point moving later as width grows. The proposal offers the truth at ranks 2, 6 and 20 on three of
  five plans at width 1,024.
- **The width is not spent on diplotypes.** Measured per layer, proposal only, on all five plans at both
  widths. Every plan's tightest layer is the same one, `ctx:4`: at width 1,024 it retains 2-4 distinct
  diplotype prefixes, with one prefix holding up to 1,018 of the 1,024 states. Four times the width does
  not buy four times the diversity -- `ctx:4` still retains 3-14 prefixes, and on `one_switch` and
  `unbridgeable` the dominant prefix simply grows to 4,080 and 4,077 of 4,096. Distinct prefixes at the
  final layer (exactly `unique_candidates_found`) go from 562-892 to 2,510-3,568 while states discarded
  go from 29-43 million to 129-180 million. What limits diplotype diversity is the selection rule, which
  ranks template-level states, not the width it is given.
- **Prefix-first selection recovers every truth, and costs more than it looks.** The layer selector now
  groups merged states by unordered prefix pair, ranks groups by aggregate mass, keeps each top group's
  best state up to the limit, then fills any spare capacity by the old state order. Parameter-free.
  Measured on the five plans at both widths, paired back to back against the state-only selector:
  - Retention: the truth survives **every** layer of **all ten** runs (`first_block_truth_not_retained`
    is `never` throughout). `every_block`, lost at both widths before (`ctx:4` at 1,024, `sv:9` at
    4,096), is now retained at full-beam ranks 711 and 783; `one_switch` at 1,024 goes from absent to
    rank 5.
  - Ranks move both ways, and this is the finding that matters for calling: `bridgeable` 6 to 4 at both
    widths and `three_switches` 20 to 12 at 1,024, but `three_switches` 8 to 12 at 4,096 and
    **`unbridgeable` 2 to 46 at both widths**. Retention is not rank: a truth retained at rank 46 sits
    outside any K the caller exactly scores, so a case that previously called well can be made worse.
  - Diversity: at `ctx:4` and `sv:5` the retained prefix count is now the full width and
    `max_states_per_prefix` is 1, against 2-17 prefixes and 417-4,080 states on the leader before.
    Final distinct candidates are exactly the width (1,024 / 4,096) against 562-892 / 2,510-3,568.
  - Cost: with every retained state carrying a distinct prefix, parent groups hold one state each and
    within-prefix merging stops. States offered per layer become exactly `limit x templates^2`
    (17,572,864 at width 1,024; 70,291,456 at 4,096) against 68,644-1,029,660 before, and states pruned
    rise from 29-43 million to 175.7 million at 1,024 and from 129-180 million to 672.3 million at
    4,096. Runtime is 3.4x to 5.3x (3.8-4.2 s to 13.1-13.6 s at 1,024; 8.1-8.2 s to 43.0-43.4 s at
    4,096) and peak RSS 1.27x to 1.41x (638-799 MB to 810-814 MB; 1.36-2.01 GB to 1.95-2.05 GB).
  - `extensions_evaluated` is **identical** between the two selectors at each width (175,745,801 and
    702,931,721), so the pre-allocation bound production checks against is unaffected.
  - **Now measured, and it was the cause.** Only each group's best state survives, so the mass on a
    prefix's other template attributions is dropped rather than carried. The grouped proposer keeps all
    of them and returns `unbridgeable`'s truth to rank 1 -- better than the rank 2 the original
    state-only selector gave and far above the 46 prefix-first gave. An earlier reading of this, that
    the truth prefix ranked about 33rd "by aggregate mass counting every attribution", was wrong: that
    aggregate was measured inside a prefix-first beam whose parents had already been truncated at every
    preceding layer, so it was never mass-complete. The loss compounds layer over layer.
  - The search is still not exhaustive anywhere, so none of this reaches block consensus.
  - **Experimental, not frozen.** `954af2f` is an experimental selector under review. No calls and no
    cohorts have been run on it, and the frozen production width and defaults are unchanged.
  - The states it drops are not duplicates. States sharing a prefix pair are different latent template
    histories spelling the same allele prefix; the current template is what the next block's transition
    reads, so keeping only the best attribution approximates that prefix's future rather than
    deduplicating it. That is the approximation the grouped-prefix design below removes.
- **Separation, not proposal, is the current difficulty**: plans that complete separate by gaps of 0.401
  and 0.358; three of five tie outright, and their tied candidates agree on all 11 block genotypes.
- **Block consensus works, via the audit, and still does not reach C4.** Demonstrated end to end: an
  exhaustive search finding 26 mosaics against a K of 16 ties on two phasings of one genotype, the audit
  scores all 36 candidates exactly (10 panel + 26 mosaic), and the result is `BLOCK_CONSENSUS` with 5 of
  5 blocks called, phase `UNRESOLVED` and no phased haplotype. `candidate_set_conditional` stays `true`
  (the prefix was truncated) while `final_candidate_scope` is `COMPLETE` (the audited set decided).
- **Why C4 is still out of reach, corrected.** The obstacle is not truncation to K. Every C4 plan prunes
  29-43 million states, so `proposal_search_exhaustive` is `false`, the audit's precondition fails, and
  no audit is attempted on any of the ten runs. The three ties remain `RESCORE_INCOMPLETE` with no calls,
  and their tied candidates are still unanimous at 11 of 11 blocks. Reaching them needs an exhaustive
  search or a different sufficiency argument -- not a wider K and not a relaxed gate.
- Two independent 2x draws of the same ACOT donor disagree at one block, and the conditional GQ flagged
  exactly that block (3.04 against thousands elsewhere). Development evidence only, not calibration.

## Trusted performance results

- Deterministic parallelism: placement 1.99x and 1.98x, totals 1.94x (C4 audit 256) and 1.52x (ACOT 2x)
  at two threads; every statistical output byte-identical across thread counts.
- Process-wide operational caps: a refusal's overshoot is bounded by `threads * piece_hit_publish_batch`.
  Measured at production scale: 2,000,003,518 piece hits against a 2e9 cap, overshoot 3,518 of 8,192.
- ACOT 30x does not run: the complete workload is 2,255,437,567 piece hits, above the cap. Work scales
  linearly with coverage.
- Read-form deduplication was measured and **rejected**: predicted piece-hit reduction 0.55% / 3.22% /
  9.24% at 2x/10x/30x, and even a perfect cache leaves 2,046,971,082 hits, still above the cap.

## The grouped-prefix proposer

Implemented **alongside** `propose_mosaic_diplotypes`, which is untouched and serves as the oracle.
Nothing here is wired into calling; prefix-first (`954af2f`) remains experimental, and no calls or
cohorts have been run on either.

**What a group carries.** Within a prefix pair every history spells the same alleles, so the diploid
emission is one scalar for the group -- and it is the only coupling between the haplotypes. The joint
therefore factorises, `W_p(i, j) = E(p1, p2) * u_p1(i) * u_p2(j)`, and a group needs two length-`T`
vectors, each shared by every pair using that haplotype prefix, plus scalars.

**Cost model, corrected.** The earlier note said `O(2LAT)` for vector construction, from scanning the
templates once per child allele. That was wrong: each template spells exactly ONE allele at a block, so
a single scan over the `T` templates produces every child at once by scattering each template into its
own allele class. All `A` children of a prefix together occupy exactly `T` entries. Construction is
therefore **`O(2LT)` per retained pair-prefix boundary**, independent of `A`.

With C4's actual per-block allele counts -- 68, 68, 100, 16, 10, 8, 118, 119, 62, 5, 126 -- at `T` =
131 and `L` = 1,024:

| term | per block | total |
|---|---|---|
| template scans, `2LT` | 268,288 | 2,682,880 (blocks 1-10) |
| pair children, `L*A_b^2` | 25,600 to 16,257,024 | 64,382,976 |
| block 0 (`T` + `A_0^2`) | -- | 4,755 |
| **grouped total** | | **67,070,611** |
| existing proposer, `T^2 + 10*L*T^2` | 17,565,164 | **175,746,801** |

So 2.6x fewer operations, and the correction matters twice: construction drops by a factor of `A`, and
the **pair term now dominates** -- 64.4M of the 67.1M. The previous estimate used `A` about 20, which
badly understates the blocks carrying 100, 118, 119 and 126 alleles.

**Memory.** Vector cells held across a layer boundary are `4LT` (parent and child layers, at most two
haplotype prefixes per retained pair): 536,576 cells at width 1,024 and 2,146,304 at 4,096. At three
arrays of 8 bytes -- sum-product, max-product, and the integer history count -- that is 12.9 MB and
51.5 MB, against measured peaks of 810 MB and 2.0 GB for the current selector. The `L * A_b^2` children
of a layer are streamed into a bounded selector and never all materialised; at block 10 that would
otherwise be 16.3 million records.

**The unordered-prefix multiplicity is part of the selection score.** An unordered pair of DISTINCT
haplotype prefixes stands for two ordered histories and carries a factor of two; identical prefixes do
not. The factor is applied to the sum-product mass and to the history count at **every** layer, not at
the final collapse, so it participates in every layer's ranking. It is not applied to the max-product
score: the swapped history scores the same, so the maximum is unchanged. A mutation omitting it is
gated.

**Quantities carried, separately.**

| quantity | definition |
|---|---|
| sum-product | `E + log T(p1) + log T(p2) + het`, where `T(p) = sum_i u_p(i)` |
| max-product | `E + max_i m_p1(i) + max_j m_p2(j)`, no `het` |
| emission component | `E`, shared by every history in the group |
| transition component | the max-product's two maxima |
| history multiplicity | `N(p1) * N(p2)`, doubled when the prefixes differ, saturating at `2^64 - 1` |

The max-product cannot use the low-rank decomposition the sum-product uses: it needs the best OTHER
template, which comes from the parent's top two values excluding the template itself.

**Deterministic tie-breaking.** Layers rank by sum-product mass, then max-product score, then prefix
ids -- total, so no tie is left to arrival order. Candidates come back by mass, then witness, then
allele paths, which is the existing rule unchanged.

**How close the two implementations are.** Candidate sets, masses, best histories, emission and
transition components and history multiplicities agree to within 1e-9, and multiplicities agree exactly.
The residual is one unit in the last place: the two sum the same histories in different association
orders, so a pair of candidates that ties bitwise in one is separated by an ULP in the other, and the
documented allele-path tie-break then engages in only one of them. Measured on the two-allele fixture at
r = 0: the existing implementation separates `011/110` and `000/110` by a single ULP, the grouped one
computes them bitwise equal. Neither is wrong, so the gates compare candidates by identity rather than
by position, and assert cross-implementation order only between candidates the data actually separates.

**Numerical scaling.** Every quantity is held in log space and no linear-space mass is ever formed, so
nothing can overflow at C4's scale of about 1e6 nats. Each reduction factors out its own maximum before
summing (`log_add_exp`), which is the only scaling required; the max-product is a plain maximum in log
space, and multiplicities are integers with saturating add and multiply.

**Refusal is transactional.** Both budgets -- the work bound above and the `4LT` vector cells -- are
computed from the exact loops the search runs and checked BEFORE anything is allocated, so a refused
search allocates nothing and leaves no partial state. Both are integer and overflow-checked.

**Measured on C4, proposal only, width 1,024, five frozen plans.** Paired back to back on the same
binary against the in-tree proposer -- which is prefix-first (`954af2f`), not the original state-only
selector. No rescoring, no calling, no cohorts.

| plan | truth rank, prefix-first | truth rank, grouped | wall (s) | peak RSS (MB) |
|---|---|---|---|---|
| one_switch | 5 | 6 | 13.14 -> 2.56 | 772 -> 383 |
| bridgeable | 4 | 6 | 13.17 -> 2.59 | 776 -> 384 |
| unbridgeable | 46 | **1** | 13.80 -> 2.85 | 759 -> 382 |
| three_switches | 12 | **8** | 13.61 -> 2.56 | 771 -> 382 |
| every_block | 711 | **87** | 13.52 -> 2.59 | 777 -> 386 |

The truth is retained at every layer of all ten runs, and both proposers return 1,024 distinct
candidates. Operations are 175,745,801 for prefix-first against 64.9-65.1 million for grouped on every
plan -- 2.7x fewer, and within the 67,070,611 bound derived above, as a bound must be. Wall time is
5.1x to 5.3x lower and peak RSS about half.

The exact proposer is therefore cheaper than the approximate one AND ranks the truth better wherever
the case is hard: it recovers the two plans prefix-first damaged most, and gives up a little on the two
easy plans (rank 5 to 6, and 4 to 6). Nothing here is a calling result: ranks are proposal ranks, and
what a caller does with them has not been measured.

## Known limitations

- Whole-locus phase is often unidentifiable from short reads; block consensus is the response.
- `every_block`-style dense switching is not recovered at either width.
- Reverse-orientation splicing is unimplemented, so only 5 of 8 frozen C4 development donors can be used
  to build synthetic mosaics (C4 is 71 forward to 60 reverse).
- `graph_sort` assigns node ids depending on the path set: node ids are **not** stable across
  independently sorted graphs. Use one canonical sorted graph for model and truth; never compare node ids
  across graphs.
- Candidate sets on large loci are `CONDITIONAL`, which by contract blocks block consensus there.
- ACOT has no hash-drawn donor manifest; HG00097 is frozen by prior use, not by draw.

## Current blockers

- ACOT 30x cannot complete under the frozen cap without a lower-complexity placement algorithm.
- Dense-switch mosaics need beam diversity work (distinct diplotype prefixes rather than many
  template-level states representing one prefix), not a larger width. Now measured rather than
  suspected: the tightest layer holds 2-4 prefixes of 1,024 retained states, and quadrupling the width
  leaves it at 3-14 while the dominant prefix's share quadruples with it. Prefix-first selection now
  does this and retains the truth on all five plans at both widths -- for 3.4x to 5.3x the runtime, and
  with `unbridgeable`'s truth falling from rank 2 to 46. Whether that trade is worth taking is an open
  decision, not a settled improvement.
- **Block consensus cannot fire on any large locus, because the search is pruned.** The truncation half
  of `candidate_set_conditional` is now handled by the full-beam audit. What remains is the other half:
  C4 discards 29-43 million merged states per plan, so the search is not exhaustive, the audit refuses to
  run, and consensus is unreachable there. This is the same wall as the beam-diversity problem, seen from
  the calling side -- a pruned search cannot certify that it omitted nothing, whatever is done downstream.

## Ordered next work

1. Sol's review of this stop point.
2. Beam diversity: select distinct diplotype prefixes so one prefix cannot crowd the layer.
3. Decided: block consensus never applies to a K-truncated prefix; the full-beam audit supplies the exact
   scores instead. What remains open is the pruned-search case, which no downstream change can fix.
4. Lower-complexity placement, required before any 30x large-locus run or cohort.
5. Reverse-orientation splicing, to widen the synthetic-mosaic donor pool.

## Provenance

Last independently reviewed commit: `2a43fa7` (harness input binding, suite 19/19). The block-consensus
implementation and the hardened input-bundle binding are under review at this stop point.

Harness: `tests/genotype_c4_mosaic_experiment.sh` runs the five frozen plans end to end, binds every link
by content signature, and asserts the curated staged model reproduces ordinary C4 exactly
(`d2b8c6c2e0b4079d`, 8,515 diplotypes, 5/6 blocks).

## Do not redo

- Do not re-run the C4 leave-zero-out coverage ladder, the ACOT/GSTM1/ANKRD36C 2x pilots, or the C4
  cohort: their results stand and are archived in the checkpoint.
- Do not re-measure read-form deduplication; it was measured and rejected.
- Do not re-derive the CYP2D6 partial-template diagnosis; it is settled and curated.
- Do not raise `max_piece_hits` to make ACOT 30x green; the linear scaling is the problem.
- Do not change state limit 1,024 or audit width 4,096 without a review decision.
- Do not treat node ids as stable across independently sorted graphs.
