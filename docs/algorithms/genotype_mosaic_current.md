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
| `COMPLETE` | one diplotype certified above every rival beyond floating-point error | genotypes, block genotypes, scores, candidates, haplotype FASTA (see below) |
| `BLOCK_CONSENSUS` | scoring finished; no whole-locus winner separable, but tied candidates agree on some blocks | genotypes and block genotypes, with disagreeing blocks as explicit `NOCALL` rows; **no** haplotype FASTA |

The FASTA a `COMPLETE` run emits holds the two members of an **unordered genotype pair**, written in a
canonical order. `MosaicDiplotype` is a bare `{first, second}` with no ordering invariant, candidates
are enumerated as unordered pairs throughout, and different phased histories can spell the same
unordered genotype at every block. So the `hap1`/`hap2` labels assert a phase only where
`phase_support` is `READS_WITHIN_CANDIDATES` or `READS_AUDITED`; under `PRIOR_ONLY` or `UNRESOLVED`
they carry no read evidence about which haplotype is which. It is therefore **not** a "phased FASTA"
in general, and is not described as one anywhere in this document.
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

State limit **256** (production) and **1,024** (diagnostic comparison width) -- counted in allele-prefix
PAIRS since the grouped proposer became production; see the unit change noted below. The 1,024 figure is
the width the ladder compares against, NOT a phase audit: `phase_audit_candidates` remains 0, nothing in
this batch measured the phase audit, and every run declared `--phase-audit 0`. K (`max_mosaic_candidates`)
**128**. `max_piece_hits` 2e9, `max_physical_origins` 5e8,
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
- **Prefix-first selection recovers every truth, and costs more than it looks.** *Historical: measured
  before the grouped proposer existed, with the state limit still counting `(prefix, template)` states
  rather than allele-prefix pairs. Kept because it is what motivated the grouped design.* The layer
  selector grouped merged states by unordered prefix pair, ranks groups by aggregate mass, keeps each top group's
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
  - **Retired.** `954af2f` is superseded and no longer exists in the tree. The grouped-prefix proposer
    below carries every template attribution exactly instead of approximating a prefix by its best
    one, and is both cheaper and better ranked, so prefix-first was removed rather than kept as an
    option: its selector, its test file and its crowding fixture are gone. The record of what it
    measured is kept because it is what motivated the grouped design.
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

**This is what production proposes with.** `propose_mosaic_diplotypes`, the original state-based
implementation, is retained unchanged as the fixture oracle it is checked against and is no longer
called by the genotyper. Prefix-first is retired entirely.

**The state limit changed units, and this is the most consequential thing in this section.** It now
counts allele-prefix PAIRS, where it counted `(prefix, template)` states. One prefix pair stands for
every template attribution of itself, so the same number buys a far wider search: on the small e2e
locus, 64 went from pruning 264 states and finding 2 mosaics to searching exhaustively and finding 26.
The then-frozen default of 1,024 changed meaning in exactly that way. It is a unit change, not a tuning
change, and no default was re-tuned to absorb it. The width has since been re-frozen at 256 by the
declared ladder rule -- see Frozen defaults above.

**Storage rolls, and is gated as three separate accounts.** Only the parent layer's and the child
layer's numerical arrays are live; the parent's are released at the boundary. What persists is compact
`(parent, allele)` ancestry, compacted after every selection to the prefixes retained pairs actually
reference and remapped. The three accounts -- rolling vector cells, persistent ancestry records, pair
records and selection storage -- are predicted from the exact loops and gated separately, because
pooling them would let one overrun hide in another's slack. Each is admitted exactly at its bound and
refused one below it, before anything is allocated, and a measured-allocation gate checks what was
actually taken against what was predicted.

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

**The production resource contract.** One fixed ceiling, `max_proposal_payload_bytes` = 256 MiB, over the
proposer's predicted PAYLOAD -- the bytes its own arrays hold -- checked before any construction, so a
refusal allocates nothing and leaves no ranking and no call. It is deliberately not a claim about
allocated memory: container capacity slack, allocator rounding and per-vector headers are excluded, and
the name says so rather than the prediction pretending to an accuracy it does not have. The status reports the ceiling, the predicted total and the observed
total, plus each account separately, so a prediction is checkable from the run's own output instead of
being trusted.

Four accounts are predicted and kept as separate internal invariants, because pooling them would let one
overrun hide in another's slack: rolling vector cells (parent and child layers, the only ones live),
persistent ancestry records, pair records and selection storage, and -- added here -- every material
proposer-owned temporary and index: the allele-class index, the per-layer child index and its remap, the
child layer held before compaction, the per-parent scatter vectors, the dense scratch, and the spelled
candidate paths. Several of those scale with the layer width times the largest allele count, so omitting
them made the "peak" a partial account of what the search takes.

**What bounds a layer.** A prefix at block `b` is an allele PATH, so a layer holds at most
`min(2L, A_0 * ... * A_b)` haplotype prefixes, and the prediction sums that per block rather than
assuming the width everywhere. An earlier version of this bound argued instead that a prefix exists only
where some template spells it, hence at most `T` prefixes per layer. That is false for a MOSAIC proposer:
whenever `r > 0` the rank-1 jump term gives a child positive mass even when no single template spells its
whole allele path, which is the entire point of the search. The measured-allocation gate caught it --
14 ancestry records against a predicted 12, exactly `2 + 4 + 8` for that fixture's allele products -- and
that is why the gate compares measurement against prediction rather than checking the formula against
itself.

A layer's PRE-COMPACTION children are a different quantity from what it retains, and the larger one: the
child layer is built in full and only then selected and compacted, so
`children_b <= min(parents_{b-1} * A_b, A_0 * ... * A_b)`, with `child_index` literally a
`[parents][A_b]` table. Bounding children by the retained count instead passes every ordinary fixture
and understates the peak exactly when the parent layer is at its width and the next block has several
reachable alleles -- which is why that case now has its own gate rather than relying on loci that happen
not to saturate.

On the small end-to-end locus the contract predicts 9,126 bytes and the run takes 8,710: admitted exactly
at the prediction, refused one byte below it, with no partial output either way.

**Measured in production, C4, width 1,024, one shared scoring pass per plan at K = 0, 8, 128.** The
standalone and production proposals were required to agree, and do, on candidate identity, rank AND
score -- the last only after the candidate table was written at full double precision, since the stream
default of six significant digits had been rounding scores near 1.7e6 to the nearest ten.

| plan | truth proposal rank | exact-score rank, K=8 | K=128 | status K=8 / K=128 | blocks exact K=8 |
|---|---|---|---|---|---|
| one_switch | 6 | 2 | 2 | COMPLETE / COMPLETE | 11 of 11 |
| bridgeable | 6 | 2 | 2 | COMPLETE / COMPLETE | 11 of 11 |
| unbridgeable | 1 | 1 | 1 | RESCORE_INCOMPLETE, 2 tied | no call |
| three_switches | 8 | 1 | 1 | RESCORE_INCOMPLETE, 2 tied | no call |
| every_block | 87 | **absent** | 39 of 8,643 | RESCORE_INCOMPLETE, 8 then 16 tied | no call |

Block edit distance is 0 for the two plans that call (11 of 11 blocks exact, against 3 and 9 at K = 0),
and undefined for the three that tie, which emit no call. Every winner set is unanimous at 11 of 11
blocks, so all three ties would call correctly if consensus were permitted -- it is not, because the
search is still not exhaustive at this width and the candidate scope is `CONDITIONAL` on all fifteen
runs. Consensus eligibility is therefore unchanged by this work.

`every_block` is the case K = 128 exists to diagnose: at K = 8 its truth is not in the scored set at
all, since proposal rank 87 falls outside the prefix. K = 128 brings it in at exact rank 39 and the run
still ties 16 ways, so a wider prefix buys the truth a place in the scored set without resolving the
tie. **K = 128 is a development diagnostic, not a new default**, and the harness keeps 0,8 as its
committed default with the wider set available only as an override.

Cost: 49-62 s per plan for the whole shared pass, with the three declared prefixes together adding
about 0.05 s -- they are one scoring pass, not three. Peak RSS 712 MB across the five-plan run. The
model-equivalence gate still holds on all five plans (`d2b8c6c2e0b4079d`, 8,515 diplotypes, 2 curated),
and the simulated reads reproduce byte-for-byte, so the input binding's earlier refusal was a signature
effect -- the probe binary is part of the bundle by design -- and not drift in the data.

**The width ladder, and why 256.** Proposal only, grouped proposer, five frozen plans, existing reads,
at 256 / 512 / 1,024 / 2,048 / 4,096. The truth is retained at **every layer of all 25 runs**, and the
proposal rank is identical at every width: 6, 6, 1, 8, 87. The width is genuinely applied -- candidates
found track it exactly (256 to 4,096), predicted bytes scale 3.6 to 57.7 MB, wall time 4.7 to 6.6 s --
so the flatness is real, not an ignored flag. A rank is a position in a mass-ordered list, and widening
appends lower-mass candidates below the truth without displacing it.

So width is not the binding constraint on these plans, and 256 is the smallest plausible width, leaving
about 3x headroom over the worst truth rank of 87. **This is five synthetic plans at one locus with one
donor pair. It is not evidence that 256 suffices anywhere else**, and dense-switch recall was exactly
where retention used to fail before the grouped proposer.

**The K ladder, and the declared rule.** At width 256, one shared exact-scoring pass per plan through
K = 256, deriving K = 8, 32, 64, 128, 256 as prefixes of that one pass; repeated at 512 for stability.

| K | truth exactly scored | calls / tied-winner consensus == K=256 | stable at 512 |
|---|---|---|---|
| 8 | **no** | yes | yes |
| 32 | **no** | yes | yes |
| 64 | **no** | yes | yes |
| 128 | yes | yes | yes |
| 256 | yes | reference | yes |

**K = 128** is the smallest K satisfying every clause. The rule is decided entirely by the first one:
`every_block`'s truth sits at proposal rank 87, so K = 8, 32 and 64 never score it at all. The second
clause never discriminates here -- block calls and tied-winner consensus are identical to K = 256 at
every K and both widths -- and the third never discriminates either. That is worth recording plainly:
one clause did the work, and a rule that looked like three tests was effectively one.

`every_block`'s winner set grows from 8 to 16 tied candidates between K = 8 and K = 32 and then holds,
unanimous at 11 of 11 blocks throughout. Scope is `CONDITIONAL` and `proposal_search_exhaustive` false
on all 50 runs, so **consensus stays refused everywhere**: the strict policy is unchanged by this batch.

Cost at the frozen settings: 62-84 s per plan for the whole shared pass, peak RSS 0.70-0.88 GB, of which
the proposer's own predicted peak is 3.6 MB at width 256 -- 1.4% of the 256 MiB ceiling.

## C4 leave-zero-out result package (2026-09-21)

Attribution: run commit `c8fe2a4`, caller binary `5c3a1525…` (byte-identical across the whole batch, so
no result moved when the harness was corrected), protocol frozen at `5bb9330` -- manifest, draw rule,
orientation and panel-allele tables and the graph are all bit-identical to that freeze, and only the
harness changed. The post-hoc evaluator is a separate read-only binary `cda6273c…`; the frozen run's
binary was never rebuilt, and all 40 experiment and 40 call signatures re-verified against it.

**Two experiments, reported separately and never pooled**, because their truths are different objects.

**Ordinary LZO at 30x** -- canonical source `tests/results/genotype_cohort_c4_r6_5937384` (round 6,
commit `5937384`). **In each arm separately**, panel and mosaic: 56/56 donors called, zero no-calls,
**56/56 fully exact**, **616/616 exact blocks**, **280/280 exact SV blocks**, and **total distance
zero** -- SV and context distance both zero. Phase in each arm: `PRIOR_ONLY` 51, `READS_AUDITED` 5.
The paired comparison in that cohort shows no winner changed and no block call changed between arms.

**Ordinary LZO across nested coverages** -- canonical source
`tests/results/genotype_coverage_c4_held_out_abb35a0` (commit `abb35a0`). Identical in both arms, and
**zero no-calls at every coverage**:

| coverage | donors fully exact | exact blocks | exact SV blocks | total distance |
|---|---|---|---|---|
| 2x | 55/56 | 615/616 | 280/280 | 1 |
| 5x | 55/56 | 615/616 | 280/280 | 1 |
| 10x | 56/56 | 616/616 | 280/280 | 0 |
| 15x | 56/56 | 616/616 | 280/280 | 0 |
| 20x | 56/56 | 616/616 | 280/280 | 0 |
| 30x | 56/56 | 616/616 | 280/280 | 0 |

Every figure above is recomputed from the per-donor rows and cross-checked against each source's own
published summary; a disagreement is a hard failure, not a warning.

**The round-5 cohort `b1684ae` is SUPERSEDED** and is excluded from every current scientific summary.
It reported 55/56 completed with one `CANDIDATE_INCOMPLETE` and a total distance of 120; those numbers
are history and must not be quoted as a current result. The directory is retained, marked with a
`SUPERSEDED.txt`, and `tests/genotype_lzo_consolidation.sh` fails if it is used as a source or cited in
a consolidated package.

**Synthetic-mosaic LZO** (10 plans x 4 coverages = 40 shared runs, 80 arm records). Arms never pooled:

| | K=0 panel arm | K=128 mosaic arm |
|---|---|---|
| completed / intended | 40/40 (100%) | 24/40 (60%) |
| no-call | 0 | 16 (40%) |
| exact over intended | 8/40 (20%) | 21/40 (52.5%) |
| exact given a call | 8/40 (20%) | 21/24 (87.5%) |
| excluding panel controls | 0/32 exact | 13/32 (40.6%), 13/16 given a call (81.2%) |
| mean exact blocks | 7.95/11 | 10.83/11 |
| total distance | 54,636 over 40 calls / 440 blocks | 111 over 24 calls / 264 blocks |

Those two distance sums have **different denominators** and must never be compared without them. The
K=0 exact calls are exactly the two panel controls, which is what a panel-only arm should recover.
Paired within plan and coverage: 13 improved, 9 unchanged, 2 worsened, 16 refused at K=128.

**Error severity.** Only three K=128 strict calls are inexact, all `one_homologue_once_1`. At 10x and
5x a single SV block (`sv:3`) is wrong by 49 edits against 488 bp of truth (normalized 0.100) while
READS-supported at conditional GQ 7,328 and 2,523 -- a confidently wrong block. At 2x two blocks are
wrong (`ctx:0` 7 edits, `sv:1` 6) under `PRIOR_ONLY` support at GQ ~67. The 10x and 5x cases are
regressions against K=0 on the same reads (distance 13 -> 49) even though they gain an exact block.

**Conditional consensus diagnostics** (`CONDITIONAL_CONSENSUS_DIAGNOSTIC`, **not calls**, excluded from
every strict rate): all 16 refused arms are unanimous on 11/11 blocks, so the ties are over
representation rather than block genotypes. Nine have consensus distance 0; the true pair is in the
tied set in eight. Winner sets are 2-16. Worst consensus is 26,379 (`three_switches_2` at 2x).

**Phase.** K=0: 40/40 `UNRESOLVED`. K=128: 3 `READS_WITHIN_CANDIDATES`, 21 `PRIOR_ONLY`, 16
`UNRESOLVED` (the refusals). No `READS_AUDITED` anywhere, because the phase audit is 0 by declaration.
`PRIOR_ONLY` is **not** read-supported phase. A haplotype FASTA is emitted for every COMPLETE run,
including all 40 K=0 runs whose phase is `UNRESOLVED`; those two records are an unordered genotype pair
in canonical order, and calling them a phased haplotype pair would assert read evidence that does not
exist. The package column is `haplotype_fasta_emitted` for that reason.

**Runtime**: 40 shared runs, 958.6 s total, peak RSS 0.65 GB (30x mean 50.5 s, 2x mean 9.7 s).

**Tables and plots**: `tests/results/genotype_c4_lzo_c8fe2a4/` -- `runs.tsv` (80 arm records),
`blocks.tsv` (880 block rows, strict calls and conditional consensus labelled separately), `tied.tsv`,
`phase.tsv`, `denominators.tsv`, `ordinary_lzo.tsv`, `ordinary_lzo_coverage.tsv`, `provenance.tsv`,
and seven `plot_*.png` each
beside the `plot_*.tsv` it was drawn from. Like every other directory under `tests/results/`, it is a
local artifact and is not tracked in git. Regenerate with `scripts/consolidate_c4_lzo.py`.

**Limitations, and they bound the claim:**

- these are unseen RECOMBINATIONS of previously seen C4 sequences, not unseen haplotypes: every donor
  was already an ordinary leave-zero-out donor;
- the search remains pruned and conditional -- `candidate_scope` is `CONDITIONAL` and
  `proposal_search_exhaustive` is false in all 80 arm records;
- the consensus diagnostics are not certified calls and may not be counted as calls;
- conditional GQ is conditional on the scored candidate set and is **not calibrated**; K=0 and K=128
  GQs are not comparable and are never pooled;
- this does not establish LOO performance;
- it does not validate any other locus.

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
  leaves it at 3-14 while the dominant prefix's share quadruples with it. **Resolved, and no longer a
  blocker**: prefix-first was the first answer to this and is retired; the grouped proposer carries every
  template attribution exactly, is cheaper, and retains every truth at every layer down to width 256.
- **Block consensus cannot fire on any large locus, because the search is pruned.** The truncation half
  of `candidate_set_conditional` is now handled by the full-beam audit. What remains is the other half:
  C4 discards 29-43 million merged states per plan (measured in the older `(prefix, template)` unit), so
  the search is not exhaustive, the audit refuses to
  run, and consensus is unreachable there. This is the same wall as the beam-diversity problem, seen from
  the calling side -- a pruned search cannot certify that it omitted nothing, whatever is done downstream.

## Ordered next work

1. Sol's review of this stop point.
2. **Resolved, not pending: beam crowding.** The grouped-prefix proposer made the allele-prefix PAIR the
   unit of retention, so a single diplotype prefix can no longer occupy a layer through many
   template-level states that all spell it. It carries every template attribution exactly, costs less,
   and retained every development truth at every layer down to width 256. Distinct-prefix (prefix-first)
   selection was the first answer to the same problem and is retired. Nothing about beam diversity
   remains open; what remains open is the pruned search, item 4.
3. Held-out synthetic-mosaic validation on C4. The protocol is frozen before any of it runs:
   `tests/genotype_c4_heldout_mosaic_manifest.tsv` fixes the donor draw, the mosaic classes, the
   simulation settings and the evaluation metrics, and `tests/genotype_c4_heldout_mosaic.sh` drives it.
   **Nothing has been simulated, called or evaluated under it.**
4. Decided: block consensus never applies to a K-truncated prefix; the full-beam audit supplies the exact
   scores instead. What remains open is the pruned-search case, which no downstream change can fix.
5. Lower-complexity placement, required before any 30x large-locus run or cohort.
6. Reverse-orientation splicing, to widen the synthetic-mosaic donor pool. Until it exists the held-out
   draw can only consider donors that are forward on both haplotypes, which is what bounds its pool.

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
- Do not change state limit 256, K 128, or the 1,024 diagnostic comparison width without a review
  decision. The 1,024 and 4,096 figures in the historical sections below are the OLD frozen width and
  its old comparison width; they are records of past measurements, not current settings.
- Do not treat node ids as stable across independently sorted graphs.
