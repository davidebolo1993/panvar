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
  leaves it at 3-14 while the dominant prefix's share quadruples with it.
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
