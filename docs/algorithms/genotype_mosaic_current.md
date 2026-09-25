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
| `MODEL_INCOMPLETE` | the structural model could not be built | status only |
| `CANDIDATE_INCOMPLETE` | proposal stopped at a declared budget | status only |
| `SPELLING_INCOMPLETE` | a candidate could not be spelled | status only |
| `RESCORE_INCOMPLETE` | a resource cap was reached, or a tie with no block agreement | status only |
| `CONDITIONAL_GENOTYPE_EQUIVALENT` | adaptive completion: the pooled top set is several candidates tied on total score that carry one unordered genotype at every block | genotypes, block genotypes, scores, candidates, equivalence set, repair trajectory; phase `UNRESOLVED`; **no** haplotype FASTA |
| `REPAIR_INCOMPLETE` | adaptive completion: the repair did not converge (round cap, a tie among improving neighbours, a cycle), the anchor did not reproduce, or a candidate scored twice did not score identically | status and repair trajectory only |
| `RUNNING` | adaptive completion in progress | status; committed rounds in `.repair_rounds.partial.tsv` |
| `TIME_LIMIT` / `RESOURCE_REFUSAL` | the watchdog stopped the process at the declared wall or physical-footprint limit | status and the committed partial rounds only |

The FASTA a `COMPLETE` run emits holds the two members of an **unordered genotype pair**, written in a
canonical order. `MosaicDiplotype` is a bare `{first, second}` with no ordering invariant, candidates
are enumerated as unordered pairs throughout, and different phased histories can spell the same
unordered genotype at every block. So the `hap1`/`hap2` labels assert a phase only where
`phase_support` is `READS_WITHIN_CANDIDATES` or `READS_AUDITED`; under `PRIOR_ONLY` or `UNRESOLVED`
they carry no read evidence about which haplotype is which. It is therefore **not** a "phased FASTA"
in general, and is not described as one anywhere in this document.

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

## Production algorithm: panel-anchored adaptive completion (2026-09-23)

The production caller is one adaptive algorithm. It is never told whether a sample is LZO, LOO or a
mosaic -- those are evaluation scenarios, which decide only the model graph a campaign hands it and the
truth an evaluation compares against afterwards. The documented command:

```text
panvar genotype --profile profiles/genotype_c4.tsv --gfa <graph.gfa> --reference-path <name> \
                --reads1 <r1.fq> --reads2 <r2.fq> --library <params.txt> --out-prefix <prefix> \
                [--panel-curation <manifest.tsv>]
```

`profiles/genotype_c4.tsv` is the committed, frozen C4 profile (`c4-panel-anchored-v1`): proposal width
256, K = 128, phase audit 0, 8 threads; the union's per-pass caps (7e8 physical origins, 6e8 state
classes, 512 MiB emission cache) from the development LOO campaign; repair with 32 local pairs, at most
11 rounds, the 384 MiB planned origin cache, a 670-candidate neighbourhood bound and cumulative caps of
4.9e9 physical origins, 4.2e9 state classes, 1.4e10 piece hits and 1.4e10 latent pairs; and a 2.5 GiB
physical-footprint limit. Unknown, repeated or malformed profile settings are refused. The status records
the profile path, name and SHA-256, every setting a flag overrode after it, and the SHA-256 of every input
(graph, reads, library, curation). Without a profile the unrepaired call from the proposal's candidate set
is made: a development view, not the production caller.

**Algorithm.**

1. Exactly score the union: every unordered panel diplotype and up to K grouped-prefix mosaics.
2. Anchor on the best exactly scored **panel** diplotype -- the first sequence candidate of the union's
   exact ranking with a panel member, its highest-path-prior panel member. A deterministic search anchor,
   not a claim about the sample.
3. Run the exact one-block coordinate repair from the anchor: each round scores every neighbour that
   replaces one block of one haplotype by one of the block's 32 locally best allele pairs (both
   assignments), and moves to the best one when it improves beyond floating-point error; visited
   incumbents are never revisited. Deterministic, 8-thread fragment-parallel window enumeration, identical
   at any thread count. The implementation is `CoordinateRepairSearch` in `genotype_repair.cpp`, the same
   code the developer probe's `repair-search` drives.
4. Pool the union and every candidate the repair scored. A candidate scored by both must score
   identically, bit for bit (`pool_score_mismatches` 0), or the run refuses.
5. Rerank the pool exactly at the sequence-pair level and decide over its **numerical top set** -- the
   candidates that cannot be separated from the best total score beyond floating-point error:
   - one candidate: the normal genotype decision, status `COMPLETE`;
   - several, all with the same unordered allele pair at every block: `CONDITIONAL_GENOTYPE_EQUIVALENT`,
     the unphased block genotype is emitted, phase is `UNRESOLVED`, `.equivalence_set.tsv` lists every
     member and no haplotype FASTA is written;
   - several that disagree at any block: refused (`RESCORE_INCOMPLETE`), the disagreeing blocks named.
   Lower-scoring candidates of the winner's genotype never turn a unique winner into a tie. No majority
   vote and no consensus across candidates that disagree anywhere.

A search that does not converge -- round cap, a tie among improving neighbours, a cycle -- emits no
genotype (`REPAIR_INCOMPLETE`), and neither does a cumulative-cap refusal (`RESCORE_INCOMPLETE`) or a
watchdog stop (`TIME_LIMIT`, `RESOURCE_REFUSAL`). A round tie that includes the incumbent is converged:
no neighbour improves on it.

**Phase support is a separate axis from the decision.** `READS_WITHIN_CANDIDATES` when no pooled rival is a
fragment contender; `PRIOR_ONLY` when some rival's fragment term is not below the winner's beyond error
and only the prior chose the unique winner; `UNRESOLVED` when total scores stay tied within the top set.
A unique `PRIOR_ONLY` winner does emit its haplotype FASTA, but as a **prior-selected haplotype
representation**, not a read-phased result: the FASTA header and the status carry
`phase_support=PRIOR_ONLY` and `representation=prior_selected` (`haplotype_representation`). The status
also reports the best same-genotype alternative below the top set -- its score, fragment and prior gaps,
and whether its fragment term is identical -- which is the evidence for `PRIOR_ONLY`, never a reason to
suppress the genotype.

**Scope and confidence.** The pool is a local search from one anchor: `candidate_scope` and
`final_candidate_scope` stay `CONDITIONAL`, `global_certified` is false and
`search_globally_exhaustive` false. Separation is reported against the best candidate carrying a
**different** genotype (`genotype_separation_gap`), beside the raw `numeric_separation_gap`.
`conditional_GQ` is conditional on the pool and labelled `uncalibrated`. The multi-start diagnostic
(9a0f4c0) found, for NA20509 30x, a second phase history at the same score outside this restricted pool;
it stays diagnostic evidence and is not added to production.

**Outputs and process contract.** Beside the standard tables: `.equivalence_set.tsv` (the complete top
set), `.repair_rounds.tsv` (the full trajectory: moves, stopping reason and every work, cache and
resource column), `.repair_memory.tsv` (per-round memory accounts at four points). While running the
status says `RUNNING` and each committed round is atomically rewritten to `.repair_rounds.partial.tsv`.
A watchdog enforces the declared wall and footprint limits: a stop removes every final output and the
staged graph, keeps the committed partial rounds, and writes `TIME_LIMIT` or `RESOURCE_REFUSAL`. Once the
run has decided, it claims itself and can no longer be stopped half-written. The staged sorted graph is
removed on every exit. The status reports `threads_peak_active_workers`,
`resources_phys_footprint_peak_bytes` and `resources_rss_peak_bytes`. Exit status is 0 when a genotype
was emitted (`COMPLETE`, `CONDITIONAL_GENOTYPE_EQUIVALENT`), 1 otherwise.

**Regressions** (`tests/genotype_production_c4_regressions.sh`, the real command on five frozen cases;
expected values from the multi-start evidence restricted to the production pool, trajectories bit-identical
to the recorded panel-anchor branch in `tests/genotype_production_c4_expected/`): NA20509 LOO 10x
(-151090.979646, distance 114, floor 91), NA20509 LOO 30x (-414918.317265, distance 113, floor 91, a
unique winner 4.40 nats from the best different genotype, `PRIOR_ONLY`: the best same-genotype alternative
has an identical fragment term and is 7.66 nats lower through the prior), NA19240 LOO 30x (distance 4 =
floor 4, with the ctx:6 repair), NA20509 LZO 30x (exact, no repair movement) and the bridgeable synthetic
mosaic (11/11 exact).

**Campaigns.** `tests/genotype_production_c4_campaign.sh plan` prepares the final held-out LZO (56
donors), eligible-donor LOO (55) and synthetic-mosaic (10 plans) campaigns at 2x, 5x, 10x and 30x -- 484
jobs -- with a (historical) runtime and disk estimate. `launch` refuses without the reviewed plan's hash. The
final campaign's binding is in "Final campaign binding" below.

## v2: exact branch-and-pool at an improving-move tie (2026-09-24)

C4 production-v1 (`c4-panel-anchored-v1`, the 484-job campaign) stays the immutable result: its profiles are
unchanged, and without `repair_tie_branching` a run behaves exactly as v1 (an improving-move tie is the
no-call `improving_move_tie`). v2 adds one step and is selected only by its own profiles,
`profiles/genotype_c4_v2.tsv` (`c4-panel-anchored-v2`) and `profiles/genotype_c4_v2_campaign_local.tsv`:
every v1 setting unchanged, plus `repair_tie_branching 1` and `repair_tie_max_branches 16`. The likelihood,
the candidate union, the repair and the pooled rule are unchanged.

**Step 3b.** When the repair stops at a round whose numerical top set excludes the incumbent, every tied
neighbour becomes a branch that continues the same exact repair from it, and a branch that stops at a tie of
its own splits again. Branches run breadth-first, within a tie in diplotype-identity order
(`tie_branch_order bfs_canonical_identity_v1`), which no input enumeration or thread schedule can change.
They share the production search's scoring context -- warm origin and emission caches. A waiting branch is a
compact descriptor only: its allele paths, the declared score components and prior, the sequence digests, its
parent tie and history identities (`tie_branch_seeding census_free_just_in_time_shared_context`). Immediately
before a branch runs, its two haplotypes are spelled from their allele paths, placed in full by the ordinary
exact path (into released slots) with their censuses retained, and the declared complete score -- every
component and the prior -- is reproduced bit for bit before any neighbour is scored; a seed that does not
spell its digests, cannot be placed, or does not reproduce fails the branching closed. Sequence digests are
memoised by allele path. Each branch has its own visited set (its ancestors' incumbents included), ranking
state and candidates. Every candidate any branch scores joins the pool; the pooled rule decides, never the
order branches ran in. Branching runs before the run is claimed, so the watchdog still bounds it.

**Branch memory.** A settle point follows the root tie's recording and every fully processed branch (its result
copied into the pool and, for a split, its children recorded). At each one the finished branch's search and
every repair slot are released; the content-keyed origin and emission caches are kept. So at every settle
point the context holds no slot, sequence, mass column or census for any branch, the frontier grows only in
descriptors, and at most the running branch's start holds a seed census; every exit (converged, capped,
refused, failed) closes with the same release. Seed placement is search work: it counts against the
cumulative caps with everything else. The frontier used before (a census-holding seed per waiting branch,
built by block delta at its tie and pinned until its branch ran) grew by one census per split and passed
4 GiB on `loo.HG01573.30`; it survives only as a gate hook, as does the older lifetime that also pinned every
parent incumbent.

**Calling contract.**

| pooled outcome | status | call class | no-call reason | FASTA |
|---|---|---|---|---|
| a unique top candidate | `COMPLETE` | `strict_call` | `none` | yes (its phase support as in v1) |
| one unordered genotype in the top set | `CONDITIONAL_GENOTYPE_EQUIVALENT` | `genotype_equivalent_call` | `none` | no; phase `UNRESOLVED` |
| different unordered genotypes | `RESCORE_INCOMPLETE` | `no_call` | `divergent_tie` | no |
| a branch cap, any cumulative cap, an exception, an unfinished branch, a split without an improving tie | `REPAIR_INCOMPLETE` | `no_call` | `tie_branching_incomplete` | no |

Every tie-branched run carries `search_conditional tie_branched`: the call is conditional on a pruned,
anchor-local search that was branched at a tie. Branch recovery improves callability; it does not
guarantee proximity to truth (in the LOO audit 18 of 37 recovered ties resolve to one genotype more than
1 kb above the floor).

**Resources.** Branches share the production context, so the cumulative repair caps
(`repair_max_physical_origins` and the others) bind production and branches together, and the status
margins are reported against that shared total. Because the caches are shared, executed work depends on the
order branches ran in; resource acceptance is defined on the canonical order. Gates that reverse the order
compare scientific output with cap headroom, never a boundary cap.

**Outputs.** `.tie_terminal_ties.tsv` (every tie: its tied neighbours, block, assignment and score
components) and `.tie_branches.tsv` (every branch: start, rounds, moves, stop, outcome, terminal, executed
physical origins, seconds). Branch-scored candidates appear in `.mosaic_candidates.tsv` and
`.equivalence_set.tsv` with `source branch` and `scored_in branch<k>_round<r>`. Status keys:
`search_conditional`, `repair_tie_branching`, `repair_tie_max_branches`, `tie_branch_order`,
`tie_branch_seeding`, `tie_branched`, `tie_ties`, `tie_branches_run`, `tie_distinct_terminal_genotypes`,
`tie_incomplete_reason`, `tie_branch_physical_origins`, `pool_branch_scored`, `tie_digest_memo_hits`,
`tie_digest_memo_misses`, `seconds_tie_branching`, and the seed accounting `tie_seeds_materialised`,
`tie_seed_haplotypes_placed`, `tie_seed_physical_origins`, `tie_seed_state_classes`, `tie_seed_piece_hits`,
`tie_seed_latent_pairs`, `seconds_tie_seed_placement`, `tie_seed_peak_live_censuses`,
`tie_peak_pending_descriptors`. `.tie_branch_memory.tsv` (`.partial.tsv` while branching runs, kept by a
refused run) has one row per settle point, per materialised seed and at close: slots, sequence, mass-column
and census bytes, origin-cache payload, footprint, pending descriptors and cumulative physical origins.

**Gates.** Engine fixtures in the command end-to-end test (unique convergence where the winner is not the
first branch, genotype-equivalent and divergent terminals, nested ties, the branch cap exact and one below,
a throwing branch, reversed order, a split without a tie); the block-delta fixture (a delta-promoted seed and
the repair from its census equal cold full seeding in every mass cell and in exposure); and
`tests/genotype_tie_branching_regressions.sh`, five frozen C4 ties (two synthetic mosaics, a low- and a
high-excess LOO tie, the 14-branch nested case) whose production outputs must equal their frozen v2/post-hoc
diagnostics: status, pool, pooled top set member for member and score for score, representative, branch
endpoints, terminal ties and the unchanged production trajectory. `tests/genotype_branch_lifetime.cpp`
replays a scripted nested branching under the census-free, eager and parent-retaining lifetimes: every seed,
neighbour and incumbent bit-identical to a fresh placement and between the census-free and eager frontiers;
at every census-free settle point no census, live slot, or sequence or mass-column growth, however large the
frontier; released slots reused.

## Representation-conflict quality guard (2026-09-24)

An annotation of calls, never a calling change: the genotype, its score, the status and call class, phase and
any branch trajectory are exactly what they are with the guard off, and a no-call is never relabelled. It is
off by default and bound explicitly by the v3 profiles (`profiles/genotype_c4_v3.tsv`,
`profiles/genotype_c4_v3_campaign_local.tsv`: each its v2 profile plus exactly `representation_guard 1`), so
the immutable v1 and v2 profiles run as they were accepted (the final campaign's binding is below). On every call of
the production (repair) path, `representation_conflict_rule_v1` examines each called block:

1. the marker (k-mer depth) channel's best available allele pair at the block -- emission descending, the
   first pair in index order on a tie -- and its margin over the called pair, `representation_marker_margin`;
2. when that margin reaches ln(1000) nats, the competitor: the called diplotype with that pair substituted at
   the block, in both assignments, scored exactly (from the pool when it is there, else in the live scoring
   context, outside the pool) and the better by total (ties: fragment term, then exposure term);
3. the block is `LOWQUAL_REPRESENTATION_CONFLICT` when min(marker margin, the fragment term's preference for
   the call, the exposure term's preference for the competitor) >= ln(1000) nats: two depth readings -- the
   marker channel and the likelihood's own exposure term -- against the call, the sequence-placement term
   alone for it. Otherwise `PASS`; a competitor that cannot be scored (a cap, an exception) gives
   `NOT_EVALUATED`, never `PASS`.

**Fail closed.** The locus filter is `LOWQUAL_REPRESENTATION_CONFLICT` if any block is flagged; otherwise
`NOT_EVALUATED` if any block did not complete the rule; `PASS` only when every block completed it and passed.
`representation_guard_complete` and `not_evaluated_blocks` keep both facts when a flag and an unevaluated
block coexist. The guard's own scoring is reported apart (`representation_guard_physical_origins`,
`representation_competitors_scored`) and the cumulative repair counters and margins are refreshed after it,
so the caps bind it with everything else.

Truth-free and locus-free: no truth, availability floor, donor or block identity enters, and the marker never
enters the ranking. **ln(1000) is a declared conservative threshold for a heuristic conflict statistic, not a
calibrated Bayes factor**: the marker depth and the exposure term are correlated evidence derived from the same
reads, and nothing here is calibrated across loci. What was observed is a broad separation on the C4
diagnostic set (below), not independence or calibration.

**Measured before integration** (`tests/results/.../v2_representation_guard/`, a ten-cell set; the truth-free
metrics table and its flags were frozen by SHA-256 before any truth was joined): the rule flags exactly the four
HG01114 sv:7 blocks (2x, and the tie-branched 5x, 10x and 30x calls; statistics 163.7, 409.3, 792.2 and 2455.3
nats; each block about 32.6 kb above the floor) and no block of the controls -- two exact LZO calls, an exact
synthetic mosaic, the floor-attaining HG00514.10. The largest unflagged statistic is 1.46 nats, so the flags are
identical for any bound in (1.46, 163.7]. The marker channel alone disagrees with exact calls by up to 63 nats
(the mosaic's sv:3), which is why a marker margin alone is not a rule. Not flagged, as the definition intends:
HG00512.2 sv:1 (2,034 bp; the marker agrees with the call, the prior decides) and HG00171.2 ctx:10 (13 bp; a
sequence-only projection difference with no exposure-scale difference).

**Representation-PASS is not high quality.** A PASS block has no certified depth-versus-sequence conflict; it
can still be wrong. HG00512.2 sv:1 stays a 2,034 bp error under PASS, because the marker channel agrees with
the call and the prior decides it at 2x -- expected under this rule, which is not a general error detector.
Downstream therefore keeps three block tables: `blocks_all_quality.tsv` (every called block),
`blocks_pass_only.tsv` (filter exactly `PASS`) and `blocks_high_quality.tsv` (`PASS` and genotype support
`READS`), each with its denominator and retained fraction. The final "high-quality" plot uses the third;
representation-PASS stays separately available. None of the three is calibrated accuracy.

**Outputs.** `.genotypes.tsv` and `.block_genotypes.tsv` gain `filter`, `filter_reason`,
`representation_marker_margin`, `representation_competitor_pair`, `representation_fragment_diff`,
`representation_exposure_diff`, `representation_likelihood_margin`, `representation_conflict_statistic` (on
the development paths, `NOT_EVALUATED`). The status gains `representation_guard`, `representation_guard_rule`,
`representation_guard_bound_nats`, `locus_filter` (`PASS`, `LOWQUAL_REPRESENTATION_CONFLICT`, `NOT_EVALUATED`;
`NA` for a no-call), `lowqual_blocks`, `not_evaluated_blocks`, `representation_guard_complete`,
`representation_competitors_scored`, `representation_guard_physical_origins` and
`seconds_representation_guard`.
An emitted FASTA header carries `locus_filter=`. A LowQual call stays in every all-quality output.

**Downstream.** The campaign classifier decides a call by status and call class, never `repair_converged`
(false for a tie-branched call), and treats a no-call carrying a locus filter as a defect. The evaluator
writes each block's `filter` into `.eval.blocks.tsv` and adds `locus_filter`, `lowqual_blocks`,
`representation_guard_complete`, `not_evaluated_blocks`, the PASS-only totals (`pass_blocks`,
`pass_blocks_exact`, `pass_unordered_distance_total`, `pass_availability_floor_total`,
`pass_excess_over_floor_total`: filter exactly PASS) and the high-quality totals (`high_quality_*`: PASS and
READS) beside the unchanged all-quality totals. The consolidator writes the three block tables and reports
`lowqual_calls`, `not_evaluated_calls` and their rates, block filter and NOT_EVALUATED rates, the three
tables' denominators and retained fractions, and the PASS-only and high-quality excess -- all beside
callability and the no-call rate, which no filter changes.

**Gates.** e2e fixtures (the conflict shape flagged, a marker-only disagreement not; the bound exact and one
ulp below; assignment order and homologue swap invariant; canonical marker ties);
`tests/genotype_representation_guard_gate.sh` (three frozen cells through the production path: guard on and
off give identical pre-filter outputs, the v1-call cells equal their frozen campaign tables, the flags and
evidence equal the frozen diagnostic, 1/4/8 threads identical, and only graph, reads and library are read);
`tests/genotype_downstream_schema_gate.sh` (CTest `panvar_genotype_downstream_schema`: NOT_EVALUATED and NA
never PASS, high quality needs READS, every rate apart); the outside-pool fixture in the command e2e (gate hooks
force every competitor into the live context: call and ranking unchanged, evidence equal to the pool's, guard
work reconciled exactly with the cumulative counters, the cumulative cap exact and one below -- the refused
block NOT_EVALUATED beside the flagged ones, the genotype untouched -- and 1/4/8 threads); and the campaign
gate's profile checks (v3 = v2 plus exactly the guard; v1 and v2 never set it).

## v4: repair round limit 32 (2026-09-25)

`profiles/genotype_c4_v4.tsv` and `profiles/genotype_c4_v4_campaign_local.tsv` are their v3 profiles key for
key except the name and `repair_rounds 32` (was 11); the campaign gate checks it. v4 is the final C4 campaign's
scientific profile (see "Final campaign binding").

Why: the three v1 round-cap no-calls (`loo.HG02282.2`, `loo.NA19384.5`, `loo.NA18989.10`) were searches still
improving (by 7.3, 11.3 and 3.7 nats) when the 11-round cap stopped them. Rerun truth-free with only the round
limit (32) and a 300 s wall limit overridden (`tests/results/.../v3_round_cap/<job>/RECORD.md`): rounds 1-11
reproduce the frozen trajectories exactly, and each converges normally (`no_improvement`) at round 12, 12 and
14 -- maximum 14, with no tie, cycle, refusal or timeout, and cumulative resource margins above 99%. 32 is frozen
with headroom; a search stops at convergence, so an easy sample never pays for unused rounds. The three become
strict calls with excess 17, 11 and 6,343 bp.

NA18989.10's 6,343 bp are 6,332 bp at sv:7 (called 65,449/59,084 bp, truth two ~65 kb alleles), the same sv:7
error NA18989's v1 calls carry at 2x, 5x and 30x. Its representation filter is PASS: the marker prefers another
pair by 232 nats, but one with the same exposure, so it is not a depth-versus-sequence conflict and the guard is
not meant to flag it -- this is not a guard failure. It has 0/11 high-quality blocks, because every block's
genotype support is `PRIOR_ONLY`: the 6.3 kb spike stays in the all-quality and representation-PASS views and is
absent from the high-quality view.

## Union placement reuse, phase 1 (2026-09-25): implemented, off, bound by no profile

The union pass places reads on haplotypes that repeat each other's blocks (58.8% of NA19434.30's window
positions repeat). Phase 1 serves the repeated local problems from a content-addressed cache without
changing a single output: likelihood, accepted origins and their multiplicity, candidates, scores, ranking,
calls, phase and filters are the reference scorer's, which stays the oracle.

- **The reuse unit.** A segment is one union block's bytes on one strand; the key is (strand, bytes), with the
  read span R and `max_gap_bp` g fixed per run. A placement starting at s depends on [s-g, s+R+g) only, so a
  segment's interior zone -- local starts [g, len-R-g] -- has placements that do not depend on its neighbours.
  Segments occurring at least twice are filled once and their interior placements translated to every
  occurrence; every start outside a served zone, including every window crossing a boundary and every short
  or empty allele, is placed as before. Replay appends before the row sort, whose (start, consumed) keys are
  unique, so rows are bit-identical. Segmentation is a hint only: a wrong one stays exact (fixture).
- **Admission.** Canonical order (occurrences desc, length asc, strand, first occurrence), admitted while
  payload + index fit `union_placement_cache_bytes`; the rest are declined and placed directly. No eviction.
- **Resource accounting.** Executed counters count executions only: the fill's, plus the direct placements.
  Replayed placements are counted by `union_cached_placements_consumed` and bounded by
  `max_union_cached_placements_consumed` (the cap accepts exactly, and one below refuses before any
  ranking). Physical origins and the other statistical caps are unchanged, and cache-served origins still
  count against them; not counting them is phase 2, which is not implemented.
- **Keys.** `union_placement_cache_bytes` and `max_union_cached_placements_consumed` (profile keys and probe
  overrides, both default 0 = off/no cap); the status records the budget, segments (distinct, eligible,
  admitted, declined), hits, misses, zone starts requested/served, placements enumerated/cached/replayed,
  payload, index and peak bytes and fill seconds. No v1-v4 profile binds them.
- **Gates.** `panvar_genotype_union_placement_cache`: substitution/insertion/deletion/empty alleles, first and
  last windows, short blocks, tandem repeats, identical content at offsets, both strands, 1/2/4/8 threads,
  permuted order, forced hash collisions, cache off/full/undersized, byte budget exact and one below,
  consumption cap exact and one below, wrong and whole-haplotype segmentations -- mass cells, exposure,
  components, totals, ranking and accepted/origin/state-class counts equal to the reference. Real:
  `loo.HG02282.2` (v4, 8 threads) cache on = off in every table, the trajectory and the status apart from the
  cache's keys and timings.
- **Benchmark** (`loo.NA19434.30`, v4 campaign-local, `max_union_physical_origins` 7e8, 300 s, 8 threads,
  512 MiB budget; `tests/results/.../v4_union/RECORD.md`). Both runs refuse at the same 700,000,001 origins
  (`RESCORE_INCOMPLETE`, no genotype artefacts); 270 haplotypes completed in both with identical accepted
  placements, origins and state classes (the one haplotype in flight at the refusal is scheduling-dependent).
  276/276 eligible segments admitted; 65% of zone starts served; 135.0M placements replayed; payload
  191.8 MB, peak 241.9 MB, fill 10.1 s. Placement-and-mass 165.3 -> 99.1 s (1.67x), seeding+emission
  worker-seconds 3.56x, wall 167.7 -> 101.5 s, CPU 1,164 -> 473 s, footprint 2.16 -> 2.02 GiB. At 2x
  (`loo.HG02282.2`) there is no overhead: wall 23.8 -> 22.2 s, footprint 0.57 -> 0.62 GiB.

## Union origin reuse, phase 2 (2026-09-25): implemented, off, bound by no profile

Phase 2 serves the physical origins wholly inside a repeated segment, on top of phase 1. Scores, ranking, calls,
phase and filters are the reference scorer's bit for bit; no v1-v4 profile binds it and no v5 exists.

- **The partition.** An origin is (fragment, arm, forward start s, insert I); its latent mass sums the forward
  placements starting at s against the reverse placements ending at e = s + I, and nothing else. When s is in
  one segment occurrence's forward interior zone [o + g, o + L - R - g] and e in the same occurrence's reverse
  interior zone [o + R + g, o + L - g], every one of those placements is a phase-1 replayed placement (both
  strands admitted), so the origin -- insert, orientation, latent mass summed in the same order -- depends on the
  segment's bytes alone. Those (s, e) pairs are skipped by the direct pairing and served; every other pair --
  across a segment or block boundary, near a locus end, in a short, declined or single-occurrence segment -- is
  paired directly. The rule is on (s, e) alone, so the two sets are disjoint and exhaustive.
- **Representation and arithmetic.** Per distinct segment, the served origins are built once from its cached
  placements with the same pairing routine the direct path uses (`pair_fragment`, shared by the reference path,
  the direct path, the construction and the self-check) and stored as sorted per-fragment StateKey
  multiplicities (`StateAggregateSet`, the repair cache's merge representation). Each occurrence adds its own
  copies, so multiplicity is preserved. Per (fragment, haplotype) the direct StateKeys are sorted and counted, the
  served aggregates merged in (`merge_sorted`), and the mass reduced with `reduce_sorted` -- ascending StateKey,
  one `log_add_exp` per key, exactly what `reduce_states` does over the reference's map. Self-check
  (`set_origin_cache_self_check`): every served fragment census is compared to a fresh enumeration.
- **Admission.** Canonical phase-1 order; admitted while the payload fits `union_origin_cache_bytes`; declined
  segments are paired directly; no eviction.
- **Resources.** `physical_origins`, `statistical_state_classes`, `latent_alignment_pairs` and `valid_fr_states`
  count executed work -- direct pairing plus the one-time construction -- and `max_physical_origins` and
  `max_state_classes` bind that. Served work is counted apart (`reused_physical_origins`,
  `reused_state_aggregates`, `reused_latent_alignment_pairs`, `reused_valid_fr_states`) and bound by
  `max_union_cached_origins_consumed` and `max_union_cached_states_consumed`, both required positive when phase 2
  is on; a haplotype's served occurrences are charged before any is merged, and crossing either cap refuses with
  no ranking. Logical totals: reference origins = executed - constructed + served (likewise pairs and valid
  states); reference classes = `merged_state_classes`. Status keys `union_origin_cache_*`,
  `union_cached_origins_consumed`, `union_cached_states_consumed`, `seconds_union_origin_cache_fill`; trace
  columns `reused_physical_origins`, `reused_state_aggregates`, `merged_state_classes`, `seconds_cached_merge`.
- **Schedule.** With phase 2 on, union haplotypes are dealt to workers round-robin instead of in contiguous
  spans (each column is still computed whole by one worker; counters are integer sums or maxima). Needed for
  the wall target: served mosaics cost a fraction of panel haplotypes and the interned order groups them -- on
  `loo.NA19434.30` the slowest worker would take 110.9 s dealt contiguously against 59.5 s round-robin.
- **Union-only diagnostic** (`diagnostic_union_only`, probe `--union-only 1`): stops after the union pass is
  scored and ranked, writes `<prefix>.union_diagnostic_ranking.tsv` (every ranked candidate's score, terms and
  prior in `%a`, plus a digest of every mass cell and exposure) and status `UNION_DIAGNOSTIC_ONLY`; never a call.
- **Gates.** `panvar_genotype_union_placement_cache` now also covers phase 2: substitution/insertion/deletion/
  empty alleles, long interiors and short blocks, inserts crossing one or several boundaries, flush locus ends,
  a tandem repeat, identical segments at offsets and twice in one haplotype, both arms, a near-duplicate segment
  whose fragments have origins of equal insert and orientation but different latent mass (1,635 such classes
  inside served interiors), forced hash collisions, origin cache off/exact/one byte below/undersized and over an
  undersized placement cache, both consumption caps exact and one below (also at 4 threads) and zero (invalid),
  the executed caps exact and one below (the reference refuses at the reused run's origin cap), 1/2/4/8
  threads, permuted haplotypes, reversed candidate requests, wrong and whole-haplotype segmentations, and the
  self-check -- mass cells, exposure, components, totals, ranking, winner set and the logical placement, origin,
  pair, valid-state and class totals equal to the reference's.
- **`loo.HG02282.2`** (full v4 call, 8 threads, phase 1+2 on vs off): every table, the 12-round repair
  trajectory and the status apart from executed work, timings and the caches' keys identical; origins 28,921,915
  = 15,252,027 executed (1,623,619 constructed once) + 15,293,507 served; classes 25,158,698 = merged; executed
  origins 1.90x fewer. Wall 24.1 -> 21.0 s, footprint 0.58 -> 0.66 GiB. The new build's reference path equals
  the phase-1 (07cd4b2) off run in every table and status key.
- **`loo.NA19434.30`** (union-only diagnostic, v4 campaign-local, unchanged 700M executed-origin cap, 512 MiB
  placement cache, 1 GiB origin cache, both served caps 2e9, 300 s): the union COMPLETES at 276,056,754 executed
  origins (39.4% of the cap). The reference (reuse off, diagnostic cap raised to 2e9) needs 700,039,154 logical
  origins -- 39,154 over the cap, which is why v4 refuses it -- and 551,131,166 classes; phase 2 reconciles to
  both exactly (executed 276.1M incl. 26.2M constructed + 450.2M served; classes 226.2M executed, 347.9M
  aggregates served, 551.1M merged). Executed origins 2.54x fewer, classes 2.44x. The certified ranking of all
  8,384 sequence candidates (TIED, winner set 16) and the digest of all 10.45M mass cells are bit-identical to
  the reference's. Union pass 72.4 s (reference 155.8 s), wall 74.8 s, CPU 500 s (1,213 s); worker-seconds
  434.5 (1,210.5). Origin cache 138 of 138 segments, 2,166 hits, 545 misses, 22.9M aggregates, payload 554.0 MB,
  fill 0.66 s; placement cache fill 10.2 s. Footprint 3.07 GiB (reference 1.93 GiB): the origin cache 0.52 GiB,
  the placement cache up to 0.23 GiB, and about 0.45 GiB of freed large blocks the allocator keeps dirty
  (`vmmap`: MALLOC_LARGE empty) -- not a leak, not yet reduced.

## Union seeding, phase 3 (2026-09-25): exact, always on

With phases 1+2 on, the union's cost is panel-haplotype seeding. Phase 3 changes how three exact steps are
computed, never what they compute; it needs no setting and changes no counter's meaning (every executed and
logical counter on `loo.NA19434.30` equals phase 2's).

- **Measurement.** One instrumented union-only run (per-haplotype block timers: seed scan, start sort,
  verification, replay, row sort; trace columns `seconds_seed_scan` ... `seconds_row_sort`) and one
  deterministic single-thread replay (`--union-seeding-profile N`, a union-only diagnostic that replays instead
  of scoring; `<prefix>.union_seeding_profile.tsv` binds the sampled haplotypes by index and FNV hash). Over 16
  evenly spaced panel haplotypes (30.4 s single-thread; 16 mosaics 1.5 s): Myers prefilter 53%, emission key and
  cache lookup 16%, start sort + unique 14%, window ACGT validation 11%, seed scan and generation 3%, DPs 2.5%
  (65k executed, 11.4 us each), storage, replay and row sort <1%.
- **Early-refusing prefilter.** The Myers last-row score is the whole read's edit distance to the target prefix
  so far, and adjacent prefixes differ by at most one edit; once score > band + remaining columns no prefix can
  come within the band, and the scan returns false -- the full scan's answer, earlier. Acceptances return where
  they did. `feasible_within_band_reference` keeps the full scan.
- **ACGT validation hoisted.** Each read form's all-ACGT flag is taken at preparation and each haplotype strand's
  once per call; when both hold, the per-window check (203 bytes per start) cannot fail and is skipped. Any
  non-ACGT byte keeps the per-window check, so a refusal fires exactly where it did.
- **Start set.** std::sort + std::unique over (form, start) replaced by an in-place counting permutation by form
  (American-flag cycle leading) and per-form insertion sorts: the same list. The replay's out-of-place version
  was 4.9x faster than the sort, but its scratch copy of each worker's proposals cost ~0.5 GiB and tripped the
  3.5 GiB footprint limit; the in-place version saves 14% of the sort in production.
- **Not done, not measured.** Batching several starts through one kernel pass, and SIMD/NEON: not attempted (an
  anchored Myers scan per start shares no state with its neighbours as written; whether a batched layout pays was
  not tested). A necessary-condition screen: not repeated (the earlier piece screen saved little net). The
  emission key and lookup (16% before phase 3) is the largest untouched share.
- **Gates.** `panvar_genotype_seeding_exactness`: 85,931 prefilter decisions over adversarial edit scripts
  (substitution/insertion/deletion/mixed, on piece boundaries, +/-4 around the true start), repeats, unrelated
  sequence, windows cut by the haplotype end, reads above the 128-bit word, bands 0-8, gaps 0-4 -- the early and
  full kernels agree, prepared and unprepared; 3,000 randomised proposal lists (duplicates, empty, one form,
  forms above the insertion threshold) equal to std::sort + std::unique; non-ACGT window or read still refuses.
  The replay compared 64 real proposal lists and 47.9M real prefilter decisions: identical. The phase 1+2 fixture
  (threads, orders, caps) and `loo.HG02282.2` (reuse off and on: every table, the trajectory and every status
  counter equal to phase 2's) pass.
- **`loo.NA19434.30`** (union-only, phases 1+2, v4 campaign-local, 700M cap, 8 threads): ranking and mass digest
  identical to the reference; union 72.4 -> 52.1 s (1.39x; 2.99x over the reference's 155.8 s), wall 74.8 ->
  54.4 s, CPU 500 -> 357 s, placement fill 10.2 -> 7.8 s, footprint 3.07 -> 3.21 GiB. Worker-seconds 423 -> 305:
  panel verification 308.8 -> 197.0, panel start sort 47.2 -> 40.6; slowest worker 57.9 -> 41.7 s.
  `loo.HG02282.2` 2x full call: 24.1 -> 21.3 s (reuse off), 21.0 -> 19.2 s (on).

## Union hot loop, final phase (2026-09-25): exact; the replay overpredicted

Measured first in the replay-only process (`--union-seeding-profile 16`; `phase4_replay_loo.NA19434.30`), then
two levers implemented, both exact and always on; the benchmark delivered 2.2%, not the ~10% predicted.

- **Replay, panel (16 haplotypes, 20.9 s single-thread after phase 3).** Prefilter (early) 11.3 s, start sort
  3.2 s, emission key 2.5 s, cache lookup 1.06 s -- of which the two clock reads 0.29, the shard lock 0.10, the
  component copy 0.26, the probe 0.43 -- seed scan 0.8 s, DPs 0.7 s.
- **Emission key: implemented.** The key hash is byte-serial, so its state after "read + separator" is a function
  of the read: each prepared form carries it (`ReadForm::key_prefix`, taken from the form's own bytes, living
  with them) and the window resumes from it -- h0 and h1 are key_view_of's exactly, so table, shard, probe order,
  hits, misses, retention and budget are unchanged, and identity is still full byte comparison. The public
  `gapped_emission` still hashes the whole read. Replay: key 2.49 -> 1.20 s, 9.28M keys checked equal. The
  lookup's clock reads were left (they feed the repair's reported `seconds_lookup`).
- **Multi-lane Myers: rejected.** Streaming 2, 4 and 8 starts of one form through the early kernel, lanes refilled
  as they decide: 14.9, 23.2 and 25.2 s against the one-start kernel's 11.3 s (143.5M decisions identical). Lane
  bookkeeping and 128-bit state in memory cost more than the overlap gains; no NEON path was built.
- **Fill balancing: implemented.** The placement-cache fill ran in batches of 8 segments, each waiting for its
  longest (7.58 s of batch makespans over 37.5 worker-s). Now bounded batches of 64, claimed longest first from a
  shared cursor; results land in canonical slots, counters merge in slot order, admission stays canonical after
  the batch. Simulated from the replay's per-segment times: 4.92 s.
- **Gates.** `panvar_genotype_seeding_exactness` adds 20,000 resumed keys equal to the full hash (and under forced
  collisions), and six scorings (full, undersized and one-byte cache budgets, with and without forced
  collisions) whose mass cells and every emission-cache counter -- prefilters, requests, hits, misses, executed
  DPs, entries, refusals to retain, payload and budget bytes -- equal those of hashing the whole read. The union
  fixture now requires the placement cache's contents, payload and admission to be identical at 1/2/4/8 threads.
  `loo.HG02282.2` (reuse off and on): every table, the trajectory and every status counter equal to phase 3's.
- **`loo.NA19434.30`** (union-only, phases 1+2): ranking and mass digest identical to the reference, every status
  counter equal to phase 3's. Union 52.1 -> 51.0 s (1.02x; 3.06x over the reference), wall 54.4 -> 53.4 s,
  footprint 3.21 -> 3.19 GiB. The fill fell 7.75 -> 6.19 s (simulated 4.92); placement did not move (41.7 -> 42.0 s
  slowest worker; panel verification 197.0 -> 193.6 worker-s, against ~15 predicted): under eight threads the
  key hashing the replay removed was not on the critical resource. Run-to-run noise on these per-phase
  worker-seconds is several seconds.
- **Wall now** (51.0 s union, 53.4 s job): placement 42.0 s (panel seeding ~83% of worker time: prefilter, then
  start sort, emission key and lookup), placement-cache fill 6.2 s, diploid combination 1.8 s, origin-cache fill
  0.7 s, fragment preparation 0.2 s; outside the union, model and markers 2.0 s.

## Structural assessment and union prefilter reuse (2026-09-25): exact, off, bound by no profile

- **Census** (`--union-verification-census <dir>`, a union-only diagnostic that replays instead of scoring;
  `census_union_verification`). Every direct verification request of `loo.NA19434.30` after phase 1 zones is
  identified by (read content, window content): read contents interned by byte equality, window contents --
  every eligible window of every haplotype strand, length min(available, read + max_gap_bp), which is exactly
  the clamp the emission's normaliser uses -- by byte equality against a representative (hashes choose slots
  only). Every repetition of an identity reproduced its decision and its column count. The stream reconciles
  exactly: 417,066,267 prefilters + the fill's 52,134,015 = production's 469,200,282.
  - 128,127,356 eligible windows hold only 743,408 distinct contents; 417.07M requests hold 8,573,016 distinct
    identities: 97.9% of requests repeat an identity (66.5% at multiplicity 65-128 -- once per panel
    haplotype; maximum 868). Mosaic requests (7.2M) are all repeats of panel identities.
  - Production kernel, 8 threads: panel accepted 37.2 worker-s, rejected 108.5; mosaic 2.2. The time spent on
    repetitions -- the exact upper bound on what remembering the first decision can remove -- is ~106 (panel
    rejected, by columns) + ~37 (accepted) + ~2 (mosaic) worker-s.
  - Accumulation: one panel haplotype holds 411,726 distinct windows; 128 hold 743,408 (fit
    411,726 + 38,728 (n-1)^0.449, within 1.3% at 64 and 128).
- **Implemented: union prefilter reuse** (`union_prefilter_cache_bytes`, `max_union_cached_prefilters_consumed`;
  off; the cap must be positive when on; one read length required). Per pass, an in-memory window-content
  index (256 sharded tables; each distinct window's bytes copied into the index's own 64 KiB slabs; identity by
  length and complete byte comparison) gives every eligible start a window ID; read forms carry their content
  ID. A fixed lock-free table keyed by (read content ID << 27 | window ID) -- exact integers, no hash decides
  equality -- holds each decision, feasible or refused, filled to half its slots and then declining (the kernel
  runs). A held decision replaces only the kernel: `gapped_emission_core` moves every logical counter as before
  and `executed_prefilters` counts kernels run; hits are published in batches against their own cap, and the
  exact merged total is what acceptance sees. Serves the placement-cache fill too.
- **Gates.** Union fixture: prefilter reuse alone, over phase 1 and over phases 1+2, compared with the same
  configuration without it -- mass cells, exposure, components, ranking, winner set, logical placement/origin/
  class totals and every logical verification counter identical, executed prefilters + served decisions equal
  to the kernels run without reuse -- at 1/2/4/8 threads, permuted haplotypes, reversed requests, forced
  window-hash collisions, an undersized table, the consumption cap exact and one below (1 and 4 threads), a
  zero cap invalid. `loo.HG02282.2`: reuse off and phases 1+2+prefilter on each equal to phase 4 in every table,
  the 12-round trajectory and every status counter (executed-prefilter and the caches' own keys apart).
- **`loo.NA19434.30`** (union-only, phases 1+2, 128 MiB decision table, cap 2e9): ranking and mass digest
  identical to the reference; every status counter equal to phase 4's apart from executed prefilters.
  Executed prefilters 469,200,282 -> 10,965,482 (458,234,800 served; 8,388,608 decisions held -- the table's
  half-load limit -- 487,756 declined). Union 51.0 -> 33.7 s (1.51x; 4.62x over the reference's 155.8 s), job
  wall 53.4 -> 35.9 s, CPU 366 -> 242 s, placement-cache fill 6.2 -> 4.4 s, footprint 3.19 -> 3.41 GiB (index
  128 MiB, table 128 MiB). Window indexing 14.8 worker-s. `loo.HG02282.2` with reuse on: 18.95 -> 19.90 s (+5.0%),
  footprint 0.66 -> 0.91 GiB; reuse off unchanged.
- **Persistent graph-side index: specified, not built.** The per-pass index costs ~15 worker-s (~2 s of wall) at
  C4, the most a persistent one could save here. Design, for larger panels: an immutable file holding only
  graph-derived content -- schema version; provenance (graph sha256, reference path, block model and spelling
  parameters, panel curation, builder version); canonical path identities (allele paths and the sha256 of each
  spelled haplotype); per window length W = read length + max_gap_bp it was built for, the distinct window
  contents (bytes, length, end flag) with their IDs and every occurrence (path, strand, offset). No read, sample
  or decision is stored. A run loads it only when the graph sha256 and every window-defining parameter (W,
  max_gap_bp, the spelling parameters) match its header; any mismatch refuses with an explicit error, and it is
  never rebuilt or updated in place by a caller. LZO uses every indexed path; LOO uses the same index with the
  held-out donor's paths masked (window IDs are content classes, so masking removes occurrences, never
  identities). A proposed or synthetic mosaic's windows are looked up by content (hash to a slot, bytes
  compared against the stored copy); a window the index lacks takes a per-run ID above the stored range, so its
  first request runs the kernel exactly as now. Another read length or gap width is another W: refused unless
  the index was built for it. With no index the per-pass in-memory index (or none) runs, and the call is
  unchanged: identical decisions either way, because IDs only ever name exact byte contents.
- **Remaining opportunities, measured (not implemented).** Out-of-place counting start set: replay 0.89 s vs
  3.19 s in-place over 16 panel haplotypes; it needs a scratch copy of each worker's proposals (~40 MB x 8), which
  puts the NA19434 footprint above the v4 profiles' 3.5 GiB watchdog (not above 4 GiB); after prefilter reuse the
  start sort and window indexing are ~58 of 198 worker-seconds. A 12 s sample of the 8-thread placement (87,805
  top-of-stack samples, ~20% idle waits): the start set (in-place permutation 15.3% + its per-form introsort 6.2%),
  the seeding body inlined in collect_form_placements 18.1%, the emission key for accepted requests 6.4%, the
  cache lookup's clock reads 4.6%, window identities and their shard locks ~5%; the pairing map does not appear
  among the top frames (pairing and reduction together are a few percent) -- replacing it is not material.
  Thread scaling of the final configuration (union-only, rankings identical at every count): 143.5 / 77.6 /
  46.2 / 33.7 s at 1 / 2 / 4 / 8 threads (1.85x / 3.11x / 4.26x). Haplotype-count scaling (projection from the census fit, not measured): distinct
  windows 0.88M / 1.05M / 1.28M at 256 / 512 / 1,024 panel haplotypes (1.18x / 1.41x / 1.73x of 128); the decision
  table and index grow with them (~226 + ~185 MiB at 1,024). Per-haplotype work outside the kernel -- seed scan,
  start sort, indexing, probes, pairing -- stays linear in the haplotype count: prefilter reuse lowers the per-
  panel-haplotype cost from 2.18 to 1.27 worker-seconds, a constant factor, not a change of order.

## v5: memory for speed, frozen (2026-09-25)

`profiles/genotype_c4_v5_campaign_local.tsv` is its v4 profile key for key except the reviewed settings below;
`panvar_genotype_v5_profile` (tests/genotype_v5_profile_gate.sh) proves that is the whole difference, that v4 sets
none of the new keys, and that the comparison catches a copy with one more changed setting. Only the campaign-local
profile has a v5; no global default moved.

| key | v5 | measured on loo.NA19434.30 (union-only, 8 threads) |
|---|---|---|
| union_placement_cache_bytes / max_union_cached_placements_consumed | 512 MiB / 2e9 | payload 191.8 MB; 135.0M consumed |
| union_origin_cache_bytes / max_union_cached_origins_consumed / max_union_cached_states_consumed | 1 GiB / 2e9 / 2e9 | payload 554.0 MB; 450.2M origins, 347.9M aggregates served |
| union_prefilter_cache_bytes / max_union_cached_prefilters_consumed | 128 MiB / 2e9 | 8.39M decisions held (half-load limit), 0.49M declined; 458.2M served |
| start_sort_scratch_bytes | 320 MiB | 40 MiB share per worker, 21.5 MB peak, 820 lists out of place, 0 fallbacks |
| footprint_limit_bytes | 4 GiB (v4: 3.5 GiB) | 3.28 GiB |
| emission_lookup_timing | 0 | status `emission_lookup_timing not_measured`; repair `seconds_lookup` NA |

repair_rounds (32), the representation guard and every calling, correctness and resource setting are v4's.

- **Out-of-place start sorter** (`start_sort_scratch_bytes`, off unless given). Each worker may hold budget /
  workers bytes, 4 per proposal (starts only; the form is the bucket), kept across haplotypes and never above
  its share; a list larger than the share is ordered by the in-place path and counted as a fallback
  (`start_sort_fallbacks`) -- the same list either way, so the contract is fallback, never refusal. Status:
  `start_sort_out_of_place`, `start_sort_in_place`, `start_sort_fallbacks`, and (thread-dependent)
  `threads_start_sort_workers`, `threads_start_sort_share_bytes`, `threads_start_sort_peak_worker_bytes`.
  std::sort + std::unique stays the fixture oracle. A first benchmark at a 128 MiB budget (16 MiB share) left 106
  of 820 lists in place (panel lists reach ~5.4M proposals); at 320 MiB none.
- **Lookup timing opt-in.** Emission-cache lookups read no clock unless `emission_lookup_timing` is on; the
  status says whether it was measured and the repair's `seconds_lookup` is NA when not. Exact; its benefit on
  NA19434 (timing on 30.1 s / 215.8 s CPU vs off 29.3 s / 209.9 s, two runs each) is within run-to-run noise.
- **Gates.** Seeding exactness: out-of-place and in-place lists equal std::sort + std::unique on 3,000
  randomised lists. Union fixture: a full budget (out of place, no fallbacks) and a one-byte budget (every list in
  place, counted as fallbacks) at 1/2/4/8 threads, permuted haplotypes and reversed requests, each identical to
  the reference, the worker peak never above its share. `loo.HG02282.2` (v4 profile, reuse off and every lever on)
  equal to phase 5 in every table, trajectory and status counter; wall off 21.09 -> 20.80 s, on 19.90 -> 19.90 s.
  With the v5 profile itself: every table and the trajectory equal v4's, 21.42 s, 0.90 GiB.
- **`loo.NA19434.30`** (every v5 setting as overrides on v4): ranking and mass digest identical to the reference;
  union 33.7 -> 28.9 s (1.16x; 5.39x over the reference's 155.8 s), job 35.9 -> 31.1 s, CPU 242 -> 207 s, footprint
  3.28 GiB; panel start sort 49.3 -> 14.8 worker-seconds, slowest worker 26.5 -> 22.1 s.
- **Persistent graph-side index: not implemented for C4.** The per-pass exact window index costs ~15
  worker-seconds (~2 s of wall) here, the most a persistent one could save. Reconsider before larger loci, with
  their real haplotype and window-content counts; the per-pass index stays the fallback and the scientific
  reference either way.

## Combined block quality (2026-09-25): quality_filter, block_quality_rule_v1

Every block row of `genotypes.tsv` and `block_genotypes.tsv`, on every call path, ends with `quality_filter` (PASS or
LOW_QUALITY) and `quality_reason`, from two truth-free facts the call already reports (`block_quality`):

| representation filter | genotype_support | quality_filter | quality_reason |
|---|---|---|---|
| PASS | READS | PASS | none |
| PASS | anything else | LOW_QUALITY | no_read_support |
| LOWQUAL_REPRESENTATION_CONFLICT | READS | LOW_QUALITY | representation_conflict |
| NOT_EVALUATED (or any other value: fail closed) | READS | LOW_QUALITY | representation_not_evaluated |

Several reasons are kept, in that canonical order, joined by ',' (e.g. `representation_conflict,no_read_support`).
The status adds `quality_rule`, `quality_pass_blocks`, `quality_lowqual_blocks`, their block-ID lists and one count
per reason. The field annotates: calls, ranking, the representation filter, genotype support and consensus are
unchanged. conditional_GQ is conditional within the scored candidate set and is not part of the rule -- it is not
sufficient for high-quality status. The evaluator carries the caller's two columns into `eval.blocks.tsv` and its
`high_quality_*` totals count `quality_filter == PASS` (a call that predates the columns is derived by the same
function, `quality_source derived_block_quality_rule_v1`); the consolidator's `blocks_high_quality.tsv` is
`quality_filter == PASS`, read, never rebuilt. Gates: `panvar_genotype_block_quality` (each case, combined reasons,
fail-closed values); the e2e test checks every block table its runs write (69 tables, 345 rows) against
`block_quality(filter, genotype_support)` of the same row; the downstream schema gate reads the column. The v5 pilot's
sealed `loo.NA19434.30` table derives 7 PASS and 4 LOW_QUALITY blocks (ctx:0, sv:1, sv:7, ctx:8, each
no_read_support), matching its high-quality view.

## v6: repair origin-cache budget 768 MiB (2026-09-25)

`profiles/genotype_c4_v6_campaign_local.tsv` is v5 key for key except `repair_origin_cache_bytes` 805306368 (768 MiB;
v5 384 MiB); `panvar_genotype_v6_profile` proves that is the whole difference. v5 stays frozen; no global default moved.
`panvar genotype --repair-origin-cache-bytes` overrides the budget and is recorded as an override.

- **Analysis from the sealed v5 pilot ledgers (projection).** With 384 MiB, rounds 2-6 of the main search each spent
  7.4-9.0 s re-enumerating windows the planned cache could not keep (43.2 s in all; 264-303 distinct windows a round,
  180-216 retained, 67-88 unselected); round 1's demand was 1.9 GB over 297 windows. The repair phase peaked at 2.57 GiB
  and the process at 3.37 GiB (union). A byte-proportional projection gave 768 MiB -11.7% total and -13.9% repair,
  with a footprint of at most 3.37 GiB; 768 MiB, the largest listed budget under 3.85 GiB, was chosen for the pilot.
  The branches have no per-round ledger, so their share was an assumption.
- **Comparison pilot (measured)**, `loo.NA19434.30`, v5 profile with only the budget overridden, clean f43710a build:
  union ranking, pooled candidates, scores and winner, the six-round main trajectory, branch starts, rounds, stops and
  terminals, genotype and block tables (apart from the new quality columns), the equivalence set and FASTA identical;
  the sealed evaluation identical (distance 26,430, floor 219, excess 26,211; high quality 7 blocks, excess 0).
  Repair 164.0 -> 127.7 s (-22.1%), total 194.8 -> 159.0 s (-18.4%), CPU 920 -> 716 s; re-enumeration 7.4-9.0 ->
  5.8-6.0 s a round; branches 67.2 -> 45.3 s; repair physical origins 1.061e9 -> 0.822e9. Footprint 3.54 GiB (repair
  phase 2.57 -> 3.43 GiB): the projection's 3.37 GiB was low -- the cache's growth cost ~2.2x its bytes -- but it
  stays under 3.85 GiB.

## Final campaign binding (2026-09-25)

The final C4 campaign runs `panvar genotype --profile profiles/genotype_c4_v6_campaign_local.tsv`. Its scientific
profile is `profiles/genotype_c4_v4.tsv` (tie branching with cap 16, the representation guard, 32 repair rounds);
v4 campaign-local carries v4's science key for key (only the name, threads and process limits differ), v5 adds
exactly the reviewed campaign-local resource and exact-acceleration settings (`panvar_genotype_v5_profile`) and v6
exactly the 768 MiB repair origin cache (`panvar_genotype_v6_profile`). The v1 profiles stay the frozen bc5dd97
campaign's and are never used by it; the v1-era tie diagnostics and production regressions keep them.

- **Plan.** `tests/genotype_production_c4_campaign.sh plan <out> <build>` refuses unless the profile chain holds, the
  build is a clean tree at the production commit `fefa6c0` (census-free tie-branch seeds, accepted with the v6
  sentinels) or at a descendant whose `src`, `include`, `CMakeLists.txt` and `external` are identical to it, and
  `<out>` is new -- never the frozen campaign, never a directory a campaign ran in -- and named
  `genotype_c4_production_<build commit, 7 hex>`. `settings.tsv` records the production and scientific profile
  paths and hashes, the intermediate profiles and gates, the production commit, the build commit, the production and
  evaluator binary hashes, the runner, harness, consolidator and mosaic-harness hashes, and the scheduler: one
  concurrent job at 8 threads unless a scheduler calibration selects more. The runner refuses a production binary,
  profile or scheduler that is not the plan's, and a resume that changes any of them.
- **Outcomes.** Completed: `COMPLETE`/`strict_call`, `CONDITIONAL_GENOTYPE_EQUIVALENT`/`genotype_equivalent_call`,
  and every declared no-call -- `REPAIR_INCOMPLETE`/`tie_branching_incomplete` (with its repair rounds and both tie
  tables) and `RESCORE_INCOMPLETE`/`divergent_tie` (the same, plus the equivalence set) among them. A no-call is never
  retried, never stops the campaign and never carries a genotype, block-genotype or FASTA output. `improving_move_tie`
  is v1-only: under v6 such a tie is branched.
- **Estimate.** `tests/genotype_production_c4_timing.tsv` was measured with the v1 profile at `b754647`; the plan's
  `estimate.tsv` is labelled historical and is not a v6 prediction.
- **Package.** The consolidator's high-quality table is the block-level `quality_filter == PASS`, never the call's
  `locus_filter`; no-call rows enter no block table or denominator; undefined-truth blocks stay in the tables flagged
  `truth_defined false`, are listed in `blocks_undefined_truth.tsv`, and are excluded from every block-accuracy
  denominator. `panvar_genotype_downstream_schema` holds an `NA20847.30`-shaped call (representation PASS, ten of eleven
  blocks LOW_QUALITY for no read support) and `panvar_genotype_production_campaign_plan` the plan's contract.

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

- **LOO can be representation-limited even above the availability floor.** The availability floor is an
  edit-distance oracle -- the panel allele pair nearest the truth in sequence -- not necessarily the
  maximum-likelihood available allele. When the truth allele is absent from the panel, the reads can prefer
  a different available allele than the sequence-nearest one, and the call then sits above the floor without
  any search or likelihood fault (HG00171 ctx:10, 13 bp at every coverage: the repair scored the floor pair
  and the reads preferred the called allele by 12.5 nats).
- **HG01114.2 sv:7 is a sequence-versus-copy-number trade-off caused by the missing truth-like short
  allele.** The per-fragment decomposition (`score-set --decompose`) shows the Poisson likelihood counting
  every module copy once: interior fragments carry exactly log(4/3), and on depth alone the likelihood
  prefers the true three copies by 18.4 nats (145.3 nats of multiplicity against 163.7 of exposure). The
  panel's only short allele (80) is 94 bp from the truth's, and the fragments covering those differences
  prefer the called long alleles by a net 44.7 nats, which reverses the call. No accounting fault, and no
  likelihood change is made; the marker score is not combined with the fragment likelihood, since both use
  the same reads.
- **Branch recovery improves callability, not proximity to truth.** A tie-branched call is one genotype of
  the pooled search, conditional on that search.
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

State as of 2026-09-24:

- **v1 is complete and immutable.** C4 production-v1 (`c4-panel-anchored-v1`, binary `bc5dd97`, evaluator
  `11627a2`): 484/484 jobs, package `tests/results/genotype_c4_production_bc5dd97/campaign/package`. Its
  profiles are never edited; nothing reruns it.
- **v2 improving-move tie branch-and-pool is integrated** at `cc765d8` (see the v2 section above) and was
  accepted in review. The 2x improving-move-tie blocker recorded below is resolved by it for v2; v1 keeps its
  `improving_move_tie` no-calls.
- **Representation-conflict filtering: accepted scientifically**, with a fail-closed and schema correction (see its
  section above). A LOO call can be the
  internally consistent likelihood optimum and still be catastrophically far from the truth, when the panel
  lacks the truth-like allele (HG01114 sv:7, above). The caller must flag such calls, truth-free, without
  changing the genotype the likelihood selects.

- **The ten non-tie no-calls of v1.** The three repair round caps are resolved (v4, `repair_rounds 32`). The seven
  union-scoring-cap refusals at 30x are safe resource refusals, addressed by acceleration, not by a higher cap.

Later work, not blocking correctness: performance engineering (the per-request merge dominates warm
branches; lower-complexity placement for 30x large loci), the final C4 campaign, and other loci.

The 30x footprint blocker below was resolved for the campaign by the 3.5 GiB local campaign profile (Sol,
2026-09-23); the global 2.5 GiB default is unchanged.

### Recorded before 2026-09-23 (historical, kept as measured)

- **The 2.5 GiB footprint limit is marginal at 30x in one process** (production campaign launch blocker).
  The union pass leaves allocator-retained pages that the repair then carries; how many varies run to
  run. NA19240 LOO 30x refused in repair round 3 at 2,684,635,872 bytes in one run and peaked at
  2,438,105,848 in another, with identical scores (repair start 1.55 GB against 1.15 GB); NA20509 LOO 30x
  peaked at 2,569,194,208. `malloc_zone_pressure_relief` returned nothing. The live accounted data is
  about 1.4 GB. Options are Sol's: a higher profile limit, the union and the repair in separate
  processes, or a smaller union emission cache.
- **At 2x a tie among improving neighbours stops the repair** (production campaign launch blocker).
  NA20509 LOO 2x: round 1's best neighbours tie with each other 36 nats above the panel anchor, the
  search cannot choose between them, and the frozen rule fails closed with `REPAIR_INCOMPLETE`. At low
  coverage few fragments link neighbouring blocks, so equal-scoring phase alternatives of the best move
  are likely common. Whether a tie among genotype-equivalent improving moves may be followed is Sol's
  decision; the rule is unchanged.
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

1. **Done:** the fail-closed and schema correction of the representation guard (`c810ef0`).
2. **Done for the three round caps:** all converge by round 14; `repair_rounds 32` is frozen in v4. The seven
   union-scoring-cap refusals (all 30x) are resource/performance refusals, handled in the acceleration phase,
   never by raising the production cap.
3. Major speed engineering, starting with the union stage.
4. One final C4 campaign (v6 campaign-local production profile, v4 science; see "Final campaign binding").
5. Final all-quality, representation-PASS and high-quality plots.
6. Other loci.

### Recorded before 2026-09-23 (historical)

1. Sol's review of this stop point: the panel-anchored production integration, its five regressions, and
   the prepared (not launched) final C4 campaigns.
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

v1: campaign binary `bc5dd97`, evaluator `11627a2`, package under
`tests/results/genotype_c4_production_bc5dd97/campaign/package`. v2: tie branch-and-pool at `cc765d8`,
accepted in review. Representation guard and downstream schema at `4b14d64`, under review. The active suite
is CTest's: 33 tests are registered, 2 of them disabled by default (`panvar_real_regressions`,
`panvar_rebuild_stats`), so 31 run -- at `4b14d64`, 31/31 (at `cc765d8`, before the downstream schema gate
was registered, 30/30). Targeted gates that
need the gitignored C4 artifacts run separately and are named with each result
(`tests/genotype_tie_branching_regressions.sh`, `tests/genotype_production_c4_regressions.sh`,
`tests/genotype_tie_branch_gate.sh`, `tests/genotype_score_set_gate.sh`).

### Recorded before 2026-09-23 (historical)

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
