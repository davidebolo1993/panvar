# M3-D: bounded, default-off integration of likelihood-guided joint-block search (experiment)

**Ancestry.** Branch `m3-d-joint-search` from `6922bc3` (M3-C). The library was unchanged since baseline `282e64e`
until this batch, and M2 stays closed.

**What this is.** An experiment, not a default promotion or a release. There is no prior change, residual editing,
dosage change, clustering, merge or push.

## 1. Integration contract (frozen)

**Option.** The profile key `repair_joint_search` is `off` (the default) or `experimental_v1`.
- `settings.tsv` records it only when a profile sets it, so a default run's output contract is unchanged.
- The `joint_*` status rows are written only when it is on.

**The rule** (M3-C, unchanged; `include/panvar/genotype_joint_search.hpp`):
1. Score every single-block alternative on each homologue, the other fixed.
2. Shortlist 4 per block and homologue by fragment + exposure, excluding the prior; ties by allele index.
3. Enumerate the shortlisted singles and every same-homologue two-block combination.
4. Rank the unordered sequence pairs, with the incumbent, by the complete score (with the unchanged sequence-level prior
   at the run's switch probability), under the repair round's own certified move rule: complete, unique best, strict
   improvement.

**Bounds.** At most 2,048 single-block alternatives and at most 2,048 catalogue candidates, checked from the model
before anything is scored. Over a bound the pass is skipped, never truncated. k, the thresholds and the bounds are never
tuned from outcomes.

**What it is not.** A bounded heuristic neighbourhood: not exhaustive search, not A*, and not a global-optimality
certificate.

**When it runs.** Once, after an ordinary trunk that reproduced, converged, and whose last round had a unique best:
- an improving tie, a round cap, a cycle or a refusal skips it, and the reason is recorded;
- so does a tie that includes the incumbent;
- so does a missing delta context, or the model bounds.

**After an improving move.** One continuation leg of ordinary single-block refinement:
- the existing leg pattern: attached context, `first_round` continuing the numbering, the remaining declared rounds,
  `start_from_slots` on the winner's census-holding slots, and the anchor and every trunk incumbent added as visited;
- joint search is never invoked again;
- the leg's rounds are appended to `.repair_rounds.tsv`;
- tie branching, rooted at the trunk, is not run after a continuation. A continuation that ends in an improving tie
  fails closed through the existing path (`improving_move_tie`), so no tied sequence is chosen.

**One context.** The refinement's own `PreparedScoringContext`; no second live context.
- Caps are cumulative and unchanged, and counters are never reset.
- A cap crossed in the pass or the leg follows the existing contract: RESCORE_INCOMPLETE (`cumulative_cap`), provisional
  where the existing rule allows, never a partial ranking.

**Candidates.** Every catalogue member, complete scores only, joins the existing pool (`source` `joint` or `joint_leg`),
so grouping, projection, phase, confidence, provenance and the FASTA come from the existing pooled decision. Shortlisting
scores are not candidates. No additional molecular phasing is claimed.

**Release.** The pass releases its slots (`release_transient`, keeping the incumbent and, on a move, the winner, whose
census is retained by delta); the leg releases its own afterwards.

**Recorded:**
- the decision, the move, the counts and the bounds;
- origins reused and recomputed, DP executions and window bases;
- full placements (expected 0);
- peak census bytes and transient payload peak and after release;
- slots;
- physical origins and state classes of the pass and the leg;
- footprint peak and timings.

The run's cumulative counters include the pass and the leg.

**Restricted** (refused before any output):
- with the library-sensitivity audit (`--library auto` with the audit on), because its views are not wired through the
  pass — unsupported in v1;
- with diagnostic tie branching, the refinement replay, `diagnostic_union_only`, or `repair_rounds 0`.

## 2. Gates (before real cells)

- **`gate.library`** (`tests/genotype_joint_search.cpp`):
  - incremental vs cold, every catalogue member (≤ 1e-9), with length changes on both homologues;
  - a two-block valley: both intermediates lose, the joint move wins and is the truth;
  - no improvement;
  - a tied winner, under the round's rule (uncertified or tie_at_best; winner set 2; no move);
  - an explicit bound skip that touches nothing;
  - a cap refusal (ResourceExhausted);
  - sequence-equivalent spellings (never a candidate, never duplicated);
  - 1 and 4 threads, plus a permuted request order: identical;
  - cleanup (transient payload and live slots return; only the winner's slots kept, with censuses);
  - 0 full placements.
  - The second-step deltas in both edit orders were gated bit-exactly in M3-C (`two-block-lik-gate`, both paths).
- **`gate.e2e`** (`tests/genotype_joint_search_e2e.cpp`, the real caller on a synthetic valley):
  - off leaves the valley uncrossed and writes no joint rows;
  - on moves, converges in the continuation, and calls the truth in the FASTA, block genotypes and final incumbent;
  - 1 and 4 threads give identical scientific outputs;
  - the bound skip (API test hook) keeps the ordinary call;
  - a cap crossed inside the pass gives RESCORE_INCOMPLETE (`refused_cumulative_cap`), with no FASTA;
  - the audit restriction is refused.
- **`gate.regress`:** the existing mosaic-core, block-delta and scoring-context exactness tests.

## 3. Bounded real validation (frozen order, inputs and comparisons)

**Cells:** the M3-C three, with the same input hashes (`INPUTS.tsv`). Controls first: `lzo_5x`, `mosaic_5x`, then
`mosaic_2x`.

**Arms**, all with 4 threads, one process at a time:

| arm | binary | profile |
|---|---|---|
| `baseline` | the previous caller (`bin/panvar_baseline`, sha `58f934ce…`, built from this branch's library before any M3-D change; library identical to `282e64e`) | v6 campaign (sha `74a3723b…`) |
| `off` | M3-D | v6 |
| `on` | M3-D | v6 plus `repair_joint_search experimental_v1` (`profiles/genotype_c4_v6_joint_experimental.tsv`, the only difference) |

**`compare`** (truth-free; writes the selection digest):
- **off vs baseline:** identical output-file sets and normalised contents (timing, thread, resource and memory lines,
  and round timing columns, excluded);
- **on vs off:** status, call, sequence and block digests, final score, moves, ties, pool, joint rows, wall time and peak
  RSS.

**`eval.<cell>.<arm>`** (truth; only after `compare`): the campaign evaluator (`genotype_probe evaluate`, scheme `mosaic`
or `lzo`) on off and on.

**Eligibility** (model only): `genotype_probe joint-bounds` on CYP2D6 (the LZO and LOO NA19650 graphs, curated) and on
ACOT. It reports candidate counts and whether the frozen bounds permit the experiment.

## 4. Budget

900 s, counting every attempt, development runs included. 4 threads, a 6 GiB envelope. A unit is not launched unless the
time spent plus its declared timeout fits.
