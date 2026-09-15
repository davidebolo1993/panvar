# Mosaic genotyper: state audit and plan for phases 4-8

Companion to `genotype_mosaic.md`, which remains the authoritative contract. This document records
what is actually built, what its own gates still owe, and the work order to reach the Phase 8
go/no-go.

## R5. Review round 5 — floor, harness, NA19240 rerun and the development ladder

The reviewer accepted the round-4 diagnosis: the scorer is internally consistent, and the substitution-only
model is inadequate for some leave-one-out indel alternatives. LZO may proceed once one blocker is
fixed, following a nine-step sequence. This section records steps 1–8. The held-out cohort (step 8) ran and surfaced three findings that need review before its result is reported (R5.8–R5.10).
Step 9, the multi-donor LOO ladder, is deferred as the reviewer ruled.

Commits:

| Commit | Content |
|---|---|
| `8ca8ff6` | One marker background floor, from depth fit to proposal scoring (step 1) |
| `203fc4e` | Development candidate-prefix harness, gated against standalone runs (step 5) |
| `28bbc3e` | Probe availability and no-call diagnostics, the `genotype` probe subcommand, the ladder driver (steps 4, 6) |
| `4488ebf` | Ladder selection script (step 7) |
| `aa76641` | Held-out cohort driver (step 8) |

Round 4's wording is corrected: the likelihood unit is the paired-end fragment, and R4 now says
"fragments" throughout.

### R5.1 One marker background floor (step 1)

The depth fit floored its background at 1e-9. The command then raised it to 1e-6 before building proposal
emissions (`genotype_command.cpp:457`, `std::max(depth.background, 1e-6)`). At an observed marker that a
candidate gives zero copies, each observation contributes `log(background)`, so that silent change moved
such scores by log(1000) = 6.9 nats per observation.

- `kMarkerBackgroundFloor = 1e-6` is now the model's only floor.
  - The fit keeps its background at or above it.
  - A declared background below it is refused; the declared-zero closed form is removed.
  - Proposal emissions refuse a background below it.
  - `proposal_marker_parameters` passes the fitted or declared background through unchanged.
- The status reports `marker_background_floor` and `marker_background_effective`.
- 1e-6 is the value round-3 proposal scoring already used. The NA19240 rerun therefore changes only the
  depth fit.
- **Gate.** The boundary fit ends on the floor, and the proposal parameters carry its background and rate
  exactly. The fixture has AAC and ACC each observed once and absent from a-rich/a-rich.
  - The mixed/mixed versus a-rich/a-rich contrast equals the Poisson kernels at the fitted background.
  - A 1000-fold background change moves it by 2·log(1000).
  - A background below the floor is refused.
  - End to end, the status reports the floor, and the effective background equals the fitted one.
- The real-C4 negative-intercept depth gate now reports background 1e-06; the other regimes are unchanged.

### R5.2 Clean worktree and NA19240 rerun (steps 2–3)

- **Build.** Isolated worktree at `8ca8ff6`, with submodules checked out at their pinned commits from the
  local checkouts, Release build. The tree is clean before and after the build, and all 16 tests pass.
- **Rerun.** `tests/genotype_probe_c4.sh`, same simulation and seed as rounds 2–4. The signature records
  `git_tree_dirty false`.

Every status reports `marker_depth_termination boundary_kkt`, rate 9.76988 and background 1e-6.

**Leave-zero-out, 1,024 merged states:**

| K | Status | Found | Pruned | Phase | Contenders | Exact | Distance | Winner |
|---|---|---|---|---|---|---|---|---|
| 8 | COMPLETE | 593 | 41,338,362 | PRIOR_ONLY | 2 | 11/11 | 0 | `285553cdf54d` |
| 16 | COMPLETE | 593 | 41,338,362 | PRIOR_ONLY | 2 | 11/11 | 0 | same |
| 32 | COMPLETE | 593 | 41,338,362 | PRIOR_ONLY | 3 | 11/11 | 0 | same |
| 64 | COMPLETE | 593 | 41,338,362 | PRIOR_ONLY | 3 | 11/11 | 0 | same |
| 128 | COMPLETE | 593 | 41,338,362 | PRIOR_ONLY | 3 | 11/11 | 0 | same |
| 256 | COMPLETE | 593 | 41,338,362 | PRIOR_ONLY | 3 | 11/11 | 0 | same |

What the corrected depth changed:

- **The candidate set.** 593 diplotypes found, where round 3 found 561, with a different pruning count.
- **Not the winner.** It is byte-identical to round 3 at every K, and exact at every block.
- **Phase at K = 8.** It is now PRIOR_ONLY (2 contenders) where round 3 had READS; all 11 block genotypes
  remain read-supported.

**Leave-one-out** (NA19240 removed):

| K | Status | Exact | Distance | Floor | Excess | Winner | Floor attained | Closest candidate: proposal rank, fragment gap, total gap |
|---|---|---|---|---|---|---|---|---|
| 64 | COMPLETE | 8/11 | 18 | 4 | 14 | `9ed92c5a` | yes | 1, 2,871.163, 2,863.035 |
| 256 | COMPLETE | 8/11 | 18 | 4 | 14 | `9ed92c5a` | yes | 1, 2,871.163, 2,863.035 |

- **K = 64 now makes the allele-11 call.** In round 3 it made the floor call (excess 0). That result
  depended on the unconverged depth fit.
- **The floor-attaining candidate is proposal rank 1.** The likelihood places it 2,871.163 fragment-term
  nats below the winner, the round-4 gap. With the proposal corrected, what remains is only the model
  misspecification R4 identified.
- **LOO is now stable across K, and consistently over the floor.**

### R5.3 LOO availability diagnostics and no-call metrics (step 4)

`genotype_probe evaluate` now reports the following. Everything that reads the truth is evaluation-only
and never enters a call.

- **Availability.** `availability_floor`, `minimum_candidate_distance`, `floor_attained_in_candidates`,
  `closest_candidate_proposal_rank`, `closest_candidate_fragment_gap_to_winner`,
  `closest_candidate_total_score_gap_to_winner` and `winner_excess_over_floor`. When several candidates
  attain the minimum, the highest-scoring one is compared.
- **Calls.** `run_no_call` means no genotype result exists. `block_unsupported` and `phase_unsupported`
  are reported separately and not counted as no-calls.
- **Accuracy.** Reported over every emitted block, and again over read-supported blocks only
  (`read_supported_blocks_exact`, `read_supported_unordered_distance_total`,
  `read_supported_excess_over_floor_total`), with SV and context distance apart.

**Checked against a known answer.** On the round-3 LOO K = 256 run, it names the K = 64 winner at proposal
rank 2 and distance 4 (the floor), with a fragment gap of 2,871.163 — the round-4 decomposition.

### R5.4 K prefixes equal standalone K (step 5)

`GenotypeCommandOptions::development_candidate_prefixes` is not a command-line option. It computes the
model, markers and proposal once. Each prefix k then passes through the same post-proposal code a
standalone run uses: refusal, spelling, exact prior, certified rescoring and every output.

- `unique_target_reached` and `candidate_set_conditional` are recomputed for k.
- Neither the beam's pruning nor its ordering reads the target.
- The production refusal rule applies to each prefix.

**Gate on the phaseable fixture.** Prefixes 2, 4, 16 and 64 must match standalone runs byte for byte:
manifest, scores, block calls, haplotypes, and status apart from stage timings.

| Merged-state limit | Prefixes that call | Prefixes refused | Identical to standalone |
|---|---|---|---|
| 64 | 2 | 2 | all four |
| exhaustive | 4 | 0 | all four |

Forcing every prefix conditional makes the gate fail.

**Real-data check.** At state limit 1,024 the ladder's extracted NA19240 K = 64 (commit `28bbc3e`) matches
the standalone K = 64 run from the clean rerun (commit `8ca8ff6`). The manifest, scores, genotypes, block
calls and haplotypes are identical byte for byte. The status files differ only in their input paths (the
`gfa` and library provenance), apart from stage timings. Both runs used the same reads.

### R5.5 The development LZO ladder (step 6)

The run:

- `tests/genotype_ladder_c4.sh`, from a clean worktree at `28bbc3e`; `git_tree_dirty false`.
- The eight manifest development donors, seed 7, 30× depth.
- Merged-state limits 256, 1,024 and 4,096.
- One proposal per limit, with K prefixes 8, 16, 32, 64, 128 and 256, and every prefix evaluated.

That is 144 prefix runs.

| Donor | Found at 256 / 1,024 / 4,096 | Winner identical wherever complete | Block calls identical | Distance, excess | K = 256 at state limit 256 |
|---|---|---|---|---|---|
| HG00733 | 138 / 661 / 2,952 | yes | yes | 0, 0 | refused (found < 256) |
| HG02666 | 224 / 819 / 2,677 | yes | yes | 0, 0 | refused |
| HG02818 | 256 / 1,024 / 4,096 | yes | yes | 0, 0 | complete |
| HG03371 | 256 / 1,024 / 4,045 | yes | yes | 0, 0 | complete |
| HG04036 | 256 / 1,024 / 4,006 | yes | yes | 0, 0 | complete |
| NA19036 | 256 / 1,024 / 3,946 | yes | yes | 0, 0 | complete |
| NA19240 | 140 / 593 / 3,016 | yes | yes | 0, 0 | refused |
| NA20509 | 224 / 872 / 3,376 | yes | yes | 0, 0 | refused |

"Found" counts unique diplotypes; no search at any limit was exhaustive.

- **Every completed setting is correct.** It is exact at every block (11/11, and 5/5 SV blocks), and each
  donor's winner is byte-identical across all 18 settings where it completed.
- **The only refusals are expected.** Four donors refused K = 256 at state limit 256, because a pruned
  search found fewer than 256 candidates. The production rule refused them and the harness recorded it.
- **Marker depth converged in every run:** `boundary_kkt` for seven donors, `interior_stationary` for
  HG03371.
- **Runtime.** A proposal plus six ranked prefixes took at most 144, 151 and 157 s at the three limits,
  with peak RSS at most 436, 424 and 442 MB. Wall time sums over donors to 730, 1,041 and 1,122 s.

**Phase support is not stable in K.** The selection rule ignores phase, but the ladder shows READS phase
at small K is not always confirmed by a larger candidate set:

- HG04036 is READS at K = 8 and 16 (state limits 1,024 and 4,096), and PRIOR_ONLY from K = 32.
- HG00733 is READS only at state limit 256 with K = 8.

A small K can drop the rival that makes phase a fragment tie. Block genotypes are unaffected: every rival
carries the same block calls.

### R5.6 Selection and freeze (step 7)

`scripts/select_genotype_ladder.py` reads the ladder table and each run's block calls, never truth distance.

- **Rule.** For each donor the reference is state limit 4,096 at K = 256. A setting is stable when every
  setting at or above it in both dimensions completed with the reference winner and unordered block calls.
- **Extension check.** No donor disagreed between 1,024 and 4,096, so none was extended to 8,192.
- **Result.** At state limit 256, the refused K = 256 prefix makes every smaller K unstable for four donors.
  At state limit 1,024 every K is stable for all eight. The selection is **state limit 1,024, K = 8**.

Frozen in `tests/genotype_frozen_settings.tsv` (commit `b1684ae`), before any held-out donor was simulated.
It is committed together with a declared **audit**: the cohort also ranks each donor's proposal at
K = 256, and reports separately whether the K = 8 winner and block calls, and any READS phase, survive
there. The audit exists because of the phase instability above, and it never changes a primary call.

The cohort driver is `tests/genotype_cohort_c4.sh`. It was smoke-tested on NA19240, a development donor,
into a scratch directory, to check its collection and aggregate; that run is not part of the cohort.

### R5.7 The held-out leave-zero-out cohort (step 8)

**Setup.**

- `tests/genotype_cohort_c4.sh`, from a clean worktree at `b1684ae`; `git_tree_dirty false`.
- The 56 held-out donors, seed 7, 30× depth.
- Frozen settings: state limit 1,024, K = 8, with the declared K = 256 audit.

**Collection.** It stopped at NA20355, because of a harness bug fixed in `5d179d2`: a refusal before the
prefix loop left no per-prefix files. The table below was re-collected from the unchanged run outputs with
the fixed collection step.

| Quantity | Value |
|---|---|
| Donors | 56 |
| Run no-calls | 1 (NA20355, R5.9) |
| Called | 55 |
| Block genotypes emitted, exact | 605, **601** |
| SV block genotypes, exact | 275, 272 |
| Read-supported blocks, exact | 604, 600 |
| Block genotypes not read-supported | 1 (NA19384 `ctx:0`, correct, PRIOR_ONLY) |
| Unordered distance (SV, context) | 120 (119, 1) |
| Availability floor, excess | 0, **120 in 4 donors** |
| Truth diplotype proposed within K = 8 | 51 of 55, at proposal and score rank 1 in all 51 |
| Globally certified | 0 (every search was pruned) |
| Audit K = 256: same winner, same block calls | 55 / 55, 55 / 55 |
| READS phase at K = 8; still READS at K = 256 | 10; **1** |
| Depth fit termination | 40 `boundary_kkt`, 15 `interior_stationary`, 1 `no_ascent` |
| Per run: wall time, peak RSS, candidates found | median 57 s (max 87); median 441 MB (max 516); median 992 (564–1,024) |

Per-run wall time covers one proposal and both ranked prefixes.

**The cohort surfaced three findings that need review** before this is reported as a result (R5.8–R5.10).
No code has been changed for any of them.

### R5.8 Finding 1 — four panel-template truths were never proposed

| Donor | Wrong block | Truth alleles (bp) | Called alleles (bp) | Distance | Certified gap (nats) | Phase at K = 8 |
|---|---|---|---|---|---|---|
| HG01573 | `ctx:4` | 374 / 374 | 374 / 374 | 1 | 8.36 | PRIOR_ONLY |
| HG02492 | `sv:3` | 244 / 244 | 244 / 244 | 1 | **0.0005** | PRIOR_ONLY |
| HG04217 | `sv:5` | 117 / 0 | 0 / 117 | 1 | 37.39 | READS |
| NA19331 | `sv:5` | 117 / 117 | 117 / 0 | 117 | **0.27** | PRIOR_ONLY |

**What the cohort run shows.** In all four, the truth is absent from the candidates at K = 8 and at
K = 256. The winner is the candidate closest to the truth: its gaps to the closest candidate are 0.

**Post-freeze diagnostics.** These are labelled diagnostics, not cohort results, run on the same reads
with the `b1684ae` build under `tests/results/genotype_cohort_c4_b1684ae_postfreeze_diagnostics/`:

- At state limit 1,024, ranking every candidate the proposal found (994–1,024): the truth is absent, and
  the winner is unchanged.
- At state limit 4,096 (3,972–4,024 found), ranking the top 256: the truth is absent, and the winner is
  unchanged.

**What this establishes.**

- In LZO the truth is a pair of panel templates, reachable with no switch. It is not among the ~1,000
  proposals at 1,024 states, nor among the top 256 of ~4,000 at 4,096. The fragment likelihood never saw it.
- Three of the errors are single bases in an equal-length allele. The fourth is a whole 117 bp SV allele,
  and it is certified by 0.27 nats.
- **The development ladder could not detect this.** Its rule measures stability, and a truth that is
  never proposed is stably absent. All eight development donors had the truth at rank 1.

**Not established:** why the proposal misses a zero-switch panel pair. One unverified possibility for the
single-base cases is that no syncmer discriminates the two alleles, so their marker emissions tie and
prior mass decides what survives the beam. That would not explain NA19331, whose 117 bp difference should
be visible to markers.

### R5.9 Finding 2 — a false depth-fit refusal (NA20355)

- **The fit.** It stopped at `no_ascent` after 4 iterations: rate 9.82334, background 0.0745, gradients
  0.00224 (background) and 0.00450 (rate), log-likelihood 1,002,762.31.
- **The test.** Interior stationarity requires |g_bg| ≤ 1e-9 · n = 2.5e-5 for 25,167 anchors. The
  background gradient is about 90× that, which is 9e-8 per anchor.
- **Likely cause.** At a log-likelihood near 1e6, the gain the Newton step still offers is below double
  precision, so backtracking finds no strict ascent. The tolerance, introduced in `2d5673d`, is tighter
  than the arithmetic can resolve.
- **Status.** Not yet verified or fixed. It is the cohort's only run no-call.

### R5.10 Finding 3 — phase support at K = 8 over-claims

- Ten donors report READS phase at K = 8. When the same proposal is ranked at K = 256, nine of them become
  PRIOR_ONLY, because a fragment contender outside the top 8 appears.
- Winners and block genotypes are identical at both K for all 55. Only the phase claim is affected.
- This matches the development ladder (HG04036), and is why the audit was declared.

### R5.11 Questions for the reviewer

1. **Finding 1.** Is the cohort reportable as it stands (601/605 blocks exact; excess in 4 of 55), or
   blocked until the proposal miss is understood? The proposed next diagnostic compares the truth pair and
   the winner on proposal emission and exact prior, block by block, then checks whether the truth path is
   pruned in the beam or never competitive.
2. **Finding 2.** Replace the absolute gradient tolerance with a criterion the arithmetic can meet, such as
   a relative gradient or a Newton decrement at the log-likelihood's floating-point resolution? The gate
   would reproduce NA20355's anchor counts. After the fix, is rerunning NA20355 part of the cohort, or a
   diagnostic?
3. **Finding 3.** Should READS phase require survival at a larger K, or should phase be reported only from
   a larger candidate set (K = 32 was stable on the development donors)?
4. **Near-zero certified gaps.** Two wrong calls are certified by 0.0005 and 0.27 nats, because
   certification only guards floating-point error. Should a call also require a minimum modelled gap, or
   report such calls as ties?

## R4. Review round 4 — diagnostic response

The reviewer made six rulings and set eight bounded diagnostic steps, closing with "Stop and report after
these checks. Do not implement an indel-aware scorer yet." All eight steps have been run. No likelihood,
emission or ranking code has changed. Commits:

- `a87ce76`: donor manifest;
- `2d5673d`: depth fit, plus `candidate_scope` and `global_certified`;
- `8cf4164`: diagnostic recording and `genotype_probe compare`.

At `8cf4164` all 16 enabled tests pass. No pilot or cohort run has been made.

**Every R3.3 ladder run predates the depth fix.** Their status files record `marker_depth_converged false`
after 200 iterations at rate 10.1043. The R3.3 LZO and LOO tables therefore have to be rerun before K or the
state limit is chosen (R4.6). The `compare` diagnostic in R4.3 is not affected: it scores fixed diplotypes
with the simulation's declared library intensity and does not use the marker fit or the proposal.

Provenance: the probe runs were built from a working tree whose only uncommitted changes are unrelated
inspect/cosigt files; no genotype source was modified. Final development runs will be built in a clean,
isolated worktree at a named commit.

### R4.1 Rulings

| Ruling | Applied |
|---|---|
| 1. Diagnose the 2,871-nat LOO difference first | R4.3–R4.4 |
| 2. True edlib distance stays primary: per-block unordered diploid distance, availability floor and excess, SV and context blocks separately, no-call rate, and phase metrics only under READS | The probe already reports all of these per donor except the no-call rate, which is a cross-donor aggregate for the pilot summariser. Its definition is question 3 in R4.7 |
| 3. Do not select K = 64 | Not selected |
| 4. A pruned search that completes is candidate-conditional, not rejected | `2d5673d`: `status COMPLETE`, `candidate_scope CONDITIONAL`, `global_certified false`. The certified gap is certified within the retained set; it is global only when the search was exhaustive, nothing was truncated and the ranking completed |
| 5. Resolve marker-depth non-convergence before choosing K or the state limit | R4.2 |
| 6. Eight development donors: the four archive-frozen ones plus four chosen by salted hash, manifest committed first | `a87ce76`, `tests/genotype_development_donors.tsv`. Pool: the 64 samples with both haplotypes in the C4 graph. The hash adds HG02666, HG03371, HG04036 and NA20509, which have not been simulated, run or examined. The other 56 are held out |

### R4.2 Marker-depth fit (step 1)

**Reproduced.** On the NA19240 reads the constrained maximum lies on the face background = floor:

- rate 9.769876;
- background gradient −5.44, so KKT holds.

Projected joint Newton steps kept clamping the background without making the rate stationary. They stopped
at the 200-iteration cap:

- rate 10.104, which is 3.4% high;
- relative rate gradient −3.3%;
- log-likelihood 284.8 nats below the optimum.

The caller used that estimate anyway.

**Fixed** with active-set Newton. When the background sits on its floor and still wants to fall, the rate is
maximised on that face and KKT is checked there.

- On the NA19240 reads it converges in 2 iterations to the boundary optimum.
- The estimate records its termination reason, both gradients and the log-likelihood, and the status file
  reports them. Termination reasons: `interior_stationary`, `boundary_kkt`, `closed_form`,
  `declared_background_root`, `no_ascent`, `iteration_cap`.
- A fit that does not converge now ends the run `CANDIDATE_INCOMPLETE`.

**Gates.** A deterministic negative-intercept fixture must stop on the boundary by KKT with a stationary
rate. A real-C4 gate uses the actual anchor multiplicities of the C4 panel:

| Regime (counts drawn on C4 anchors) | Rate | Background | Termination | Iterations |
|---|---|---|---|---|
| Interior: Poisson(1.5 + 9.32·2m) | 9.19987 | 1.73563 | interior_stationary | 3 |
| No background: Poisson(9.77·2m) | 9.69415 | 0.13901 | interior_stationary | 3 |
| Exact negative intercept: 22m − 4 | 9.02185 | 1e-09 (floor) | boundary_kkt | 2 |

**Reconciliation.** Rounds 1 and 2 said the fit "converged in three iterations". That was measured on the
HG00096/HG00268 simulation, whose optimum is interior (background 1.48). It was generalised without checking
NA19240, which never converged. The claim is retracted.

### R4.3 The 2,871-nat difference (steps 2–8)

NA19240 LOO, `ctx:10`:

- **A** is the K = 64 winner `6529c792`, with allele 8 on haplotype 1.
- **B** is the K = 256 winner `9ed92c5a`, with allele 11 on haplotype 1.
- Haplotype 2 is identical in A, B and the truth (252,695 bp), and it is the only other haplotype in both
  candidates. Below it is "the shared homologue".

| Check | Result |
|---|---|
| 2. A and B scored alone and together | `candidate_set_independent true` |
| Fragment terms | A −396,285.630, B −393,414.467; gap B − A = **2,871.16306226** |
| Exposure terms | A −35,877.837, B −35,877.612 |
| 3. Per-fragment deltas summed | 2,871.16306226; residual −1.36e-11 against a tolerance of 6.25e-6 |
| Fragments favouring B / A | 74 / 167 (all others identical) |
| 4. Signal under B only | 14 fragments, **+3,815.5 nats (133% of the gap)** |
| Signal under A only | 3 fragments, −817.5 nats |
| Smallest set explaining 90% | 10 fragments, 2,742.7 nats; all have signal under B only, each worth 268–276 nats |
| 6. Accelerated against index-free origins, dominant 10 | 10 origins; 0 identity mismatches, 0 mass mismatches |
| 7. Calibration, reads simulated from A | 35,882 fragments; prefers A by 3,098.3 nats |
| Calibration, reads simulated from B | 36,133 fragments; prefers B by 9,879.5 nats |
| 8. Band 3 / 4 / 5 / 6 / 7 (diagnostic only) | gap 2,856.1 / 3,103.0 / 3,103.0 / 2,871.2 / 2,871.2; the sign never changes |

**Step 7 passes in both directions.** By the reviewer's rule, this is not a scorer defect; it is
substitution-model misspecification.

**Step 5, edlib traceback.** Coordinates are block-relative, in the truth frame. Allele lengths: truth
63,448 bp, allele 8 63,450 bp, allele 11 63,447 bp. Edit distances: truth–8 is 2, truth–11 is 16, 8–11 is 16.

| Truth position | Truth → allele 8 | Truth → allele 11 | Local sequence |
|---|---|---|---|
| 6,834 / 6,854 / 6,856 | — | I4, I1, I3 | (AAAG)n repeat, 8 bp |
| 17,808 / 17,810 | — | D1, D3 | (AAAG)n repeat, 4 bp |
| 24,150 | — | X | substitution |
| 33,792 | D1 | D1 | same event in both alleles |
| 39,637 | — | D1 | end of a 25 bp poly-A, then CAGTCATTCTCC |
| 40,769 | D1 | — | end of a 25 bp poly-A, then TCTCTGCAAAAC |
| 45,808 | — | D1 | 11 bp poly-G |

In edlib's convention here, D1 means the allele carries one base the truth lacks. Both alleles carry the
same event at 33,792. Beyond it, allele 8 differs from the truth by one extra A in a poly-A at 40,769. Allele
11 differs by repeat-unit indels, one substitution and two homopolymer indels. The net length difference
says nothing about this structure.

**Dominant fragments.** All ten place on allele 11's haplotype with one mismatch, footprints between 40,401
and 40,999 in allele 11's frame, each straddling the poly-A site. None of them places on allele 8's
haplotype or on the shared homologue.

### R4.4 Why the substitution model prefers allele 11

For each site where the two alleles differ, the table follows the truth-haplotype-1 fragments that span the
site on the allele that matches the truth. It then asks where those fragments are placed in the other
candidate. The window columns search for the truth-haplotype-1 bytes within ±30 bp of the site.

| Site (truth frame) | Allele that differs | Spanning | Still placed on the differing allele | Only on shared homologue | No signal | Truth window in allele 8 hap / allele 11 hap / shared homologue |
|---|---|---|---|---|---|---|
| 6,834 (AAAG)n | 11 | 14 | 12 | 0 | 2 | 1 / 0 / 0 |
| 17,808 (AAAG)n | 11 | 22 | 22 | 0 | 0 | 1 / 1 / 1 * |
| 24,150 substitution | 11 | 28 | 28 | 0 | 0 | 1 / 0 / 1 |
| 39,637 poly-A | 11 | 16 | 6 | 9 | 1 | 1 / 0 / 1 |
| **40,769 poly-A** | **8** | **22** | **8** | **0** | **14** | **0 / 1 / 0** |
| 45,808 poly-G | 11 | 27 | 5 | 22 | 0 | 1 / 0 / 1 |

\* A repeat-unit indel can be placed anywhere in the repeat. This ±30 bp window happens to fall where every
haplotype spells the same bytes.

- **Where allele 11 differs, fragments survive.** Only 3 of 107 spanning fragments lose all signal. Hamming
  placement tolerates the (AAAG)n shifts and the substitution within the band. At 24,150, 39,637 and 45,808
  the shared homologue also spells the truth bytes and absorbs the fragments that the differing allele cannot
  hold (9 of 16, and 22 of 27).
- **Where allele 8 differs, fragments are lost.** 14 of 22 spanning fragments lose all signal. The shared
  homologue does not carry the truth bytes at 40,769, so nothing absorbs them. They fall to the background
  arm at about 270 nats each; 14 × 270 ≈ 3,800 matches the 3,815.5 nats measured.
- **The contrast between the two poly-A sites is the homologue, not homopolymer mechanics.** Both are 1-bp
  slips at the end of a 25 bp poly-A, yet they lose 1 of 16 and 14 of 22 fragments. At 39,637 the donor's other
  haplotype spells the local truth; at 40,769 it does not.
- **Not measured:** why 8 of the 22 fragments at 40,769 keep a placement on allele 8. Presumably the site
  falls close enough to a mate's end for the shifted tail to stay within the band.

**Conclusion.** The substitution-only likelihood counts fragments lost to background; edlib counts bases.
Allele 11 is 16 bases from the truth, but its differences are either Hamming-tolerable or covered by the
homologue. Allele 8 is 2 bases away, but its one uncovered indel costs 14 fragments. The two rankings
therefore disagree. The effect is specific to LOO, where the truth allele is absent and the nearest
available alleles differ from it by indels. In LZO the truth allele is present and places every fragment.

### R4.5 Corrections to round 3

- R3.3 said an indel-versus-substitution explanation "is contradicted by the lengths". This was wrong:
  net allele length does not constrain where indels lie. It is retracted; R4.3 shows the structure.
- The "three iterations" convergence claim is retracted (R4.2).
- The R3.3 ladder was run with a non-converged depth fit and must be rerun before it informs K (above).

### R4.6 Proposed next steps, pending review

1. From a clean isolated worktree at the reviewed commit, rerun NA19240 LZO and LOO with the fixed fit.
   Confirm `marker_depth_termination boundary_kkt` and record whether the proposal or the winner changes.
2. Run the two-dimensional ladder on the eight development donors, LZO and LOO. For each state limit, run
   the largest K once and evaluate its prefixes. Report the ruled metric set, with `candidate_scope` and
   `global_certified`.
3. Choose K and the state limit only if LZO is stable across settings on the development donors; LOO stays
   provisional. Then run the LZO cohort on the 56 held-out donors.

No indel-aware emission is planned unless the reviewer rules otherwise.

### R4.7 Questions for the reviewer

1. **Misspecification.** Does R4.3–R4.4 settle step 7 as substitution-model misspecification, so the LZO
   ladder can proceed with LOO reported but provisional? Or should anything else be checked first?
2. **An LOO flag.** Should the probe report, per LOO donor, whether an availability-floor diplotype was in
   the candidate set, and its fragment-term gap to the winner? That is evaluation-only, since it uses the
   truth, and it would expose this failure wherever it recurs.
3. **No-call rate.** Is a no-call a run whose status is not `COMPLETE`, a block whose `genotype_support` is
   not read-supported, or both reported separately?
4. **Ladder grid.** Proposed state limits: 256, 1,024 and 4,096 merged states, with K prefixes 8–256 from a
   single K = 256 run at each limit. Acceptable, or should the grid differ?

## R3. Review round 3 — response

The reviewer accepted the round-2 scorer and single-donor measurements, then identified one proposal
blocker and one likelihood-normalisation issue, and set a ten-step sequence. Steps 1–7 are done, each
committed with its gate. Steps 8–10 — choosing K, the development pilot and the held-out cohort — have not
been started: the rerun surfaced a leave-one-out instability that bears directly on K, and the archive does
not identify eight development donors. Section 0 and below are the round-2 and round-1 record.

Branch `genotype-mosaic`, 16 enabled tests pass. The unrelated inspect/cosigt work remains uncommitted.

### R3.1 Item by item

| Review step | Commit | Resolution | Gate |
|---|---|---|---|
| 1. Signal/background normalisation | `c8998ed` | Exposure term is `-(1-eta) * lambda * (E_a + E_b)`, the candidate-dependent part of the integrated superposed intensity; the contract formula is updated | The complete signal-plus-background intensity summed over all 4,096 observations of a 14 bp locus equals `(1-eta) * lambda * E + eta` (2.35 = 2.35 at eta = 0.25); the pre-review `lambda * E` visibly fails it |
| 2. Per-layer equivalence merging | `8ba6942` | State key is the unordered pair of (allele-prefix id, current template); log-summed mass propagated, best-history witness kept, physical multiplicity counted, pruning only after merging; prefixes interned as `intern(previous, allele)` | Exhaustive enumeration at r = 0, 0.25/0.6 and r = 1 on three- and four-template panels: log-summed mass, witness score and components, multiplicity, homologue swap; a fixture where template–prefix association changes the future; identical templates collapse into one mosaic with every history |
| 3. Proposal status reported separately | `8ba6942` | `proposal_states_pruned`, `proposal_search_exhaustive`, `proposal_target_reached`, `unique_candidates(_found)`, `candidate_set_conditional` | Beam-level: exhaustive below K is unconditional; exhaustive above K, and pruned reaching K, are conditional; pruned missing K reaches neither |
| 4. Invert the wrong-call test | `8ba6942` | A pruned search that misses K ends CANDIDATE_INCOMPLETE with a diagnostic manifest and no ranking or call | End-to-end: a starved proposal refuses with no ranking or call; an exhaustive search below K completes unconditionally |
| 5. Crowding fixture under the same budget | `8ba6942` | See R3.2 | Under 64 merged states the diagnostic manifest contains the truth; at an exhaustive budget it is recovered and read-phased |
| 6. Honest phase support | `8ba6942` | `phase_support` = READS / PRIOR_ONLY / UNRESOLVED over fragment contenders (rivals whose fragment term is not below the winner's beyond numerical error); per-block `genotype_support` | READS on the phaseable truth; UNRESOLVED on the unphaseable one; a fixture with a second carrier of one haplotype makes a prior-decided call reported PRIOR_ONLY with all five block genotypes read-supported |
| 7. Rerun NA19240 LZO and LOO | `1546df8` | Probe reports the new status; phased distance and phase switches are NA unless phase support is READS | R3.3 |

### R3.2 The crowding fixture

Tiny locus, phaseable truth, neither haplotype a panel row. It has 36 distinct diplotypes.

| Merged-state limit | 16 | 64 | 256 | unbounded |
|---|---|---|---|---|
| Unique diplotypes found | 2 | 4 | 16 | 36 (0 pruned) |
| Truth proposed | no | **yes** | yes | yes |
| K = 16 | refused | refused | complete | complete |

The same 64-slot budget that held two wrong phasings before merging now proposes the truth. With K = 16 the
run still refuses, correctly: 64 merged states hold four of the 36 diplotypes. One consequence for K: a
pruned search that reaches a small K can still omit the truth in its top-K truncation — at 64 states and
K = 2 the returned pair did not contain it. Under the reviewer's rule that call is allowed and marked
conditional.

### R3.3 NA19240 rerun

Same simulation and seed as round 2; checkpoint signature under `tests/results/genotype_probe_c4/`.
Grouped-state limit 1,024 throughout.

**Leave-zero-out** — 561 diplotypes found, 43.1M merged states pruned, target reached at every K:

| K | Truth rank (proposal / score) | Blocks exact | Distance | Phase support | Contenders | Certified gap | Proposal s | Total s | Peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 1 / 1 | 11/11 | 0 | READS | 0 | 34.7 | 2.75 | 10.7 | 347 MB |
| 16 | 1 / 1 | 11/11 | 0 | PRIOR_ONLY | 1 | 16.9 | 2.78 | 12.8 | 422 MB |
| 32 | 1 / 1 | 11/11 | 0 | PRIOR_ONLY | 1 | 16.9 | 2.69 | 15.7 | 346 MB |
| 64 | 1 / 1 | 11/11 | 0 | PRIOR_ONLY | 1 | 16.9 | 2.62 | 21.8 | 437 MB |
| 128 | 1 / 1 | 11/11 | 0 | PRIOR_ONLY | 1 | 16.9 | 2.62 | 29.4 | 453 MB |
| 256 | 1 / 1 | 11/11 | 0 | PRIOR_ONLY | 1 | 16.9 | 2.69 | 46.1 | 490 MB |

- **The round-2 blocker is gone.** At K = 8 the truth was never proposed and a wrong call was certified;
  now it is rank 1 and exact.
- The winner is byte-identical at every K. From K = 16 one fragment contender — a rival with an identical
  fragment term — is retained, so phase is honestly PRIOR_ONLY while all 11 block genotypes stay
  read-supported.
- Total time grows with K through rescoring, not proposal.

**Leave-one-out** — NA19240 removed; 129 templates. Haplotype 1 is not spellable from the retained panel:

| K | Blocks exact | SV exact | Distance | Floor | Excess | Phase support | Contenders | Total s | Peak RSS | Winner |
|---|---|---|---|---|---|---|---|---|---|---|
| 64 | 8/11 | 5/5 | 4 | 4 | **0** | PRIOR_ONLY | 3 | 23.1 | 438 MB | `6529c792` |
| 256 | 8/11 | 5/5 | 18 | 4 | **14** | PRIOR_ONLY | 3 | 49.7 | 424 MB | `9ed92c5a` |

**The winner is not stable between K = 64 and K = 256, and the larger set is worse.** What is established:

- The two calls differ only at `ctx:10`, in one haplotype: allele 8 at K = 64 (63,450 bp; edit distance 2,
  equal to the availability floor) against allele 11 at K = 256 (63,447 bp; distance 16). The truth allele
  is 63,448 bp and absent from the retained panel.
- **The reads choose it.** In the K = 256 ranking, the K = 64 winner is score rank 183. Its fragment term is
  2,871 nats worse than the new winner's; the exposure terms differ by 0.23 nats, and the prior favours the
  K = 64 winner by 8.4 nats.
- The K = 64 winner was proposal rank 2; the K = 256 winner is proposal rank 76. K = 64 was accurate only
  because the candidate the likelihood prefers was outside its top 64. That accuracy came from truncation,
  not from evidence.
- Every block genotype is reported read-supported at both K, which is correct: the reads did choose.

What is **not** established is why the fragment likelihood prefers an allele 16 edits away over one 2
edits away. A first hypothesis — that one allele differs from the truth by substitutions and keeps reads in
register while the other shifts them by an indel — is contradicted by the lengths: both alleles differ in
length from the truth (+2 and −1 bp). No explanation is offered here.

### R3.4 Questions before choosing K or starting the pilot

1. **The LOO instability.** Proposed diagnostic, not yet run: localise the 2,871-nat difference to individual
   fragments — which fragments carry it, where they place on each allele, and whether they fall outside the
   mismatch band (and so into the outlier arm) on one allele but not the other. Should that precede any K
   choice?
2. **Model versus metric.** The substitution-only truncated likelihood and the edlib accuracy metric can
   rank alleles differently, and here they rank them in opposite order. Is edit distance still the primary
   metric for the pilot, and is indel-aware emission in scope before the cohort, given it is a stated
   version-one non-goal?
3. **K.** LZO is exact and stable at every K ≥ 8. LOO is worse at the larger K. Is K = 64 still the
   provisional choice, knowing its LOO result depended on truncation?
4. **The development donors.** Step 9 names eight. The archive froze four as development donors — HG00733,
   HG02818, NA19036, NA19240 — and lists no eight. Which donors form the pilot?
5. **Small-K conditional calls.** A pruned search reaching a small K can omit the truth in its top-K
   truncation and still call, marked conditional. Acceptable, or should accuracy-eligible runs also require
   the truth-agnostic condition that the returned set is not truncated?

## 0. Review round 2 — response

The reviewer reproduced the suite and the C4 LZO result, then rejected the certification as resting on
a candidate-dependent shortcut and listed seven corrections. Every item is addressed below, each in its
own commit with its gate. Sections 1–7 are the round-1 record; where round 2 supersedes them, that is
marked in place rather than rewritten, so the reviewer can see what changed.

Branch `genotype-mosaic`, 16 enabled tests pass (two slow suites disabled as before). The working tree
also carries unrelated uncommitted inspect/cosigt work from another session; none of it is in any of
these commits.

### 0.1 Item by item

| Review item | Commit | Resolution | Gate |
|---|---|---|---|
| Blocking: candidate-dependent certification | `9c3e305` | Declared-outlier policy, omitted-tail bound and band-deepening loop removed, not patched | The reviewer's fixture: failed before (A's upper bound −42.2 alone, −186.9 with B); now A's score, error bound and disposition are bit-identical alone, alongside B, and ranked with B |
| Choose one band semantics | `9c3e305` | Truncated likelihood: mate emission conditioned on at most d substitutions, zero beyond; in-band mass is exact relative to the declared model | Truncated emission sums to one over every mate sequence at every band |
| Derive strand and exposure normalisation | `9c3e305` | 1/2 orientation prior in the mass; lambda summed over both orientations; library-normalised insert probability; exposure averaged over the observed physical-floor distribution | Summing the signal mass over all 4,096 observations on a 14 bp locus equals the exposure exactly, at a truncating and a non-truncating band; each orientation arm carries exactly half |
| Marker depth from invariant anchors | `ff7a8cc` | Poisson fit over anchors of known diploid copy number; background declared or fitted jointly; refuses when unidentifiable | Recovers known rate and background from 6,000 simulated anchors; closed form and 1-D root checked; multiplicity-blind estimator shown wrong on the same data |
| Emit all 11 blocks | `3f1c140` | `block_genotypes.tsv` carries every block with a `reportable` flag; `genotypes.tsv` stays SV-facing | End-to-end: block table row counts equal the status file's variant and context block counts |
| Aggregate sequence-equivalent mosaics | `7d60303` | Exact sequence-level prior by forward recursion over (block, offset, template); log(2) decided on bytes; candidates ranked as unordered sequence pairs | Brute force over every allele path and ordered pair grouped by spelled bytes, at r = 0, 0.3, 1; the path-based prior shown wrong on a sequence homozygote |
| Resource refusal | `3f1c140` | Declared caps on piece hits, proposed starts, physical origins, state classes, mass cells and beam extensions; `ResourceExhausted` returns nothing | Each cap throws and certification returns no ranking; end-to-end runs end CANDIDATE_INCOMPLETE or RESCORE_INCOMPLETE with no ranking and no call |
| Bypass fail-closed invariance gate | `3f1c140` | Blocks whose boundary a bypass makes representation-dependent are flagged and the model refused | A fixture whose bypassing template still round-trips is flagged at `ctx:2, sv:3, ctx:4` and refused; the control is trimmed and accepted |
| Off-panel synthetic truth | `3f1c140` | Recombinant truth, neither haplotype a panel row, checked block by block through the model frame | See 0.3 |
| Real-donor LOO probe with the LZO ladder | `7ebc9ed` | Developer probe: seeded simulator under the model's own generative process, LOO graph filter, fingerprint-gated truth evaluator | See 0.4 |

Also in round 2: a run removes earlier outputs under its prefix before starting; a structural failure
writes MODEL_INCOMPLETE instead of escaping with no status; and the evaluation projection is gated —
all 131 C4 templates project onto exactly their own model rows.

### 0.2 Corrections to round-1 claims

- **The 17.5-nat certified gap was the prior, not the reads.** Ranking by sequence pair exposed it:
  the runner-up is the same two haplotypes phased the other way across the last block, with a
  bit-identical fragment term. The status file now reports `winner_decided_by` and the fragment-term
  gap, and that run reads `prior`, gap 0.
- **"LZO cannot exercise recombination" was too absolute.** LZO exercises the implementation, but its
  truth does not require recombination, so it cannot validate off-panel recovery.
- **The candidate-set narrowness was not LZO having nothing to explore.** Round 1 first flagged
  8 unique candidates at K=64 as a possible defect, then retracted that as expected behaviour. The
  retraction was wrong: see 0.3.
- **The oracle is an independent search path, not a wholly independent derivation.** It shares the
  primitive formulas — emission, insert probability, FR geometry — with the accelerated scorer.

### 0.3 What the off-panel truth found

Both truth haplotypes recombine block alleles present in the panel; neither is a panel row.

**The proposal's raw-history width is consumed by equivalent template histories.** The marker
emission is phase-blind, so the truth and the zero-switch template phasings sharing its genotype
receive identical emissions, and among them only the prior decides; the truth needs a switch on both
haplotypes. The raw width then fills with many template histories of the same two zero-switch allele
paths before the truth is reached:

| Tiny locus, phaseable | Width 64 | 256 | 1,024 | 4,096 |
|---|---|---|---|---|
| Unique candidates | 2 | 2 | 4 | 4 |
| Truth offered | no | no | **yes, rank 1** | yes, rank 1 |

The same crowding holds at C4, which round 1 misread:

| C4 LZO, 30x | Width 256 | 1,024 | 4,096 |
|---|---|---|---|
| Unique candidates | 8 | 20 | 26 |
| Proposal seconds | 0.73 | 2.85 | 11.46 |
| Total seconds | 7.7 | 13.7 | 24.3 |

The unique target was never reached at any width. The raw-history limit is now an explicit option,
and the status file reports it, whether it was reached, and `candidate_set_truncated`. No default was
changed after seeing these numbers; one end-to-end arm pins the default-width limitation.

**Where the reads can phase it and the proposal is exhaustive, the truth is recovered.** With context
spans shorter than the insert and a width of (4×4)^5, derived from the locus's 4 templates and
5 blocks, every block is called and phased correctly and the winner is decided by fragments.

**Where the reads cannot phase it, phase is not claimed.** With every context span wider than every
insert, no fragment links neighbouring variants. The two zero-switch phasings outrank the truth by
prior alone and tie with each other, so the run ends RESCORE_INCOMPLETE; both carry the truth genotype
at every block, and the truth phasing is retained with an identical fragment term.

### 0.4 C4 LZO ladder and LOO probe

`tests/genotype_probe_c4.sh`, donor NA19240 — one of the four donors frozen as development donors in
`archive/genotype-factor-2026-09-11:tests/genotype_baseline_cohort.sh`, so the held-out cohort is
untouched. Reads simulated under the model's own generative process at 30x (35,628 fragments, 100 bp
mates, insert 350±40, 0.2% substitution error, seed 7). No cohort run has been made.

Checkpoint signature: commit `7ebc9ed`; tree dirty only with the unrelated inspect/cosigt work;
`panvar` sha256 `b2f50a2b…`, probe `1e1d40e2…`, script `c3d1d6ee…`, graph `b13af390…`. Every run below
was evaluated by rebuilding its model, reproducing its fingerprint, and projecting the donor's truth
paths onto that model's blocks.

**Leave-zero-out** (NA19240 in the panel):

| K | Raw width | Unique | Truth offered | Truth rank (proposal / score) | Blocks exact | SV exact | Distance | Decided by | Certified gap | Seconds | Peak RSS |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 64 | 2 | no | absent | 10/11 | 4/5 | 117 | fragments | 37.2 | 6.2 | 343 MB |
| 16 | 64 | 2 | no | absent | 10/11 | 4/5 | 117 | fragments | 37.2 | 6.0 | 398 MB |
| 32 | 128 | 3 | yes | 1 / 1 | 11/11 | 5/5 | 0 | fragments | 74.2 | 6.2 | 395 MB |
| 64 | 256 | 4 | yes | 1 / 1 | 11/11 | 5/5 | 0 | prior | 16.9 | 7.1 | 399 MB |
| 128 | 512 | 10 | yes | 1 / 1 | 11/11 | 5/5 | 0 | prior | 16.9 | 9.8 | 423 MB |
| 256 | 1,024 | 10 | yes | 1 / 1 | 11/11 | 5/5 | 0 | prior | 16.9 | 10.9 | 406 MB |

- **At K = 8 and 16 the truth was never proposed, and the run still certified a wrong call.** The
  winner called `sv:5` as a 117 bp / 0 bp heterozygote where the truth is homozygous 117 bp —
  117 edits — decided by fragments with a certified gap of 37.2 nats. Certification is correct relative
  to the retained set; the only warning is `candidate_set_truncated = true`.
- From K = 32 the truth is rank 1 by both proposal and score, every block is exact, and the winner's
  haplotypes are byte-identical at every K from 32 to 256.
- From K = 64 the winner is decided by the prior: a rival with a bit-identical fragment term enters the
  set, and the whole 16.9-nat gap is prior.
- The contract's stop condition — truth absent at K = 256 — is not triggered. The runtime and memory
  gates pass at every K.

**Leave-one-out** (NA19240's two paths removed from every proposal source; 129 templates, 697 alleles):

| K | Raw width | Unique | Truth spellable (hap1 / hap2) | Blocks exact | SV exact | Distance | Availability floor | Excess | Phase switches | Decided by | Seconds | Peak RSS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 64 | 256 | 9 | no / yes | 8/11 | 5/5 | 4 | 4 | **0** | 0 | fragments | 8.5 | 400 MB |
| 256 | 1,024 | 9 | no / yes | 8/11 | 5/5 | 4 | 4 | **0** | 0 | fragments | 10.5 | 427 MB |
| 1,024 | 4,096 | 9 | no / yes | 8/11 | 5/5 | 4 | 4 | **0** | 0 | fragments | 20.1 | 447 MB |

- The call reaches the availability floor exactly. Its only errors are in three context blocks — 1, 1
  and 2 edits at `ctx:2`, `ctx:6`, `ctx:10` — where haplotype 1's allele is absent from the remaining
  panel, and in each the called allele is the closest one available.
- This LOO is easy on the SV side and should be read that way: all five of NA19240's SV alleles remain
  in the retained 129 templates, so no SV allele had to be inferred off-panel. It does test mosaic
  assembly — haplotype 1 is not spellable from any single template — but only at small context
  differences.
- The winner is identical at every width, decided by fragments with a 103-nat fragment gap.
- The fitted marker background sits at its floor (1e-9): the anchors prefer no background at all.

**What this establishes, and what it does not.** One donor, one seed, 30x only, reads simulated under
exactly the declared model. It shows the pipeline recovers a panel truth whenever the proposal offers
it, and reaches the availability floor for one held-out donor. It does not show robustness to a
misspecified or real library, behaviour at lower coverage, an off-panel SV allele, or any cohort-level
rate.

### 0.5 Questions for the reviewer before any cohort run

> **Answered in review round 3.** Crowding is fixed structurally (R3.1); truncated sets now refuse; phase support is reported; the bypass refusal and declared library parameters are accepted for version one.


1. **Proposal crowding.** The raw width fills with equivalent template histories of the same allele
   paths, so unique candidates stay far below K (C4: 26 at width 4,096). A structural fix is to merge
   beam histories that share the current template pair and the allele-path prefix, by log-sum, during
   expansion: exact for the retained mass, and it frees width for distinct allele paths. That changes
   Phase 5's certified beam and needs its exhaustive oracle re-run. The alternative is to leave the beam
   as specified and fix the width from a ladder. Which?
2. **Truncated candidate sets.** LZO K = 8 and 16 produced COMPLETE, fragment-decided, certified — and
   wrong — calls because the truth was never proposed. Should `candidate_set_truncated` runs end
   CANDIDATE_INCOMPLETE and be excluded from accuracy tables, rather than COMPLETE with a flag?
3. **Prior-decided phase.** From K = 64 the LZO winner is decided by the prior with a fragment gap of 0.
   Block genotypes stay fragment-supported, but phase across blocks the inserts cannot link is chosen by
   the prior. Should the block table carry phase sets, or an unresolved-phase flag, when that happens?
4. **Representation-dependence refusal.** Version one refuses any model in which a bypass allele
   prevents the divergence trim. With internally called bubbles this did not arise at C4; with supplied
   bubbles it can. Acceptable scope for version one?
5. **Library parameters** are still declared, not estimated. Adequate for simulated donors; a real gap
   before real sequencing.
6. **Cohort settings.** Given item 1, the pilot's K and raw width decide truth-in-candidate rates. Which
   K and width should be frozen for the development pilot, and should the pilot add nested coverages now?


## 1. Verified state

Commit `d842deb` plus an uncommitted working tree. Build: Release, `-DPANVAR_BUILD_TESTING=ON`.
Both genotype test binaries pass.

| Phase | Component | Code | Gates |
|---|---|---|---|
| 1 model foundation | `genotype_model`, `ls_hmm` | committed `fefd4f9` | pass |
| 2 immutable panel/block | `genotype_model` | committed `715e8ed` | pass |
| 3 strict spelling | `mosaic_spelling` | committed `d842deb` | pass; graph checkpoint met at `51b6a49` |
| 4 marker proposal | `marker_model` | committed `c3b8c38` | pass; gate 3 met on real C4 at `51b6a49` |
| 5 deterministic top-K | `mosaic_beam` | committed `c3b8c38` | pass against exhaustive oracle |
| — graph adapter | `genotype_locus` | committed `51b6a49` | pass on real C4 |
| 6 fragment likelihood | `fragment_likelihood` | committed `98324c9`..`6bb382c` | pass; oracle-certified |
| 7 command integration | `genotype_command` | committed `7c5a42b` | pass end-to-end on C4 |

`panvar genotype` is now registered in `main.cpp`, which the contract permits from Phase 7.

The k-way merge over already-sorted allele lists is present in `build_marker_panel`
(`src/marker_model.cpp:198`), replacing the second N log N sort.

`tests/genotype_mosaic_proposal.cpp` is substantive, not decorative: it carries an independent
exhaustive enumerator over ordered template-pair histories and compares the beam against it at
r=0, r=0.25/0.6, and r=1, asserting best score, retained log-sum, emission and transition
components, and history counts separately. Non-vacuity is asserted explicitly
(`summed_histories_are_live`), and chunk-size variation stands in for scheduling nondeterminism.

## 2. The C4 block structure, measured

Panvar is an SV caller, so `panvar bubble`'s default `--min-variant-bp 50` is the correct block
alphabet: it keeps bubbles whose interior span is SV-sized. On `tests/real_data/c4.gfa.gz` that is
**5 SV bubbles** out of 1005 snarls in the complete decomposition.

The other 1000 snarls are all polymorphic, and every one of them lies in one of the six spans
*between* the SV bubbles — none is nested inside one:

```text
gap  0: 101 snarls    gap  6: 313 snarls
gap  2: 132 snarls    gap  8: 100 snarls
gap  4:  10 snarls    gap 10: 344 snarls
```

So at C4 no inter-block span is invariant, and `invariant_segments` — one shared string per gap —
cannot respell 131 templates byte-for-byte. Phase 3's round-trip gate forces the adapter to **tile**
the locus instead: the six context spans become blocks in their own right, alternating with the five
SV blocks. That yields **11 blocks**, which is the figure the contract's performance paragraph
assumes.

Measured per-block census after canonicalising orientation:

| Block | Kind | Alleles | Min bp | Max bp |
|---|---|---|---|---|
| 0 | context | 68 | 12291 | 12313 |
| 1 | SV | 68 | 277 | 7947 |
| 2 | context | 100 | 20784 | 20808 |
| 3 | SV | 16 | 203 | 506 |
| 4 | context | 10 | 274 | 324 |
| 5 | SV | 8 | 88 | 205 |
| 6 | context | 118 | 47756 | 47808 |
| 7 | SV | 119 | 26667 | 98508 |
| 8 | context | 62 | 13857 | 13862 |
| 9 | SV | 5 | 387 | 456 |
| 10 | context | 126 | 63407 | 63478 |

700 alleles over 11 blocks. The contract's "C4 contains blocks with 118--126 alleles" matches blocks
6, 7 and 10 exactly, which confirms the tiled decomposition is the one it was written against.

Two further structural facts the adapter must handle, both already anticipated by the contract:

- **0 partial frames.** All 131 panel paths traverse every one of the ten SV boundary nodes exactly
  once. No template is refused.
- **60 of 131 templates are antiparallel** to the reference chain. `mosaic_spelling`'s reverse-frame
  handling is load-bearing on real data, not a defensive branch.

130 of the 131 templates spell a distinct whole-locus sequence; the locus is 194,089--265,891 bp.

## 3. Defect found and fixed

`run_beam` stored each retained layer by moving `best`, which after `flush()`'s `best.swap(pool)`
owns the expansion pool's buffer — capacity `expansion_chunk_size`, default 262144 nodes. Layers are
retained for every block, so one chunk-sized allocation stayed alive per block. Storing exactly the
retained nodes instead, measured at C4's shape (131 templates, 11 blocks, r=0.01):

| K | Extensions | RSS before | RSS after | Wall |
|---|---|---|---|---|
| 64 | 11,000,201 | 121 MB | 26 MB | 0.08 s |
| 256 | 43,949,321 | 121 MB | 26 MB | 0.20 s |

Both test binaries remain green after the fix.

Two things worth recording. First, the contract's estimate — "131 squared ordered template states,
K=256, and ten transitions imply roughly 44 million extensions" — is exact: the measured figure is
43,949,321. The proposal stage is not a performance risk, and the stay/switch factorization the
contract permits is not needed for C4.

Second, the waste is proportional to block count, so the fix matters more elsewhere than here. On a
locus with a thousand blocks the same code reached 4.58 GB at K=64, over the Phase 8 memory gate on
its own. At C4's eleven blocks it is 121 MB, which no gate would have caught. It is fixed either way.

## 4. Steps 1 and 2: complete

Both landed at `51b6a49`, gated by `tests/genotype_locus_adapter.cpp` against the real C4 panel.

### What the adapter does

It is the only place where graph structure becomes statistical input: every later stage reads
`StructuralLocusInput` and never touches nodes, edges or walks. It tiles the locus — SV bubbles at
the production `--min-variant-bp`, interleaved with context blocks covering the spans between and
beyond them — and emits one source allele per distinct walk, leaving `prepare_genotype_model` to
collapse sequence-identical walks and retain their aliases.

Four decisions were carried over from the archived block chain
(`archive/genotype-factor-2026-09-11:src/genotype_blocks.cpp`), each of which the archive records as
a defect it had already paid for:

- bubbles contained in another bubble are dropped, because a chain is a partition and a contained
  bubble puts the same sequence in two blocks, after which marker confinement strips both bare;
- a bubble's `source`/`sink` are ordered by node position rather than taken at face value, because a
  reverse-oriented bubble stores them swapped and the neighbouring context spans would then run
  straight across its interior;
- abutting bubbles sharing one boundary node trim it from the later one, or the chain counts those
  bases twice;
- a template that does not walk a block **bypasses** it, with an explicit empty allele. Treating
  that as missing data makes a template whose deletion removed a block's own boundaries
  unrepresentable at that block, so it could never be called and never be scored.

### The canonicalisation defect, found by the gate

The node-split gate failed on the first run. Splitting every node in half left the concatenated
locus byte-identical and the census identical — 11 blocks, 700 alleles, 71 forward / 60 reverse —
but every block's sequence shifted, some by hundreds of bases, and the fingerprint with it.

The cause is that a bubble's source and sink are traversed by every path through it, so their bases
sit inside the variant block for every allele. Where the graph cuts those boundary nodes is
representation: a split moves bases between a variant block and its context neighbour.

The fix is to trim each variant block down to the span over which its alleles actually differ,
pushing the shared prefix and suffix outward into the context neighbours. The maximal trim is the
same whichever way the boundary node was cut, so both representations converge; the shared flanks
are identical for every allele, so no template's spelled locus changes; and the result is also the
honest definition of an SV allele.

The two representations move different amounts and still converge, which is the evidence that the
trim is doing the work: the original C4 graph moves 1,159 bases outward, the node-split graph moves
660 — because its context blocks already absorbed half of each cut boundary node — and both arrive
at the same 11 blocks, the same 700 allele sequences and the same fingerprint. Blocks carrying a
bypass allele have no common flank and are left alone; that is a remaining limitation, not exercised
at C4.

### Measured C4 census, after canonicalisation

| Block | Kind | Alleles | Min bp | Max bp |
|---|---|---|---|---|
| ctx:0 | context | 68 | 12326 | 12348 |
| sv:1 | SV | 68 | 0 | 7670 |
| ctx:2 | context | 100 | 20811 | 20835 |
| sv:3 | SV | 16 | 117 | 420 |
| ctx:4 | context | 10 | 325 | 375 |
| sv:5 | SV | 8 | 0 | 117 |
| ctx:6 | context | 118 | 47951 | 48003 |
| sv:7 | SV | 119 | 26346 | 98187 |
| ctx:8 | context | 62 | 14267 | 14272 |
| sv:9 | SV | 5 | 0 | 69 |
| ctx:10 | context | 126 | 63396 | 63467 |

736 source walks collapse to 700 canonical alleles. 71 forward and 60 antiparallel templates, zero
partial frames, zero bypasses, zero contained bubbles, zero abutting boundaries. The three SV blocks
with a zero-length allele are deletions: after trimming, the shortest allele of a deletion bubble is
empty, which is what it should be.

### Gate results

1. **Round trip.** All 131 templates respell byte-for-byte in their own frame, with the 60
   antiparallel ones canonicalised. Orientation is derived structurally, from the order in which a
   template visits the chain's boundary nodes — never by comparing bytes against the authoritative
   sequence, which would make this gate tautological.
2. **Fingerprint invariance.** Splitting every graph node of length >= 2 in half — changing every node
   identifier, every walk signature and every edge — leaves the fingerprint at `a58e2656e6209e0c`.
   The split model also passes its own round trip. Synthetic arms cover walk renaming, record
   reordering and sequence-identical aliases as invariants, and allele-sequence and block-kind
   changes as non-invariants.
3. **Block intervals tile** the spelled haplotype with no gap or overlap.
4. **Partial frames refuse**, proven synthetically since C4 has none.
5. **Phase 4 gate 3.** The marker panel covers exactly the 700 canonical alleles.

### Marker results on real C4

At k=31, s=15 over the prepared model: marker universe 47,960, of which 25,167 are invariant
anchors; **114 markers dropped for varying in more than one block** (0.24% of the universe); no block
is left without a discriminating marker; 694 of 700 alleles carry at least one.

That answers the `dropped_multiblock_markers` question with a measurement rather than a guess: the
C4A/C4B duplication does not cost the proposal model anything material, so the current semantics
stand.

The six alleles without a marker are three empty deletion alleles and three non-empty alleles
shorter than k, which carry no syncmer at all. Those three are invisible to the marker proposal and
can only be separated by the fragment likelihood. This is the same mechanism that makes GSTM1 hard,
now quantified at C4 and small.

One honesty note for the reviewer: `scored_alleles == total_alleles` holds by construction, because
version one has no allele cap anywhere in the code. On its own that equality proves nothing. The
load-bearing assertions are that the panel covers exactly the model's canonical allele set and the
per-allele marker coverage census above.

## 5. Phase 6 and Phase 7: complete

### The fragment likelihood

Four commits, in the order the contract names them.

`fragment_likelihood` owns the single definition of read emission, insert probability, valid-FR
geometry, exposure, fragment contribution and log-addition. `log_add_exp` had three copies in the
tree before this and now has one. The exhaustive `enumerate_placements_exhaustive` is certified
against hand arithmetic, not against a second implementation of the same idea, and then certifies
the accelerated scorer: every `(fragment, haplotype)` mass agrees to 1e-9 on planted substitutions,
tandem repeats, mixed read lengths and twelve random seeds.

The accelerated path inverts the index — pieces of the reads are indexed and each haplotype is
scanned once, rather than a position hash built per candidate. Piece keys are prefixes of the
pigeonhole chunks, which preserves the band guarantee because an exact chunk has an exact prefix.
Reads are bucketed by piece length so mixed read lengths do not lose placements. Unique haplotypes
are interned by bytes. Statistical classes are compressed by `(insert, edits, edits, orientation)`
with an exact log-multiplicity, which is the sum rather than an approximation of it; physical
origins and stored classes are counted separately.

There is no occurrence cap, placement top-K, start bin, midpoint deduplication or best-placement
mode anywhere in this path.

### Two findings the gates forced

> **Superseded in round 2 (`9c3e305`).** The declared-outlier policy described here was rejected by review as candidate-dependent and has been removed; the band is now part of an exact truncated model. See 0.1.


**Exactness is unreachable through the band.** The scorer is exact only when the band is at least
as wide as the longest mate, and at that band the pigeonhole index has no pieces left to key on. So
the upper bound is what certifies, and the reference it is certified against is the oracle on an
independent code path — not the scorer in a different mode, which would have been circular.

**One unplaceable fragment can swallow every real candidate separation.** A fragment that places on
no candidate at the working band contributes the outlier arm to the lower bound but the entire
omitted-tail bound to the upper one. Measured on C4 at 30x: exactly **one fragment in 34,415**
placed nowhere, and that single fragment produced a **244-nat interval against an 18-nat candidate
separation** — identical to four decimals across all five candidates, so unable to reorder anything,
but wide enough that no band the index can support would ever close it.

The resolution is a declared model statement, not a relaxed proof: a fragment that places on no
candidate is an outlier, so its contribution is the outlier arm exactly and its omitted tail is
zero. That closes the interval from 244 nats to 0.52. It is named, counted in the status file as
`fragments_declared_outlier`, and can be switched off for a fully unconditional bound. The
containment gate tests the bound with the policy off; a separate gate proves the policy touches only
fragments that place nowhere, and never the score itself — only the upper bound.

### The command

`panvar genotype` with the eight declared options and no others. It runs the declared execution
order, including the byte-for-byte panel round trip on real input before any statistic is computed
from those blocks. The exact Li-Stephens prior is passed *into* certification rather than added
after it, because pruning on the fragment term alone would discard candidates the prior would have
rescued. The proposal score is not passed and may never be. Context blocks are modelled but never
reported as calls. The fingerprint is keyed on the graph's reference frame, not the output prefix.

### An LZO run on real C4

This is a **leave-zero-out run**, not a precursor to one: both truth haplotypes are panel templates,
which is what LZO means. An earlier draft of this document called Phase 8 "a held-out LZO donor",
which conflated two senses of held-out. Phase 9's held-out cohort is held out from *development
inspection*, not removed from the panel; only Phase 10 removes a donor's haplotypes from the panel.
LZO always has the donor in the panel.

Simulated at 30x from `HG00096#1` and `HG00268#2`, the second one of the 60 antiparallel templates:
34,415 read pairs, 100 bp mates, insert 350±40, 0.2% substitution error, declared library
parameters.

| K | Status | Certified gap | Total | Peak RSS |
|---|---|---|---|---|
| 16 | COMPLETE | 17.02 nats | 5.0 s | — |
| 64 | COMPLETE | 17.02 nats | **6.3 s** | **544 MB** |

Against the Phase 8 gate of 90 s and 4 GB: 14x under time, 7.5x under memory. Stage split at K=64:
model 0.96 s, markers 0.95 s, proposal 0.37 s, spelling 0.00 s, rescoring 4.05 s. Rescoring
dominates, as designed for.

#### Accuracy, verified at block level

Both called haplotypes are byte-exact matches to the two truth haplotypes — 226,137 and 232,743 bp,
the second in reverse complement, correct because that template is antiparallel and the caller emits
in canonical chain orientation.

Byte equality does not by itself prove the per-block calls are right: a different combination of
block alleles could in principle concatenate to the same string. So the allele paths were compared
directly against the panel rows:

```text
truth HG00096#1  25,53,61,13,6,0,2,97,22,1,51     called hap1  25,53,61,13,6,0,2,97,22,1,51
truth HG00268#2  63,65,4,7,4,6,42,64,12,1,32      called hap2  63,65,4,7,4,6,42,64,12,1,32
```

All **22 block-allele assignments** are correct (11 blocks x 2 haplotypes), including the six context
blocks, and correctly **phased**: the alleles sit on the right haplotypes, not merely in the right
unordered pair. The five reported SV calls are `sv:1` 53/65 HET, `sv:3` 13/7 HET, `sv:5` 0/6 HET,
`sv:7` 97/64 HET, `sv:9` 1/1 HOM, each matching truth.

#### What LZO can and cannot establish

> **Corrected in round 2.** The narrowness reading below is wrong — it is raw-width crowding — and "cannot exercise recombination" is too absolute. See 0.2 and 0.3.


**LZO structurally cannot exercise recombination.** When the donor's haplotypes are panel templates,
the truth is reachable by copying two panel rows, and the Li-Stephens prior actively prefers that,
because every switch costs. The winning mosaic here is exactly two panel rows with zero switches. So
the mosaic machinery — the capability the whole design exists for — was idle throughout this run.

That reframes the candidate-set narrowness. The beam returned only 8 unique candidates at K=64,
which an earlier draft flagged as a possible defect. In LZO there is little for it to explore,
because the winner needs no switches; narrowness is the expected behaviour rather than evidence of a
problem. Whether it is a problem becomes visible only at LOO, where the truth is not copyable.

What this run does establish: the plumbing is correct end to end — adapter, markers, beam, spelling,
scorer, certification, projection — the runtime gate passes with an order of magnitude to spare, and
the proposal retains and the ranking certifies a truth that is present in the panel.

Two further caveats. The reads were simulated under exactly the declared library model, so this says
nothing about robustness to a real library. And the "donor" is a chimera: `HG00096#1` and
`HG00268#2` are haplotypes of two *different* samples, so this diploid is not one anybody carries.

## 6. Remaining work order

> **Superseded in round 2.** The Phase 8 protocol has now been run for one development donor; see 0.4 and the questions in 0.5.


### Step 1 — Phase 8 one-donor C4 go/no-go

The scenario is already exercised by the LZO run above; what is missing is the protocol, not the
setting:

- a **real donor**: one sample's two haplotypes. The run above used haplotypes of two different
  samples, which is a chimera rather than a diploid anyone carries;
- the **K ladder** 8..256, with truth-in-candidate rank at each K and winner stability between K and
  2K;
- per-block true edlib distance, score gap and candidate-set probability;
- a frozen checkpoint signature covering binary, manifest, parameters, scripts, simulator and seed.

Each run takes about six seconds, so the ladder is minutes of compute.

### Step 2 — Phases 9 and 10

The development pilot, the held-out LZO cohort, then LOO. Unchanged from the contract.

### Known limitations, carried forward

- **No LOO evidence yet.** Everything measured so far is LZO, which by construction cannot test
  whether the caller can build a haplotype absent from the panel. Until a LOO run exists, no claim
  about the mosaic capability is supported by measurement.

- **Library parameters are declared, not estimated.** `--library` is required. The contract asks for
  estimation with declaration as an option; version one has only the declaration. For a simulated
  donor this is the right input anyway, but it is a real gap before any real sequencing run.
- **`--threads` is recorded and not used.** The contract forbids scorer multithreading before the
  serial algorithm meets its performance gate. It now does, so threading is unblocked but unbuilt.
- **Blocks carrying a bypass allele are not divergence-trimmed**, so node-split invariance is not
  proven for them. C4 has no bypasses.
- **The beam's bypass-allele arm is exercised only through the adapter**, not in the beam's own
  exhaustive-oracle fixture.

## 7. Reviewer questions

> **Round 1's questions, answered by review.** The open questions are now in 0.5.


1. **Tiling.** Confirm that promoting the six polymorphic context spans to blocks is the intended
   reading of "ordered blocks and their invariant contexts". At C4 no gap is invariant, so the
   alternative — one shared context string per gap — cannot pass Phase 3's byte-for-byte gate, and
   would inject systematic false bases into every spelled candidate that the fragment likelihood
   would then have to absorb. Related: confirm context blocks are modelled but excluded from
   per-block genotype output, since they are not SV calls.
2. **Long blocks and switch resolution.** SV block 7 reaches 98.5 kb and context block 10 reaches
   63.5 kb. Inside such a block the mosaic cannot recombine at all, so the caller is effectively
   choosing a whole template over those spans. Is that the intended resolution for C4, or should
   long blocks be sub-divided at a declared coordinate? This changes the legal switch set and so is
   a model change needing its own gate.
3. **`dropped_multiblock_markers` — answered, confirmation requested.** 114 of 47,960 markers are
   dropped for varying in more than one block, and no block loses its discriminating markers. The
   C4A/C4B duplication costs the proposal model nothing material, so the current semantics should
   stand. Confirm rather than decide.
4. **Divergence trim.** The adapter canonicalises block boundaries by trimming each variant block to
   the span over which its alleles differ. This is a block-boundary definition, so by the contract's
   own rule it is a model choice and not merely an implementation detail. It is what makes the
   node-split fingerprint gate pass. Confirm it, and confirm that leaving blocks with a bypass allele
   untrimmed is acceptable for version one.
5. **Oracle provenance.** Confirm that re-deriving the Phase 6 oracle independently, rather than
   porting the archived exhaustive scorer, is the intended reading of "the archived exhaustive
   substitution-only fragment model is the oracle". It is now certified against hand arithmetic and
   in turn certifies the accelerated path on twelve random seeds.
6. **The declared-outlier policy.** This is the one place where a model statement was added to make
   certification reachable, and it is the question most worth a reviewer's attention. A fragment
   that places on no candidate at the working band is declared an outlier, closing a 244-nat
   interval to 0.52. The alternative is that every real run returns `RESCORE_INCOMPLETE`, because
   the unconditional bound cannot be closed by any band the pigeonhole index supports. Confirm the
   policy, its default, and that reporting `fragments_declared_outlier` is sufficient auditing.
7. **Candidate-set narrowness, and what LZO can measure.** The beam returns 8 unique candidates at
   K=64 on the LZO run. Because LZO structurally cannot exercise recombination — the truth is two
   panel rows and the switch prior prefers exactly that — narrowness there is expected rather than
   diagnostic. The consequence is that Phase 8's go/no-go can establish the plumbing, the runtime
   gate, and truth retention, but not the caller's central claim: that it can assemble a haplotype
   the panel does not contain. Confirm that this is the intended reading of the phase ordering, and
   whether a single LOO probe should be run alongside Phase 8 rather than waiting for Phase 10 —
   the ordering as written defers the first test of the mosaic a long way.
