# M3-G result: the `three_switches_1` 2× ctx:0 error is a model preference, not a proposal or search failure

**Run.** Branch `m3-g-proposal-trace`: freeze `453774d`, M3-F qualifications `eb4aa1b`. No caller change. Results in
`tests/results/m3g_proposal_trace_20261019/`.

**Budget.** 4.5 of 600 s. Four units, all exit 0 and first-time, with no amendments.

**Checks, all passed:**
- **C1:** the probe reproduces all 6 recorded move ranks.
- **C2:** the ctx:0 pairs scored in each of the 7 saved rounds equal the probe's top 32.
- **C3:** all 5 saved totals reproduce, which proves the probe scores as the caller did.

## The route, stage by stage (saved evidence first)

| stage | finding | source |
|---|---|---|
| Where the call came from | The call (4,52) is the **union** pool's winner: a 2-member phase-equivalence set W1/W2 at −31870.9194. The repair converged after 7 rounds to R, ctx:0 = (40,52), at −31871.3999, and lost the pooled decision by 0.48. | saved status, equivalence set, rounds |
| Distance to truth | The truth path is R with ctx:0 set to (14,56) in assignment 0. That is **one ordinary single-block move from R**. | truth.tsv vs saved R |
| 1. Proposal rank of (14,56) as an unordered pair | **90th of 2,346** at ctx:0: 9.60 below the top emission and 6.22 below the 32nd. (4,52) is exactly 32nd, tied in emission with 33rd (6,52) and kept by index order. The repair's own pair (40,52) is 723rd; it came in from the anchor. | `rank` (repair rule; C1/C2) |
| 2. In the local-pair shortlist? | **No.** The list is the fixed top 32 by emission, built once. | `rank`, code |
| 3. Generated / suppressed / rejected / scored / retained? | **Never generated** in any of the 7 rounds; a pair outside the list cannot be. Each round scored all 32 listed pairs, both assignments (64 neighbours). (4,52) was scored every round and lost, by 0.613 in the final round. The only pool members carrying (14,56) at ctx:0 are 2 panel diplotypes, at −33605 and −33059. | `saved.trace`; C2 |
| 4. If scored, why would it lose? | Under the unchanged model it **loses in every context**, below. | `score` |

## Diagnostic scoring (truth-guided, post hoc; never fed back)

The ctx:0 alternatives inside the repair incumbent R. The same differences hold inside W1 and W2.

| ctx:0 pair | Δfragment | Δexposure | Δlikelihood | Δprior | Δtotal vs R (40,52) |
|---|---|---|---|---|---|
| truth (14,56), either assignment | +4.889 | +0.020 | **+4.909** | **−17.018** | **−12.109** |
| called (4,52), assignment 0 | +6.221 | +0.020 | +6.241 | −6.855 | −0.613 |
| called (4,52), assignment 1 | +6.221 | +0.020 | +6.241 | −15.226 | −8.985 |

Against the call (W1 or W2, ctx:0 (4,52)), the truth has Δfragment **−1.332**, Δexposure 0, Δprior −10.164 and Δtotal
−11.496.

| ordering at ctx:0 | first | second | third |
|---|---|---|---|
| likelihood (F + E) | (4,52) | (14,56), −1.33 | (40,52), −4.91 more |
| prior | (40,52) | (4,52) | (14,56) |
| total | (40,52) | (4,52) | (14,56), about −11.5 to −12.1 |

- **Sequence identity:** the truth-substituted candidates are distinct sequences. R_t0 is the truth pair.
- **Reachability:** every R_* candidate is one existing single-block move from R.
- **Nothing beats R or the call** among the truth-substituted candidates.

## Decision

**Frozen class: PRIOR_REVERSES**, with no context dependence; the ordering is identical in R, W1 and W2.
- **Against the refinement incumbent R:** the likelihood prefers the truth by +4.91, and the sequence prior reverses
  that by −17.02.
- **Against the reported call:** the truth already loses on likelihood, by 1.33. This is weak 2× fragment evidence for
  the wrong pair (4,52).
- **So the error does not depend on the prior alone.** With the prior removed, (4,52) would still beat the truth.
- **The shortlist omission is real but does not bind.** (14,56) ranks 90th and is never generated. Had it been
  proposed, it would lose by 12.1 to R and by 11.5 to the call.
- **Coverage-limited.** At 5×, the twin calls ctx:0 = (14,56) DATA_SUPPORTED with gap 6.68. At 2×, support v2 already
  labels this block PRIOR_DEPENDENT / LOW_QUALITY.

**Not the cause:**
- search bookkeeping (no suppression involved);
- phase or context dependence;
- the joint catalogue's shape.

**Separate observation** (diagnostic; it changes no prediction). The union winner is never refined.
- Moving ctx:0 from (4,52) to (40,52) inside W1/W2 (W*_r0) scores **+0.613 above the reported call**.
- The pooled decision therefore sits off a single-block optimum.
- That move goes to (40,52), which is also wrong, so it would not fix this cell.

## Recommended next experiment (one; not implemented)

**Per-fragment and prior decomposition of the single comparison truth (14,56) vs called (4,52) at ctx:0, at 2× and at
the 5× twin.** Use the existing `score-set --contrast` per-fragment deltas and the M3-A prior trace (switches and
carriers behind the −10.2 to −17.0 prior).

| outcome | conclusion |
|---|---|
| The 1.33 likelihood deficit comes from a few weak reads and disappears at 5×, and the prior gap is the simulated mosaic's switches | Expected behaviour of a copying prior at low depth. The LOW_QUALITY label is the correct output; no algorithm change. |
| Reads that should be neutral or favour the truth (for example paralogous placements) favour allele 4 | A fragment-model question. Only then is a likelihood experiment justified. |

**Do not change the prior first.** Removing it would still leave (4,52) ahead.

**This does not bear on C4 sv:7.** The NA19650 sv:7 allocation error is a different mechanism (M1-B/M1-C):
- panel-supported 2+2 alternatives were evaluated (about 5–6k);
- the **fragment term** decides for 3+1 (271 at 10×), on content the LOO panel cannot spell;
- the prior and the proposal are not the cause.

Solving this small context error does not address that allocation problem.

## Qualifications appended at M3-H (2026-10-19; Sol's review; the record above is not rewritten)

1. **No global optimum is established.** The truth correction loses in every *tested* context (R, W1, W2): fixed
   substitutions in three contexts. That does not show a global optimum, or that no other context exists in which
   (14,56) wins.
2. **"Not binding" is specific to what was tested.** The shortlist omission (rank 90) does not bind for *these tested
   corrections*. It is not a general statement about the shortlist.
3. **The unrefined union winner is a separate search limitation.** Convergence of the refinement trunk does not
   establish local optimality of the final pooled call: a one-block move inside W1/W2 scores +0.613 above it. That move
   is still wrong at ctx:0, so it is recorded separately and does not explain this error.
4. **One coverage pair does not establish calibration.** The 5× twin calling (14,56) DATA_SUPPORTED is consistent with
   low-depth uncertainty at 2×.
