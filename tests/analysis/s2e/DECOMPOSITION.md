# S2-E: why fragment-blocked union scoring reproduces bulk scoring bit for bit

**The target:** the implemented combination, `combine_candidate_scores` / `combine_with_audits`. For a candidate {a, b}:

```
exposure_term   = main ? -(main.exposure_scale * (E_a + E_b)) : -expected_signal_fragments(E_a + E_b)
absolute        = |exposure_term|
fragment_term   = 0
for f = 0 .. N-1 (in order):
    combined    = log_add_exp(M[a][f], M[b][f])
    if combined != -inf: ++fragments_with_mass
    t           = main ? prepared_contribution(main, combined, bg[f]) : log_fragment_contribution(combined, bg[f], lib)
    fragment_term += t ; absolute += |t|
    (each audit view k: the same with its own constants, its own exposure_term_k and absolute_k)
log_likelihood  = exposure_term + fragment_term
```

## 1. Which quantities are fragment-local

**M[h][f]** depends only on fragment f's reads, haplotype h, the library, the band and the scorer options. In the indexed scorer, f's forms are proposed from piece hits and verified against unit windows. Its paired origins are then formed per fragment (`pair_fragment`) and reduced by `reduce_states` or `reduce_sorted`, which are equal by I1.2's gates whether or not zones serve it.

**What a batch adds is reuse, not content.** A batch's piece index, prefilter cache, emission cache and origin aggregates are reuse structures over its own fragments. They change executed work, never a placement's emission or a fragment's state multiset. So a block of fragments yields the same M[·][f] as the whole run.

**This is checked, not assumed:** per-fragment masses are compared on fixtures.

## 2. Which quantities are whole-run, and are computed once from the whole run

| quantity | why it is whole-run | how the blocked route computes it |
|---|---|---|
| E_h | `exposure_for_length(|h|, N, fragments_by_floor, insert)` weights the observed floor distribution over ALL N fragments. A block's batch would weight its own floors. | once, from every fragment's floor; the batch's per-block exposure is ignored |
| bg[f] | `log_background_probability` | per fragment |
| the prefilter-cache single-read-length condition | the bulk adapter refuses mixed read lengths under that cache. A block could be uniform when the run is not. | checked over the whole run before any block |
| the mass-cell cap | bulk refuses when unique × N > max_mass_cells | checked with the whole run's N before any block, so splitting never admits what bulk refuses |
| the sequence prior | applied by the ranking outside the union scorer | untouched |

## 3. Association

Each candidate's `fragment_term`, `absolute` and `fragments_with_mass` (and the audit views' term and absolute) are carried across blocks. The next block's fragments are added to the same accumulator, in the original fragment order. That is the bulk loop's exact operation sequence; the scores are never block subtotals added afterwards.

The exposure term is computed once and seeds `absolute`, exactly as in bulk. `log_likelihood = exposure_term + fragment_term` is formed once, at the end. The prepared/scalar choice (`main`) and the audit constants are the bulk route's.

**Threads.** Accumulators are per candidate and fragments are walked in order inside a candidate, so the thread count changes only which worker owns which candidates.

## 4. Caps

- **Cumulative.** Executed-work caps and consumption caps are seeded into each block's batch from the running totals: piece hits in flight; physical origins in flight; acceptance on the merged executed counters; cached placements, origins, states and prefilters consumed. A cap therefore refuses exactly when the cumulative total crosses it.
- **Per block.** Memory budgets — the placement store, the origin-aggregate budget, the emission and prefilter caches — govern simultaneously held storage. They apply to each block's batch, which is destroyed before the next exists.
- **On refusal.** A refusal in any block ends the union with no ranking: no accumulated partial score is returned.
