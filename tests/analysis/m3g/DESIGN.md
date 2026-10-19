# M3-G: trace the ordinary single-block proposal failure at `three_switches_1` 2×, ctx:0 (diagnostic only)

**Ancestry.** Branch `m3-g-proposal-trace` from `30be3fd` (M3-F), plus the M3-F qualifications (`eb4aa1b`).

**No caller change.** The binary is the M3-F one, `a728eb71…`, from `f31b9e6` sources. Joint search stays default-off
and is not used here.

**Case.** The saved M3-F **off** run (`run.ts1_2.off`).
- **Call:** CONDITIONAL_GENOTYPE_EQUIVALENT, wrong only at ctx:0 (called (4,52), truth (14,56), excess 8).
- **Equivalence set:** two members (W1, W2), both from the union pool, scoring −31870.9193977016.
- **Repair:** converged after 7 rounds (6 moves) to R = `40,30,40,4,4,6,84,47,35,1,67 | 52,46,19,7,4,5,52,67,3,1,0`,
  scoring −31871.3998675. R loses the pooled decision to W by 0.48.
- **Truth path:** `14,30,40,4,4,6,84,47,35,1,67 | 56,46,19,7,4,5,52,67,3,1,0`. It equals R with ctx:0 set to (14,56)
  in assignment 0, so it is one ordinary single-block move from R.

## Stage 1: saved evidence (`saved.trace`, no rescoring)

For each repair round, from the saved pool (`source repair`, `scored_in round<r>`):
- every scored neighbour that differs from the round's start incumbent at ctx:0 only, with its unordered pair, assignment
  and saved total;
- whether (14,56) or (4,52) was scored;
- the round's aggregate counters: generated, scored, suppressed_visited, duplicate, self_skipped, rejected_infinite_prior.

Also, every pool candidate carrying (14,56) at ctx:0.

The pool records only scored identities. A pair that was generated but suppressed, a duplicate, or prior-rejected is
visible only through the aggregate counters. Stage 2 resolves this from the ranking rule; it is not inferred.

## Stage 2: proposal rank (`rank`, truth-free)

`score-set --rank-blocks` ranks every unordered allele pair at a block exactly as `CoordinateRepairSearch` builds its
fixed local-pair list: proposal emission descending, then indices, top `repair_local_pairs` = 32. Blocks ctx:0, ctx:2,
ctx:4, sv:3 and sv:5 are ranked; the last four are cross-checks.

**Cross-checks** (run in `analyse`; if either fails, the class is UNRESOLVED):
- **C1:** every saved move's `local_pair_rank` equals the probe's rank of that pair at that block.
- **C2:** the set of ctx:0 pairs scored in every saved round equals the probe's top 32 at ctx:0.

## Stage 3: diagnostic scoring (`score`, truth-guided post hoc; never fed back)

`score-set` scores the frozen list `candidates.tsv`: 15 path candidates, at most 30 unique haplotypes (preflight), with
the unchanged fragment model, library, band, exposure and sequence prior (r = 0.01).

| context | candidates |
|---|---|
| **R** (repair final incumbent) | R; ctx:0 → (14,56) in both assignments (R_t0 is the truth path); ctx:0 → the called (4,52) in both assignments (saved round-7 neighbours) |
| **W1**, **W2** (the saved equivalence set, which is the call; complete, 2 members) | each member; ctx:0 → (14,56) in both assignments; ctx:0 → the repair's (40,52) in both assignments |

**Identity.** Five saved totals are declared as expectations (R, R_c0, R_c1, W1, W2). The probe refuses the whole report
unless they reproduce within 1e-6 (check C3).

**Scope limitation.** These are fixed-context substitutions in three contexts. This is not a search, and it makes no
exhaustive claim.

## Stage 4: decision (`analyse_m3g.py`)

**Reported per comparison:**
- fragment, exposure, prior and total differences;
- unordered block genotype identity and sequence identity (interned haplotype ids);
- whether it beats R, or the call;
- whether it is one existing single-block move from R.

**Primary class, in precedence:**
1. PROPOSAL_SHORTLIST_OMISSION: rank above 32, never scored, and the best target assignment beats R on total.
2. SEARCH_BOOKKEEPING: rank 32 or better, yet never scored.
3. PRIOR_REVERSES: wins on F+E, loses on total.
4. EVALUATED_LOSES.
5. UNRESOLVED.

**CONTEXT_DEPENDENT** is added when the preference for (14,56), or for (4,52) over (40,52), flips between the R and W
contexts. The analysis also reports whether the best target candidate would outscore the saved pooled winner.

## Bounds

- 600 s counting every attempt, 4 threads, existing resource settings.
- **Declared timeouts:** `saved.trace` 30 s, `rank` 120 s, `score` 240 s, `analyse` 30 s.
- No caller change, simulation, cap or support-rule change, residual editing or clustering, and no full joint catalogue.
- Inputs are hashed in INPUTS.tsv and verified before every unit.
