# S2-G: complete-reference validation and narrow caller integration (experimental)

Branch `s2-g-blocked-integration` from `2062867` (S2-F). Runner `s2g_run.py`; Part A was frozen before any execution, Part B
after the gates were developed and before its execution.

## Part A — one complete ACOT 30× reference (done before integration)

Ordinary blocking (S2-E route) with the two validation-only extensions declared by Sol: physical-origin cap 1 B and
executed state-class cap 1 B; everything else unchanged (4 GiB guard, 300 s timeout, 4 threads). All 108,939 candidate
score components and the numerical top set compared bit for bit with the S2-F combined run (`blocked_fine.acot30`,
`3ab2a11`; the kernel source was unchanged at the time of the reference run). The ACOT inputs carry no audit libraries, so
there are no audit views to compare at 30× (audit views are covered by the S2-F kernel gates).

## The integration

- Profile key `union_fragment_blocks` = `off` (default) | `fine_core_4096_experimental`. It selects the existing indexed
  scorer through `score_diplotypes_indexed_blocked(..., fine_origin_aggregates = true)`: fine views, exact store,
  4,096-pair blocks, fine-core and whole-segment aggregates under the existing half/half split of
  `union_origin_cache_bytes`, every cap cumulative. Its scores and audit views replace the bulk route's in the unchanged
  certification, pooling and downstream decision code. The S2-A early pass is unchanged (bulk, one pair).
- Contract (`union_fragment_blocks_refusal`, checked in `run_genotype_body` before the prefix and in the CLI validation):
  needs the indexed backend, the exact store, `fine_views_experimental`, origin reuse and both origin consumption caps;
  refused with `union_progressive_schedule`, `adaptive_emission_handoff`, `diagnostic_union_only`, the evidence dump, the
  work trace, the seeding profile and the verification census (all need the full mass matrix). The kernel dump (inputs
  only) stays allowed. Progressive scheduling, persistent emissions and the emission handoff stay off.
- Refusal: no union score is published; the existing `ResourceExhausted` path runs (the S2-A checkpoint, when on, stays
  PROVISIONAL). The status reports completed blocks (`union_block<k>_*`, `union_completed_blocks_*`), the refused block
  (`union_refused_block_*`; direct/served breakdown `unavailable`) and the whole attempt (`union_attempt_*`).
- Provenance: `union_route`, `union_fragment_block_pairs`, `union_origin_aggregate_policy`, the budget split, every cap
  actually in force (`union_cap_*`, cumulative across blocks), the settings row `union_fragment_blocks`.
- Reporting completeness (kernel, reporting only): the blocked driver now fills the run-level origin-cache, prefilter and
  start-sort reports (counts summed over blocks; held bytes = the largest block), instead of leaving misleading zeros.
- Footprint checkpoints (route on only; existing `physical_footprint()`): union completion (or union refusal), candidate
  pool constructed, refinement context constructed, initial incumbent placed, first committed refinement round, status
  written (completion or refusal); also in a stopped run's minimal status.
- Test hook (API only): `test_union_fragment_block_pairs` (0 = 4,096), recorded as `union_fragment_block_pairs`.

## Gates

`panvar_genotype_blocked_union_route` (focused e2e; also in the full e2e): ordinary vs combined complete runs at 4,096 and
37-pair blocks; 1 vs 4 threads; CLI profile route; prepared-locus route (cold/warm/ordinary); default off; a later-block
refusal keeping the early checkpoint PROVISIONAL with distinct completed/refused/attempt work; nine incompatible
configurations refused with the prefix untouched; prefix reuse. Comparison uses explicit allowances only (named keys or
named families: route-only provenance, timing, process resources, executed/served work; thread rows only when threads
differ; profile identity only on the CLI comparison). The S2-F kernel gates and the full e2e are re-run.

## Caller runs (Part B)

C4 and CYP2D6 5× control/integrated with `s2g_compare.py` (explicit allowances, calibrated on one counted C4 development
pair before freezing); one integrated ACOT 30× attempt with the original main-arm caps and unchanged refinement (before
ACOT 5×); ACOT 5× if budget permits. Outcome classes: computational completion, scientific ambiguity (completed search
with a tie or set-level call), resource refusal (cap or guard), implementation failure (exception/FAILED).
No accuracy claim; no further caps or variants.
