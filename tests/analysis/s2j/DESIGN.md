# S2-J: reuse-aware origin-cache admission in refinement round 1 (experimental, default off)

Branch `s2-j-origin-admission` from `b47e60e`.

## The policy
Profile key `repair_round1_origin_admission`: `existing` (the default) or `reuse_only_experimental`.

**Rule.** In the production search's round 1, an origin-window entry is retained only when the round's complete plan requests its exact identity for **two or more distinct new haplotypes**. The identity is the full key: bytes, both edges and relative coordinates. The count comes from `plan_origin_cache`'s existing multiplicity (distinct haplotype hashes; held haplotypes excluded), never from raw requests, reads or aliases.

**Where the restriction applies** (`PreparedScoringContext::set_reuse_only_admission`):
- planned selection: a single-requester known window is never a candidate;
- the admission of a previously unseen window: the "size unknown until this miss" route cannot bypass it.

**What stays unchanged:**
- records are still built and used exactly;
- group-holder reuse of adjacent identical windows, hits, request order, candidates, emissions, origin reduction and arithmetic;
- the cache byte budget and every work and consumption cap.

**Scope.** `RepairSearchOptions::round1_reuse_only_admission`: an RAII scope in `run()` turns the restriction on for round index 1 only, covering the planned batch and the promotion. It restores the previous state when the round ends, however it ends. Later rounds use the existing policy. Tie branches and legs have the flag cleared, along with S2-H's live-file path.

**Reporting:**
- context: `OriginAdmissionReport` (restricted plans; planned, singleton and eligible windows; misses while restricted; singletons not retained; retained; hits and group reuses while restricted);
- status: rows `repair_admission_*` when the policy is on;
- live file: appended `admission_*` columns.

## Replay identity
`repair_setup_text` appends `;round1_admission=reuse_only_experimental` only when the policy is on, so default setups are unchanged. The replay compares the recorded setup against the effective setup with **only this field** reset to the recorded value. It records `declared_setup_override`, `repair_setup_recorded` and `repair_setup_effective`. Every other identity check (model, inputs, emissions, every other setup field) is unchanged, and the gate checks that any other difference is still refused.

## Gates
- **`panvar_genotype_round1_admission`** (block_delta focused mode):
  - long-repeat fixture: a shared window (two distinct haplotypes) is retained and reused once enumerated; single-use windows are not retained; a duplicated request does not raise multiplicity;
  - masses and scores are identical to the existing policy;
  - forced key-hash collisions, permuted requests, and 1 vs 4 threads;
  - a two-round multi-block fixture: round-1 single-use windows are revisited and enumerated again in round 2; the extra executed origins equal the re-enumerated windows' origins exactly; round 2 re-admits singletons.
- **`test_round1_admission_caller`** (focused e2e):
  - outputs and trajectory are identical to the existing policy, with explicit allowances: the policy rows, six refinement work and margin rows, timing and resources;
  - exactly one restricted plan;
  - the replay accepts only the declared override and reproduces round 1.

## Runs (`s2j_run.py`; manifest frozen before execution)
1. The gates.
2. One ACOT 30× first-round replay with the policy.
3. C4 and ACOT 5× existing vs reuse-only full callers, compared with `s2j_compare.py` (explicit allowances).
4. **Conditional:** one integrated ACOT 30× full caller with the policy (S2-A, blocked/fine union route, unchanged caps, 4 GiB guard), with live refinement memory for round boundaries.
5. **Conditional, declared by Sol before execution:** if that caller stops first on the 4 GiB footprint guard, one 6 GiB-envelope follow-up (`footprint_limit_bytes 6442450944`), everything else unchanged. It is never reported as a 4 GiB result.
