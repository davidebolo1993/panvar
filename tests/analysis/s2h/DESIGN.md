# S2-H: attribute refinement memory; a first-round replay (measurement only)

Branch `s2-h-refinement-memory` from `81a2345`. No change to the refinement's candidates, order, batching, caps, scores or decisions.

## Checkpoint correction

The S2-G checkpoint after `search->start()` is renamed `repair_search_started`. `start()` places nothing. `run()` reports the seed's actual placement and census through a new `RepairHooks::on_checkpoint`: `initial_incumbent_placed` and `initial_incumbent_census_retained`. The S2-G record carries an appended qualification.

## Live first-round memory (profile key `diagnostic_refinement_memory 1`)

The production search writes `<prefix>.refinement_memory_live.tsv`. Each row is opened, appended and closed, so it survives a watchdog `_Exit`. Rows are written only from safe execution points on the search's own thread; the watchdog traverses nothing.

**Points:**
- in `run()`:
  - round entry;
  - before and after seed placement (`add_haplotypes`);
  - before and after census retention;
  - before and after neighbourhood construction;
  - before and after candidate combination;
  - promotion, release, round commit;
- through `PreparedScoringContext::set_progress_probe`:
  - in `add_repaired_haplotypes`: before and after repair-window planning, after origin-cache planning, neighbour materialisation after requests 1, n/4, n/2 and 3n/4, and after the last;
  - in `score_candidates`: combination start, combination scores done, flattened mass matrix built.

**Columns:**
- process footprint (current and peak) and RSS (current and peak);
- the context's `ContextMemoryAccount`, live versus capacity:
  - slots, sequences, segments, mass columns;
  - censuses;
  - origin cache, merge buffers, origin demand;
  - window and origin scratch;
  - emission table and arena;
  - prepared state;
  - the last scoring result;
- emission-cache occupancy: entries, entry metadata, content arena;
- requests done out of total;
- haplotypes created and unique;
- cumulative physical origins, state classes and piece hits.

## First-round replay (profile key `diagnostic_refinement_replay <inputs>`)

**The full run records the inputs.** With `diagnostic_refinement_memory` on, the run writes `<prefix>.refinement_replay_inputs.tsv` when the anchor is selected:
- the model fingerprint;
- every input SHA-256;
- a digest of the proposal emissions;
- the refinement setup as text, built by the shared `make_repair_setup`;
- the anchor's labels, its exact score (hex float) and its union rank.

**The replay** runs the production preparation in a fresh process (model, markers, proposal emissions, union candidate spelling).
- It checks every identity in the inputs file, and refuses on any mismatch.
- It then runs the production `CoordinateRepairSearch` from the recorded anchor, with the same setup and `max_rounds = 1`.
- It does not run the union pass or the early pass.
- **Outputs:** `<prefix>.refinement_replay{,_rounds,_scored}.tsv` and the status `REFINEMENT_REPLAY_DIAGNOSTIC_ONLY`. It never produces a genotype.
- **Contract:** `output_schema legacy`, `early_incumbent off`, repair rounds > 0.

The production repair setup now comes from the same `make_repair_setup`. That refactor is mechanical, and the same assignments are made.

## Gates and runs (`s2h_run.py`; manifest frozen before execution)

**Focused caller gate** (`panvar_genotype_blocked_union_route`, which now also runs `test_refinement_diagnostics`):
- instrumentation leaves every output and the trajectory unchanged, with explicit allowances;
- every round-1 point is recorded;
- the rows survive an induced stop inside round 1;
- the replay, run in a separate process, reproduces the anchor's exact score and round 1's candidates, scores and decision;
- a tampered identity is refused;
- the contract is checked before the prefix is touched.

**Regressions:** the context and repair ctests, and the full e2e.

**Real data:** one instrumented ACOT 30× caller (the S2-G integrated arm with the diagnostic key), then one first-round replay from that run's recorded anchor. The original caps and the 4 GiB guard apply to both, and neither is repeated.
