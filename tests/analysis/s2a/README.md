# S2-A: early evaluated incumbent and safe budget-stop reporting

**What it is:** opt-in experimental reporting, enabled with the profile key `early_incumbent experimental_v1`. It requires `sequence_first_v2` and the adaptive repair. The likelihood, prior, proposal rules, default caps, index format and legacy output are unchanged. The contract is documented in `docs/modules/genotype.md`, under "Early evaluated incumbent".

## Selection and evaluation

**Which candidate is evaluated first.** After the markers, every panel pair is examined in the proposal emission table: `sum_b at(a1_b, a2_b)` plus the exact path prior with log 2 multiplicity. The pairs are not fragment-scored. The argmax is taken, with ties going to panel order.

**How it is evaluated.** The chosen pair is evaluated once, through the union's own route (indexed or reference) and caps, with the ranking's sequence prior.

**What it costs.** The work and time before the first prediction are reported as:
- `early_incumbent_seconds_preparation`: model, reads, library, markers and emissions;
- `early_incumbent_seconds_selection`;
- `early_incumbent_seconds_evaluation`.

**How it is checked.** The union pass must reproduce the early candidate's score bit for bit.

**Reuse.** The early evidence is reused at the score level only: its placements are not handed to the union pass. Its cost is reported beside the union's and never subtracted.

## Checkpoints

| stage | when |
|---|---|
| `early_panel_incumbent` | after the early evaluation |
| `union_best` | after the union pass, with the group's tied groups, projection ambiguity and explored labellings |
| `repair_pool_best` | at a round commit, when a repair-scored candidate beats the incumbent |

**How a commit is written.** Each commit is a v2 PROVISIONAL report, written under the watch mutex: content files first, then the locus record.

## Stops

**Watchdog.** After a checkpoint, the watchdog rewrites only the pre-rendered locus record (completion `budget_stop`) and exits 3. Before a checkpoint, it writes a typed NO_PREDICTION and exits 124 or 125.

**Declared caps.** These publish the incumbent.

**Failures.** These are explicit: a union mismatch, a failed checkpoint publication (exit 2), or an exception, which removes the checkpoint.

**Test hook.** `test_interrupt_point` takes `before_first_evaluation`, `after_checkpoint`, `union_scoring` or `adaptive_refinement`.

## Tools

- **`s2a_fixtures.py`:** F1 (equivalent labellings), F5 (genuine tie), F7full and F7, each run off, on and stopped.
- **`s2a_run.py`:** the frozen manifest runner.
- **`s2a_summary.py`:** the real-data summary table.
- **e2e:** `test_early_incumbent` in `tests/genotype_command_e2e.cpp`.
