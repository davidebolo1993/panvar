# S2-H: refinement memory attribution and a first-round replay (diagnostic only)

- **Live memory file:** profile key `diagnostic_refinement_memory 1`. It writes `<prefix>.refinement_memory_live.tsv`, with rows from safe execution points, and records the replay inputs.
- **Replay:** profile key `diagnostic_refinement_replay <inputs>`. A fresh-process first-round replay from the recorded anchor, with every identity checked; diagnostic only.
- **Checkpoint fix:** `repair_search_started` replaces the misleading S2-G `initial_incumbent_placed`. The actual seed placement and census are reported from `run()`.
- **Gate:** `test_refinement_diagnostics`. **Runner:** `s2h_run.py`.

**Outcome (ACOT 30×):** the full caller and the replay both cross 4 GiB inside round 1. The replay crosses it from a 1.25 GB start, so refinement-owned growth is sufficient on its own.
- Owners: origin cache at its 768 MiB budget, emission cache 462 MB, mass columns, scratch, censuses.
- About 1.08 GB of growth is outside every account (602 MB inside seed placement).

Results: `tests/results/s2h_refinement_memory_20261018/README.md`.
