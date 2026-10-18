# S2-C: demand-driven execution of the existing indexed scorer (experimental, measurement only)

**What it is:** a selective route on `ContentIndexBatch`, `discover_for(slots)` and `score_selected(candidates)`, plus the driver `score_diplotypes_indexed_progressive`. It uses the existing fine discovery, byte-verified sharing, pairing, reduction and caps, charged cumulatively across calls.

**What it never does:**
- discover units the requested slots don't name;
- build their columns or entity aggregates;
- export the mass matrix;
- score a slot whose evidence is incomplete — that is refused, because unknown evidence is never zero mass.

**The bulk path** (`discover`/`score`) is the control. Its code is shared and unchanged; a batch is one or the other.

**Tools:**
- the gate: `tests/genotype_progressive_kernel.cpp` (ctest `panvar_genotype_progressive_kernel`);
- the measurement: `tests/s2c_progressive_harness.cpp`, with bulk, progressive (under a footprint guard) and compare modes;
- the inputs: the diagnostic profile key `diagnostic_union_kernel_dump`, which writes the caller's exact union inputs;
- the runner: `s2c_run.py`.

**The schedule** (declared in the harness):
1. the S2-A early pair;
2. the next 31 pairs by the truth-free proposal score (S2-A's formula, maximised over each pair's labellings);
3. the rest of the pool.

Results: `tests/results/s2c_demand_driven_20261018/README.md`.
