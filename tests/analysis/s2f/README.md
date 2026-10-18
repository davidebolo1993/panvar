# S2-F: fragment blocking with the I4-E fine-core aggregates (experimental kernel, measurement only)

- **What it is:** `score_diplotypes_indexed_blocked(..., fine_origin_aggregates = true)`. Each 4,096-pair block's batch builds and owns its own whole-segment aggregates and fine-core aggregates. The unchanged `union_origin_cache_bytes` is split half/half between them, and nothing is reused across blocks. A refusal records the refused block's observed work before its batch is destroyed: `completed_totals`, `refused_block` and `attempt_totals`.
- **Design:** `DESIGN.md`, frozen before execution.
- **Gate:** `tests/genotype_blocked_fine.cpp` (ctest `panvar_genotype_blocked_fine`).
- **Harness:** `s2c_progressive_harness blocked_fine`; `blocked` takes an optional declared origin cap.
- **Runner:** `s2f_run.py`.

**Outcome:**
- Bit-identical at 5× on C4, CYP2D6 and ACOT. Executed origins fall 3.8–4.2×, runtime is flat, and the footprint peak rises 5–7 %.
- ACOT 30× with unchanged caps completes in 117.7 s at a 1.26 GB peak, with 164.5 M executed origins of the 500 M cap.
- The extended-resource control (ordinary blocking, 1 B origin limit) refuses in block 9 on the 600 M state-class cap.

Results: `tests/results/s2f_blocked_fine_aggregates_20261018/README.md`.
