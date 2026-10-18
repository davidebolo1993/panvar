# S2-F: fragment blocking with the I4-E fine-core aggregates (experimental kernel, measurement only)

Frozen before execution. Starting point: `f1ba915` (S2-E). This is a kernel experiment: no caller integration and no refinement change.

## What is ported (minimum, from I4-E `17f0e8f`, geometry from I4-D.1 `4715eab`)

- **Geometry.**
  - `CoreWindowGeometry`, `core_window_geometry_violation`, the per-segment core builder, and `core_window_geometry_of_segment`.
  - Not ported: the haplotype geometry, the window key, the legacy construction, and the census.
- **Aggregation.** In `ContentIndexBatch::score_impl`, bulk batches only (a selective batch refuses):
  - fine windows on every segment not served whole, keyed exactly (an FNV hash selects, bytes decide);
  - eligible with two or more occurrences;
  - canonical order: occurrences descending, then shorter window, then entity and core;
  - built from fine views by `build_origin_entry`, in batches of 16;
  - admitted as a prefix;
  - fine zones translated to each occurrence;
  - a per-fragment counting-sort merge index.
- **Accounting.** Construction is executed work, charged to the physical-origin cap. Served origins and states charge the consumption caps. The I4-E report fields are kept unchanged.
- **Test hook.** `set_fine_window_hash_collapse`.
- **Bulk adapter.** `score_diplotypes_indexed_union(..., fine_origin_aggregates)`, used as the gate's bulk-fine route.

## Allocation policy

The total origin budget is `union_origin_cache_bytes`, unchanged; the inputs give 1 GiB. It is split once per block's batch:
- whole-segment aggregates get `floor(B/2)`, under the control's own admission rule;
- fine-core aggregates get the rest, charged as payload plus lookup.

Admission is deterministic. Each origin is owned exactly once: a segment with a whole zone never gets fine zones. Junction, edge, empty, tiny and declined pairs fall back to direct pairing.

Each block's batch builds and owns its own read-dependent whole and fine aggregates. Nothing is reused across blocks: the batch, with its aggregates, is destroyed before the next block's batch exists.

## Kept from S2-E

- 4,096 read pairs per block;
- whole-run exposure and validation;
- original fragment-order arithmetic;
- all candidates and all fragments;
- every executed-work and consumption cap, cumulative across blocks;
- the harness's 4 GiB footprint guard;
- the main arm keeps the 500 M executed physical-origin cap.

## Refusal reporting

`BlockedUnionResult` reports a refusal in three parts:
- `completed_totals`: the completed blocks;
- `refused_block`: the refused block's observed work, captured before its batch is destroyed;
- `attempt_totals`: the whole attempt, completed plus refused.

`refused_block.stage` is `discovery` or `columns`. The observed counters come from `ContentIndexBatch::observed_charges()`:
- piece hits and physical origins, read from the caps' own counters, so aborted construction or pairing is included;
- the consumption counters.

Every observed counter starts at the block's baseline, so a refusal before pairing observes zero of its own pairing work. No score is returned on refusal.

## Gates

The gate is `tests/genotype_blocked_fine.cpp` (ctest `panvar_genotype_blocked_fine`). It compares combined blocking against ordinary blocking, the bulk indexed control, bulk fine, and the reference scorer.

**Fixtures:**
- the I4-E fixtures, under two libraries;
- the S2-E shapes, with the prefilter cache.

**Checks:**
- **Identity.** Bit-identical masses (per block), score components, audit views and certified top sets.
- **Origin accounting.** Exact-once, per block and cumulatively: executed − constructed (whole + fine) + served equals the reference's logical origins, and logical state classes are equal. Executed and served work are reported apart.
- **Admission.** Partial admission, complete fallback, and the exact fine-partition boundary inside a block (at C and at C − 1 byte).
- **Robustness.** Forced fine-key hash collisions; 1 vs 4 threads giving identical science, work and admission.
- **Refusal.** A cumulative refusal in block 3 at four points: the physical-origin cap during construction and during pairing, the served-origin consumption cap, and the piece-hit cap during discovery. Each must report completed, refused and attempt totals, with no scores.
- **Invalid configurations.** Fine aggregates without fine views or without origin reuse are refused.

## Measurements (the runner `s2f_run.py`; manifest in the runner)

1. **The gate and the indexed regressions:**
   - blocked union;
   - progressive kernel;
   - union placement cache;
   - full e2e.
2. **5× cells** (C4, CYP2D6, ACOT; S2-C's union kernel inputs): ordinary blocked and blocked+fine at this revision. Each pair is compared bit for bit, and blocked+fine is also compared with S2-E's bulk scores.
3. **ACOT 30× blocked+fine**, with the actual 30× pool (605 haplotypes, 108,939 pairs, 40,021 fragments) and unchanged caps. Reported:
   - completion, or the first binding limit;
   - runtime and footprint;
   - executed, constructed and served origins;
   - state classes;
   - consumption-cap margins;
   - whole and fine admission and storage per block.
4. **Conditional control, extended-resource validation.** One ordinary-blocked ACOT 30× run with a physical-origin limit of **1,000,000,000**, declared here before execution. It runs only if the combined union in step 3 completes and budget remains.

## Budget and stop rule

- **Budget:** 1,200 s cumulative, development runs included; one process; at most 4 threads; at most 300 s per unit. A unit runs only if the time spent plus its declared timeout fits.
- **Stop rule:** no further variants and no cap sweeps.
  - If the combination completes ACOT 30×, recommend the smallest caller integration, with refinement unchanged at first.
  - Otherwise, report the measured reason.
