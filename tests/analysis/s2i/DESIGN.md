# S2-I: close the identified memory-accounting gaps (diagnostic only)

Branch `s2-i-memory-accounting` from `5d1d309`. Read-only accounting; no change to placement, threading, census capture, scoring, caps or decisions.

## 1. Allocator in-use bytes
- I4-G's `allocator_bytes_in_use()` (`malloc_zone_statistics` over every zone, `size_in_use`) is ported, read-only, and added as the live column `allocator_in_use`.
- Unsupported platforms write `unavailable`.
- It counts malloc'd bytes only: a different quantity from the physical footprint. Their difference is not reported as free pages.

## 2. Worker storage outside the S2-H ledger
`ContextMemoryAccount` gains new fields:
- `partition_scratch_bytes` / `_capacity`: each enumeration partition's ordinary `WindowScratch`, meaning rows, touched lists, proposed starts, groups, offsets, accumulators, accepted flags and origins. Its `WindowFineStore` is already in `window_scratch_capacity` and is not counted again.
- `partition_piece_index_bytes` / `_capacity`: each partition's forward and reverse `PieceIndex` vectors, plus the partition vector itself.
- `partitions`, and `capacity_total_extended` = `capacity_total` (its S2-H definition, unchanged) + the two new capacities.

The live file appends these columns and `footprint_minus_capacity_extended`. Whatever remains is outside the current ledger; it is not evidence of a leak.

**Gate.** `panvar_genotype_partition_accounts` (block_delta in a focused mode) compares the new capacities with `recount_partition_storage_for_testing()`, an independent recount written apart from `memory_account()`. It does so with reference and fine window discovery, at 4 threads.

## 3. The seed placement's own report
When the diagnostic is on, the seed `add_haplotypes` call receives an `IncrementalBatchReport*`. That parameter is output-only and does not change the call. The report goes to `<prefix>.refinement_seed_report.tsv`:
- `ScorerMemoryAccounts`: placement rows and their capacity peak, proposed starts and their capacity peak, piece index, forms, sequences, emission storage;
- the batch's work.

The `*_peak` values are per-worker transient peaks summed over workers. They are reported apart from the safe-point accounts and never added into a live total.

## 4. Runs
- Gates: the partition recount; the scoring-context gates; the caller diagnostics gate (the focused e2e, which now also checks the allocator column, the full row width and the seed report, and that the seed report survives the induced stop).
- One fresh-process ACOT 30× first-round replay, from S2-H's recorded anchor: original caps, 4 threads, 4 GiB guard. No union rerun and no full caller.
