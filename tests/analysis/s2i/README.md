# S2-I: closing the refinement memory-accounting gaps (diagnostic only)

**What it adds:**
- the allocator's in-use bytes (read-only) in the live refinement rows;
- the enumeration partitions' ordinary scratch and piece-index storage as separate size/capacity accounts, plus an extended total, checked against an independent recount (`panvar_genotype_partition_accounts`);
- the seed placement's own batch report (transient peaks) in `<prefix>.refinement_seed_report.tsv`.

**Outcome (ACOT 30× replay):** the ledger now explains allocator growth during refinement (seed +569 vs +573 MB; materialisation 94 %). The largest refinement-owned live item is the origin cache: 747 MB at request 323 of 646, nearly every window stored.

Proposed (not implemented): limit first-round origin-cache admission to windows the round requests more than once. Expected saving about 0.6–0.7 GB; no round-1 recomputation.

Still unknown: about 300 MB of the seed-placement footprint jump.

Results: `tests/results/s2i_memory_accounting_20261018/README.md`.
