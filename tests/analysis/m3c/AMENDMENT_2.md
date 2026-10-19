# M3-C Amendment 2: the caller's refinement-context configuration

Declared after `lik.lzo_5x` (74.1 of 600 s spent, exit 1) and before any further execution.

## Failure

`lik.lzo_5x` stopped after 66.2 s (peak RSS 1.52 GB) with `fragment scoring reached its statistical-class cap`.

- **Cause: the tool's configuration.** It copied **`score-set`'s caps**, which are the caller's union-scoring caps
  (`max_state_classes` 600M), and drove every replacement through the single-request reference path.
- **How caps count.** In a context they are cumulative, and they count executed work.
- **What the caller actually does.** Its refinement builds its context with `make_repair_setup`:
  - the **repair** caps: `repair_max_state_classes` 4.2e9, physical origins 4.9e9, piece hits and latent pairs 1.4e10;
  - the origin-cache path with the run's `repair_origin_cache_bytes` (805,306,368; planned policy);
  - planned, batched `add_repaired_haplotypes`, so requests sharing a window enumerate it once.
- **What was missing.** The single-request path re-enumerated each window for every allele. No progress was recorded
  before the refusal, so how far the singles got is unknown.

## Change (configuration only; no library change, no scoring mathematics changed)

1. **`two-block-lik --config repair`** (now the default) uses the caller's refinement configuration:
   - the repair caps from the run's status;
   - `enable_origin_cache_path(true)` with `set_origin_cache_limits(repair_origin_cache_bytes, 0)`;
   - one `add_repaired_haplotypes` batch for all singles and one for all joints.
   - Promotion stays `add_repaired_haplotype(..., retain_census = true)`.
   - Cache hits change work, not scores. The planned cache is exact by construction and its scores must reproduce,
     which the gate and `compare.a2` check.
2. **`gate.a2`** re-runs the frozen fixtures with every second-step replacement through
   `add_repaired_haplotypes` on the origin-cache path, under the same requirements: 1e-9, byte identity, same slot for
   both orders, 0 full placements.

The rule, top-4, both 2,048 bounds, the run order (LZO, then 5×, then 2×), tolerances, selection and evaluation are
unchanged. `lik.lzo_5x` stays recorded with exit 1.

## Units (manifest A2)

| unit | declared timeout |
|---|---|
| `gate.a2` | 60 s |
| `lik.a2.lzo_5x` | 180 s |
| `lik.a2.mosaic_5x` | 180 s |
| `lik.a2.mosaic_2x` | 180 s |
| `rank.a2` | 30 s |
| `compare.a2` | 30 s |
| `eval.a2.<cell>`, ×3 | 60 s each |
