# M15-A result: short-read release candidate frozen; output/filter contract and final evaluation protocol proposed; stopped for approval

**Run.** Branch `m15-a-release-candidate`; verification runner `6b43eaa`; amendment A1 `5003586`; documents committed after.
- **Results:** `tests/results/m15a_release_candidate_20261021/`.
- **Budget:** 693.0 of 900 s, including 60 s of disclosed development.

**Not done:**
- no genotyping of reserved families and no reserved truth opened;
- no merge, push, default promotion or cleanup.

**Deliverables** (in `tests/analysis/m15a/`):
- `RELEASE_CANDIDATE.md`: source, binary, dependencies, index format, effective configuration, feature classification,
  validated-together scope and consolidation requirements;
- `OUTPUT_AND_FILTER_CONTRACT.md`: six questions mapped to fields; validity checks, annotations, **no quality threshold**
  (an experimental release with uncalibrated support); decisions for Davide;
- `FINAL_EVALUATION_PROTOCOL.md`: evidence classes, the reserved-ACOT manifest and input recipe, analysis rules, and what can
  and cannot be established;
- `RELEASE_BLOCKERS.md`.

**Verification:**

| unit | result |
|---|---|
| `build` | every target, 23.3 s |
| `ctest` (frozen focused selection, −j1) | 58 passed, 2 disabled; `genotype_command_e2e` hit the 120 s per-test cap (stale cost data); unit watchdog at 420 s while `smoke_real_c4` ran |
| `ctest.a1` (A1: exactly those two, 280 s per test) | **both passed** (`command_e2e` 175 s, `smoke_real_c4` 10.7 s), so **60 of 60** enabled selected tests pass |
| `reproduce` | **15 of 15** M14-A full-panel runs (current sources, re-linked `panvar` `47042a5c…`) identical to M11-D's archived calls, with every RC output present |

**Source check:** `src/` and `include/` are unchanged since `ed10d53` (M10-A + M11-C), so M12, M13 and M14 cannot affect
predictions.

**Closures recorded:** M14-A (NOT USEFUL accepted, with qualifications and block/cell/donor counts) and its manifest
binary-hash correction.

**Stopped for approval before the final validation.**
