# M8-C: paired ACOT genotyping challenge, constant r = 0.01 vs r = 0.1 (frozen before reads or outcomes)

**Ancestry.** Branch `m8-c-copying-challenge` from `06971e3`. M8-B is accepted as exploratory masked-score evidence.

**Scope.** A mechanistic development test on two donors, not release validation.
- **Excluded:** default promotion, fitting, grid extension, C4/CYP2D6 campaigns, distance mode, clustering, residual
  editing, output redesign, merge, push and cleanup.
- **Budget:** 1,200 s counting preparation and failures. Admission: spent + timeout + 60 s (contrast + summary) ≤ 1,200.

## Samples and panels (`prepare`)

- **Donors:** from M8-B's 30 ACOT outer-evaluation families, the first two by sha256("M8-C:" + family) whose donor is a
  complete diploid (templates #1 and #2). Nothing about scores, complexity or expected outcome enters. The selected
  identities are kept whatever the results.
- **Calling panels:** curation manifests on the **fixed full graph** `acot.gfa`, with reference
  `GRCh38#0#chr14:73416151-73668669`.
  - **LOO:** M8-B's frozen training bank, i.e. every outer-evaluation template excluded (the queried family included).
  - **LZO:** the same bank plus the queried donor's two haplotypes.
  - No other outer-evaluation donor enters any panel.
  - The frame (graph, bubbles, blocks) is derived from the full 467-template panel; this is disclosed.
- **Simulator:** the current build. Its lineage is re-checked by regenerating the archived C4 NA19650 30× source byte for
  byte. Reads use `simulate --depth 30 --seed 20260916 --seed-salt <donor>` (the ACOT pilot rule), then
  `subsample --coverage 5 --name-salt <donor>`. **Each donor's four runs use byte-identical reads and library.**
- **Profiles:** the M4-B/P1 composition of `perf/acot_5.tsv` with the full performance overlay. The validated combined
  route; 4 threads; the base profile's 6 GiB footprint and 900 s guards are unchanged. Each arm adds only
  `copying_transition constant_v1` and `switch_probability r`.
  - **Gate:** the two arms differ only in the profile name and `switch_probability`.
  - The caller command is otherwise identical in both arms: `--dosage-mode repair`, the curation, threads.
  - Joint search stays off (default) in both.

## Runs (frozen order)

| step | units |
|---|---|
| 1 | `prepare` |
| 2 | 8 calls, each 120 s: donor 1 LZO r0.01, LZO r0.1, LOO r0.01, LOO r0.1, then donor 2 the same |
| 3 | `checks`: every recorded `switch_probability` equals the arm's r (no legacy/experimental mixture); any transition digests agree; inputs byte-identical per donor |
| 4 | 8 `evaluate` runs, each 150 s, after all predictions are frozen. Order: LOO (donor 1 r0.01, r0.1; donor 2 r0.01, r0.1), then LZO in the same pattern. A unit the budget refuses is recorded NOT RUN, and no setting changes |
| 5 | `contrast`: for a donor × scheme whose reported unordered sequence pair differs between arms, whether each winner was scored by the other run (search reach vs ranking), and fragment / exposure / prior of both winners under both runs, from the existing `score-set` |
| 6 | `summary`: every cell in the denominator (ties, no-calls and failures included) |

**`summary` reports:**
- call class and phase support;
- evaluation floors, distances and excess;
- phased versus unordered results;
- exact blocks;
- sequence-based panel membership of each called haplotype: truth template, other panel member, or off-panel
  mosaic/reconstruction (SHA-256 against every template spelling and its reverse complement);
- search work, runtime and memory;
- paired verdicts.

## Decision rules (declared)

- **Error:** excess over the availability floor > 0. In LOO, a call at the floor is not necessarily the exact truth.
- **No errors in either arm:** the comparison is **non-discriminating**.
- **r = 0.1 better in LOO but worse in LZO:** reported as a **trade-off**, with no tuning.
- **LZO inspection:** whether r = 0.1 introduces an off-panel or mosaic call with excess > 0. A different donor label alone
  is not an error.

## Amendment 1 (declared after `prepare` exited 1; 3.4 s spent)

**What failed.** `subsample` does not create its output directory: it reported success and wrote nothing. The simulator
lineage re-check passed before that step.

**Fix.** `prepare` creates the directory first. The unit re-runs as `prepare.a1`, and every later unit reads its outputs.
Donor selection, seeds, panels and profiles are unchanged.
