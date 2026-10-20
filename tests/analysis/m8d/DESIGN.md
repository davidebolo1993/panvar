# M8-D: continuation of the interrupted HG01975 LOO comparison (r = 0.01 vs r = 0.1)

A separate continuation record for M8-C. **M8-C, its timed-out runs, its outputs and its tables are not changed.**
Results: `tests/results/m8d_hg01975_loo_20261020/`.

## Why

M8-C's runner killed both HG01975 LOO calls at its declared 120-s unit bound. The caller's own profile guard is 900 s,
so the caller never had its declared runtime. M8-D completes only this pair.

## What is rerun (exactly)

| unit | command |
|---|---|
| `call.r0.01` | M8-C's `call.HG01975.LOO.r0.01` command, with `--out-prefix` in a new directory |
| `call.r0.1` | M8-C's `call.HG01975.LOO.r0.1` command, with `--out-prefix` in a new directory |

**Inputs reused unchanged** (all checked by `hashcheck` before any call; a mismatch stops the batch):
- panvar `98f339e8…` and genotype_probe `af2ed281…`, full SHA-256 equal to M8-C `prepare.tsv`;
- `acot.gfa` `05a1a82a…`;
- reads and library from M8-C `prepare.a1/reads/HG01975/cov5`, prefixes equal to `prepare.tsv`;
- the curation `curation.LOO.HG01975.txt` and profiles `profile.r0.01.tsv` and `profile.r0.1.tsv`. M8-C did not hash
  these, so `hashcheck` re-derives them with M8-C's frozen prepare rule (the split file, by its recorded SHA-256, and the
  composed `acot_5.tsv` profile) and requires byte equality;
- M8-C's recorded command lines, which must reference exactly these files.

**Unchanged:** reads (not regenerated), every scientific setting, the search bounds, the 900-s caller guard, the 6-GiB
footprint guard and the filters. No new donors, coverages, parameter values or fitting.

## Execution

- **Sequential, one process at a time**, nice, 4 threads, `p1_run.run_timed`: own session/process group, whole-group
  SIGTERM then SIGKILL, descendants reaped.
- **External watchdog: 930 s per call** (900 s plus 30 s for shutdown and output finalisation).
- **Budget: 2,100 s** of new execution, including evaluation and analysis, counted as wall time spent plus the declared
  timeout of the next unit.
- **Both arms are reserved before starting.** A unit is admitted only if: spent + its timeout + 930 s for each call still
  pending + 200 s (two evaluations at 60 s, contrast 60 s, summary 20 s) ≤ 2,100.
  - `hashcheck`: 0 + 30 + 1,860 + 200 = 2,090.
- **Binding guard or watchdog.** If the caller guard or the watchdog binds on an arm, that outcome is retained as final for
  the arm and is never rerun; no limit is increased.
  - Interpretation, declared: the other arm of the frozen pair still runs within its already reserved budget, so that
    both outcomes are recorded.
  - With an incomplete arm, the pair is reported as not assessable.

## Manifest P (frozen)

| unit | kind | timeout (s) |
|---|---|---|
| `hashcheck` | frozen-input check | 30 |
| `call.r0.01` | caller | 930 |
| `call.r0.1` | caller | 930 |
| `freeze` | prediction freeze | 10 |
| `eval.r0.01` | truth evaluation (only if COMPLETE) | 60 |
| `eval.r0.1` | truth evaluation (only if COMPLETE) | 60 |
| `contrast` | search and score comparison | 60 |
| `summary` | corrected paired summary | 20 |

## Prediction freeze and evaluation

**`freeze` runs before any truth is read.** For each run it records:
- the outcome: **COMPLETE** only if `call.locus_genotype.tsv` exists with `run_status = COMPLETE`, otherwise INCOMPLETE;
- `run_status`, and `switch_probability` from the settings (the arm check);
- the SHA-256 of `locus_genotype`, `equivalence_set` and `mosaic_haplotypes.fa`.

**Evaluation** uses the original evaluator (`genotype_probe evaluate --scheme loo`, with the same model, truth graph and
curation as M8-C) and runs only on COMPLETE outputs.
- An INCOMPLETE run is recorded as not evaluated. Partial refinement state (`*.partial.tsv`) is never evaluated as a
  genotype.
- Reported: unordered sequence prediction (haplotype SHA-256 and sequence-based panel membership via M8-C's
  `analyse.py summary` method), projection ambiguity, floor/distance/excess, exact blocks and SV blocks where assessable,
  phase separately, and runtime and resources (wall, peak RSS, exit, timeout).

## Search comparison (`contrast`)

- **The M8-C partial traces already differ:**
  - the first three repair moves coincide (ctx:4, ctx:6, ctx:8), with only the prior differing;
  - later, r = 0.01 moves ctx:0 where r = 0.1 moves sv:13.
- **"Search work identical" in M8-C applies only to its completed pairs.** `summary` writes both the partial and the
  completed move sequences.
- **If the final unordered sequence pairs differ,** both winners are rescored under each run with the existing
  `score-set` probe: same model, curation and reads, `--switch-probability` set to that run's r, and the caller's total
  as the expected value for its own winner. This gives fragment, exposure and prior for the same contrast.
- **Reachability is reported separately from ranking.** A winner is *reached* by a run if its allele pair appears among
  that run's scored mosaic candidates or among its repair-round states.
  - Changed ranking: the other winner was reached but scored lower.
  - Different reachability: it was never reached.
- **Defects corrected before freezing.** The contrast code M8-C inherited never ran, because no call changed, and had
  two latent defects, now fixed:
  - `score-set` defaults to r = 0.01, so r is now passed explicitly;
  - reachability now includes the repair-round states.

## Corrected missing-result semantics (`summary`)

`same_call(a, b)`:
- **both complete:** compare unordered haplotype-sequence pairs (SHA-256);
- **either incomplete or missing:** **UNKNOWN** (not assessable).

**Fixtures,** which must pass or `summary` exits 1:
- both missing → UNKNOWN;
- one missing → UNKNOWN;
- one incomplete → UNKNOWN;
- both complete with swapped order → True;
- both complete and different → False.

The corrected table is a new file (`paired_corrected.tsv`) in this package. **M8-C's `summary.a2/paired.tsv` (where the
timed-out pair reads `same_call = True`) is not overwritten.**

**Combined table:**
- M8-C's six completed cells;
- M8-C's two HG01975 LOO TIMEOUT rows, kept and labelled as historical;
- the two M8-D outcomes.

## Closure rule (frozen)

**No default promotion follows.**
- **Predictions unchanged:** r = 0.1 is parked as a masked-score improvement without demonstrated genotyping benefit.
- **Predictions differ:** the benefit or deterioration (excess over floor, distance, exact blocks) is reported without
  retuning.
- **Either arm incomplete:** the pair is not assessable, the outcome is retained, and the experiment closes on the six
  M8-C cells, which give the non-discriminating verdict.

## Pre-freeze check (disclosed)

`analyse.py hashcheck` was dry-run once into the scratchpad (about 3 s, no caller execution) to validate the check. It
passed, and its time is counted in RESULT.
