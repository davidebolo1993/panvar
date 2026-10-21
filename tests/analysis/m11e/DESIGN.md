# M11-E: bounded attribution of the five M11-D catastrophic blocks (frozen before execution)

**Status of earlier batches:** M11-D is closed as a negative development comparison. Its frozen results are preserved and
read-only. The 18 ACOT validation families are not consumed.
- **Branch:** `m11-e-catastrophic-attribution`, from `m11-d-dev-pilot` `3039462`.
- **Results:** `tests/results/m11e_catastrophic_attribution_20261021/`.
- **Budget:** 900 s including gates, of which 120 s is disclosed development:
  - inspecting the saved traces and the score-set interface;
  - the control hash selection;
  - scratchpad dry runs of `fixes`, `prepare` and the runner regression (no scoring).
- **Resources:** one scoring process at a time; 4 threads; 6 GiB envelope.

**Not done:**
- no caller rerun, new reads, new inference algorithm or confidence fit;
- no cutoff, new metric or floor search;
- no validation, clustering, prior change, default promotion, merge, push or cleanup.

## 1. Reporting and eligibility corrections (`fixes`; nothing in M11-D is rewritten)

**1. Interval artefact totals.**
- **Defect:** `report.a4/blocks_summary.tsv` summed each interval's artefact bp over every flagged block in it.
- **Correction:** `artefact_correction.tsv` and `artefact_intervals.tsv` count the bp once per interval (the summed block
  distance minus the concatenated-interval distance, from the saved M11-D alignments). They also re-count the report.a4
  figures, to show where the original number came from.
- **Record:** the original report stays as it is. The correction is appended to M11-D's RESULT.

**2. Calibration-eligibility annotations** (`eligibility.tsv`, one row per M11-D block). Original labels are copied
unchanged with `relabelled = no`, and nothing is refitted.

| eligibility | rule | fit-informative |
|---|---|---|
| `UNDEFINED_TRUTH` / `NO_PREDICTION` | as in M11-D | no |
| `ALL_EMPTY_COUNTED` | both truth and both called alleles 0 bp. Counted in coverage | no |
| `REPRESENTATION_DEPENDENT_INTERVAL_EXACT` | the block lies in an evaluated degeneracy interval with artefact bp > 0, **and** the concatenated predicted interval equals the truth interval exactly (distance 0) | interval level only. The block label is **not** relabelled |
| `REPRESENTATION_DEPENDENT_UNRESOLVED` | as above, but the concatenated distance is > 0. The interval's availability floor was never computed, and computing it would be a new floor search | no |
| `BLOCK_OUTCOME` | everything else | yes, if assessable |

**3. The QVfrac per-homologue breakdown** of the catastrophic blocks (`qvfrac_catastrophic.tsv`). It shows how summed QVs let
one accurately predicted homologue carry a block past 0.8. QVfrac stays descriptive; excess and normalized excess stay primary.

**4. The runner** (`runner.py`):
- **Pass rule:** a unit passes only when it exits 0 **and** the watchdog verdict is `ok`.
- **Dependents:** a dependent runs only if every unit it needs has passed.
- **`next`** returns non-zero on a failure or NOT RUN, and 2 when everything is recorded.
- **Regression `gate.runner`** (`test_runner.py`, on throwaway units):
  - a failing analysis under verdict `ok` stops its dependent, and so does a timeout;
  - an independent unit still runs;
  - read-only: the sealed M11-D ledger shows the old precondition admitted `fit` after `outcomes` exited 1.

## 2. Cases (frozen)

| case | cell | block | role |
|---|---|---|---|
| K1 | c4.NA24385 | sv:7 | catastrophic (excess 39,155) |
| K2 | acot.NA19909 | sv:13 | catastrophic (27,270) |
| K3 | c4.NA19434 | sv:7 | catastrophic (26,210) |
| K4 | cyp.NA19650 | sv:9 | catastrophic (12,159) |
| K5 | cyp.NA19983 | sv:9 | catastrophic (12,075) |
| R1 | cyp.HG02953 | sv:9 | control |
| R2 | c4.HG00358 | sv:7 | control |
| R3 | acot.HG03209 | sv:13 | control |

**Control rule:** among each locus's floor-attaining `repeated_unit` blocks in M11-D, the lowest
sha256("M11-E:control:" + cell + ":" + block).
- CYP2D6 pool: HG02769, HG02953.
- C4 pool: NA19705, NA19317, HG03065, HG00358, HG03732.
- ACOT pool: HG03209, HG00639, HG02074.

## 3. Method (saved inputs and existing probes only)

**Incumbent:** the reported winner, the representative rank-1 row of the run's `equivalence_set.tsv`. Its allele vectors
and full-precision total are taken from there.
- **Check:** its block alleles must equal the evaluator's called pair.

**Replacements:** only the affected block changes, to the evaluator's best-available pair (`truth1_closest_allele`,
`truth2_closest_allele`). Both homologue assignments are built:
- `assign_A`: hap1 gets the truth1-best allele, hap2 gets the truth2-best allele;
- `assign_B`: the reverse.

**Verification (`prepare`, before scoring):** each haplotype is spelled from the block-allele FASTAs of the run's own model,
and then:
- the incumbent's SHA-256 must equal the run's recorded haplotype SHA-256;
- every replacement keeps all other blocks byte-identical;
- the replaced block is exactly the declared allele.

**Scoring:** `genotype_probe score-set` with the run's model, curation, reads, library and status-recorded scoring
settings. The transition specification is checked against the run's.
- **Reproduction:** the incumbent is declared with its expected total, and score-set refuses unless it reproduces within
  1e-6.
- **Byte identity:** the probe's scored bytes (FNV-64 and length) must equal the verified spellings.
- **`--rank-blocks <block>`:** the emission rank of the best pair and whether it is in the repair's local-pair list
  (`repair_local_pairs`).
- **`--decompose`:** each replacement against the incumbent, fragment by fragment, with placement classes. The origin
  re-summation must equal the scorer's fragment term within 1e-6.

**Trace (`trace`), from saved files only:**
- **Run record:** whether each replacement diplotype is in the run's scored-candidate record (`mosaic_candidates`, union +
  repair). Where it is, the run's recorded total must equal score-set's within 1e-5 (scores are printed to about 12
  significant digits).
- **Proposed:** whether the best pair at this block appears in any recorded candidate, in any context.
- **Repair rounds.**
- **Dosage audit row:** reported, never treated as proof.

## 4. Classification (frozen order; per catastrophic case; controls labelled `CONTROL:` and expected to reproduce)

1. **`UNRESOLVED_ATTRIBUTION`** (harness): the incumbent is not reproduced, or the scored bytes are not verified.
2. **`NUMERICAL_OR_IMPLEMENTATION_INCONSISTENCY`:** any of the following.
   - The run's recorded total for a replacement differs from score-set's.
   - The origin re-summation disagrees with the scorer.
   - The run itself scored a replacement with a higher total yet reported the incumbent.
3. **`SEARCH_MISSED_BETTER_ALTERNATIVE`:** the better replacement's Δtotal exceeds +1e-6 and the run never scored it.
4. **`OBJECTIVE_PREFERS_FARTHER_PREDICTION`:** both replacements score lower by more than 1e-6. The term differences say
   which component carries the preference.
5. Otherwise **`UNRESOLVED_ATTRIBUTION`** (a tie).

**Interpretation limits:**
- A losing local replacement is evidence about the objective **in the winner's context only**. It does not establish a
  global optimum, and it is not by itself a likelihood bug.
- A case is never forced into a known explanation. Anything not covered by the rules above stays unresolved.

## 5. Units (`m11e_run.py`; dependencies gate execution)

| unit | needs | declared timeout (s) |
|---|---|---|
| `build` | | 300 |
| `gate.runner` | build | 60 |
| `fixes` | gate.runner | 120 |
| `prepare` | gate.runner | 120 |
| `score.K1..K5`, `score.R1..R3` | build, prepare | C4 / CYP2D6 120, ACOT 240 |
| `trace` | prepare | 120 |

**Expected time:** about 150–300 s (score-set loads one model, scores three diplotypes and decomposes two).

**Admission:** a unit runs only if spent + its timeout ≤ 900 s; otherwise it is recorded NOT RUN.

**Failures:**
- a score unit's failure (for example, a reproduction refusal) is that case's recorded result (UNRESOLVED), and the other
  cases continue;
- an analysis-code failure needs a declared amendment before anything is re-run.
