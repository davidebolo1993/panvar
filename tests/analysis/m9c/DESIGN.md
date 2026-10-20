# M9-C: bounded ACOT reliability pilot under the current output contract

A descriptive check of whether the v2 projection-aware block error masses track exact block-allele agreement on ACOT
families never previously genotyped or confidence-inspected.

**Not done:** confidence fitting, GQ remapping, PASS rule, prior tuning, clustering, default promotion, caller or schema
change, merge, push or cleanup.

**Results:** `tests/results/m9c_acot_reliability_20261020/`.

**First step, before any execution:** M9-B's qualification is appended to `tests/analysis/m9b/RESULT.md`. It covers:
- legacy no-call versus v2 prediction;
- the existing `error_mass_lower`/`error_mass_upper`;
- the v2 FASTA policy;
- the narrower C4 allocation evidence;
- conservative enclosing intervals.

## Selection (`prepare.py`)

**Pool.** M8-B's 30 complete-diploid ACOT outer-evaluation families (`split.ACOT.tsv`, SHA-256 recorded), minus the two
genotyped in M8-C/M8-D: SL17 = HG03209 and PEL25 = HG01975. That leaves **28 families**.

**Rule.** The first **six** families by sha256("M9-C:" + family). Each family's single complete diploid donor is its
member; every one of the 30 is single-donor with #1 and #2, which is asserted.

**Cells:**
- **six LOO cells at 5×;**
- **two LZO controls:** the first two selected donors, with byte-identical reads to their LOO cells.

**Disclosed.** These assemblies took part in M8-B's score-level outer evaluation of r (masked copying-prior scores). This is
new genotype-confidence evidence, **not a pristine holdout for every modelling decision.**

**Panels** follow M8-C's recipe on the **fixed full-panel graph `acot.gfa`**.
- **Graph frame:** blocks, labels and bubbles are derived from the full 467-template panel, eval families included.
  Disclosed.
- **LOO:** every outer-evaluation template is excluded. The list is asserted identical to M8-C's.
- **LZO:** adds only the queried donor's two haplotypes.

**Reads** follow M8-C's recipe:
- `simulate --depth 30 --read-length 100 --seed 20260916 --seed-salt <donor>`, then `subsample` to 5× with
  `--name-salt <donor>`;
- the simulator lineage check is M8-C's (`prepare.a1` PASS). It carries over because `genotype_probe` is asserted
  byte-identical to M8-C's.

## Route and settings (frozen before calls)

**Executables:** panvar `98f339e8…` and genotype_probe `af2ed281…`, asserted equal to M8-C's. The relinked test binaries
are not used for calls.

**Profile:** M8-C's r = 0.01 rule is re-derived and asserted byte-identical to M8-C's `profile.r0.01.tsv`. Then:
- the profile name is changed;
- `output_schema sequence_first_v2` is added;
- `repair_joint_search off` is added.

**What that gives:**
- the **validated combined computational route**: P1 composition of `acot_5.tsv` with the full
  `genotype_dev_performance_v1` overlay;
- `copying_transition constant_v1`, `switch_probability 0.01`;
- four threads;
- `footprint_limit_bytes 6442450944` (6 GiB) and `time_limit_seconds 900`;
- every other scientific setting and caller guard unchanged.

**On the command line:** `--dosage-mode audit` (declared).

**This is not "untouched production defaults".**

**Check after calls, before evaluation:** `freeze` reads each run's settings and status and compares `output_schema`,
`switch_probability`, `copying_transition`, `repair_joint_search`, `threads`, `time_limit_seconds` and
`footprint_limit_bytes` (plus `dosage_mode` from the sequence locus) against these values. Mismatches and unverified keys
are recorded.

## Contract fixtures (gates; calls are refused if either fails)

The existing focused tests are relinked against the current library (`build`); the binaries predate the M4-B library
change. Then they run.

**`gate.sequence_report`** runs `genotype_sequence_report_test`:
- `test_coupled`: a legacy label-only distinct tie still yields a v2 `PREDICTION`/`RESOLVED`. Sequence-identical label
  ambiguity does not erase the prediction.
- `test_bounds_count_groups_once`: mixed non-winning projections give `error_mass_lower 0` and `upper` = the group mass
  counted once.
- `test_panel_membership`: panel pair gives no FASTA with extraction references; off-panel in either haplotype gives
  FASTA.

**`gate.e2e_sequence_first`** runs `genotype_command_e2e_test` with `PANVAR_E2E_ONLY=sequence_first`: the caller-level v1/v2
schema, lifecycle and FASTA contract.

**If the contract or the fixtures disagree, stop before real-data calls.**

## Execution

- **One caller process at a time.**
- **Process-group watchdog: 930 s per call.** The caller's own guard is 900 s; the earlier 120-s bound is not repeated.
- **Budget: 2,400 s in total,** including build, preparation, gates, failures and evaluation.
- **Admission rule:** spent + the unit's timeout + the declared timeouts of every pending freeze/eval/analyse unit (550 s)
  ≤ 2,400.
- **Refused units** are recorded NOT RUN. Donors are never substituted.
- **Order:** build, prepare, the two gates, `call.LOO.d1` … `call.LOO.d6`, `call.LZO.d1`, `call.LZO.d2`, `freeze`, the
  eight evaluations, `analyse`.

## Freeze before evaluation

`freeze` runs only after every call has a recorded outcome, and reads no truth. For each cell it records:
- **outcome:** resource failure, no sequence report, `NO_PREDICTION`, `PROVISIONAL`, or `PREDICTION` with its
  `sequence_status`;
- **the v2 prediction:** representative SHA-256s, haplotype origins, FASTA decision;
- block-status counts and representation intervals;
- legacy status and no-call reason;
- SHA-256 of the authoritative `.sequence_*` outputs;
- the settings check.

**Every selected cell stays in the denominator,** including resource failures, provisional predictions and unresolved
genotypes. No replacement is chosen while calling.

## Evaluation and descriptive analysis

**Evaluation** uses the archived evaluator, `genotype_probe evaluate`, with the same model, curation and truth graph. Only
its per-block truth columns are used for v2 agreement: truth-closest alleles, distances and availability.
- **Labels** are distinct sequences within a block. So exact block-allele agreement means both truth alleles are present
  in the calling panel (distance 0, i.e. floor zero) and the IDENTIFIED unordered label pair equals them.

**Unit of analysis.** Reportable (SV) blocks with `block_status IDENTIFIED`. All blocks of every cell enter the overall
accounting, by kind and status.

**Confidence value:**
- `error_mass_upper` is used directly; `error_mass_lower` is reported beside it.
- **Never** 1 − rounded posterior.

**Predeclared bins** on `error_mass_upper`:
- printed 0;
- (0, 1e-6);
- [1e-6, 1e-3);
- [1e-3, 1e-2);
- [1e-2, 1e-1);
- [0.1, 1].

**Strata:**
- **floor-zero:** both truth alleles are in the panel;
- **truth allele unavailable in the calling panel:** exact agreement is impossible;
- **truth undefined.**

**Reported per scheme (LOO / LZO):** counts, exact, not exact, and families.

**Family-level table.** One row per cell: outcome, reportable, identified, exact, errors, high-confidence errors, phase
support, resources. Blocks are not independent samples.

**Legacy GQ comparison.** Compare `error_from_GQ` = 10^(−GQ/10) with `error_mass_upper`, only where the assumptions
coincide: `lower == upper` (no mixed groups at the block) and a single decision-set group. Otherwise "assumptions differ".

**High-confidence errors** (`error_mass_upper` < 0.01, floor zero, not exact) are inspected **after** evaluation and
classified with this precedence:
1. **representation ambiguity:** the block is inside a representation interval, or lower ≠ upper;
2. **correct block pair absent from the evaluated candidates;**
3. **correct pair present but losing as whole hypotheses:** reported as the best carrier's whole-hypothesis
   total / likelihood / prior gap.
4. **Unsupported attribution:** none of the above can be assigned.

The third category is never called a block-local likelihood defect, because there is no isolated comparison.

**Kept separate:**
- **exact block-sequence correctness** versus **structural correctness**, which is not assessed (no structural target);
- **phase**, reported from the legacy phase support per cell.

**No calibration curve, fitting, GQ remapping, PASS rule, prior tuning, clustering or promotion.**
