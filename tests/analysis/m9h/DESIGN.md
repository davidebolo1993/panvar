# M9-H: bounded paired validation of JP1.1 on additional ACOT families

**This validates the frozen mechanism.** It is not optimisation or tuning. `repair_junction_proposals experimental_v1`
stays default off, and the binaries are M9-G's, asserted by hash.

**Results:** `tests/results/m9h_jp11_validation_20261020/`.

## 1. M9-G closure (`correction`, before anything else)

**Tie-class wording.** JP1.1 compares each score with its **class head**: a new class starts when head − score > 1e-9.
Scores are not chained through adjacent scores. M9-G's DESIGN and RESULT said "chained". The code is unchanged; the
correction is appended to M9-G's RESULT.

**`same_sequence_pair`.** It now compares **unordered** pairs of haplotype-sequence digests: homozygous as a pair of
equal digests. It returns **UNKNOWN** if either prediction is absent or incomplete (not `PREDICTION`, or a missing digest).
- **Fixtures:** swapped haplotypes give True; a different pair gives False; homozygous versus heterozygous gives False;
  missing, provisional and NO_PREDICTION give UNKNOWN.
- **Re-applied** to M9-G's saved outputs, with the conclusion check recorded. M9-G's tables are untouched.

## 2. Validation set (`prepare`, frozen before reads or outcomes)

**Pool.** M8-B's 30 ACOT outer-evaluation families (`split.ACOT.tsv`, SHA-256 recorded) minus every family already used
for genotyping or confidence inspection:
- SL17 (HG03209) and PEL25 (HG01975) in M8-C/D;
- 2427, NA18747, HG02178, PR03, NA21106 and BD27 in M9-C to M9-G.

The remaining count is recorded and asserted to be **22**.

**Selection:**
- **Families:** the first **4** by ascending sha256("M9-H:" + family).
- **Donor:** per family, among complete diploid donors (exactly #1 and #2), the smallest sha256("M9-H:donor:" + donor).
  The eligible set is recorded.

**Cells:**

| cell | scheme | family |
|---|---|---|
| L1–L4 | LOO | families 1–4 |
| Z1, Z2 | LZO | families 1 and 2, fixed now |

**Disclosed.** These families are new to read-based genotyping validation, but they are not untouched holdouts: M8-B used
their assemblies for masked prior scoring.

**Panels** follow M8-C's recipe on the full-panel graph frame `acot.gfa`; the frame is derived from all templates
(disclosed).
- **LOO:** excludes every outer-evaluation template (asserted equal to M9-C's list).
- **LZO:** adds back only the selected donor's own two haplotypes.

**Reads** follow M8-C's recipe: 30× (seed 20260916, salt = donor), subsampled to 5× (name-salt = donor).
- **Lineage re-check:** genotype_probe changed in M9-G, so regenerating M9-C's NA19909 reads must reproduce M9-C's
  recorded SHA-256s byte for byte.
- **Within each pair:** reads and library are identical.

**Profiles:**
- **off:** M9-C's `profile.tsv`, the validated M9-G baseline;
- **on:** M9-G's `profile.on.tsv`. It is asserted to differ only in the profile name and `repair_junction_proposals
  experimental_v1`.

**Common to both arms:**
- the validated combined route;
- r = 0.01, `sequence_first_v2`, joint search off;
- `--dosage-mode audit`;
- four threads;
- the 6 GiB and 900 s caller guards, existing work caps;
- a 930-s process-group watchdog.

**`models`.** The relinked dump tool writes the block alleles of the LOO model and of both new LZO models, used for
interval spelling. Their fingerprints are checked against each run.

**`gate.unit`.** The JP1.1 unit gates are rerun as a smoke check; the binaries are unchanged since M9-G.

## 3. Run order (12 caller runs, one at a time)

**Cells:** L1, L2, L3, L4, Z1, Z2. **Arms alternate deterministically** by cell position: off then on for L1, L3 and Z1;
on then off for L2, L4 and Z2.

**Then:**
1. `freeze`: every prediction of both arms, before any truth-based step;
2. 12 evaluations (`genotype_probe evaluate`);
3. `truth.<cell>`: `block-sequences sv:5,ctx:6` with the truth projection;
4. `intervals`;
5. `analyse`.

## 4. Evaluation (`m9h.py`)

**Primary outcomes, per cell and arm:**
1. **Paired sequence distance over the coupled sv:5..ctx:6 interval.**
   - **Predicted bytes:** the v2 representative labelling's sv:5 + ctx:6 alleles, spelled from the dumped model.
   - **Truth bytes:** the truth projection's sv:5 + ctx:6.
   - **Distance:** the exact edit distance (`genotype_probe edit-distances`), minimised over the homologue exchange.
2. **Exact reportable-block agreement.** For v2 IDENTIFIED SV blocks with floor-zero truth, the unordered labels must
   equal the truth-closest alleles. Not-exact, unattainable (truth not in the calling panel) and not-identified blocks are
   counted separately.
3. **LZO control stability.**

**Also reported:**
- whole-locus distance, availability floor and excess (non-normalised);
- every block whose v2 genotype changed;
- eligible blocks, proposal invocations, cache hits and cutoff ties;
- trunk moves and tie branches;
- fragment, exposure, prior and total;
- wall time, peak memory, refusals and timeouts;
- error-mass bounds at changed blocks, reported as pool-conditional and **not** calibrated probabilities;
- phase support, which is not phase accuracy.

**Independent units:** families, not blocks.

## 5. Predeclared interpretation

**Per family (off versus on).** Three metrics: interval distance (lower is better), exact reportable blocks (higher is
better) and locus excess (lower is better).

| outcome | rule |
|---|---|
| UNCHANGED | the same unordered sequence pair |
| IMPROVED | changed, at least one metric better, none worse |
| DETERIORATED | changed, at least one metric worse, none better |
| MIXED_WITHIN_FAMILY | changed, both better and worse metrics |
| CHANGED_NEUTRAL | changed, no metric moves |
| INCOMPLETE_PAIR | either arm lacks a complete prediction; never counted as unchanged |

**Batch verdict:**

| verdict | rule |
|---|---|
| **CONTROL_HARM** | a correct LZO control deteriorates |
| **NON_DISCRIMINATING** | no changed sequence prediction anywhere |
| **SUPPORTIVE** | at least one LOO family IMPROVED, none deteriorated, both LZO controls UNCHANGED |
| **MIXED** | improvements and deteriorations (every deterioration reported) |
| otherwise | stated as found |

**Incomplete pairs limit the conclusion.**

**None of these verdicts authorises default promotion.** There is no retuning of k, thresholds, the prior, activation or
the tie rules.

## 6. Budget

**Total: 3,600 s,** including 2 s of disclosed pre-freeze inspection.

**Admission:** spent + the unit's timeout + the declared timeouts of every pending
freeze/eval/truth/intervals/analyse unit (for preparation, gates and calls) ≤ 3,600. Refused units are recorded NOT RUN,
never replaced.
