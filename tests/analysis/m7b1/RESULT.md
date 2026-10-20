# M7-B.1 result: with an N-free, validated calibration, Locityper's objective improves neither failure and breaks the C4 control's copy number (fixed pool, synthetic calibration)

**Run.** Branch `m7-b1-resume`: freeze `6abc72d`, Amendment 1 (attempt verification). Results in
`tests/results/m7b1_fixed_pool_20261020/`.

**Budget.**
- **93.9 s of new execution, of 800.** M7-B's 66.0 s is recorded separately in its own package, which is unchanged.
- Two recorded failures: both control units, on the corrupt `sol_ext` attempt check (Amendment 1); both were re-run as
  `.a1`.
- 4 threads, one process at a time.

**Scope.** No downloads, production changes, scoring reimplementation, merge, push or cleanup.

## 1. What changed and what was reused

**Background (remedy (a), recorded before execution).** `m7bg` is the direct concatenation of GRCh38 ACOT (252,518 bp,
SHA-256 `c3f8bd5f…`) and GSTM1 (240,715 bp, `4d41ff16…`): 493,233 bp, **only A/C/G/T**, background SHA-256 `06a677be…`.
- **Artificial junction** at positions 252,518 | 252,519.
- The same homozygous diploid representation.
- Background reads were regenerated with the M7-B-verified lineages and unchanged depths, libraries and salts: 0
  fragments skipped; λ and template z-scores pass.
- References, indexes, jellyfish counts and `preproc` were rebuilt. No old calibration or reads were reused.

**Reused after binding checks (`reuse`: all pass).** From M7-B: provenance (both lineages PASS), canonical pools and
mapping, Panvar D/X scores (all expectations reproduced) and truth distances. 92 files are bound by SHA-256, and every
executable's SHA-256 is unchanged.

## 2. Calibration (explicit validation: all four cells PASS)

| cell | background reads | insert mean / var (NB) | mismatch rate | windows (size 67, neighbourhood 300) | GC bins with ≥ 100 windows |
|---|---|---|---|---|---|
| C4 HG00171 LZO | 12,302 | 350.2 / 1,601 | e^−6.233 ≈ 0.0020 | 4,429 | 33–52 (20 bins) |
| CYP HG02769 LZO | 12,479 | 349.9 / 1,629 | e^−6.226 | 4,418 | 33–51 (19) |
| C4 NA19650 LOO | 24,735 | 350.2 / 1,594 | e^−6.230 | 4,429 | 33–52 (20) |
| CYP NA19650 LOO | 12,332 | 350.1 / 1,580 | e^−6.214 | 4,418 | 33–51 (19) |

- Every fitted value is finite and valid: insert NB, error profile, and all 101 depth-NB GC bins.
- Locityper recovers the simulator's library (insert 350 ± 40, substitution 0.002, no indels).

**Weakly supported GC (reported, as declared).** Many target windows (weight ≥ 0.001) fall in GC bins with fewer than
100 background windows, where Locityper's depth model uses its extrapolated tail values:

| locus | target-window fraction relying on extrapolated GC | target GC range | supported range |
|---|---|---|---|
| C4 | **56.3–56.5 %** | 16–82 % | 33–52 % |
| CYP2D6 | **35.0–35.3 %** | 20–75 % | 33–51 % |

**This is a material limitation for the depth term, especially at C4.** The background is a small synthetic calibration
set made from real sequence. It provides neither genome-wide k-mer uniqueness nor real-WGS calibration.

## 3. Run validity (every unit used passes)

**Checks passing in every run:**
- the database equals the pool exactly (nothing discarded);
- the `res.json` options equal the priors set;
- every canonical genotype has a stage-2 row in `sol.csv.br`;
- Stage II is logged as `a=20`, and every `lik_sd` is finite.

**Read accounting.**
- Recruited: 99.84–99.87 % of reads.
- Used read pairs: C4 LZO 5,699 (1 poorly mapped); CYP LZO 5,396; C4 LOO 11,692 (4); CYP LOO 5,423 (5).

**Disclosed.** Locityper's per-attempt `sol_ext..csv.br` is corrupt with 4 threads: brotli decodes only 3–4 genotypes.

**Noise scale.** `res.json`'s `lik_sd` is `Ln::to_log10(variance)` in the pinned source, not a standard deviation.
- Values of 100–570 imply per-attempt SDs of about 7–16 log₁₀ units, and standard errors of the 20-attempt mean of about
  1.5–3.6 log₁₀.
- Gaps of a few log₁₀ units are therefore within search noise.

## 4. Results (within each method only; truth labels are evaluation only)

**Locityper `lik_mean` (log₁₀), D pools:**

| cell | Locityper top (s1) | correct structure | stability |
|---|---|---|---|
| **C4 HG00171 LZO (control)** | **saved 1+2: wrong total (3)**, excess 32,748 | best 2+2 ranks 2nd, **93.8 behind**; the caller's 2+2 (excess 0) ranks 4th, 97.0 behind | the same-seed repeat is identical (STABLE) |
| **CYP HG02769 LZO (control)** | a 2+2 with excess 7 (count-selected) | the caller's 2+2 (excess 0) ranks 3rd, **2.8 behind**, within noise (SE ≈ 2.8) | no seed repeat by design |
| **C4 NA19650 LOO** | reported 3+1 | true 2+2 best at rank **5** under both seeds: **44.8 (s1) / 48.3 (s2)** behind | the top within the wrong-allocation class changes between seeds (margin 1.95 < seed difference 2.58: INCONCLUSIVE for the top); the allocation verdict is stable |
| **CYP NA19650 LOO** | reported 2+2: wrong total (4) | best 2+3 (excess 0) ranks 3rd, **62.0 behind**. A 1+3 sits between them, also total 4 | top and margin STABLE (margin 54.0 ≫ seed difference 1.25); rank changes only among near-tied 2+3s |

**Panvar on the same pools** (M7-B, nats; for within-method contrast only):

| cell | F + E | total |
|---|---|---|
| C4 LZO | caller first | caller first |
| CYP LZO | caller first | caller first |
| C4 LOO | 3+1 first; true 2+2 rank 5, −252.5 | 3+1 first; true 2+2 rank 5, −271.4 |
| CYP LOO | **2+3 first by 3.75** | 2+2 first by 3.64 |

**The two questions.**
- **CYP: does the alternative pipeline favour the correct total dosage?** **No.** Locityper ranks the wrong-dosage 2+2
  first, stably, ahead of every total-5 candidate by ≥ 62 log₁₀ units, which is far outside its noise. Panvar's fragment
  + exposure term favoured the correct total, and only its copying prior flipped it.
- **C4: does it improve allocation ranking despite missing content?** **No.** The true 2+2 stays at rank 5, as in Panvar.
  It is behind every 3+1 / 1+3 candidate by ≥ 44.8 log₁₀ units under both seeds.

**Correct-control deterioration.**
- **C4 LZO: yes, at the copy-number level.** A 1+2 (total 3) ranks first, far outside noise, and reproducibly. Panvar
  ranks it 9th, 733 nats behind.
- **CYP LZO:** the correct structure ranks first, but a sequence-different 2+2 leads the caller by 2.8 log₁₀ units, which
  is within noise. Deterioration is nominal, not established.

**Pool membership (diagnostics; never deployable claims).**
- **CYP D+X:** adding one haplotype pair shifted D genotypes by up to 9.0 log₁₀ units and reordered near-ties (true-class
  ranks 3, 4, 5, 7…). The top (reported 2+2) was unchanged.
- **C4 D+O:** D genotypes shifted by up to 26.4 log₁₀ units. **Oracle constructions ranked 1–3**: `m1e_e13_all`,
  `m1e_e13_T1`, `m1e_e13_U2`, the truth-patched content in the **1+3** arrangement. The true-class D candidates fell to
  ranks 21+.
- So Locityper's scores depend substantially on database membership. That membership changes mapping caps and
  best-relative alignment thresholds.

## 5. Interpretation and limitations

**On this fixed pool and this calibration,** Locityper's alignment + local-depth objective:
- did not improve the CYP dosage failure or the C4 allocation failure;
- introduced a copy-number error in the C4 control.

**Limitations that keep this from being a general verdict on the objective:**
1. **Synthetic, reduced calibration.** About 0.49 Mb of real sequence with an artificial junction at 252,518 | 252,519;
   no genome-wide uniqueness; and **56 % (C4) / 35 % (CYP) of target windows in extrapolated GC ranges**.
2. **Simulator-matched data.** No GC bias, no indel errors and no off-target reads.
3. **Binary identity.** The installed Locityper's source identity is unproven (checkout `7a6f2993`, same version string).
4. **Pool membership.** Membership measurably moves Locityper's scores. A larger or different pool could reorder D
   candidates.
5. **Small evaluation.** Four inspected development cells; one seed pair per failure cell; one same-seed control repeat.
6. **Not attributed.** These are whole-pipeline differences: placement retention, assignment versus marginalisation, the
   depth term, and calibration. Which component produced the control's copy-number error is not attributed.

No favourable result arose, so nothing here motivates attribution of a benefit. Nothing here concerns truth-free
candidate discovery or general genotyping.

## Addendum at M7-C (2026-10-20; the record above and its sealed package are unchanged)

M7-C's read-only audit (`tests/analysis/m7c/RESULT.md`) corrects four points:

1. **The attempt files are intact.** `sol_ext..csv.br` holds 4 concatenated Brotli streams; they are not corrupt. The
   direct gate (attempts 1–20 for every genotype) passes in all 11 runs.
2. **The "56 % / 35 % extrapolated GC" figure counts sparse individual GC bins.** Locityper's actual blurred tails
   (outside < 29 / > 67) cover 3.8 % (C4) and 2.1–2.2 % (CYP) of target windows.
3. **`lik_sd` is the variance in natural-log units divided by ln 10,** not "log₁₀ of the variance". The derived noise
   figures stand.
4. **Score terms added.** The depth term carries the wrong copy-number preference in the C4 control and the CYP failure.

The scoped conclusion is unchanged.
