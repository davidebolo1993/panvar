# M7-C: read-only audit and closure of the fixed-pool comparison (addendum to M7-B.1; sealed packages unchanged)

**Run.** Branch `m7-c-audit`, audit script `731b61e`. Results in `tests/results/m7c_audit_20261020/`.

**Execution.** 1.9 s of new execution (of 300). No genotyping, weight or calibration change, Panvar change, or campaign
expansion. The M7-B and M7-B.1 packages are read only.

**Reuse of Sol's script.** Its decoding approach, `audit_streams.py` (libbrotlidec through ctypes, large-window parameter,
one decoder per member, complete-consumption assertion), was **reviewed before reuse** and re-implemented in
`brotli_members.py`. The gates, negative controls, GC reproduction and contrasts are new code in `audit.py`.

## Recovered facts (from the archived run outputs)

**F1. The attempt files are intact.**
- Every `sol_ext..csv.br` of all 11 M7-B.1 Locityper runs consists of **4 concatenated Brotli streams**, decoded with
  large-window support, with **all input bytes consumed**.
- **The original direct gate is restored and passes for all 11 runs:** every declared genotype has exactly attempts 1–20 at
  stage 2 (11, 11, 29, 11, 11, 11, 13, 13, 14, 13 and 13 genotypes).
- **Negative controls, all detected in every run:**
  - truncating the last 7 bytes: decode error;
  - truncating mid-way through stream 2: decode error;
  - removing one stage-2 attempt: the gate fails;
  - duplicating one attempt: the gate fails.
- **Reconstructions match `res.json.gz`:**
  - the attempt-mean `lik` equals `lik_mean` (max error 3.0e-8);
  - `lik_sd` equals the attempt **sample** variance × ln 10 (max error 1.4e-6; the population variance is off by 17–30);
  - every attempt satisfies `lik = 0.15·aln_lik + 1.85·depth_lik` (max error 1.45e-7; prior 0).
- **The M7-B.1 statement "sol_ext is corrupt with 4 threads" was wrong.** The `brotli` CLI decoded only the first stream.
  Amendment 1's substitute check (logged `a=20` plus finite `lik_sd`) is superseded by the direct gate, which passes.
  The two control units were not re-run.

**F2. GC support: three distinct categories.** Target windows have weight ≥ 0.001, the same inclusion rule as M7-B.1.
- **Locityper's cumulative-tail boundaries** (`bg/depth.rs`, `blur_boundary_values`; `min_tail_obs` = 100): reproduced
  from the per-GC background window counts as **< 29 or > 67**, identical to the boundaries logged by all four `preproc`
  runs.
- **Observed background GC range:** 23–79 %.

| runs | target windows | (1) in sparse individual GC bins (< 100 windows) | (2) outside the cumulative-tail boundaries (blurred) | (3) outside the observed background range |
|---|---|---|---|---|
| C4 LOO D s1 / s2 | 62,609 | 35,256 (56.3 %) | **2,395 (3.83 %)** | 260 (0.42 %) |
| C4 LOO D+O | 144,306 | 81,375 (56.4 %) | 5,582 (3.87 %) | 626 (0.43 %) |
| C4 LZO (s1, .a1, s1r) | 37,661 | 21,259 (56.4 %) | **1,417 (3.76 %)** | 144 (0.38 %) |
| CYP LOO D s1 / s2 | 48,586 | 17,114 (35.2 %) | **1,085 (2.23 %)** | 8 (0.02 %) |
| CYP LOO D+X | 55,133 | 19,436 (35.3 %) | 1,230 (2.23 %) | 9 (0.02 %) |
| CYP LZO (s1, .a1) | 45,008 | 15,739 (35.0 %) | **961 (2.14 %)** | 8 (0.02 %) |

**F3. Score terms at the recovered attempts.** Each value is the mean over 20 attempts; all values are in `lik` units
(log₁₀). The weighted decomposition reconstructs each Δ to ≤ 2.4e-8.

| run | contrast (candidate − Locityper top) | Δ total | 0.15 · Δ alignment | 1.85 · Δ depth |
|---|---|---|---|---|
| C4 HG00171 LZO s1 (= .a1 = s1r) | best correct-class 2+2 (`m6a_g2_m1`, rank 2) − top 1+2 | −93.81 | **+47.77** | **−141.59** |
| | **exact caller** 2+2 (rank 4) − top 1+2 | −97.03 | **+47.75** | **−144.78** |
| CYP HG02769 LZO s1 | best correct-class = top (rank 1; a 2+2 with excess 7) | 0 | 0 | 0 |
| | **exact caller** 2+2 (excess 0, rank 3) − top | −2.82 | +2.14 | −4.95 |
| C4 NA19650 LOO s1 | best correct-class 2+2 (`m6a_g4_m1`, excess 0, rank 5) − top 3+1 (= the caller) | −44.82 | −27.40 | −17.41 |
| C4 NA19650 LOO s2 | same candidate (rank 5) − top (`m6a_g3_m1`, 1+3) | −48.29 | −24.35 | −23.94 |
| | exact caller (rank 2) − top | −0.92 | +3.62 | −4.54 |
| CYP NA19650 LOO s1 | best correct-class 2+3 (`sub_truth_b`/`m6a_g3_m1`, excess 0, rank 3) − top 2+2 (= the caller) | −62.01 | **+0.68** | **−62.69** |
| CYP NA19650 LOO s2 | same candidate (rank 3) − top (the caller) | −60.64 | +1.12 | −61.76 |

**F4. Variance semantics.** For every audited genotype, `lik_sd` equals the recovered attempt variance (log₁₀ units²) ×
ln 10. Equivalently, it is the variance in natural-log score units divided by ln 10. Among the audited genotypes (the top 3 per run), per-attempt SDs are 7.9–12.4 log₁₀
units, and standard errors of the 20-attempt means are 1.8–2.8.

## Source-based interpretation

**I1. Variance semantics.** In the inspected source, `lik_sd: Ln::to_log10(pred.lik_var)`, with `Ln::to_log10(x) = x / ln 10`
applied to a natural-log **variance**. The exported value is therefore *variance in natural-log units divided by
ln(10)*, matching F4.
- **Correction:** M7-B.1 described it as "`Ln::to_log10(variance)`" and as "log₁₀ of the variance". The phrase is
  misleading and is withdrawn. M7-B.1's derived noise figures (SDs about 7–16, SEs about 1.5–3.6 log₁₀ units) used the
  correct arithmetic and stand.
- **This spread is variability across stochastic optimisation attempts.** It is not calibrated genotype uncertainty, and
  not a posterior or sampling error of the data.

**I2. GC.**
- **Category 1 was mislabelled.** The 56 % / 35 % figure M7-B.1 reported as "relying on extrapolated GC" counts windows
  whose **individual** GC bin is sparse. Inside 29–67 % those bins take LOESS-fitted values, and Locityper does not blur
  them.
- **Category 2 is Locityper's actual extrapolation:** constant means and inflated variances beyond its logged boundaries,
  at about 3.8 % (C4) and 2.1–2.2 % (CYP) of target windows.
- **Category 3** windows lie beyond any observed background GC.
- **None of this establishes adequate calibration.** The background is still a ~0.49-Mb synthetic set from real sequence,
  with an artificial junction at 252,518 | 252,519 and no genome-wide uniqueness.

**I3. Score terms.**
- These are **decompositions at jointly optimised read assignments.** Each genotype's assignment maximises the combined
  objective, so the terms are not counterfactual alignment-only or depth-only results.
- **Within that qualification:**
  - **C4 control and CYP failure:** the depth term supplies the preference for the wrong copy number (−141.6 and
    −62.7 / −61.8). The alignment term slightly favours the correct candidate (+47.8 and +0.7 / +1.1). In both cases the
    wrong candidate has **fewer** copies (total 3 against 4; 4 against 5).
  - **C4 allocation failure:** both terms favour the wrong allocation under both seeds.
  - **CYP control:** the caller's deficit (−2.8) is within the attempt noise (SE about 2.8).
- **Which part of the depth model is responsible is not established.** Candidates include the CN-posterior
  normalisation, GC-dependent NB parameters, window weights, the calibration set, or the interaction with assignment.

## Corrected conclusions (closure)

**Scoped conclusion, unchanged:** in this fixed-pool experiment, under the stated synthetic calibration, Locityper's
objective
- **did not improve** the CYP NA19650 dosage failure or the C4 NA19650 allocation failure;
- **deteriorated** the C4 HG00171 control's copy-number prediction.

**Corrections to M7-B.1:**
- The attempt files were never corrupt, and the direct 1–20 attempt gate passes.
- The GC "extrapolation" figure is replaced by the three categories in F2 / I2.
- The `lik_sd` wording is corrected as in I1.
- Score-term attribution (F3 / I3) is added.

**Limitations retained:**
- the reduced, synthetic, simulator-matched background;
- sensitivity of Locityper's scores to pool membership (D+X shift up to 9.0, D+O up to 26.4 log₁₀ units);
- the installed Locityper binary's source identity (checkout `7a6f2993`) remains unestablished;
- four inspected development cells.

**Not claimed:** nothing here generalises to Locityper overall, or to all local-depth models. Nothing here bears on
truth-free discovery or general genotyping.
