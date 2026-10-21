# M12-A result: the rates and the copy-number arithmetic are consistent; with the true block content the objective prefers the true copy number in K1, K4 and K5; the errors come from the available alleles' content

**Run.** Branch `m12-a-rate-exposure-audit`; freeze `823a09d`; amendment A1 `7e36171`; result after.
- **Results:** `tests/results/m12a_rate_exposure_audit_20261021/`.
- **Budget:** 215.3 of 900 s: 150 s of disclosed development plus 65.3 s of ledger wall.
- **Units:** 18 recorded. `gate.fixtures` failed on a fixture-code defect, declared as A1 and re-run as `gate.fixtures.a1`,
  which passed. All other units passed first time. No watchdog or envelope event; peak 1.44 GiB.
- **Context read:** after execution I read the marker rates of all 30 sealed M11-D runs (§1). It is context only; no
  rule uses it.

**Not done:**
- no caller rerun, simulation, tuning, changed prediction, dosage-repair application, prior change, fit, validation family,
  merge, push or cleanup;
- the dosage audit was not repaired and KMAX was not raised.

**Truth was used only diagnostically** (the truth-in-context candidates, the noiseless count profile, provenance).

**M11-E closure:** appended to its RESULT with the local scope explicit:
- the tested replacements lose only in the winner's context;
- multi-block alternatives and other correct-copy-number alternatives are not excluded;
- K2 and K3 show that the tested alleles are inconsistent with the reads, not that no read-consistent allele exists.

## 1. The rate and exposure contract (`runs/contract/contract.tsv`; derivation in DESIGN §3)

**λ, the fragment-start intensity.** It is in templates per base per chromosome copy, both orientations summed; the signal
part is (1−η)λ.
- **Simulator:** declares λ_sim = depth / (4R) / (1−η). (1−η)λ_sim = 0.0750000 exactly reproduces each `simulation.tsv`
  "expected templates" (0.075 × (L − 349)).
- **5× subsample:** declares λ = λ_sim × 1/6. **declared / (sim × fraction) = 1.0000000000007 in all six cases.**
- **Read length:** 100 in the reads and in the simulation.
- **Error and insert:** e = 0.002, η = 0.001, insert 350 ± 40 (4σ), the same in the library and the simulation.

**Locus fragment counts.** The observed count against the simulator's own expectation (Σ expected templates × 1/6):
z = −0.70, −0.68, +0.74, 0.00, −0.58, −0.70.
- **The exposure term matches it exactly.** −exposure_term for the truth-in-context candidates equals that expectation
  (K1: 5,809.0 against 5,809.1).

**Marker rate (`marker_depth_rate_per_copy`).**
- **Expected:** 2 mates × (1−η)λ × (R − k + 1) × (1−e)^k = 1.6443 per copy (canonical closed syncmer 31-mers, exact, both
  mates).
- **Status ÷ expected:**

  | case | ratio |
  |---|---|
  | K4 | 0.999 |
  | K5 | 1.003 |
  | R1 | 0.987 |
  | R3 | 0.991 |
  | **K1** | **0.938** |
  | **R2** | **0.946** |

- **The two C4 cases fall outside the frozen 5%.**
- **Context, read after execution:** all 12 M11-D C4 runs average 0.923 (range 0.79–1.02). CYP2D6 averages 0.989 (6
  runs) and ACOT 0.968 (12 runs).

**Count depth (count-attrib, stable-copy 31-mers / 2).**

| case | depth | expected | ratio |
|---|---|---|---|
| K1 | 1.6036 | 1.6443 | 0.975 |
| K4 (M6-A) | 1.6522 | | 1.005 |
| K5 | 1.6592 | | 1.009 |

Every ratio is within 5%. The stable k-mers' truth multiplicity is 1.9995 (K1) and 1.9998 (K5), against the ideal 2.

**Reading:**
- **The intensity the objective uses is exact, and the exposure arithmetic is exact.**
- **The C4 marker rate is 5–8% low on average and variable.** It is an estimate from anchors assumed to be invariant.
- **It is not an input to the scored objective.** score_diplotypes uses only the declared library λ. The marker rate feeds
  proposal emissions, the junction rate and the dosage audit.
- **Its cause is unresolved.** One candidate is anchor multiplicities that do not hold for LOO truths in the C4 region. I
  did not test this, and I do not offer it as an explanation.

## 2. Analytic fixtures (`tests/m12a_copy_balance.cpp`; `runs/gate.fixtures.a1`): **3,976 of 3,976 checks pass**

| gate | result |
|---|---|
| G1 | 1,926 comparisons of origins on H_0..H_5: oracle = brute force = analytic (interior k, unit-junction k−1, flank/array boundary 1, flank-only 1) |
| G2 | scorer masses = oracle = independent n·½·P(i) (1,926 comparisons). Exposure = analytic Σ P(i)(L−i+1); one unit adds exactly U = 80. The fragment and exposure terms of all 21 candidate pairs = the independent Poisson-process formula |
| G3 | **noise-free tiling at the declared intensity:** the score peaks at the true total of 6 copies; (3,3) = (2,4) = (1,5) exactly. One copy more: fragment +72.7 against exposure −80.0; one copy fewer: −89.7 against +80.0 |
| G4 | two allele paths with the same bytes intern to one haplotype and score exactly as the homozygote |
| G5 | gapped model: an indel-ambiguous mate is **one** origin with up to 15 latent alignments (summed); two tandem copies give twice the origins (8 → 16); scorer = oracle on 63 gapped interior/junction masses |

**Rate sensitivity** (G3, reported, not gated):
- λ × 0.8 or × 0.9 → the argmax is **7** copies;
- × 1.1 → 6;
- × 1.25 → 5.

**The copy-number balance is correct and narrow.** In the fixture the true count wins by 7.3 nats per copy, and a 10% λ
underestimate overcalls. **In these cases λ is exact,** so this is not the cause here. It matters for real data, where λ is
estimated.

**A1:** the first run failed 48 checks. In every one, the oracle and the brute-force count agreed (1 origin) and only my
analytic class was wrong: flank-only templates on H_0. The fixture was corrected, not the implementation.

## 3. Real cases (`runs/analyse/audit.tsv`)

**Every likelihood run reproduced the incumbent within 1e-6,** and every decomposition reconciled with the scorer.

**Truth-in-context** = the winner's haplotypes with only the block's bytes replaced by the truth's (the better homologue
assignment is shown). Its likelihood (fragment + exposure) is defined; its prior is −inf.

| | K1 C4 NA24385 (+1 unit) | K4 CYP2D6 NA19650 (−1 unit) | K5 CYP2D6 NA19983 (+1 unit) |
|---|---|---|---|
| **truth-in-context − incumbent: Δlikelihood** | **+398.8** (fragment +69.1, exposure +329.6) | **+146.7** (+298.5, −151.8) | **+74.4** (−96.9, +171.2) |
| its multiplicity-scale δ (\|δ\| ≤ ln 2) + exposure | −159.0 + 329.6 = **+170.6** | +196.1 − 151.8 = **+44.3** | −123.4 + 171.2 = **+47.8** |
| its emission-scale δ | +228.1 (32 fragments) | +102.3 (23) | +26.6 (5) |
| **best-available (assign_A) − incumbent: Δlikelihood** | −55.7 | −24.8 | −132.7 |
| its multiplicity-scale δ + exposure | −149.9 + 329.6 = **+179.7** | +179.9 − 151.9 = **+28.0** | −214.0 + 151.7 = **−62.3** |
| its emission-scale δ | **−235.4 (54 fragments)** | **−52.8 (25)** | **−70.4 (14)** |
| block fragments: observed / expected under truth (z) | 1,680 / 1,644.9 (+0.86) | 775 / 847.1 (**−2.48**) | 675 / 676.3 (−0.05) |
| … / expected under incumbent (z) | / 1,974.6 (**−6.63**) | / 695.3 (+3.02) | / 847.5 (**−5.93**) |
| locus: observed / expected under incumbent (z) | 5,756 / 6,138.6 (−4.88) | 5,435 / 5,333.4 (+1.39) | 5,372 / 5,489.3 (−1.58) |
| count loss, real reads: best_A − incumbent / truth − incumbent | −19,131 / −21,551 | −4,057 (M6-A; best group = best_A; caller rank 164) / not measured | −7,019 / −12,159 |
| count loss, noiseless truth profile: best_A − incumbent | −26,808 | not measured | −9,642 |
| truth k-mers outside the panel dictionary | 244 | not measured | 41 |
| **frozen classification** | **INCONSISTENT_RATE_OR_CONVENTION** (C4 marker rate 0.938; see below) | **SEQUENCE_CONTENT_OR_PLACEMENT_MODEL_MISMATCH** | **SEQUENCE_CONTENT_OR_PLACEMENT_MODEL_MISMATCH** |

**About K1's label.**
- **Why it applies:** the frozen rule order tests every rate before anything else, and the C4 marker rate is outside 5%.
  The label stands.
- **Why it is not the cause:** that rate is not used by the objective; the objective's λ is exact; the locus count fits
  the truth (z = −0.70); and the count depth is within 2.5%.
- **What the rest of K1 shows:** every other K1 measurement matches K4 and K5. With truth content, the objective prefers
  the true copy number by +399, and the block count rejects the incumbent's extra unit (z = −6.6).
- **The frozen rule's limitation:** it did not separate rates that enter the objective from rates that do not. I record
  the limitation and do not relabel the case.

**Controls** (R1 CYP2D6, R2 C4, R3 ACOT): every incumbent is already the best-available pair.
- **Truth-in-context gains** +105.8, +100.3 and +35.2. The floor-level content differences cost likelihood even when the
  copy number is right.
- **Count losses agree:** the truth beats the incumbent by −178, −1,085 and −246.
- **R2** carries the same C4 marker-rate flag (0.946), so its frozen label is INCONSISTENT_RATE_OR_CONVENTION.

## 4. Diagnosis

**Implementation and normalization: consistent.**
- The fixtures pass and the scorer equals the oracle and the independent calculation.
- Equivalent representations collapse correctly.
- Latent alignments are summed, not multiplied.
- Exposure is linear per unit.
- Every real decomposition reconciles.

**The rates and conventions the objective uses: consistent.**
- λ equals the simulator's × the subsample fraction exactly.
- The exposure term equals the simulator's own expected fragment count.
- The count depth agrees within 2.5%.
- **The one deviation is the anchor-based marker rate in C4** (K1, R2 and the M11-D C4 runs). It is outside the scored
  objective. Its cause is unresolved and is not offered as an explanation of these errors.

**The copy-number balance itself: in the true copy number's favour in all three.**
- Changing the block to the truth's content, the multiplicity-scale gain plus exposure favours the true copy number:
  K1 +171, K4 +44, K5 +48.
- **In K1 and K5 the block fragment counts clearly reject the incumbent** (z = −6.6 and −5.9) **and fit the truth**
  (+0.9, −0.05).
- **K4 is the case where sampling contributes.**
  - Its block has 775 fragments: 2.5 σ below the truth's expectation and 3.0 σ above the incumbent's.
  - Even so, the truth content wins by +147.
  - The best-available pair's balance is only +28, so it loses to its own 53-nat content penalty.

**What decides the calls: the available alleles' content** (sequence-content mismatch).
- In K1 and K4 the best-available pair's emission-scale penalty is larger than its copy-number advantage.
  - K1: 54 fragments cost −235 against a balance of +180.
  - K4: 25 fragments cost −53 against +28.
- In K5 the best-available allele (15,327 bp, 1,621 edits from truth2) is wrong in content and length.
  - Its multiplicity balance is negative too (−62), whereas the truth's is +48.
- The incumbent's extra or missing unit lets reads that carry truth content place somewhere without a mismatch.
  - That is legitimate under the model: reads covering sequence the best allele lacks are explained by the incumbent's
    other unit copies.
  - The truth's own content removes the trade-off.
- **This is the classification `SEQUENCE_CONTENT_OR_PLACEMENT_MODEL_MISMATCH`, and specifically the content part.** The
  placement model passes its fixtures and the true-content comparison.

**Count evidence** (reported beside the likelihood, never added to it) **prefers the best-available pair over the
incumbent in all three cases.**
- Real reads: −19,131 / −4,057 / −7,019.
- The noiseless truth profile gives the same.
- Squared k-mer residuals weigh a unit's dosage excess, about 28,000 interior k-mers in K1, far above a few dozen
  content-mismatched reads.
- So the two models rank the same candidates differently because they weigh dosage and content differently. Neither
  computation is inconsistent.

**Unresolved:** the C4 marker-rate deficit (outside the objective), and whether panel content can come close enough to the
truth's.

## 5. Recommendation

**No rate, exposure, normalization or depth-recalibration correction is supported.** The rates agree, the analytic
fixtures agree, and the objective prefers the true copy number when given the true content. I do not propose a depth
recalibration.

**One correction experiment is supported:** test whether correct-copy-number candidates whose block content is closer to
the reads can be built from panel material, and score them with the unchanged objective.
- **What to build:** unit-level recombinations of the panel's unit alleles inside the `repeated_unit` block.
- **The gap to close:** the emission-scale penalties of the current best-available pairs — K1 235 nats over 54 fragments,
  K4 53 over 25, K5 70 over 14.
- **What success requires:** the truth-content margins (+399 / +147 / +74) show that a candidate with content close
  enough would win under the current objective. That is no change to the likelihood.
- **Caution:** M1-D found panel-only unit composition recovered only 3–6% of the truth-block gain for C4 NA19650 sv:7.
  The experiment may show the content is simply unavailable in the panel, which would put these cases with K2 and K3 as
  availability limits.

**Separate, lower-priority item:** the C4 anchor marker-rate deficit (0.79–1.02). It affects proposals and the dosage
audit, not the score.
