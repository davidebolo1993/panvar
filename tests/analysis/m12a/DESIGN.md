# M12-A: repeated-unit placement/exposure consistency audit (frozen before execution)

**Status of earlier batches:** M11-E is closed with its local scope appended to its RESULT. M10-A plus M11-C remains the
caller checkpoint.
- **Branch:** `m12-a-rate-exposure-audit`, from `m11-e-catastrophic-attribution` `05dd39b`.
- **Results:** `tests/results/m12a_rate_exposure_audit_20261021/`.
- **Budget:** 900 s, including 150 s of disclosed development:
  - code tracing and provenance reads;
  - the probe and fixture builds, about 45 s of compute;
  - read-length scans.
- **Resources:** 6 GiB; one scoring process at a time; 4 threads.

**Not done:**
- no caller rerun, simulation, parameter tuning, changed prediction, production sequence or dosage-repair application;
- no prior change, clustering, confidence fit, validation family, merge, push or cleanup;
- the whole dosage audit is not repaired, KMAX is not raised, and the enumeration is not expanded.

**Truth is used only for diagnostic checks.** It never fits a rate, never selects a feature, and never enters a prediction.

## 1. Code changes (diagnostic only)

- **`genotype_probe count-attrib --skip-catalogue 1`:** scores only the frozen candidates. The catalogue table records that
  it was skipped; nothing else changes.
- **`tests/m12a_copy_balance.cpp`:** a new diagnostic fixture executable. It is not part of ctest.
- **Caller and scorer:** untouched.

## 2. Frozen cases

| case | cell | block | role | count evidence |
|---|---|---|---|---|
| K1 | c4.NA24385 | sv:7 | catastrophic, one unit too many | count-attrib, this batch |
| K4 | cyp.NA19650 | sv:9 | catastrophic, one unit too few | **M6-A count-fit, reused** (same model fingerprint `e41f77fd…`, same 5,435 read pairs, same winner context; its fragment, exposure and prior deltas equal M11-E's) |
| K5 | cyp.NA19983 | sv:9 | catastrophic, one unit too many | count-attrib, this batch |
| R1 / R2 / R3 | cyp.HG02953 sv:9 / c4.HG00358 sv:7 / acot.HG03209 sv:13 | | M11-E controls | count-attrib, this batch |

**K2 and K3 stay out of scope.** They are a separate content-mismatch problem and are not attempted here.

## 3. The rate and exposure contract (derived from the code; checked in `contract`)

**Where each quantity comes from:**
- **λ (`fragment_start_intensity`)** = templates per base **per chromosome copy**, both library orientations summed
  (`genotype_model.hpp`, `fragment_likelihood.hpp`). The signal intensity is (1−η)λ.
  - **Simulator:** `genotype_probe simulate` declares λ_sim = depth / (4R) / (1−η), so that (1−η)λ_sim × 2 copies × 2 mates
    × R = depth. The template count is 0.075 × E per copy at 30×, matching `simulation.tsv`.
  - **Subsample:** declares λ = λ_sim × requested fraction (fnv-nested selector).
- **Exposure E_h** = Σ_{i ≥ floor} P(i)·(L_h − i + 1)⁺, averaged over the fragments' physical floors (the longer mate,
  100 bp). That floor never binds when the insert support is 190–510 bp.
  - For L ≫ insert, E_h = L_h + 1 − μ, so **one added unit of U bp adds exactly U to E**.
  - The score's exposure term is −(1−η)λ(E_a + E_b): **the expected number of signal fragments** for that candidate.
- **Fragment term:** Σ_f log((1−η)λ(M_a,f + M_b,f) + η·P_bg(f)).
  - M_h,f = Σ over physical origins of ½ × P(insert) × the two mates' normalised emissions.
  - **Distinct repeat copies are distinct origins.** Latent alignments of one origin are summed into that origin, never
    counted as extra origins.
- **`marker_depth_rate_per_copy`:**
  - the Poisson MLE of the expected exact observations of one marker per chromosome copy, from invariant anchors (count ~
    Poisson(bg + rate × 2 × m_k));
  - markers are canonical closed syncmer 31-mers (s = 15), counted in **both mates** with exact matching.
- **Expected relationship (interior positions):**
  rate_per_copy = 2 mates × (1−η)λ × (R − k + 1) × (1−e)^k, where k = 31, R = 100 and e = 0.002.
  - (R − k + 1) is the number of read starts that cover a k-mer.
  - (1−e)^k is the probability the k-mer is error-free.
  - Strand does not matter, because codes are canonical.
  - Boundaries reduce coverage only within one insert of a haplotype end, so anchors and stable k-mers are interior.
  - **At 5×:** 2 × 0.0125 × 70 × 0.9399 = 1.645.
- **M6-A count depth:** `depth_haploid_kmer_mean` = the mean count of stable-copy 31-mers / 2. The same expectation
  applies, scaled by (stable k-mers' truth multiplicity) / 2. Marker depth is never compared with λ without this conversion.

**Contract checks:**
- λ_declared / (λ_sim × fraction) = 1 within 1e-6;
- read length equals the simulated one;
- marker rate / expected within 5%;
- count depth / expected within 5%;
- locus fragment count against the simulator's own expectation (Σ expected templates × fraction): |z| ≤ 4.

## 4. Fixtures (`gate.fixtures`; deterministic; frozen in `tests/m12a_copy_balance.cpp`)

**Setup:** identical tandem units (80 bp), random flanks (60 bp), 20 bp mates, inserts of 40 ± 3, and (1−η)λ = 1 template per
base per copy.

**What the fixtures are compared against** — an independent calculation:
- brute-force origin counts;
- the analytic per-origin mass, ½·P(i), since band 0 gives emission 1;
- the analytic exposure;
- the Poisson-process score.

| gate | tests |
|---|---|
| G1 | origins on H_0..H_5 for interior / unit-junction / boundary fragments: oracle = brute force = analytic (k, k−1, 1) — adding and removing copies |
| G2 | scorer masses = oracle = independent; exposure = analytic; ΔE per unit = U; fragment and exposure terms = independent, for every candidate pair |
| G3 | noise-free tiling of the truth (H_3, H_3) at the declared intensity: the score peaks at the true total copies, and equal-total allocations tie exactly. Rate sensitivity is reported, not gated |
| G4 | equivalent representations: two allele paths with the same bytes are one haplotype and score as the homozygote |
| G5 | gapped model: an indel-ambiguous mate is one origin with several latent alignments; two tandem copies give two origins; scorer = oracle masses on gapped interior and junction fragments |

**Fixture failures are not repaired after the fact.** A failure is either an implementation inconsistency or, if it is a
defect in the fixture code, a declared amendment.

## 5. Real-case comparisons and classification (frozen)

**Likelihood (`lik.<case>`, `genotype_probe score-set`)** — in the winner's context:
- the incumbent, with its expected total, which must reproduce;
- assign_A (M11-E's best-available pair);
- **truth-in-context**: the incumbent's haplotypes with only the block's bytes replaced by the truth's, in both assignments
  (`--explicit`). Its prior is −inf, so only its likelihood (fragment + exposure) is used.

Each candidate is decomposed against the incumbent, and every decomposition must reconcile with the scorer.

**Per-fragment scale split** (descriptive):
- Changing the copy count k → k±1 (k ≥ 1) moves a fragment's log contribution by at most ln 2. |δ| ≤ ln 2 is called the
  multiplicity scale.
- |δ| > ln 2 is the emission scale; one substitution is ≈ 7.3 nats at e = 0.002.

**Observed vs expected fragments:**
- **Locus:** the observed count against −exposure_term for each candidate.
- **Block:** fragments with a block placement on a candidate (decomposition classes) against
  (1−η)λ × (Σ_h block bp + 2(μ − 1)). That is the expected count of templates overlapping the block.

**Count evidence (`count.<case>`; `genotype_probe count-attrib --skip-catalogue 1`):**
- the loss of incumbent, best_A and the truth block in both assignments (oracle) — under real reads with the stable-copy
  depth, and under the noiseless truth profile (sensitivity to missing candidate content);
- the contrast decomposition by k-mer location; k-mers differing (marker availability); truth k-mers outside the dictionary;
  the stable k-mers' truth multiplicity;
- K4 reuses M6-A's count-fit (depth, best group, caller rank).

**Count and likelihood scores are never added**, because their evidence overlaps.

**Classification (frozen order):**
1. **`IMPLEMENTATION_OR_NORMALIZATION_INCONSISTENCY`:** a fixture gate fails, the incumbent does not reproduce, or a
   decomposition does not reconcile.
2. **`INCONSISTENT_RATE_OR_CONVENTION`:** a §3 contract check fails.
3. **Controls:** `CONTROL`. The numbers are reported.
4. **`SEQUENCE_CONTENT_OR_PLACEMENT_MODEL_MISMATCH`:** the truth-in-context likelihood exceeds the incumbent's. With the
   true block content, the placement/exposure balance prefers the true copy number, so the available alleles' content
   decides.
5. If even the truth-in-context likelihood loses, the block's observed fragment count is compared against the truth
   length's expectation:
   - it deviates toward the incumbent's length with |z| ≤ 3: **`SAMPLING_COMPATIBLE_DISAGREEMENT`**;
   - it favours the true length: **`SEQUENCE_CONTENT_OR_PLACEMENT_MODEL_MISMATCH`** (the placement model);
   - otherwise: **`UNRESOLVED_CAUSE`**.

**A depth-recalibration explanation** is offered only if the rate contract fails. If the rates and fixtures agree, the
result says so.

## 6. Units (`m12a_run.py`; dependency-gated; `next` stops on any failure)

| unit | needs | declared timeout (s) |
|---|---|---|
| `build` | | 300 |
| `gate.fixtures` | build | 120 |
| `contract` | | 120 |
| `prepare` | | 60 |
| `lik.K1 K4 K5 R1 R2 R3` | build, prepare | C4 / CYP2D6 120, ACOT 240 |
| `count.K1 K5 R1 R2` | build, prepare | 180 |
| `count.R3` | build, prepare | 240 |
| `analyse` | contract, prepare | 120 |

**Expected time:** about 150–300 s.

**Admission:** spent + timeout ≤ 900 s.

**Failures:**
- a failed case unit is recorded as a limitation, and the other units continue;
- an analysis-code failure needs a declared amendment.
