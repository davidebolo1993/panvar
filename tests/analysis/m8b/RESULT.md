# M8-B result: on the masked observed-block score, every locus's training profile favours ~10× more donor redraw than legacy; ACOT held-out families confirm it for the constant mode; C4 and CYP2D6 have no clean evaluation population

**Run.** Branch `m8-b-copying-calibration`: freeze `90ac3c5`, Amendment 1 (gate fixture naming). Results in
`tests/results/m8b_copying_calibration_20261020/`.

**Execution.** 151.7 of 300 s. One recorded failure: `gates`, on a fixture naming defect, re-run as `gates.a1`.

**Scope.** No caller integration, genotyping, default change, mutation model, clustering, merge or push. M8-A is kept as
the historical record.

## 1. Corrected inventory (`inventory`; frozen before any parameter scoring)

**Sources.** Authoritative sample-level records under `tests/results` (per-donor directories, LOO model files, mosaic
`stage.log` donors, simulations, campaign signatures) plus M8-A's list, giving **69 development donors**.
- NA19705 is confirmed: CYP2D6 and C4 LOO models, cohort signatures.
- HG03520 is confirmed: C4 LOO models, coverage and cohort signatures.
- Status is expanded to whole IGSR pedigree families.

| locus | families | development | reference genomes | uninspected | outer evaluation |
|---|---|---|---|---|---|
| C4 | 61 | **59** | 2 | **0** | **INFEASIBLE** (< 4 uninspected families) |
| CYP2D6 | 61 | **59** | 2 | **0** | **INFEASIBLE** |
| ACOT | 235 | 13 | 3 | 219 | **30 families (60 haplotypes) by sha256 rule** |

**M8-A's "uninspected" C4 and CYP2D6 families were all previously used** (HG02018, HG03520, HG04036, NA19705, NA19836,
NA20355). Disclosed: M8-A had already examined support and minimum-switch summaries over that population.

## 2. Gates (`gates.a1`: all pass)

| gate | result |
|---|---|
| G1 masked forward vs exhaustive tiny-history oracle (200 fixtures, including duplicate rows) | 3.7e-15 |
| G2 no masks = `log_haploid_mosaic_prior` | 0 |
| G3 all masked = log 1 | 7.8e-16 |
| G4 masked interior block = removed block with composed transition 1 − (1 − r₁)(1 − r₂) | 8.9e-16 |
| G5 donor-row permutation; duplicating every row unchanged (declared uniform row weighting) | 4.4e-16; 1.8e-15 |
| **G6** AB \| C control: whole-sequence prior finite (log 0) while the labelled path (AB, C) is −∞ | PASS |

## 3. Whole-sequence support (`support`, 25 s per locus)

`log_haploid_sequence_prior` was run on same-graph models with the eval templates and the queried family excluded:

| locus | role | haplotypes | tested | whole-sequence supported | unsupported | untested (time bound) |
|---|---|---|---|---|---|---|
| C4 | training queries | 131 | 69 | **0** | 69 | 62 |
| CYP2D6 | training queries | 126 | 77 | **0** | 77 | 49 |
| ACOT | eval | 60 | 14 | **0** | 14 | **46** |
| ACOT | training queries | 407 | 0 | — | — | 407 |

- **Every tested haplotype is whole-sequence unsupported** under the copying-only prior, including cross-block spellings.
- Untested cases are reported, not extrapolated.
- The masked score below is therefore **not** a predictive probability of these haplotypes.

## 4. Training profiles

The training objective is the masked labelled-path score: inner leave-one-family-out, every training haplotype, masks
recomputed per bank, no −∞ cases. Values are differences from each mode's best, in nats.

| r | C4 | CYP2D6 | ACOT |
|---|---|---|---|
| 0.0003 | −420.5 | −1,486.9 | −3,929.9 |
| 0.001 | −319.8 | −1,150.8 | −3,016.9 |
| 0.003 | −228.7 | −845.9 | −2,190.1 |
| **0.01 (legacy)** | **−132.0** | **−518.4** | **−1,306.7** |
| 0.03 | −52.6 | −238.3 | −566.5 |
| 0.1 | **0** | −0.5 | **0** |
| 0.3 | −65.2 | **0 (BOUNDARY)** | −269.8 |

| λ / λ_HM | C4 (λ_HM 1.87 Mb) | CYP2D6 (1.13 Mb) | ACOT (1.22 Mb) |
|---|---|---|---|
| 1/30 | −24.6 | **0 (BOUNDARY)** | −23.5 |
| 1/10 | **0** (187 kb) | −106.0 | **0** (122 kb) |
| 1/3 | −57.2 | −360.4 | −595.1 |
| 1 | −136.9 | −643.4 | −1,339.2 |
| 3 | −225.0 | −942.7 | −2,145.2 |
| 10 | −324.6 | −1,276.9 | −3,051.1 |
| 30 | −416.5 | −1,583.6 | −3,884.3 |

**Profile shape:**
- Every profile is steep below the optimum, so it is far from flat.
- **C4 and ACOT** have interior optima: r = 0.1 and λ = λ_HM/10.
- **CYP2D6** selects grid edges in both modes. Its constant profile is nearly flat between 0.1 and 0.3 (0.5 nats), and the
  selection is **family-sensitive**: removing any one of 21 of the 61 training families changes it. Per the stopping rule
  the CYP2D6 grid is **not extended**, and its selection is reported as BOUNDARY / unresolved.
- Selections at C4 and ACOT, and the CYP2D6 distance selection, are unchanged by removing any single family.

**Minimum-switch counts (M8-A)** were a feasibility signal only. These profiles show the score responds to the parameter,
and C4, despite few required switches, still prefers r ≈ 0.1.

**Isolation.** Training profiles and per-family contributions are **byte-identical** with every eval row scrambled.
- For C4 and CYP2D6 this is vacuous, since there are no eval rows.
- **Fixed and evaluation-derived:** the graph frame (blocks, allele labels, bubbles from the full panel, eval families
  included) and the reference-core geometry used for distances.

## 5. Outer evaluation (ACOT only; 30 families, 60 haplotypes, frozen training bank of 205 families)

All 30 families are included.
- None is all-masked: mask fractions 0.18–0.45 (median 0.26).
- Every haplotype has observed blocks.

| comparison | Δ total | per haplotype | families positive | largest positive share | meaningful (≥ 1 nat/haplotype) | consistent |
|---|---|---|---|---|---|---|
| fitted constant (r = 0.1) − legacy (r = 0.01) | +203.8 | **+3.40** | **29/30** | 0.08 | **yes** | **yes** |
| fitted distance (λ = 122 kb) − fitted constant | −19.6 | −0.33 | 10/30 | 0.24 | no | no |

Per-family Δ for constant vs legacy ranges from −2.9 to +16.8 (median +6.8); only family VN066 is negative.

## 6. Interpretation and stopping

**On this masked observed-block score,** the panels prefer a donor-redraw parameter about 10× the legacy value (r ≈ 0.1
per boundary). This is supported in training at all three loci and on held-out ACOT families. Distance-based redraw does
not improve on a fitted constant at ACOT.

**Limits, all binding:**
- **The score is conditional on the masking procedure.** About 26–41 % of blocks are masked, and every tested whole
  sequence has probability 0 under the copying-only prior. A higher r may partly reflect the copying model accommodating
  panel incompleteness, rather than a property to transfer to genotyping. r stays a copying-model parameter, not a
  crossover rate.
- **C4 and CYP2D6 have no uninspected families.** Their selections are training-only (CYP2D6 BOUNDARY and
  family-sensitive).
- **The graph frame and geometry are derived from the full panel.**
- **ACOT is the only locus with an outer evaluation:** 30 families, one panel.

**What this authorises.** Only **consideration** of a later genotyping validation of r ≈ 0.1 (constant mode), with
controls, at ACOT first. No parameter or default is promoted.
- Nothing here addresses C4's missing sequence.
- Nothing was tuned to the CYP2D6 NA19650 failure, which remains a later diagnostic only.

**Stopping.** No grid extension, no re-fitting. Results returned.
