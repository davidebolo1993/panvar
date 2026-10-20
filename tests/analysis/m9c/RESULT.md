# M9-C result: on six new ACOT families every prediction completed under sequence_first_v2; outside sv:5 the block error masses track exact agreement, but sv:5 is wrong in all six LOO cells at extreme confidence

**Run.** Branch `m9-c-acot-reliability`.
- **Commits:** freeze `c3df6eb`, which includes the M9-B qualification; Amendment 1 `ca317cd`, an analysis-only defect.
- **Results:** `tests/results/m9c_acot_reliability_20261020/`.

**Execution: 604.3 of 2,400 s.**
- No caller guard or watchdog bound.
- No unit was refused or substituted.

**Recorded defect:** the first `analyse` misclassified its high-confidence error list (Amendment 1). It is kept;
`analyse.a1` is authoritative.

**Scope:** descriptive only. No fitting, GQ remapping, PASS rule, prior tuning, clustering or promotion.

## 0. M9-B corrected first

The qualification is appended to `tests/analysis/m9b/RESULT.md`. It covers:
- legacy no-calls versus sequence_first_v2 predictions;
- the existing projection-aware `error_mass_lower`/`error_mass_upper`;
- the v2 FASTA policy;
- the narrower C4 allocation evidence: weak, non-zero discrimination among the examined alternatives, with no
  impossibility proof;
- representation intervals, which are conservative enclosing intervals, not minimal ones.

## 1. Pilot as frozen

**Families.** The first six by sha256("M9-C:" + family), from the 28 ACOT outer-evaluation families never genotyped
before. The donors are:

| cell | donor | family |
|---|---|---|
| d1 | NA19909 | 2427 |
| d2 | NA18747 | NA18747 |
| d3 | HG02178 | HG02178 |
| d4 | HG00639 | PR03 |
| d5 | NA21106 | NA21106 |
| d6 | HG03831 | BD27 |

- **Disclosed:** these assemblies took part in M8-B's score-level outer evaluation. This is new genotype-confidence
  evidence, not a pristine holdout.

**Cells:**
- six LOO at 5×;
- LZO controls for d1 and d2, with byte-identical reads.

**Panels:**
- M8-C's recipe on the full-panel graph frame (`acot.gfa`; blocks derived from all 467 templates; disclosed);
- **LOO** excludes all 60 outer-evaluation templates (asserted identical to M8-C's list);
- **LZO** adds only the donor's two haplotypes.

**Route.**
- **Executables:** panvar `98f339e8…` and genotype_probe `af2ed281…`, both asserted equal to M8-C's.
- **Profile:** M8-C's r = 0.01 profile (validated combined route, re-derived byte-identically) plus
  `output_schema sequence_first_v2` and `repair_joint_search off`.
- **Command line:** `--dosage-mode audit`, four threads.
- **Guards:** 6 GiB and 900 s, with a 930-s process-group watchdog.
- **Settings check:** every call's settings match, and none is unverified.
- **Status:** this is not "untouched production defaults".

**Contract gates:**
- `genotype_sequence_report_test` passed: the coupled label ambiguity yields a v2 PREDICTION, the mixed-bound fixture
  holds, and panel/off-panel FASTA behave as specified.
- The e2e `sequence_first` mode passed, after both binaries were relinked against the current library.

## 2. Callability (frozen before evaluation)

| cell | v2 outcome | legacy | phase support | wall / peak |
|---|---|---|---|---|
| LOO d1, d2, d3, d4, d6 | PREDICTION, RESOLVED, 19/19 blocks IDENTIFIED | COMPLETE | PRIOR_ONLY | 40–74 s / ≤ 1.93 GiB |
| **LOO d5** | **PREDICTION, SEQUENCE_TIED** (2 decision-set groups): 18 IDENTIFIED, 1 context block GENOTYPE_UNRESOLVED | **RESCORE_INCOMPLETE, `divergent_tie` no-call** | UNRESOLVED | 47.7 s / 1.98 GiB |
| LZO d1, d2 | PREDICTION, RESOLVED, 19/19 IDENTIFIED | COMPLETE | PRIOR_ONLY | 27–29 s / ≤ 2.13 GiB |

**Prediction rates:**
- 8/8 cells predicted.
- **72/72 reportable blocks IDENTIFIED.**
- No representation intervals, and `error_mass_lower = error_mass_upper` at every block, so no MIXED groups occurred.

**d5** is the case M9-B's qualification describes: a legacy no-call that still carries a complete v2 prediction for all
reportable blocks.

**FASTA policy as specified:**
- **LOO:** off-panel pairs (d4 is panel + off-panel), so `.sequence_haplotypes.fa` was emitted.
- **LZO:** panel pairs, so no FASTA, with extraction references instead.
- **Legacy FASTA:** absent in every cell.

## 3. Reliability at reportable blocks with IDENTIFIED genotypes (`error_mass_upper`, predeclared bins)

| scheme | stratum | bin | blocks | exact | not exact | families |
|---|---|---|---|---|---|---|
| LOO | floor-zero (both truth alleles in panel) | (0, 1e-6) | 34 | 28 | **6** | 6 |
| LOO | floor-zero | [1e-6, 1e-3) | 7 | 7 | 0 | 4 |
| LOO | floor-zero | [1e-3, 1e-2) | 0 | — | — | — |
| LOO | floor-zero | [1e-2, 1e-1) | 2 | 0 | 2 | 2 |
| LOO | floor-zero | [0.1, 1] | 1 | 0 | 1 | 1 |
| LOO | truth allele unavailable in calling panel | (0, 1e-6) | 4 | 0 | 4 | 4 |
| LOO | truth unavailable | [1e-6, 1e-3) | 3 | 0 | 3 | 2 |
| LOO | truth unavailable | [1e-2, 1e-1) | 3 | 0 | 3 | 3 |
| LZO | floor-zero | (0, 1e-6) | 18 | 18 | 0 | 2 |

No block printed an error mass of exactly 0.

**Accounting.** All 152 blocks are covered:
- LOO: 54 SV IDENTIFIED, 59 context IDENTIFIED, 1 context GENOTYPE_UNRESOLVED;
- LZO: 18 SV and 20 context, all IDENTIFIED.

**LOO floor-zero overall:** 35 of 44 reportable blocks exact. Errors are at both extremes of error mass:
- **Moderate-mass blocks:** all three are errors (d2 sv:9 at 0.022; d6 sv:9 at 0.013; d6 sv:1 at 0.254).
- **The lowest bin:** 6 errors in 34 blocks. **Five of those six are sv:5.**
- **Without sv:5,** the lowest bin is 28 exact out of 29; the only error is d3 sv:9 at 1.7e-7.

**Truth-unavailable blocks:** 10. Exact agreement is impossible by construction, and 7 of them carry error mass below
0.01 (for example d1 sv:13 at 1e-38, with floor 620). They show again that the mass is conditional on the pool, not
correctness.

**Legacy GQ comparison.** Where the assumptions coincide (lower = upper, single decision-set group: 63 of 72 blocks),
10^(−GQ/10) equals `error_mass_upper` within a log10 ratio of 5e-4. The 9 d5 blocks are not compared, because d5 has a
decision set of 2.

**Family-level** (floor-zero errors / of which high-confidence):

| family | floor-zero errors | high-confidence | also not exact, truth unavailable |
|---|---|---|---|
| d1 | 1 | 1 | 1 |
| d2 | 2 | 1 | 1 |
| d3 | 2 | 2 | 2 |
| d4 | 1 | 1 | 1 |
| d5 | 0 | 0 | 3 |
| d6 | 3 | 1 | 2 |
| LZO d1, LZO d2 | 0 | 0 | 0 |

Every LOO family has at least one not-exact reportable block. Five of six have a high-confidence floor-zero error.

## 4. High-confidence errors (floor-zero, `error_mass_upper` < 0.01; inspected after evaluation)

| cell | block | called / truth | error_mass_upper | category | evaluated hypotheses with truth pair | best carrier's whole-hypothesis gap (total / likelihood / prior) |
|---|---|---|---|---|---|---|
| d1 | sv:5 | 5/71 vs 6/77 | 1.1e-38 | present, losing as whole hypotheses | 36 | 589 / 638 / −49 |
| d2 | sv:5 | 76/89 vs 76/81 | 5.1e-28 | present, losing | 40 | 1,656 / 1,712 / −56 |
| d3 | sv:5 | 45/107 vs 34/94 | 1.5e-24 | present, losing | 5 | 5,280 / 5,329 / −49 |
| d3 | sv:9 | 2/2 vs 2/10 | 1.7e-7 | present, losing | 9,000 | 191 / 240 / −49 |
| d4 | sv:5 | 59/75 vs 59/78 | 1.3e-68 | present, losing | 8 | 6,413 / 6,479 / −66 |
| d6 | sv:5 | 0/21 vs 13/21 | 9.1e-7 | present, losing | 1 | 13,819 / 13,914 / −95 |

**Category counts:**
- **Representation ambiguity:** none.
- **Correct pair absent from the evaluated candidates:** none.
- **Correct pair present but losing as whole hypotheses:** all six.

**Not block-local.** These carriers differ from the winner at other blocks too: the likelihood gaps of hundreds to
thousands of nats are whole-hypothesis contrasts, not isolated block comparisons.

**sv:5 is a systematic finding:**
- not exact in **all six LOO cells**: five floor-zero, and d5 with truth unavailable at 3.5e-29;
- exact in **both LZO controls**;
- not exact in M8-D's HG01975 LOO too (D = 1.0, 10 bp);
- the legacy distances are 1–9 bp except d6 (43 bp), so the errors are small exact-sequence differences in 20–30-bp
  alleles of a bubble that is "SV" only because some panel path is ≥ 50 bp (M9-B §2).

**Not concluded:**
- that sv:5 has a block-local likelihood defect;
- that the error is structural;
- what causes it.

Answering these needs an isolated comparison: the winner with only its sv:5 alleles exchanged. That was not run.

## 5. Kept separate

- **Exact block-sequence correctness** is the only accuracy assessed. **Structural correctness is not assessed:** ACOT
  has no structural target defined.
- **Phase:** PRIOR_ONLY in 7 cells and UNRESOLVED in d5 (the sequence tie). There was no phase evaluation, because truth
  is not spellable under LOO.

## 6. Limitations

- **Sample size:** six families and 44 floor-zero LOO reportable blocks, which are not independent within a family. The
  middle bins hold 0–7 blocks. No calibration curve is warranted or fitted.
- **One locus, one coverage (5×), one route and one panel frame.** The selected assemblies were seen by M8-B's
  score-level evaluation.
- **Bins and stratum** are defined on conditional masses within each run's explored, non-exhaustive pool.
  "Truth unavailable" means unavailable in the calling panel, not biologically absent.
- **LZO controls are two families.**

## 7. What this supports

**The v2 contract behaved as documented on real data:**
- all cells predicted;
- a legacy no-call was carried as a v2 sequence-tied prediction;
- the projection-aware bounds agree with the legacy GQ where their assumptions coincide;
- the FASTA policy held.

**Descriptively,** the block error mass separates correct from incorrect reportable blocks well **except at sv:5**, where
it is confidently wrong in every LOO cell.

**Next (not run):** an isolated sv:5 comparison per affected cell, as the attribution step, before any calibration or
contract claim covers that block.
