# Automatic library estimation: diagnostic and measurements (2026-09-28)

**Status:** off by default and diagnostic only. The estimate is written beside a run; no call ever uses it.
- **Code:** `src/library_estimate.cpp`.
- **Entry points:** `panvar genotype --library-diagnostic`, which writes `<prefix>.library_estimate.tsv`, and
  `genotype_probe library-estimate`.
- **Enabling it:** `--library auto` (`genotype_cli.md` §4) is not enabled. This document is the evidence for that
  decision.

## 1. What a paired FASTQ reveals, and what it does not

**Read directly from the reads:**
- read-length distribution per mate;
- base-quality distribution;
- the fraction of ambiguous (non-ACGT) bases;
- the number of pairs.

**Not readable from the reads.** Headers carry no insert size. The insert distribution, the error rate and the
fragment-start intensity need placement against the model.

## 2. Method (truth-free; the model graph and the reads only)

**Anchors:**
- A canonical 31-mer is an anchor if it occurs exactly once, on one strand, in every panel template's spelled haplotype
  (canonical block-chain frame).
- A deterministic quarter of k-mers is kept, those whose hash is 0 mod 4, to bound the index.
- An anchor's position may differ between templates.

**Placement of a mate:** every anchor hit in a mate must imply the same orientation and the same start on every
template, or the mate is inconsistent.

**Decisive, concordant pairs.** A pair is accepted when all of these hold:
- one mate is forward and one reverse, facing inward (FR);
- the implied insert (outer distance) is positive;
- the insert is identical on every template.

The insert then depends on no candidate genotype, so it cannot be circular with the call.

**Rejected pairs are counted, not used:** unanchored, mate-inconsistent, orientation-discordant (FF/RR/outward) and
not decisive.

**Sampling:** the pairs with the smallest `fnv1a64(salt/name)` are taken, up to `max_pairs` (20,000), in hash order.
This does not depend on file order or thread count.

**Estimators:**

| quantity | estimator |
|---|---|
| insert centre | median, and the mean within 6 robust sd of it |
| insert spread | 1.4826 × MAD, and the sd within the same bound |
| outliers | the fraction beyond 6 robust sd |
| multimodality | best two-group (Otsu) split of the non-outliers. Flagged when the groups are ≥ 4 within-group sd apart and the smaller holds ≥ 5%. A single Gaussian gives about 2.7 sd |
| error | mismatch rate at invariant sites of accepted pairs: sites where every template carries the same base, first 2,000 pairs |
| fragment-start intensity, anchor route | 5%-trimmed mean anchor count ÷ 2 copies ÷ (2 (L − k + 1)(1 − e)^k) |
| fragment-start intensity, marker-fit route | the caller's own marker depth fit (invariant anchors, Poisson MLE with background), same conversion |

**Failure conditions:**
- **FALLBACK** on fewer than `min_accepted_pairs` accepted (default 500, from §3), or on multimodality.
- **WARN** when the outlier fraction is above 5%.
- A fallback value comes from the declared library. Without one the value is `NA` and the source is `none_available`:
  nothing is invented.

**Table columns:** `parameter`, `value`, `method`, `pairs_examined`, `pairs_accepted`, `spread`, `status`, `reason`,
`source` (`reads`, `model`, `placement_inference`, `explicit_or_preset_fallback`, `none_available`), and
`declared_value`.

## 3. Measurements (`tests/results/library_estimate_26401b4`, `measure.py`, no genotyping)

**Scope:** 18 read-set runs: 8 at 2×, 2 at 5×, 1 at 10× and 7 at 30×. Four thread runs and six sample-size runs are
added to these.
- LZO 2/5/10/30×, one simulation's nested subsamples.
- Identical reads under the LZO model (`c4.gfa`) and under the donor's LOO model: HG03456 and NA19434 at 2× and 30×.
- Synthetic mosaics (dense switching at 2× and 30×; panel control at 30×).
- Three fresh 2×/5× read sets.

**Comparison:** the frozen declared library, insert 350 ± 40, error 0.002, and fragment-start intensity as recorded.
The declared intensity is the realised fragment count ÷ 2 × locus length.

**Reads:**
- Read length is 100/100 in every set.
- Base quality is Q40 on every base: the simulator writes constant qualities.
  - The Phred-implied error is therefore 1e-4, a factor of 20 below the true 0.002.
  - Base quality cannot calibrate the error model on these data.

| quantity | 2× (8 runs) | 30× (7 runs) | declared |
|---|---|---|---|
| accepted pairs | 1,145–1,227 (46–56% of pairs) | 9,148–11,135 of 20,000 sampled | — |
| insert median | 344–352 | 348–349 | 350 |
| insert mean (within 6 sd) | 346.1–351.8 | 348.0–349.2 | 350 |
| insert sd | 38.8–41.8 | 39.6–40.2 | 40 |
| outliers; multimodal | 0; none | 0; none | — |
| invariant-site error | 0.00185–0.00203 | 0.00192–0.00199 | 0.002 |
| intensity, anchor route | −3.8% to −8.2% | −0.1% to −2.9% | — |
| intensity, marker-fit route | −0.02% to −11.5% | −0.1% to −4.1% | — |

**Other findings:**
- **Scenario independence.** The same reads give the same insert under LZO and LOO models: HG03456 30× 348.682 against
  348.678, NA19434 2× 347.318 against 347.378. The accepted counts differ by ≤ 1.2%.
- **Anchor route under LOO.** At 30× it runs about 1.3% lower under the LOO model (HG03456: 0.0739 against 0.0729). A
  panel anchor need not be present in a donor outside the panel.
- **Caller agreement.** The marker-fit route reproduces the caller's recorded fit exactly: `loo.HG03456.2` 0.6500
  against the campaign status's 0.650007.
- **Threads.** At 1/2/4/8 threads the table and every accepted insert are byte-identical (timing row excepted).
- **Sample size.** On `lzo.HG00096.30`, nested hash samples:

  | accepted pairs | 57 | 109 | 277 | 532 | 1,088 | 2,792 | 11,135 |
  |---|---|---|---|---|---|---|---|
  | insert centre (median / mean) | 332 / 338.1 | 339 / 344.2 | 343 / 345.9 | 345 / 347.2 | 345 / 346.2 | 346 / 347.2 | 348 / 348.0 |
  | insert sd | 40.2 | 39.8 | 40.7 | 39.6 | 40.5 | 40.9 | 40.0 |

  The sd is stable from about 100 pairs; the centre from about 500.
- **Residual bias.** The centre keeps about −2 bp (−0.6%) at every size. That is about 5% of the 40 bp spread, and
  immaterial to the model. The likely cause is selection: short fragments are more often decisive.
- **Cost.** The estimate takes 1.3–1.6 s of wall, about half of it the anchor index over 131 templates. A full probe
  run including read and model loading takes about 2.4 s. Peak footprint is 0.27–0.35 GiB, dominated by the model and
  reads.
- **No genotype effect.** `panvar genotype --library-diagnostic` on `lzo.HG00096.2`: every output is identical to the
  run without it (paths and timings normalised). The only addition is the table.
- **Failure conditions.** The fixtures (`tests/library_estimate.cpp`) check each one on a synthetic locus:
  - FALLBACK on 60 pairs;
  - FALLBACK multimodal on a 300/700 bp mixture;
  - outward pairs are never accepted;
  - N(350, 40) is recovered within 6 bp in centre and 5 bp in spread.
- **Not observed in C4 data:** contamination, multimodality and insufficient placement. The warning and fallback paths
  are exercised by fixtures only.

## 4. Recommendations (from these measurements)

1. **Infer read lengths.** They are measured exactly.
2. **Infer the insert centre and spread.** Use the median and the MAD sd, falling back on the within-bound mean and sd.
   - **Minimum accepted pairs: 500.** The centre is within about 3 bp of the large-sample value from there; below about
     300 it drifts 5–16 bp low. Every read set here has at least about 1,150.
   - **Sample size: 20,000 pairs.** That gives about 10,000 accepted at 30×, and every pair below about 20,000.
   - **Failure conditions:** fewer than 500 accepted; multimodal (≥ 4 sd, ≥ 5%); outlier fraction above 5% (warn).
     Also warn on an orientation-discordant fraction above 5% (at most 0.02% observed here).
   - **Fallback:** refuse and name `--library FILE`. No generic insert centre is substituted, because a wrong centre
     moves every fragment likelihood.
3. **Keep `substitution_error` a technology-preset parameter (0.002).** The invariant-site estimate matches it within
   10%, but:
   - base qualities cannot calibrate it here;
   - the model uses it as a fixed emission parameter across loci.
   - The estimate is reported as DIAGNOSTIC.
4. **Fragment-start intensity: inferable, not yet safe to enable.** It is a sample parameter: the declared value itself
   was the realised fragment count.
   - Both routes are low by up to about 8–11% at 2× and within about 1–4% at 30×.
   - It scales the exposure term, so before it replaces a declared value, a genotype sensitivity check is needed:
     intensity ×0.9 and ×1.1 on a few cells.
   - If that check shows insensitivity, use the marker-fit route: it is the caller's own fit and scenario-stable.
5. **Keep as preset or model parameters, never inferred:** `support_sigmas`, `outlier_weight`, the gap model, marker
   k/s, and the dosage policy.
6. **`--library auto` stays off.** It needs item 4's check. It also needs the insert refusal path tested on real data
   whose placements are insufficient or multimodal; none exists in the C4 set.

## 5. `--library auto` (implemented `c758ae2`; not the default)

**What `auto` uses:**
- **Insert centre and spread:** the estimator of §2, taking the mean and sd within 6 robust sd.
- **Validation:** at most 20,000 pairs are examined, and at least 500 must be accepted. The run is refused
  (`LIBRARY_REFUSED`, naming `--library FILE`) when the estimate is:
  - insufficient;
  - multimodal;
  - contaminated in the insert, with more than 5% beyond 6 robust sd;
  - contaminated in orientation, with more than 5% of anchored pairs not inward-facing FR.
- **Every other value from `--preset short-read`:** support 4, substitution 0.002, outlier 0.001, the built-in gap
  model, and fragment-start intensity 0.05. That intensity is a declared convention, λ = depth / (4 L) at 30× with
  2 × 150 bp pairs. It is not inferred.
- **Records:** `<prefix>.library_estimate.tsv` holds the estimate and a `library.<parameter>` row with the value and
  source of every used value (`placement_estimate`, `short_read_preset`, `reads`). The status carries `library_mode
  auto` and the provenance.
- **`--library FILE`** is unchanged.
- **Fixtures:** acceptance, each of the four refusals and its message, and the preset's scope are fixtures in
  `tests/library_estimate.cpp`.

**Sentinel gates** (`tests/results/library_auto_gates_c758ae2`):
- **Setup:** frozen C4 profile, dosage off. Each sentinel runs with `--library FILE` and with `--library auto`.
- **`--library FILE`** reproduces the sealed campaign byte for byte in 8 of 8.
- **Estimate determinism:** identical at 1/2/4/8 threads.
- **LZO against LOO, identical reads (HG03456 2×):** 346.59 ± 41.72 bp against 346.74 ± 41.77 bp. The anchor sets
  differ, so the two are close but not identical.

**`--library auto` changed scientific outputs in 6 of 8 sentinels.** Unchanged: `lzo.HG00096.30` and
`mosaic.dense_switching_1.30`.

| sentinel | changed |
|---|---|
| `loo.HG03456.2` | 4 block genotypes, trajectory, tie branches, locus filter |
| `lzo.HG03456.2` | block genotypes, trajectory |
| `mosaic.dense_switching_1.2` | block genotypes, trajectory |
| `loo.HG03456.30` | filters, trajectory |
| `lzo.HG00096.2` | one block's read support |
| `loo.NA19650.30` | the equivalence set's winner rows (same unordered genotype) |

The winner's exposure term moves by about −18,700 nats at 2× and +10,500 to +12,100 at 30×.

**Attribution:** arm A is the estimated insert with every declared value; arm B is the declared insert with the preset
intensity.
- **The preset intensity (B) reproduces every genotype, trajectory and tie-branch change,** and the 30× filter and
  winner changes. Against the declared intensity it is 10× too high at 2× and 1.5× too low at 30×.
- **The estimated insert (A) changes no genotype, trajectory or winner in any of the 6 cells tested.** Its only effect
  is the `lzo.HG00096.2` ctx:0 read-support flip: a fragment-contender rival crosses the error bound after a 3.5 bp
  shift in the insert centre.

**Decision: auto stays off by omission** (the batch's stop condition).
- A fixed preset intensity cannot represent samples of different depth.
- The intensity has to be declared, or inferred after the deferred sensitivity study.

## 6. `--library auto` with the marker-fit intensity (2026-09-29; `af4c2ce`; auto stays opt-in)

Frozen before any call: `tests/results/library_auto_lambda_study/SPEC.md`, seal `889f4fc6…`.

**What changed:**
- `--library auto` takes the fragment-start intensity from the invariant-marker depth fit.
- The anchor-depth estimate and the ratio of the two are recorded as diagnostics.
- The run refuses (`LIBRARY_REFUSED`, naming `--library FILE`) when the fit did not converge or the intensity is not
  finite and positive.
- The fixed preset intensity is no longer used by auto.

**Experiment** (8 sentinels, frozen reads, C4 profile, dosage off). Its outputs were sealed (`seal_outputs.sha256`,
431 files) before the frozen evaluator was run.

| arm | what it is | scientific outputs against A |
|---|---|---|
| A | the sealed explicit-library campaign output (not rerun) | — |
| B | declared library, marker-fit intensity | unchanged in 8 of 8 |
| C | full auto | primary call unchanged in 8 of 8; secondary differences in 3, below |

**The secondary differences in C:**
- `lzo.HG00096.2`: ctx:0 read support READS → PRIOR_ONLY. This is the insert-estimate knife-edge of §5.
- `mosaic.dense_switching_1.30` and `loo.NA19650.30`: the genotype-equivalent sets have the same members and the same
  unordered genotype. Only row order among tied members, or the phase-unresolved representative, changes.

**The intensity estimates:** marker-fit against declared ranges from 0.885 to 1.000, and marker-fit against anchor
from 0.92 to 1.07 (`comparison.tsv`).

**Other gates:**
- **Threads:** at 1/2/4/8 threads every output of the auto candidate is identical.
- **`--library FILE`:** byte-identical to the sealed campaign (2 cells rerun with the new binary).
- **Truth** (post-seal): no new catastrophic error, no LZO or mosaic call made inexact, and no LOO excess regression in
  any arm.

**Sensitivity** (auto insert fixed, λ at 0.9, 1.0 and 1.1 × the marker-fit estimate; 4 cells):
- Stable on 3 of 4: `lzo.HG00096.2`, `lzo.HG00096.30` and `loo.NA19650.30`.
- **Changed on `loo.HG03456.2` at 0.9×:** 5 block genotypes changed. ctx:0 4/9 → 4/6, sv:1 6/16 → 6/17, ctx:2 74/45 →
  74/44, ctx:6 90/110 → 90/100 and sv:7 6/91 → 6/100 are 1–2 bp allele swaps, and ctx:10 lost read support.
  - Against truth it is no worse: excess 4 → 2 bp, sv:1 truth undefined.
  - But it is a primary genotype change.

**Score decomposition (`hg03456_2_start_decomposition.tsv`).** The two best panel diplotypes compete for the repair's
start.
- The fragment term prefers the 0.9× start by 120.58 nats at every λ.
- The exposure term prefers the 1.0× start by 116.68, 129.64 and 142.61 nats at 0.9, 1.0 and 1.1 × λ, which is linear
  in λ.
- The prior is flat (0.0005 nats).
- The start therefore flips at about 0.93 × the estimate, and the repair then follows each start's lineage.
- §3 measured the marker-fit estimate as up to 11% low at 2×, so the flip lies inside the estimator's own uncertainty.

**Decision (the frozen rule):** a primary genotype changed inside ±10%, so `--library auto` stays opt-in and an omitted
`--library` does not mean auto. `--library FILE` remains the reproducibility path.

## 7. The library-intensity sensitivity audit (`library_sensitivity_audit_v1`; 2026-09-29; `f64733c`; auto stays opt-in)

The contract was frozen before implementation: `tests/results/library_sensitivity_audit/SPEC.md`, seal `a100d8bb…`.

**What it does (under `--library auto` only):**
- The reported call is also followed at 0.9× and 1.1× the marker-fit intensity, with every other library value
  identical. The central value alone decides the reported genotype.
- **One combination pass, three views.** Every diploid combination (union and prepared context) also scores the two
  endpoint libraries from the same mass columns. This is gated bit-identical to a run declaring each library, and the
  run's own scores are unchanged. Nothing is placed again because λ changed.
- **Shared trajectories.**
  - Each repair round's move rule is evaluated at every view from the already-scored neighbourhood.
  - While an endpoint decides as the call did, it shares the call's trajectory and scores nothing.
  - At its first different decision, the endpoint either stops (that round's items are already scored) or forks from
    its winner into the shared context. The same applies to its anchor and to every central tie branch.
  - Tie branching runs per view under the production rule.
- **Verdict.**
  - The common pool (the union plus everything any trajectory scored) is reranked at each view.
  - A block is library-sensitive when a view's top-set representatives do not carry exactly the reported unordered
    pair. Row order, phase, representatives and scores are not grounds.
  - Any incomplete view makes the verdict NOT_EVALUATED.
- **Filters.**
  - SENSITIVE puts `library_sensitive` into `quality_reason` (LOW_QUALITY) of the affected blocks only, and adds
    `LOWQUAL_LIBRARY_SENSITIVE` to the locus filter. The call itself never changes.
  - The record is in `<prefix>.library_sensitivity.tsv`, written on every run (NOT_APPLICABLE for `--library FILE`),
    and in the status under auto.
- **Recorded beside the verdict (diagnostic only):** each endpoint's independent call. That is its own pool (union set,
  trunk and branches) decided at its own view, which is what a run declaring that intensity reports.

**Acceptance on the 8 sentinels** (`gates.py`, frozen binary `f64733c`; outputs sealed at `7e408a2c…`, 346 files,
before truth):

| gate | result |
|---|---|
| central call identical to the `af4c2ce` auto run | 8 of 8, every deterministic file |
| audit on against audit off (same binary) | 8 of 8 identical but for the new record |
| `--library FILE` | byte-identical to the sealed campaign (2 cells); NOT_APPLICABLE; no new status line |
| threads 1/2/4/8 (`loo.HG03456.2`) | identical |
| stable cells not flagged for row order or phase | 8 of 8 STABLE; no `phase_or_representative_only` |
| **`loo.HG03456.2` detected SENSITIVE** | **fails: STABLE** |
| cells whose trajectories never fork: added placement and DP work | 0 on 6 of 6 (4 of them reuse or stop-fork central tie branches) |
| those six cells: wall, audit on over off (median) | 0.99–1.04 (audit section 0.15–0.32 s) |
| anchor-fork cells, reported separately | `loo.HG03456.2` +5.0 s (1.28×; 1,992 haplotypes placed, 0 new DPs); `loo.HG03456.30` +21.7 s (1.24×; 2,077 placed, 0 new DPs) |
| unit fixtures (`library_sensitivity_test`) | pass |

**Why `loo.HG03456.2` is STABLE:**
- At 0.9× the repair anchor flips (§6). The endpoint forks at its anchor and converges.
- Its independent call reproduces the study's 0.9× genotype exactly: ctx:0 4/6, sv:1 6/17, ctx:2 44/74, ctx:6 90/100
  and sv:7 6/100.
- But at 0.9× itself the reported call scores −33623.216 and the independent call −33624.095. The reported genotype is
  still the 0.9× likelihood's best candidate among everything scored, by 0.88 nats.
- The §6 change is therefore a search-path effect: a different anchor reaches a lower local optimum. It is not a
  genotype the 0.9× likelihood prefers. `loo.HG03456.30` forks at its anchor the same way, and its independent call
  agrees with the reported one.

**Truth (post-seal, `truth_eval.py`):**
- The frozen rule flags nothing. All 88 blocks are retained (75 exact); no block changed against `af4c2ce` auto; there
  is no new catastrophic error; the LOO excess totals are identical.
- The recorded independent-call diagnostic would flag the 5 `loo.HG03456.2` blocks. Of the 4 with defined truth, 3
  have an inexact reported call (ctx:0, ctx:6 and sv:7; excess 1, 1 and 2 bp) and ctx:2 is exact. sv:1's truth is
  undefined.

**Decision (the frozen rule):** the known sensitive cell was not identified, so an omitted `--library` does not mean
`--library auto`. The audit stays on under `--library auto`. Whether the verdict should be likelihood-rank sensitivity
(the frozen common-pool rule) or pipeline sensitivity (the recorded independent call) is open.
