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
