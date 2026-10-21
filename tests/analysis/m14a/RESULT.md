# M14-A result: NOT USEFUL — whole-haplotype representative panels give at most 1.22× pooled call speedup, and they introduce and worsen structural catastrophes

**Run.** Branch `m14-a-representative-panel`; freeze `f714375`; result after.
- **Results:** `tests/results/m14a_representative_panel_20261021/`.
- **Budget:** 1,701.9 of 2,400 s: 60 s of disclosed development plus 1,641.9 s of ledger wall.
- **Units:** 108 of 108 passed first time, with no amendment and no watchdog or envelope event.
- **Peaks:** 2.52 GiB per call; 1.64 GiB per evaluation.

**Not done:**
- no default, integration, FASTA, merge, push or cleanup;
- the 18 reserved ACOT families are untouched.

**This is a speed–accuracy measurement,** not an equivalence or calibration result.

## 1. Mechanism and clustering (frozen; panel sequences only)

**Mechanism:** whole panel haplotypes are removed from the active panel through the existing `--panel-curation`.
- **Unchanged:** the graph frame and blocks.
- **Uncallable:** alleles carried only by removed haplotypes.
- **Copying-prior donors:** the retained templates, each counted once, with no cluster-size reweighting.
- **Representatives:** actual panel templates. All member-to-representative mappings are in `runs/cluster/membership.tsv`.

**Clustering:** a deterministic farthest-point selection on Mash distance (21-mers, bottom-2,000), seeded at the locus
reference template. Each member is assigned directly to its nearest representative.

| locus | active templates | r2 representatives | r4 representatives | member→representative distance, median / p90 / max (r2) | (r4) |
|---|---|---|---|---|---|
| CYP2D6 | 124 | 62 | 31 | 4.8e-5 / 7.2–8.4e-5 / 9.6e-5 | 1.6–1.9e-4 / 2.5–2.8e-4 / 3.4–3.9e-4 |
| C4 | 129 | 65 | 33 | 4.8–6.0e-5 / 0.8–1.0e-4 / 1.1–1.2e-4 | 1.3–1.4e-4 / 2.3–2.6e-4 / 2.9–3.1e-4 |
| ACOT | 407 | 204 | 102 | 6.0e-5 / 1.3e-4 / 1.7e-4 | 2.8e-4 / 3.5e-4 / 4.3e-4 |

The distances are Mash per-base divergence estimates. The largest clusters hold 6–7 members (r2) and 23–37 (r4).

**Preparation cost** (from the ledger):
- **Sketches:** 15 units, 35.7 s in total — about 2–3 s per C4/CYP2D6 model and 5 s for ACOT. One sketch serves every
  sample of a model, so this amortises.
- **Clustering:** 39.9 s. It is a pure-Python all-pairs Mash computation, dominated by ACOT's 407 templates × 3 cells; it
  could be computed once per model.
- **Sketch files:** 96 MB as text in total.
- **In total,** preparation is about 4% of the 1,035 s of call wall time.

## 2. Speed and work (call wall time, 15 cells, matched full-panel runs in this batch)

| locus | full (s) | r2 (s) | speedup | r4 (s) | speedup |
|---|---|---|---|---|---|
| C4 (6 cells) | 73.0 | 63.4 | 1.15× | 74.6 | **0.98×** |
| CYP2D6 (6) | 140.4 | 131.2 | 1.07× | 145.7 | **0.96×** |
| ACOT (3) | 175.8 | 125.3 | 1.40× | 106.0 | 1.66× |
| **pooled** | 389.2 | 319.9 | **1.22×** | 326.3 | **1.19×** |

**Executed work fell as expected,** summed over cells (full → r2 → r4):
- **panel diplotypes:** ÷4 and ÷16;
- **unique haplotypes:** C4 1,350 → 1,057 → 883; CYP2D6 1,313 → 909 → 733; ACOT 1,556 → 975 → 675;
- **physical origins:** −30 to −42% (r2) and −45 to −60% (r4);
- **marker time:** ÷2–6; **rescoring time:** −20 to −54%;
- **ACOT peak footprint:** about 2.6 → 1.5 → 1.3 GiB per run.

**Wall time barely moved:** placement and rescoring are a minority of call time. The adaptive search, repair and dosage stages
do not shrink with the panel, and in several reduced runs they took longer:
- CYP2D6 NA19983 r4: 38.3 s against 15.5 s;
- C4 HG00513 r4: 20.3 s against 12.5 s;
- CYP2D6 NA19650 r2: 32.7 s against 18.4 s.

## 3. Accuracy (per block; primary = excess against the ORIGINAL full-panel floor)

| | full | r2 | r4 |
|---|---|---|---|
| summed excess vs the full floor: all 236 blocks | 53,678 | **78,887** (+25,209) | **90,567** (+36,889) |
| C4 | 125 | 264 | 6,736 |
| CYP2D6 | 26,226 | 50,790 | 55,438 |
| ACOT | 27,327 | 27,833 | 28,393 |
| catastrophe introduced (< 1,000 → ≥ 1,000) | — | 1 | **4** |
| catastrophe worsened | — | 3 | 3 |
| catastrophe reduced / removed | — | 0 / 0 | 0 / 0 |
| control damaged (≤ 10 → > 10) | — | **17** | **28** |
| other blocks worse / better / unchanged | — | 52 / 24 / 139 | 63 / 21 / 117 |

**Catastrophic transitions** (excess vs the full floor; in brackets, the full floor → the reduced floor):

| setting | block | excess (bp): full → reduced | floors | transition |
|---|---|---|---|---|
| r2 | CYP2D6 HG02953 sv:9 | 0 → 12,176 | 7 → 198 | introduced |
| r2 | CYP2D6 NA19983 sv:9 | 12,075 → 24,251 | 1,625 → 1,625 | worsened |
| r2 | CYP2D6 NA19650 sv:9 | 12,159 → 12,162 | | worsened |
| r2 | ACOT NA19909 sv:13 | 27,270 → 27,276 | | worsened |
| r4 | CYP2D6 HG02953 sv:9 | 0 → 1,820 | 7 → 1,602 | introduced |
| r4 | CYP2D6 HG03456 sv:9 | 2 → 13,681 | 29 → 57 | introduced |
| r4 | CYP2D6 NA19705 sv:9 | 27 → 12,192 | 41 → 1,764 | introduced |
| r4 | C4 HG00513 sv:7 | 9 → 6,427 | 8 → 23 | introduced |
| r4 | CYP2D6 NA19983 sv:9 | 12,075 → 13,910 | | worsened |
| r4 | CYP2D6 NA19650 sv:9 | 12,159 → 12,216 | | worsened |
| r4 | ACOT NA19909 sv:13 | 27,270 → 27,282 | | worsened |

**The reduced-panel floor rises:**
- r2: 61 of 236 blocks rise, by 714 bp in total (max 191);
- r4: 99 of 236 rise, by 7,376 bp in total (max 1,748).

**Most new error is not the lost availability.** In several introduced catastrophes the reduced floor rises by far less than
the error: HG02953 r2's floor went 7 → 198 while its excess went to 12,176, and HG03456 r4's went 29 → 57 while its excess
went to 13,681. With fewer haplotypes, the panel's copy-number and allocation choices go wrong. Subtracting only the reduced
floor would hide this.

**Completion and ambiguity:**
- no no-call and no incomplete search in any setting;
- **CONDITIONAL_GENOTYPE_EQUIVALENT** runs: 1 at full panel (CYP2D6 NA19705); 3 at r2 (CYP2D6 NA19650 and NA19983, C4 NA12329), with NA19705 becoming resolved; 2 at r4 (NA19705, NA19983).

**Locus-level whole-region distance** (complementary diagnostic, `runs/analyse/runs.tsv`):
- **It is dominated by region-end differences** between the called and simulated regions (up to about 220 kb in some CYP2D6
  cells), so only changes against the full panel are meaningful, and even those can mislead. For example, CYP2D6 HG02953 r2:
  the locus distance fell by 2.3 kb while that block's excess rose by 12.2 kb.
- **Not computed for 5 of 45 runs:** both haplotypes were panel paths, for which the caller writes no FASTA by design.
- **It does not change any conclusion here.**

## 4. Decision (frozen): **NOT USEFUL**

- **Speed:** pooled call speedup is below 1.25× at both settings (1.22× at r2, 1.19× at r4). C4 and CYP2D6 are no faster at
  r4.
- **Structure:** there are introduced catastrophes (1 at r2, 4 at r4), and catastrophes worsened or introduced (4 and 7)
  outnumber those reduced or removed (0 and 0).
- **Controls:** 17 (r2) and 28 (r4) near-floor controls were damaged.

**Only ACOT accelerates materially** (1.40× and 1.66×, with a lower memory footprint), at an accuracy cost of +506 / +1,066
bp of summed excess on 3 cells. That is too small a sample to support an ACOT-only route, and nothing here was tuned toward it.

**What the result shows:** representative selection reduces the work it can reach (placement, origins, markers), but that
work is not where end-to-end time goes. Removing haplotypes also removes the copy-number and allocation context the search
relies on.

These are inspected development findings on 15 cells. They are not validation and they are not a confidence claim.

## Closure (2026-10-21): NOT USEFUL accepted for the tested route; qualifications appended

**Accepted:** NOT USEFUL for the tested **whole-haplotype representative route with Mash similarity.** It does not enter
the release, and no further clustering, profiling campaign or ACOT-only tuning follows.

**Qualifications:**
1. **No conclusion is drawn about structure-aware clustering methods in general.** The broader idea stays a future research
   option.
2. **Three changes were not isolated:** removing templates changed allele availability, copying-prior composition and
   search behaviour together. The experiment establishes harmful outcomes, but not which change caused each error.
3. **Damaged-control counts refer to near-floor blocks, not independent samples.** Block and cell/donor counts:

   | setting | transition | blocks | cells | donors | cells by locus |
   |---|---|---|---|---|---|
   | r2 | control damaged | 17 | 11 | 10 | ACOT 3, C4 3, CYP2D6 5 |
   | r4 | control damaged | 28 | 11 | 10 | ACOT 2, C4 4, CYP2D6 5 |
   | r2 | catastrophe introduced | 1 | 1 | 1 | |
   | r4 | catastrophe introduced | 4 | 4 | 4 | |
   | r2 / r4 | catastrophe worsened | 3 / 3 | 3 / 3 | 3 / 3 | |

   NA19705 contributes a C4 and a CYP2D6 cell, so there are 10 distinct donors in 11 cells.
4. **Small Mash distances do not certify similar repeat counts or arrangements,** so tight clusters were not necessarily
   safe structural substitutes.
