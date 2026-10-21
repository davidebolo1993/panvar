# M14-A: bounded representative-panel experiment (frozen before execution)

**Question:** can a smaller representative panel substantially reduce computation while preserving useful structural
predictions?

**What this is:** a deliberately approximate, default-off speed–accuracy experiment.

**What this is not:** an equivalence test or confidence calibration. Exact haplotype identity is not required, and
additional structural error is measured, not hidden.

**M13-A** is closed (INCONCLUSIVE), with its five qualifications appended. Its diagnostic code is preserved.
- **Branch:** `m14-a-representative-panel`, from `m13-a-variation-tolerant` `5cb91d5`.
- **Results:** `tests/results/m14a_representative_panel_20261021/`.
- **Budget:** 2,400 s including 60 s of disclosed development:
  - the probe build (16 s);
  - three development sketches (9 s);
  - one distance timing (1.5 s).
- **Resources:** 4 threads; one process at a time; 6 GiB.

**Not done:**
- no default, integration, FASTA, merge, push or cleanup;
- the 18 reserved ACOT families are untouched.

## 1. Mechanism (one only)

**Representatives replace whole panel haplotypes.** They are selected through the caller's existing `--panel-curation`:
every non-representative template is added to the cell's base curation.

**Graph frame:** unchanged. Curation never changes blocks (docs/algorithms/genotype.md): the removed paths still define
block structure.

**What is removed:**
- **Candidate generation:**
  - mosaic generation and the panel baseline run over the retained templates only (panel diplotypes: N(N+1)/2);
  - alleles carried only by removed templates are uncallable, so no stage proposes, repairs or reports them.
- **Placement:** fewer distinct block alleles means fewer unique haplotypes to place reads on.
- **Preparation:** model and marker building, with no persisted index in this configuration.

**What changes:**
- **Copying prior:** the donors are the retained templates, **each counted once**. There is no cluster-size weighting:
  donor multiplicities are not reweighted, and that consequence is declared.
- **Unchanged:** every scientific and resource setting (the M10-A unified D configuration and M11-D's per-locus watchdogs:
  C4 / CYP2D6 150 s, ACOT 300 s; time guard = watchdog − 20).

**Representatives are actual panel sequences.** Each is a template spelled through the model, and every member's mapping is
recorded in `runs/cluster/membership.tsv`.

## 2. Clustering (frozen; panel sequences only)

**Sketches.** `genotype_probe panel-sketch` spells every active template through its own allele path. It keeps a bottom-2,000
MinHash of canonical 21-mers (splitmix64 hash), at about 2–5 s per model; there is no all-pairs alignment matrix.

**Distance:** the Mash distance −ln(2J/(1+J))/21 from the union bottom-2,000.

**Selection:** deterministic farthest-point (k-centre).
- **Seed:** the locus reference template, active in every cell.
- **Next pick:** the template farthest from the current representatives (ties go to the lexicographically smallest name).
- **Stop:** at K = ceil(N/2) (setting **r2**) and K = ceil(N/4) (setting **r4**).
- **Assignment:** each member goes to its nearest representative (ties by name).
- **No chaining:** every member-to-representative distance is computed directly and reported, with distribution and
  maximum.

**Status of the settings:** experimental operating points, not biological-equivalence thresholds.

**Not used:** held-out assemblies, truth distances, or which candidate would fix a failure.

## 3. Cells (frozen manifest: the M11-D development prefix per locus)

| locus | cells |
|---|---|
| CYP2D6 | all 6 |
| C4 | the first 6 (NA19705, NA19317, HG03065, HG02018, HG00513, NA12329) |
| ACOT | the first 3 (HG03209, HG01975, NA19909) |

**Every cell is run at full, r2 and r4 in this batch,** so each reduced run has a matched full-panel control. Each call is
followed by its evaluation.

## 4. Evaluation (`analyse`; truth opened only here)

**Two floors, kept distinct:**
- **The full-panel floor** (from the matched full-panel evaluation): **excess_vs_full_floor = distance − full floor.** This
  is the primary measure, of degradation against the existing task. It is never offset by the reduced floor.
- **The reduced-panel floor** (the evaluator run on the reduced active model): a diagnostic of approximation cost (the floor
  increase).

**Per block:** distance, both floors, both excesses, normalized excess (against the full floor) and the change against the
matched full-panel run.

**Transitions per block** against the full-panel run:
- catastrophe introduced / worsened / reduced / removed / unchanged (≥ 1,000 bp);
- control damaged (≤ 10 → > 10);
- otherwise worse / better / unchanged.

**Per run:**
- status, completion, ambiguity and no-call;
- call wall time and peak memory;
- the stage seconds (model, markers, proposal, rescoring);
- panel diplotypes, mosaic and sequence candidates, unique haplotypes, physical origins and accepted mate placements;
- the **locus-level whole-region distance** — called haplotypes against the simulated truth haplotypes, the minimum over
  pairings, by exact edit distance. It is a complementary diagnostic and **not** the M10-D chain-interior contract; the
  region ends differ slightly, so it is reported as a change against the full panel.

**Preparation:** sketch and cluster costs are reported apart from the calls. **Amortization:** sketches depend only on the
model, so one sketch serves every sample of that model.

## 5. Decision (frozen)

- **PROMISING TRADE-OFF:** at r2 or r4, a call-wall speedup of at least 1.5× (pooled) with a clearly described accuracy
  cost, and no structural collapse.
- **NOT USEFUL:** pooled speedup below 1.25× at both settings, or structural degradation that would damage the task: any
  introduced catastrophe **and** more catastrophes worsened or introduced than reduced or removed.
- **INCONCLUSIVE:** preparation cost, incomplete comparisons or resource bounds prevent a conclusion.

**No outcome-driven threshold change** and no wider campaign.

## 6. Units (`m14a_run.py`; dependency-gated; `next` stops on any failure)

| unit | needs | declared timeout (s) |
|---|---|---|
| `build` | | 300 |
| `sketch.<cell>` ×15 | build | 60 |
| `cluster` | all sketches | 120 |
| `call.<cell>.<setting>` ×45 | build, cluster | per-locus watchdog |
| `eval.<cell>.<setting>` ×45 | the call | 120 |
| `analyse` | cluster | 600 |

**Admission:** call and eval units are admitted only if spent + timeout + 600 (the analyse reserve) ≤ 2,400.

**Expected time:** about 1,300–1,600 s.
