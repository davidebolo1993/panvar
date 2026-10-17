# S1 implementation proposal: sequence-first reporting, schema v1 (for review; not started)

**What it implements.** Contract section 7 of `tests/genotype_sequence_label_reporting_design.md`, as a **versioned output-schema change**.

**What it leaves alone.** No score, prior, search, candidate-set, cap or default-behaviour change. The current schema stays the default until the new one is validated.

**Base.** Branch `s1-a-sequence-reporting`, from `genotype-cyp2d6-pilot` `85ae189`. Its `src/` and `include/` are identical to `genotype-mosaic` `cee58de`, plus the docs-only contract-v2 commits. No experimental performance backend, and no I4/I5 diagnostic machinery, is carried.

## 1. Switch

**New setting** (profile key and CLI flag): `output_schema` = `legacy` (default; today's outputs byte for byte) or `sequence_first_v1`. It is recorded in the settings table and the run status.

## 2. Affected functions

1. **`decide_pooled_genotype`** (`src/genotype_repair.cpp`) and `PooledDecision` (`include/panvar/genotype_repair.hpp`).
   - **Unchanged:** grouping, certification, the existing `TopSetDecision` and every current field (history and the legacy schema use them).
   - **Added, computed from the same groups:**
     - `group_block_genotypes[g][b]`: the set of unordered genotypes over group g's member labellings, i.e. the explored projections;
     - `sequence_status`, as in section 7.3;
     - per-block `block_status`: `IDENTIFIED`, `REPRESENTATION_AMBIGUOUS` (one decision group, several genotypes) or `GENOTYPE_UNRESOLVED` (distinct decision groups disagree);
     - the block's alternatives;
     - representation intervals (maximal runs of ambiguous blocks), with each label tuple of the representative group;
     - the canonical representative group: highest certified score, then canonical digest pair (already the group order, so `ranking.candidates.front()`);
     - per identified block, error-mass bounds over the pool, one weight per group: `lower` = Σ p(g) for groups whose every explored genotype differs, `upper` = Σ p(g) for groups with any differing genotype;
     - the alternatives: the decision groups plus the top 5 by score, compared with the representative block by block.
   - **Cost:** O(groups × members × blocks). No pairwise sequence distances.
2. **`run_genotype`, pooled path** (`src/genotype_command.cpp`), under `sequence_first_v1` only:
   - **The `DistinctTie` stop** (`stop_with(RescoreIncomplete, …)` after the equivalence set) no longer ends the run without outputs. `run_status`, `call_class` and `no_call_reason` keep their legacy values. The prediction outputs below are written with the representative group, and the process exit code is unchanged.
   - **The non-converged stop** (`!report.converged && !report.tie_branched`) with stop reason `round_cap` or a cumulative cap writes a `PROVISIONAL` prediction:
     - it comes from the last committed incumbent (exactly scored; already in the report as `final_incumbent_*`), with `sequence_status=NOT_DETERMINED`;
     - it has no block confidence;
     - it carries the stop reason in the completion fields.

     The repair trajectory stays as it is today.
   - **Unchanged:** watchdog stops (`stop_run`), refusals before a completed ranking, and hard failures. They write no prediction.
3. **Output writers on the pooled path**, under `sequence_first_v1`:
   - `.locus_genotype.tsv`: new columns `schema_version`, `prediction` (PREDICTION / NO_PREDICTION), `sequence_status`, `provisional`, `completion_reason`, `representative_hap1_sha256`/`_bp`, `representative_hap2_sha256`/`_bp`, `representation_intervals`, `pool_conditional_posterior` (advanced; labelled candidate-conditional), `structure_confidence` (`NA`), `structure_confidence_reason` (`UNCALIBRATED`), `projection_coverage` (`EXPLORED`). Existing columns are kept, with their legacy definitions.
   - `.block_genotypes.tsv` / `.genotypes.tsv`: new columns `block_status`, `block_alternatives`, `error_mass_lower`, `error_mass_upper`. `conditional_GQ` is kept and documented as `representative_projection_GQ`. Allele columns are `NA` on `REPRESENTATION_AMBIGUOUS` and `GENOTYPE_UNRESOLVED` blocks. Headline filters are unchanged, and computed only on identified blocks.
   - `.mosaic_haplotypes.fa`: written for every `PREDICTION`. The header carries `sequence_status`, `provisional` and the representation intervals. For `SEQUENCE_TIED` and `GENOTYPE_RESOLVED_SEQUENCE_SET` it says the record is the canonical representative of a set.
   - **New** `.representation_intervals.tsv` (interval, blocks, every label tuple per haplotype, the interval sequence digest per haplotype) and **new** `.sequence_alternatives.tsv` (the alternatives table).
   - `write_status`: the new locus fields, mirrored.
4. **Not changed in S1:**
   - the panel-only and union-only reporting paths (`rank_set` consumers), which are listed as legacy-only under v1;
   - dosage repair and dosage reporting on non-identified blocks;
   - phase reporting beyond the current `phase_support`;
   - the representation-conflict guard, which still runs only on identified blocks.

## 3. Acceptance tests

1. **Fixtures as oracle.** `tests/analysis/s1a/s1a_fixtures.py` gains a mode comparing the caller's `sequence_first_v1` outputs with `s1a_project.py`'s projection of the same run's legacy tables, on F1–F7, F7full and the rev/ren/sam variants. Required:
   - identical `prediction`, `sequence_status`, provisional flag, representative digests, block statuses and alternatives, intervals and error-mass bounds;
   - order independence across variants.
2. **Legacy identity.** With `output_schema legacy`, every e2e test and the review-regression expectations (`tests/genotype_review_regressions.expected`) are byte-identical.
3. **No prediction where it must not exist.** The watchdog test hooks (time limit, footprint) and an invalid-input case still write no prediction under v1.
4. **The adversarial contender and permutation cases** stay deterministic: the canonical representative is independent of candidate, template and alias order.
5. **Sealed CYP2D6 pilot.** Explained only, never rerun: its reporting projection is updated to section 7's classes. Prediction yield, high-confidence yield (not available: uncalibrated) and accuracy are reported separately (C5). No accuracy gain is claimed from no-calls that become predictions.

## 4. Out of scope (later, each separately reviewed)

- calibration and any PASS threshold;
- τ and structural classes;
- MBR;
- projection weighting;
- complete projection sets from model-level representation intervals;
- phase and dosage handling of intervals;
- the panel-only path;
- any performance-backend integration.
