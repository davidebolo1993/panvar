# M8-A: copying-prior calibration feasibility and study contract (frozen before any unit)

**Ancestry.** Branch `m8-a-copying-calibration` from `5294d29`. M7 is closed: M7-C is accepted, and there are no further
Locityper experiments.

**Scope.**
- Preparation, plus a bounded panel-only feasibility check: block structure and sequence support only.
- **Excluded:** fitting, genotyping, caller changes, default changes, clustering, residual editing, new dosage models,
  merge, push and cleanup.
- Budget: 300 s.

**Reused, not reimplemented:**
- `build_copying_transitions` / `CopyingTransitionSpec` (`constant_v1`, `distance_ref_core_v1_experimental`);
- `LocusGeometry` and the M4-A geometry outputs (block ids);
- `log_haploid_sequence_prior` / `log_haploid_mosaic_prior` (`ls_hmm.hpp`);
- the M4-A gates G1–G10 and the M4-B ctests;
- `genotype_probe block-sequences`.

## Units

| unit | does |
|---|---|
| `audit_gates` | M7-C maintenance: explicit numerical gates (mean ≤ 1e-6, lik_sd ≤ 1e-4, weighted ≤ 1e-6, attempts 1–20) and six perturbation controls, recomputed from the archived M7-B.1 outputs |
| `support.full` | `block-sequences` on the full models (C4 `c4.gfa`; CYP2D6 `cyp2d6.gfa` + curation; ACOT `acot.gfa`), all blocks from the M4-A geometry. Inventory: templates, donors, IGSR-pedigree families, status (UNINSPECTED / INSPECTED_DEVELOPMENT per `DEVELOPMENT_DONORS.tsv`, expanded to families / REFERENCE_GENOME). Per (family, haplotype, block): label support and sequence support from templates outside the family. Also: identical complete training paths, duplicate allele sequences and incomplete templates |
| `support.loo` | frame-leakage check on at most 8 UNINSPECTED donors per locus that already have a LOO model: their true haplotypes projected onto the LOO block chain, with support tested against training alleles (relatives excluded) |

**Development donors (`DEVELOPMENT_DONORS.tsv`, 63).** Donors named in memory notes or in `tests/analysis`
designs/results as of M8-A, campaign mosaic-plan donors, and the donors of M6-F's evaluation families (whose panel
sequences were scored).

**Not used for fitting:** the known CYP failure, which is a later diagnostic only.
