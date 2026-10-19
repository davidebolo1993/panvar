# M4-A: integration surface and a minimal default-off proposal (not implemented)

## 1. Inventory: every active construction or consumption of transition probabilities

Today every site builds its own `std::vector<double>(B−1, options.switch_probability)`, or receives the scalar. The
lines below are at `f31b9e6` lineage sources.

### Production construction and consumption

| stage | site | use |
|---|---|---|
| Settings | `genotype_command.cpp:2274/2327/2359/2452/7888`, `genotype_command.hpp:50` | parse, range (0, 0.5], default 0.01, CLI `--switch-probability`; status `:1365` |
| Cache keys | `genotype_command.cpp:223,249` | repair/tie settings digest includes `switch=<hex>` |
| Early incumbent (S2-A, experimental) | `genotype_command.cpp:3370–3376, 3457–3458` | path and sequence priors |
| **Proposals** (mosaic beam) | `genotype_command.cpp:3500–3523` → `mosaic_beam.cpp:315–512`, `mosaic_beam_grouped.cpp:318–544` | the beam's transition at every block |
| **Ranking** (union pool) | `genotype_command.cpp:3625` | `log_haploid_sequence_prior` per interned haplotype, giving the sequence-pair prior in the total score |
| **Path priors** (candidate tables, consensus and phase audit) | `genotype_command.cpp:3788–3793, 3993–4007, 4293` | `log_haploid_mosaic_prior` (labelled path; reporting and audit) |
| **Refinement** | `genotype_repair.cpp:258` (`switches_`), `:709` (neighbour prior), `:957` (pooled combination) | sequence-pair prior |
| Repair seeds and tie branches | `genotype_command.cpp:5332, 5463`; start check `:6650–6697` | sequence-pair prior |
| Joint search (parked) | `genotype_command.cpp:5122` → `genotype_joint_search.cpp:102, 262` | sequence prior |
| **Projection / representation guard** | `genotype_command.cpp:6528, 6565` | sequence-pair prior of guard competitors |
| **Dosage** (audit, shadow, repair, integrate) | `genotype_command.cpp:8370–8382, 8412` → `genotype_dosage_shadow.cpp:94–331`, `genotype_dosage_repair.cpp:295/403/1484`, `genotype_dosage_integrate.cpp:256/282` | sequence-pair prior |

`genotype_dosage_integrate.cpp:221/844` (`template_switches`, `switch_basis panel_templates_min_switches_left_greedy`)
is a switch **count**, not a probability. It is unaffected, but its name must not be confused with the transition.

### Diagnostics and supported audits

| tool | use |
|---|---|
| `genotype_probe score-set` (`--switch-probability`) | its own constant vector |
| `prior-pool` | the M3-A fixed model set |
| `two-block-*`, `tie-branch`, `repair-search`, `dosage-replay` | each builds its own constant vector |
| tests (`m3a_prior_representation`, `m4a_distance_prior`, `genotype_mosaic_core`) | their own vectors |

### Persisted and reusable data

- The prepared-locus file (`prepared_locus_file.cpp`, I3-C/D v1) stores model, panel, markers and digests. It holds
  **no prior and no geometry**.
- The model fingerprint and model-content digest cover blocks, alleles, invariant segments and panel rows. They do
  **not** cover the reference path, so a distance geometry would not be identified by any current digest.

## 2. One shared construction and provenance contract

```
struct CopyingTransitionSpec {            // from settings, validated once
    enum class Mode { constant_v1, distance_ref_core_v1 } mode = Mode::constant_v1;
    double constant_r = 0.01;              // constant_v1 only
    double copying_length_bp = 0;          // distance mode: lambda (declared) ...
    bool hazard_matched_control = false;   // ... or the diagnostic control lambda_HM (frozen per model fingerprint)
    std::string redraw_measure = "uniform_template_rows_v1";
};
struct LocusGeometry {                     // NEW adapter output (build_locus_input), never derived from reads/truth/candidates
    std::string reference_path_name, reference_path_sha256;
    std::vector<std::pair<double, double>> core_interval;   // per block, reference axis, cores trimmed against the FRAME
    std::string geometry_digest;           // over the above + the frame's block ids
};
struct CopyingTransitions {
    std::vector<double> r;                 // one per boundary
    double lambda_bp = 0;                  // the value actually used
    std::string spec_digest, geometry_digest, vector_digest;
};
CopyingTransitions build_copying_transitions(const PreparedGenotypeModel&, const LocusGeometry*, const CopyingTransitionSpec&);
```

**Rules.**
1. **Built once per locus**, before proposals. Every consumer in section 1 takes `const CopyingTransitions&`; none reads
   `options.switch_probability` directly.
2. **No silent mixing.** Each stage records the `vector_digest` it used. The status writes one digest, and the run
   refuses if any stage reports another.
3. **Provenance.** Cache keys (repair/tie settings digest, union caches) replace `switch=<hex>` with `vector_digest`.
   Settings and status record mode, λ (and whether it is the matched control, with the model fingerprint it was frozen
   on), U, `geometry_digest` and `vector_digest`.
4. **Fail closed.** The distance mode refuses when:
   - the geometry is unavailable: the reference does not tile, or is Uncovered or inconsistent;
   - the run uses a v1 prepared-locus index, which has no geometry;
   - `switch_probability` is also set explicitly (ambiguous);
   - representation-dependent blocks are present (already refused).

   There is no fallback to constant r.
5. **Index compatibility.** A future index version would persist `LocusGeometry` with its digest. Until then, distance
   mode prepares from the GFA. A v1 index loaded under distance mode is `incompatible_preparation: geometry absent`.
6. **Leave-one-out.** The model built from a LOO graph uses that graph's own frame geometry. The digest change is
   reported, and comparisons across scenarios state it.

## 3. Minimal default-off integration (for a later batch; not done here)

| step | change |
|---|---|
| 1 | Adapter: compute `LocusGeometry` from the reference path, with cores against the frame including curation-inactive alleles. Gate: tiling plus agreement with the reference template, as the probe already checks. |
| 2 | Library: `build_copying_transitions`. Gates: port G1–G10 from the test prototype unchanged; the oracle stays the brute-force history sum. |
| 3 | Caller: thread `CopyingTransitions` through the section-1 consumers with **no behaviour change in `constant_v1`**. Gate: every output byte-identical to the current caller on the M3-F/M3-D cells under the existing volatile set. |
| 4 | Profile key `copying_transition` = `constant_v1` (default) \| `distance_ref_core_v1_experimental`, plus `copying_length_bp` or `copying_length_control hazard_matched`. Settings and status rows as above. |
| 5 | Probes: `score-set` and `prior-pool` accept the same spec, so diagnostics use the identical vector. |

## 4. Proposed bounded fresh-search validation set (next batch; runs nothing now)

**Purpose.** Measure what the **geometry** changes at equal total hazard, with fresh search, not fixed pools.

**Arms.** `constant_v1` against `distance_ref_core_v1` with the **hazard-matched control λ_HM**, frozen per model
fingerprint. No physical λ is used until the independent study (CONTRACT section 3b) has chosen one.

| class | cells (frozen campaign inputs) | why |
|---|---|---|
| Ordinary LZO controls | C4 HG00171 2× and 5× (the first sample by ID other than HG00096 with exact calls at both; verified in the campaign); the next such sample by ID, HG00268, at 2× and 5× | ordinary panel behaviour, low and moderate depth |
| Mosaic panel controls | C4 `panel_control_1` 2× and 5× | the paths are real panel rows: the prior should not invent switches |
| Mosaic recombinants | C4 `both_homologues_1` 5×; `one_homologue_once_1` 2× (the M3-D case); `three_switches_1` 2× and 5× (the ctx:0 case, included but **not** a selection criterion); `dense_switching_1` 10× | switch-rich truths, where the hazard placement matters most |
| Structural control (expected unchanged) | C4 LOO NA19650 10× | sv:7 is content (M1-B/C); a change here would be a warning, not a success |
| Cross-locus ordinary | CYP2D6 LZO 5× (pilot profile); ACOT one ordinary 5× cell | different block counts and λ_HM (1.13 and 1.22 Mb vs C4 1.87 Mb) |

**Report.** Matched off/on outcomes, harmful changes as prominently as recoveries, refusals in every denominator,
prior-vs-likelihood decomposition of any changed call, and costs.

**Not in scope there.** Choosing λ, promotion, or campaigns.

**Size.** About 14 cells × 2 arms, an estimated 700–1,000 s at 4 threads, so 1,200 s is a reasonable bound.
