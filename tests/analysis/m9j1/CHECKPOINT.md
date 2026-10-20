# Development checkpoint (M9-J.1, 2026-10-20)

**What this is.** A frozen development baseline: the source, the experimental settings that were actually validated, and
the artifacts needed to reproduce them. It is **not a release**. Nothing is merged, pushed or promoted to a default. Every
switch listed under "experimental" stays opt-in.

## Source

| item | value |
|---|---|
| repository worktree | `/Users/davide.bolognini/panvar-worktrees/s1-d-panel-references` |
| branch | `m9-j1-resolution-fix` (from `m9-j-compound-integration` `3481bca`) |
| library and caller source | `e61c09c` (M9-J.1 resolution fix). `59c9919` changes the runner only |
| build | `build/` (CMake Release). Targets: `panvar`, `genotype_probe`, `compound_annotation_test`, `junction_proposals_test`, `genotype_sequence_report_test`, `genotype_command_e2e_test`, `cyp2d6_block_alleles` |
| submodules | `external/edlib`, `external/abpoa`, `external/minigraph`, at the worktree's recorded commits |
| preserved, untouched | primary worktree WIP `tests/genotype_probe.cpp` (sha256 `311a41adc1eddc0c…`), stash `dbb170c5`, every earlier package under `tests/results/` |

## Effective experimental settings (as validated; read from each run's `settings.tsv`)

**Common to every validated configuration:**
- `output_schema sequence_first_v2`;
- `switch_probability 0.01`, `kmer_size 31`, `repair_rounds 32`, `repair_tie_branching 1` (short-read preset values);
- dosage audit only (`dosage_audit shadow`, `dosage_policy_v2`): it never changes the prediction;
- `repair_joint_search off`;
- the performance backends of `profiles/genotype_dev_performance_v1.tsv`, all exact. They are listed in the profiles below.

**Scientific switches that differ between the validated configurations:**

| configuration | used for | JP1.1 (`repair_junction_proposals`) | M1-B overlay (`early_incumbent`, `union_fragment_blocks`, `repair_round1_origin_admission`) | `compound_annotation` | other |
|---|---|---|---|---|---|
| ACOT (M9-H on-arm) | `m9h…/runs/prepare/profile.on.tsv` + `--dosage-mode audit` | `experimental_v1` | not set | `contract_v2` | `copying_transition constant_v1`, `repair_max_neighbours 1004`, 900 s limit |
| CYP2D6 (M1-B) | `--preset short-read`, `m1b…/call.cyp_30/profile.{cyp_30,m1b}.tsv` | off | set | `contract_v2` | `repair_max_neighbours 958`, 280 s limit |
| C4 (M1-B) | `--preset short-read`, `m1b…/call.c4_10/profile.{c4_10,m1b}.tsv` | off | set | `contract_v2` | 100 s limit |

**Not validated:** a single combined profile (JP1.1 + M1-B overlay + compound annotation) on any locus. The checkpoint is
the source plus these three configurations, not one unified profile.

## Declaration artifacts (`compound_interval_declarations_v1`)

All are under `tests/results/m9j_compound_integration_20261020/runs/declare/`.

| model | declarations sha256 | profile sha256 | content |
|---|---|---|---|
| ACOT `12cbe6be…` | `ac0610b3d472a2d8…` | `21188529d141a5f3…` | sv:5..ctx:6 |
| CYP2D6 `e41f77fd…` | `bf7efd8eb1eed191…` | `aeda5f67acef2f31…` | sv:3..sv:5 |
| C4 `14851130…` | `4155f42cb96842b0…` | `8c1a580a5b9d456f…` | none (a control) |

**Model binding.** Each declaration is bound to a model fingerprint and an ordered list of block ids. A run on any other
model reports `INVALID_FINGERPRINT` or `INVALID_BLOCK_IDS` and ignores the declarations. A new model needs a new artifact
(`m9j.py declare`, from that model's inventory).

## Required inputs and indexes

| locus | graph | reference path | reads, library, curation |
|---|---|---|---|
| ACOT | `genotype_perf_bounded_20261015/inputs/acot/acot.gfa` | `GRCh38#0#chr14:73416151-73668669` | `m9h…/runs/prepare/reads/HG02083/cov5/`, `curation.LOO.txt` |
| CYP2D6 | `cyp2d6_shortread_pilot_v1/data/graphs/loo.NA19650.gfa` | `grch38#1#chr22:42031864-42245566` | `cyp2d6_shortread_pilot_v1/data/reads/NA19650/cov30/`, `tests/genotype_cyp2d6_curation.tsv` |
| C4 | `genotype_c4_production_11fe681/graphs/NA19650.loo.gfa` | `HG00096#1#haplotype1-0000024:31848049-32074186` | `genotype_c4_production_11fe681/reads/NA19650/cov10/` |

All paths are under `/Users/davide.bolognini/panvar/tests/results/`.

- **Persisted index.** None is required: every validated run prepares the locus in memory, with no prepared-index option.
- **Model inventories.** `m9h…/runs/models/LOO.block_alleles.tsv`, the M1-B and M9-J C4 inventories. They are needed only
  to regenerate declarations and to replay saved annotations.
- **Kept regardless.** The hard-linked long-read indexes (`h3e` / `h3g_c21`) are not used here, but must never be deleted.

## Validation state

| what | evidence |
|---|---|
| compound helper | `compound_annotation_test` 26/26, including Sol's two cases (R1, R2) and R3–R5 |
| reviewer's reproducer | `sol_repro.a1`: both cases `INCOMPLETE_DECISION_SET`; no `the_resolved_target` |
| replay against the Python reference | 34/34, and byte-identical to the saved M9-J replay in every saved column |
| saved caller annotations | ACOT (4 and 1 threads) and CYP2D6 unchanged (`RESOLVED`, 2 decision-set groups, 0 missing); the C4 control still has no rows |
| caller | JP1.1 (M9-G/H), sequence-report and e2e `sequence_first` gates pass; M9-J smoke pairs: pre-existing outputs identical with the annotation off and on |

## Limitations (to carry forward)

1. **Concurrent calls in one process are not validated.** The validated contract is one caller invocation per process. The
   M9-J helper-concurrency fixture (4 concurrent `compound_assess` calls) and within-sample threads (1 against 4 threads,
   identical outputs) are **not** multi-sample caller validation.
2. **No speed claim from M9-J timings.** M9-J's smoke timings overlapped storage compaction. Do not read the off/on
   wall-time differences (for example ACOT 170.3 against 161.1 s) as an annotation speed effect. The annotation's own
   measured cost is under 0.1 s.
3. **Compound masses are uncalibrated.** They are pool-conditional: explored candidates and recorded labellings only. There
   is no PASS threshold and no user-facing quality decision.
4. **The compound target is exact sequence.** It includes within-interval phase and is not structural equivalence. It is
   therefore stricter than the production target proposed in `QUALITY_TARGET.md`.
5. **The decision set is the caller's `top_groups`.** An `INCOMPLETE_DECISION_SET` row is reachable only if a decision-set
   group has no recorded labelling. No saved run has one.
6. **Known open scientific limits** (unchanged by this batch): C4 sv:7 allocation (3+1 against 2+2) is weakly identifiable
   from short reads under LOO; ACOT high-confidence errors (M9-D) are attributed but not removed by default; CYP2D6
   label-ambiguity no-calls.
7. **Configurations differ by locus** (table above). No combined profile has been validated.
