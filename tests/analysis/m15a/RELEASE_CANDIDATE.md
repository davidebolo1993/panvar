# Short-read release candidate (M15-A freeze)

**One reproducible short-read release candidate (RC).** It establishes neither calibrated confidence nor production
readiness. Model and performance research has stopped:
- M13 is INCONCLUSIVE;
- M14's tested route is rejected;
- M12's dosage-first rule is parked.

None of them enters the RC.

## 1. Source, binary and dependencies

**Source:**
- **Repository:** worktree `s1-d-panel-references`, branch `m15-a-release-candidate`.
- **Caller sources:** `src/` and `include/` are **byte-identical to the M10-A + M11-C checkpoint** (`ed10d53`; a zero-line
  `git diff ed10d53..HEAD -- src include`).
- **Changes after that checkpoint:** only `tests/` (probes, fixtures, analysis scripts), test targets in `CMakeLists.txt`,
  and records.

**Binaries:**

| item | value |
|---|---|
| `panvar` | sha256 `47042a5cbae5b6418d4c02f147c27772ec1a85791685dcf7c8f11bd8fa2af568` |
| `genotype_probe` (evaluation and diagnostics only) | sha256 recorded in `tests/results/m15a_release_candidate_20261021/MANIFEST.tsv` |

- **Why `panvar` has a new hash:** it was re-linked from the unchanged sources at M14-A's build unit. M11-D's binary
  (`c03cbba4…`) came from the same sources.
- **Behavioural identity:** M14-A's full-panel runs, made with this binary and M11-D's exact command, reproduce M11-D's
  archived calls on **15 of 15 cells** (rank-1 score, block paths, haplotype SHA-256, status, model fingerprint), with
  every release-candidate output present. This is `runs/reproduce`.
- **Build and regression** (`runs/build`, `runs/ctest`, `runs/ctest.a1`):
  - every target built;
  - **60 of 60** enabled tests passed in the frozen focused selection (2 are disabled in the suite);
  - three heavy tests (`genotype_report`, `genotype_branch_lifetime`, `genotype_block_delta`) were excluded for budget and
    remain a pre-tag blocker.

**Toolchain:** Apple clang 21.0.0 and CMake 4.3.3.

**Submodules:**

| submodule | revision |
|---|---|
| abpoa | `1b7cb00` |
| edlib | `ec2310e` |
| eigen | 3.4.0 `3147391` |
| minigraph | `2f569eb` |
| minimap2 | `e5066c7` |

zlib and threads are also linked.

**Index format.** The RC computes every prediction from the graph. **No persisted prepared-locus index is used** anywhere
in the archived evidence (M11-D, M12, M14-A).
- The persisted index (`.panel.pvpl`, I3-C v1 with the I3-D panel extension) is an **optional accelerator**.
- M10-A showed cold and index-loaded runs identical on all six cells.
- It is **not part of the RC contract**.

## 2. Effective configuration (`genotype-dev-m10a-v1`, M10-A arm D)

**Command line:**

`panvar genotype --preset short-read --profile <composed> --threads 4 --graph <gfa> --reference-path <ref> [--panel-curation <file>] --reads1 … --reads2 … --library <library.txt> --dosage-mode audit --out <prefix>`

**The composed profile** (`tests/analysis/m10a/compose_m10a.py`, arm D) has four layers, in order:
1. **The locus base profile** (`genotype_perf_bounded_20261015/inputs/profiles/<locus>_5.tsv`): `footprint_limit_bytes`
   6 GiB, `max_state_classes` 6e8, and `repair_max_neighbours` (CYP2D6 958, ACOT 1004; C4 at its default).
2. **The performance overlay** (`profiles/genotype_dev_performance_v1.tsv`): computational backends, caches and capacities.
3. **The M10-A overlay** (`profiles/genotype_dev_m10a_v1.tsv`):
   - `early_incumbent experimental_v1`;
   - `union_fragment_blocks fine_core_4096_experimental`;
   - `repair_round1_origin_admission reuse_only_experimental`;
   - `repair_junction_proposals experimental_v1`;
   - `output_schema sequence_first_v2`;
   - `compound_annotation contract_v2`;
   - `repair_joint_search off`.
4. **Operational resource overrides only:** `time_limit_seconds` = the declared watchdog − 20 s.

**Every effective setting and its source** is written to `<prefix>.settings.tsv`.

## 3. Classification of features

| class | features | evidence |
|---|---|---|
| **Enabled in the RC (affect computation)** | sequence-first v2 output (panel membership, extraction references, conditional FASTA); early evaluated incumbent; blocked union `fine_core_4096`; round-1 reuse-only origin admission; JP1.1 junction proposals; compound annotation contract v2; the dosage **audit** (diagnostic; never changes a prediction); the performance backends | M10-A: backend, route, cold/loaded and determinism identity (with the gaps below). JP1.1 changes science slightly: one sequence candidate removed from the pool, and one block's conditional GQ in CYP2D6/ACOT. Used unchanged in M11-D (30 cells), M12 and M14-A |
| **Optional experimental (not in the RC)** | persisted prepared-locus index; the distance copying prior (M4-B, default off); full-data block support v1x/v2 (M3-E/E.1, default off); `repair_joint_search` (M3-D, off); `--dosage-mode repair` (refused under v2; legacy only) | each default-off, gated in its batch |
| **Diagnostics only (never change predictions)** | `block_diagnostics diagnostic_v1` (M11-C, off/on identity verified); `timing_detail`; `genotype_probe` modes (evaluate, block-sequences, score-set, count-attrib incl. `--features`/`--skip-catalogue`, vt-score, panel-sketch); fixtures (`m12a_copy_balance`, `m13a_vt_gates`) | test/probe code only; no `src/` change |
| **Parked or rejected** | M12 dosage-first selection (parked); M13 variation-tolerant scoring (inconclusive; diagnostic code only); M14 representative panels (rejected); confidence candidates C1/C2 (negative against baseline, M11-D); I5 pruning and S2-B pre-placement (parked) | **none of these is reachable from `panvar genotype`:** they live in `tests/`, and `src/`/`include/` are unchanged |

## 4. Validated together, and what is not

**Validated together** (M10-A, computational identity):
- C4, CYP2D6 and ACOT at 5×;
- C4 and CYP2D6 at 30× (route and JP);
- ACOT 30× (JP, cold/loaded, repeat).

**Gaps kept explicit:**
- **ACOT 30×** has no reference or ordinary-route comparison: the ordinary route refuses at 500 M physical origins under
  6 GiB.
- **Backend equivalence** was tested at 5× only.
- **The M11-D accuracy evidence** is development evidence at 5× on simulated reads.

**Remaining consolidation requirements** (see `RELEASE_BLOCKERS.md`):
1. **The RC is "defaults + a declared profile",** and several keys carry `experimental` names. Promoting them to defaults,
   or shipping a named profile, is a release decision that has not been made.
2. **`repair_max_neighbours` is locus-specific** (CYP2D6 958, ACOT 1004). A release needs either per-locus declared
   profiles or a general rule.
3. **The time guard** is operational (watchdog − 20 s). A release default must be chosen.
4. **Record correction:** M14-A's `MANIFEST.tsv` holds the `panvar` hash taken before its build unit re-linked the binary.
   The calls ran on `47042a5c…` (identical sources). A correction is appended to M14-A.
