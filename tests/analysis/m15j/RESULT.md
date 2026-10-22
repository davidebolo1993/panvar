# M15-J result: one public reporting contract, audit-only dosage, removed legacy spellings refused by name; predictions, sequences, scores, settings, assessments and the index unchanged on C4, CYP2D6 and ACOT; full suite 64/64

**Run.**
- **Branch:** `m15-j-interface`, from `m15-i-release` (`31d99d7`).
- **Commits:** design frozen `8cd7dbb`; implementation and harness frozen `8b71661`; amendment A1 `b73e2d8` (harness only);
  result after.
- **Results:** `tests/results/m15j_interface_20261022/` (`LEDGER.tsv`, `MANIFEST.tsv`, `runs/`, and `baseline/panvar.m15i`, the
  M15-I binary copied before any change).
- **Budget:** 1,636.5 of 2,400 s (240 s of disclosed development).
- **Units:** 11 of 12 passed on the first pass; `ex.index_c4` failed on a harness expectation (§4) and passed as
  `a1.ex.index_c4`.
- **Not done:** no inference change, output removal, rename or restructuring; no donors, reserved data, tuning, merge,
  push or tag. WIP, stash, indexes and earlier results untouched.

| binary | SHA-256 |
|---|---|
| `panvar` (`8b71661`) | `57d1ffb2fba7505caeaf0dfbaf2d3f9e4a680f5ea19cb0300208a3d86a28c2a7` |
| baseline `panvar.m15i` (M15-H/M15-I) | `3f7e03a41c9d507d8e4600df819140202aa069b8feb717cad67abc96b17b462b` |

## 1. Completed changes

**One reporting contract.**
- `--output-schema` is removed: refused by name, exit 2.
- Ordinary genotyping reports `sequence_first_v2`, recorded as before in `.settings.tsv` (`output_schema`, source
  `preset short-read (default)`). The option inventory carries it as a fixed row ("reporting contract").
- A profile may only restate `output_schema sequence_first_v2`; `legacy`, `sequence_first_v1` or anything else is refused.
  Restating is allowed because the three evidence profiles and one development profile carry that row, and their bytes are
  cited by the sealed packages.
- The legacy and `sequence_first_v1` contracts remain reachable through the internal API only (`GenotypeCommandOptions`),
  where e2e and the curation test use them.

**Dosage: audit only.**
- `--dosage-mode` and `--dosage-audit` are removed (refused, pointing to `--dosage-check`).
- `dosage_mode`, `dosage_audit` and `dosage_check` are refused as profile keys by name.
- `--dosage-check` is unchanged: it audits the published alleles and never replaces the prediction.
- Dosage repair is reachable only through the internal API under the legacy contract. Adaptive refinement is untouched.

**Removed spellings, refused by name with what replaces them:** `--phase-audit-candidates`, `--repair-origin-cache-bytes`,
`--repair-tie-seed-cache-bytes` (their profile keys remain), `--dosage-allocation-test`, `--library-diagnostic`.
- `--max-candidates` stays refused.
- Removed options are checked before every other argument (even before `--help`), so nothing is read or written.

**`--profile`.**
- Documented in `--help-advanced`, the guide (§2) and the module page:
  - optional overrides, never required;
  - literal `key<TAB>value` lines, with `#` comments and blank lines ignored;
  - the required `profile<TAB>NAME` row;
  - a minimal working example;
  - precedence default < preset < profile < flag;
  - a key set twice, an unknown key and a malformed value are refused.
- **One behaviour change, declared in the design:** a repeated `--profile` is refused. Before, files applied in order and
  the status named only the last file.

**Documentation and tests.**
- **Docs:** help, advanced help, the generated inventory (guide §3), the guide (§2, §5, §8, §8a, §10), the module and
  algorithm pages, and the pending release notes.
- **Rewritten:** the CLI contract test (reporting contract, removed options, profile rules, help and inventory).
- **Moved to the internal API:** the curation test and one e2e command (the removed flags).
- **Annotated:** two unregistered historical gates (`genotype_dosage_shadow_identity.sh`, `genotype_dosage_mode_gate.py`).
  They now reproduce only with a pre-M15-J binary.

## 2. Verification

| unit | result |
|---|---|
| `build` | all targets |
| `cli_contract` | CLI and curation contract tests pass |
| `ex.ordinary_c4` / `_cyp2d6` / `_acot` | 22/22, 22/22, 21/21. The M15-G `a3.ex.complete_*` command with automatic library estimation: **every output identical** to M15-G apart from names, timing and the exact caches' scheduling-dependent counters (the M15-H rule). This covers predictions, sequences, candidate scores, every effective setting (`.settings.tsv`), the assessments and the old headline columns. `output_schema sequence_first_v2` recorded |
| `ex.dosage_check_cyp2d6` | 27/27: the M15-G `ex.dosage_cyp2d6` command's outputs identical to M15-G; prediction, sequences and scores identical to the same command without the flag; `dosage_check` UNAVAILABLE with it, NOT_REQUESTED without |
| `ex.dosage_check_c4` | 8/8: prediction, sequences and scores identical to `ex.ordinary_c4`; `dosage_check` AGREES |
| `ex.refusals` | 32/32: 16 cases (11 removed spellings, 4 forbidden profile settings, a repeated profile); each exits 2 with its message, and the output prefix is untouched (sentinel files byte-identical, no new file) |
| `ex.index_c4` | 37/38; failed only the index byte-identity expectation (§4) |
| `a1.ex.index_c4` | 39/39: the index equals M15-E's apart from the builder-revision field and its header digest; cold = indexed; both equal M15-E's runs (prediction and science files) |
| `ex.unrefined` | 13/13: `--adaptive-repair off` on C4 gives the same outputs from the M15-J and M15-I binaries; legacy report, recorded as the preset's v2 made inactive with its reason |
| `ctest.full` | **64 of 64 passed** (1,008.8 s); disabled, as before: `panvar_real_regressions`, `panvar_rebuild_stats` |

## 3. Unresolved dependencies

**The unrepaired view still reports `legacy`** (`--adaptive-repair off`, `--repair-rounds 0`). It is unchanged and documented
as the one exception, not a silent fallback. Findings (DESIGN.md §2.2):
- **The scores allow v2.** `decide_pooled_genotype` is a pure function of exactly scored pool candidates, and a pool of the
  union candidates alone needs **no new scoring**. The refined route builds its pool from exactly those already-computed
  quantities.
- **The decision rules differ.** The unrepaired route decides by its own rules over the same scores:
  - block consensus on ties;
  - the optional K-prefix phase audit (`phase_audit_candidates`, `READS_AUDITED`).

  v2 would publish the top-set decision (`SEQUENCE_TIED`, `GENOTYPE_RESOLVED_SEQUENCE_SET`, `GENOTYPE_UNRESOLVED` blocks) and
  the pooled phase support instead. Tie and phase outcomes would change; scores, candidates and the union ranking would not.

**The assessment reads compatibility reports.** `assessment_v1` reads `.block_genotypes.tsv` for:
- the verification that the legacy alleles equal the v2 genotype;
- the representation-guard `filter` (`genotype_check`);
- `genotype_support` (`read_support`);
- the dosage-stage fields (`dosage_check`).

It also reads the legacy status `.mosaic_status.tsv` (call class, library sensitivity, phase) and the dosage tables. The
full inventory and a retirement plan are in DESIGN.md §3:
1. the existing producers also write these already-computed facts to a v2 or internal facts file;
2. `read_facts` reads from it;
3. gate: `assessment_v1` byte-identical on the 606 saved prefixes and the three controls;
4. only then retire the writers, in a separately approved step.

Nothing was retired.

**Observation (pre-existing, not changed):** the index records the source revision at which the build directory was
**configured** (`e1b900d` here), not the revision of the sources built (`8b71661`). Recorded provenance can therefore be
stale.

## 4. Amendment (harness only; the failed unit stays recorded)

**A1.** `ex.index_c4` expected the index to be byte-identical to M15-E's.
- The header records the builder's source revision (40 hex characters, bytes 372–411) and a 32-byte header digest covering
  it (bytes 24–55), so 71 bytes differ by design: `e1b900d…` against M15-E's `42cc16f…`.
- Every behavioural check had passed.
- The check now allows exactly those two fields and requires every other byte to be identical. It passed, 39/39.

## 5. Decisions for Davide

1. **The unrepaired view's report** (§3):
   - **A (recommended):** v2 from the zero-round pool. No new scoring; tie and phase outcomes follow v2's rules.
   - **B:** keep legacy as the documented development-view exception.
   - **C:** remove the public refinement controls. Not proposed.
2. **The repeated-profile refusal**, introduced here (declared in the design): confirm, or restore applying several files in
   order.
3. **`phase_audit_candidates`:** it is kept as a profile key and only affects the unrepaired view, so it goes with decision 1.
4. **Retiring the compatibility reports:** approve the §3 plan as a separate batch, or not.

**Stopped for review.**
