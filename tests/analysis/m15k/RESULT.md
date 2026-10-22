# M15-K result: unrefined runs report v2 from the already-scored union pool; the profile rules kept; the index records the source actually built; refined outputs unchanged; full suite 65/65

**Run.**
- **Branch:** `m15-k-unrefined-v2`, from `m15-j-interface` (`fa88c00`).
- **Commits:**
  - design frozen `3f77834`;
  - implementation and harness frozen `6d1102d`;
  - `38ab8f5`, the empty commit made by `ex.provenance_rebuild`;
  - amendments A1 `ead726c` and A2 `1d0c1c6`, harness only;
  - result after.
- **Results:** `tests/results/m15k_unrefined_v2_20261022/`, including `baseline/panvar.m15j`, the M15-J binary.
- **Budget:** 1,807.5 of 2,400 s (300 s of disclosed development).
- **Units:** 10 of 11 passed on the first pass. `ex.index_c4` failed on two harness expectations; A1 then also failed, on
  two others; it passed as `a2.ex.index_c4` (§4).
- **Not done:** no scoring, search, prior, filter, assessment or compatibility-rule change; no output file deleted, renamed
  or restructured; no tuning, new families, merge, push or tag.

| binary | SHA-256 |
|---|---|
| `panvar`, the final build (`1d0c1c6`; C++ identical to `6d1102d`) | `5c4707a575f23ed4eab3b392143ad9ed1dbef8d9a3feb1e7ee78f12a9cc182e6` |
| baseline `panvar.m15j` | `57d1ffb2fba7505caeaf0dfbaf2d3f9e4a680f5ea19cb0300208a3d86a28c2a7` |

## 1. Completed changes

**Unrefined runs report the v2 contract** (`--adaptive-repair off` / `--repair-rounds 0`).
- **How:** the unrefined tail builds the report from the union candidates it already ranked. Each is converted exactly as
  the refined route seeds its pool: labels, sequence digests, the union `DiplotypeScore`, the haploid sequence priors plus
  ln 2 for a heterozygous pair, and the path prior. Then come `decide_pooled_genotype`, `build_sequence_report`,
  `annotate_v2` and the compound annotation.
- **Nothing is added:** no anchor, seeding, census, repair search, tie branching or representation guard.
- **Phase:** the status's phase support is the pooled decision's.
- **Internal, legacy contract only:** the K-prefix phase audit and development prefix ladders are refused beside v2.
- **Historical behaviour:** reachable through the internal API (`output_schema legacy`).

**Profiles.**
- A repeated `--profile` stays refused; a flag overrides the profile.
- `output_schema` may only restate `sequence_first_v2`, so the evidence profiles still parse.
- `phase_audit_candidates` may only restate 0, so the historical C4 profiles still parse; a positive value is refused.
- No new public reporting exception exists.

**Index builder provenance.**
- `cmake/source_provenance.cmake` runs on every build. It writes `generated/panvar_source_provenance.h`, rewritten only
  when its content changes:
  - **`source_commit`:** `git rev-parse HEAD`, or `unavailable`;
  - **`source_tree`:** `clean`, `dirty` or `unknown`.
- The index stores them in the provenance section's existing two string fields. The layout, versions, identity and load
  rules are unchanged.
- The build record now reports `producer_source_commit`, `producer_source_tree` and `producer_build_type`, separately from
  the compatibility identity. Provenance is reported, never compared on load.

**Documentation:**
- **Guide:** §2 on unrefined runs and the outcome mapping, §5, §8, and §8a on producer provenance against identity.
- **Module page:** option table and the v2 section.
- **Also updated:** release notes (two new rows), help and the generated inventory.

## 2. Unrefined reporting: expected changes, and any unexpected difference

**The development controls** (`ex.unrefined_*`: `--adaptive-repair off`, M15-K binary against M15-J binary; every
difference classified against DESIGN.md §1):

| cell | M15-J | M15-K | `.mosaic_candidates.tsv`, `.mosaic_scores.tsv` | expected differences | **unexpected** |
|---|---|---|---|---|---|
| C4 | COMPLETE, exit 0 | PREDICTION RESOLVED, exit 0 | identical | 56 | **0** |
| CYP2D6 | COMPLETE, exit 0 | PREDICTION RESOLVED, exit 0 | identical | 54 | **0** |
| ACOT | COMPLETE, exit 0 | PREDICTION RESOLVED, exit 0 | identical | 54 | **0** |

**The expected differences on all three:**
- the v2 files are new (`.sequence_locus.tsv`, `.sequence_blocks.tsv`, `.sequence_haplotypes.tsv`/`.fa`,
  `.sequence_alternatives.tsv`, `.representation_intervals.tsv`, `.sequence_compound.tsv`);
- `.mosaic_haplotypes.fa` is no longer written (the v2 FASTA rule), and `main_fasta` in the locus table now links the v2
  FASTA;
- the `assessment_v1` columns and rows are added;
- `.settings.tsv`: `output_schema` legacy → `sequence_first_v2`; `compound_annotation` off → `contract_v2`; `early_incumbent`
  stays off, its inactive reason now "needs the adaptive repair";
- `sequence_*` and `compound_*` rows are added to the status;
- the library estimate differs in its timing row only.

**Also checked:**
- On every cell the representative's score equals the best union score.
- Phase is `PRIOR_ONLY` before and after; the pooled phase equals the historical one on these cells.
- The assessment reads `GENOTYPE_NOT_ASSESSED` (no representation guard on this route), `SUPPORT_PRIOR_DEPENDENT`,
  `PHASE_PRIOR_ONLY`.

**Fixtures** (end-to-end `test_unrefined_v2`, legacy against v2 on one binary; candidates and scores identical in every case):

| case | historical (legacy) | v2 |
|---|---|---|
| unique | COMPLETE, exit 0 | PREDICTION RESOLVED, exit 0, phase READS_WITHIN_CANDIDATES |
| genotype tie (two phasings of one genotype) | BLOCK_CONSENSUS, exit 1 | PREDICTION GENOTYPE_RESOLVED_SEQUENCE_SET, exit 0, phase UNRESOLVED |
| split tie (three heterozygous pairs disagreeing at a block) | BLOCK_CONSENSUS, exit 1 | PREDICTION SEQUENCE_TIED, exit 0, the block GENOTYPE_UNRESOLVED |
| label-equivalent (one sequence, several labellings) | COMPLETE, exit 0 | PREDICTION RESOLVED, both blocks REPRESENTATION_AMBIGUOUS |
| prior-only phase | COMPLETE, exit 0, PRIOR_ONLY | PREDICTION RESOLVED, PRIOR_ONLY |

The two tie rows are the documented change of outcome: block consensus (exit 1) becomes a v2 prediction (exit 0) with the
tie reported in the sequence status and the blocks.

**Unexpected differences: none.**

## 3. Verification

| unit | result |
|---|---|
| `build` | all targets; the generated header named the commit built |
| `contract` | CLI and curation contract tests, `panvar_source_provenance` (clean, dirty, new commit without reconfigure, unchanged revision not rewritten, outside git), and the e2e unrefined-v2 gate |
| `ex.refined_{c4,cyp2d6,acot}` | 20/20, 20/20, 19/19: **every output identical to M15-J** (predictions, sequences, scores, settings, assessments), apart from names, timing and scheduling-dependent cache counters |
| `ex.unrefined_{c4,cyp2d6,acot}` | 6/6 each (§2) |
| `ex.provenance_rebuild` | 4/4: an empty commit, `cmake --build` with no reconfigure, and the next index records the new HEAD, clean |
| `a2.ex.index_c4` | 46/46: the index records the HEAD built, clean. It equals M15-J's apart from the provenance strings, the header digest and two length fields grown by exactly 4 bytes (`1` → `clean`). The compatibility identity is unchanged; cold = indexed; the indexed run equals M15-J's |
| `ctest.full` | **65 of 65 passed** (1,019.2 s), including the new `panvar_source_provenance`; disabled, as before: `panvar_real_regressions`, `panvar_rebuild_stats`. Run at `38ab8f5`: the C++ is identical to the final revision, and A1/A2 changed only `tests/analysis/m15k` |

## 4. Amendments (harness only; failed units stay recorded)

**A1:** `ex.index_c4` failed 44/46.
- It assumed M15-J's stored tree string was `0`; M15-J's build directory had been configured with a dirty tree, so it is
  `1`. The tree is now read from its own length-prefixed field.
- The indexed run's status differs from M15-J's in exactly two rows by design: the index file's SHA-256 and its path.
- A1 itself failed 43/45, for the two reasons below.

**A2:**
- **Stale binary:** the binary had last been built at the provenance-rebuild commit, so the index recorded that commit,
  correctly the source actually built, while the check compared with the later HEAD. The unit now rebuilds first.
- **Length fields:** the provenance section grew by 4 bytes, so the payload length and the section's length field differ.
  They are allowed only when they grow by exactly the tree string's change.
- A2 passed 46/46.

## 5. Notes for review

**Unrefined runs and the representation guard.** Unrefined runs never ran the representation guard. Under v2 they report
`genotype_check NOT_ASSESSED`, with the reason per block. Running the guard there would score substituted candidates, so
it was not added.

**Indexes written before M15-K** describe their provenance with the configure-time values (`source_commit` holds that commit,
`source_tree` `0`/`1`). They load unchanged.

**Stopped for review.** The output reorganization (M15-L) follows as a separate batch.
