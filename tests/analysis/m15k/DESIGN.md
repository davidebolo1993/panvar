# M15-K design: v2 reporting for unrefined runs, the profile rules kept, honest index-builder provenance (frozen before implementation)

**Branch:** `m15-k-unrefined-v2`, from `m15-j-interface` (`fa88c00`).

**Scope:** reporting, profile acceptance and index provenance only.
- **Not done:** no change to scoring, candidates, the search, priors, refinement, filter or assessment rules, or the index's
  compatibility rules; no output file deleted, renamed or restructured; no tuning, no new evaluation families; no merge, push,
  tag or publication.

**Budget:** 2,400 s for the recorded units.

## 1. Unrefined runs report the v2 contract

**Route.** `--adaptive-repair off` / `--repair-rounds 0` end in the unrefined tail of `finish` (after `rank_set` and the
full-beam audit). That tail is kept as it is: its candidates, the union scores (objective components) and the union ranking,
and its legacy tables and status. Under `output_schema sequence_first_v2` (the preset), the tail additionally builds the v2
report from the **already-scored union pool**:
1. **The pool:** the union candidates of the ranking the tail used (`chosen`, the audited full beam when the audit ran),
   converted to pool candidates exactly as the refined route seeds its pool. Each carries its allele labels, the interned
   sequence digests, the union's `DiplotypeScore`, the sum of the haploid sequence priors (+ ln 2 for a heterozygous pair)
   and the exact path prior. Duplicate identities keep the first entry, as the refined route does.
2. **The decision:** `decide_pooled_genotype` over that pool. It ranks existing scores; it scores nothing.
3. **The report:** `build_sequence_report`, `annotate_v2` (panel membership, extraction references, the FASTA rule) and,
   when enabled, `compound_pooled`. All are reporting over the decision.

**No additional scoring or refinement.** The tail runs no anchor, seeding, census, repair search, tie branching or
representation guard; the unrefined tail never ran the guard and keeps reporting it `NOT_EVALUATED`.

**Phase follows v2.** The status's `phase_support` (which the assessment reads) is the pooled decision's phase support, the
same quantity the refined route reports, set before the tail writes its tables. The historical K-prefix phase audit
(`phase_audit_candidates`) is internal (§2) and is refused beside the v2 contract.

**Expected differences for unrefined runs** (relative to M15-J, which reported legacy there), and nothing else:
- the v2 files are written: `.sequence_locus.tsv`, `.sequence_blocks.tsv`, `.sequence_haplotypes.tsv`, and
  `.sequence_haplotypes.fa` only when the v2 FASTA rule writes it; `.sequence_alternatives.tsv`,
  `.representation_intervals.tsv`; `.sequence_compound.tsv` when enabled;
- the assessment columns and rows (`assessment_v1`), which are written only under v2;
- the exit status follows v2 (0 PREDICTION, 3 PROVISIONAL, 1 NO_PREDICTION);
- the legacy `.mosaic_haplotypes.fa` is removed by the existing v2 FASTA rule;
- `.settings.tsv`: `output_schema sequence_first_v2` with source `preset short-read (default)`, no longer inactive; and the
  `compound_annotation` row loses its inactive marker if it carried one;
- status `phase_support`, and the legacy block tables' phase-derived fields, can change where the pooled phase support
  differs from the union ranking's;
- tie outcomes follow v2. A tie whose candidates spell one genotype is a PREDICTION (`GENOTYPE_RESOLVED_SEQUENCE_SET`); one
  that disagrees at a block is reported `SEQUENCE_TIED` with that block `GENOTYPE_UNRESOLVED`. The legacy block-consensus
  tables are still written beside it.

**Unchanged:** `.mosaic_candidates.tsv`, `.mosaic_scores.tsv` (candidates, objective components, union ranking), and every
other legacy table apart from the phase-derived fields above.

**Historical behaviour** stays reachable internally: `output_schema legacy` through `GenotypeCommandOptions` (e2e, the
curation test). `sequence_first_v1` is still refused without repair rounds (internal, unchanged).

## 2. Profiles (decisions confirmed)

- **One optional profile:** a repeated `--profile` stays refused; an explicit flag overrides the profile.
- **The fixed contract:** a profile may restate `output_schema sequence_first_v2` (the evidence profiles do); any other value
  is refused.
- **The phase audit is internal:** a profile may restate `phase_audit_candidates 0` (the historical C4 profiles do); a
  positive value is refused. Internally, `phase_audit_candidates > 0` is refused under `sequence_first_v2` (it would be a
  second, public reporting exception) and stays available with the internal legacy contract.

## 3. Index builder provenance

**Defect:** CMake captures `git rev-parse HEAD` and the dirty state once, at configure time, so an index built after a later
commit records an earlier revision (M15-J: `e1b900d` recorded while `8b71661` was built).

**Fix:**
- A build-time step (a custom target run on every build, `cmake -P cmake/source_provenance.cmake`) writes
  `generated/panvar_source_provenance.h`, replacing it only when its content changes. `prepared_locus_file.cpp` includes the
  header, so a new revision recompiles that one file and relinks; no manual reconfigure is needed.
- **What is recorded:**
  - **`source_commit`:** `git rev-parse HEAD` at build time, or `unavailable` (not a git checkout, or git missing);
  - **`source_tree`:** `clean`, `dirty` (tracked changes, untracked files ignored as before) or `unknown`.
- **Format:** the two strings occupy the provenance section's existing fields. Layout, versions and digests are unchanged,
  so compatibility is unchanged.
- **Labels:** the description text (the build record's and status's provenance value) changes from
  `configure_commit=…;configure_dirty=…` to `source_commit=…;source_tree=…`. Indexes written before M15-K describe with
  their stored values (a configure-time commit and `0`/`1`).

**Producer provenance against compatibility identity:**
- **Provenance:** writer, `source_commit`, `source_tree`, compiler and build type describe who produced an index. It is
  reported and never compared.
- **Compatibility identity:** schema, preparation, panel-extension and marker-panel versions, graph SHA-256, reference,
  bubbles, `min_variant_bp`, curation, k and s, verified by digests. It decides whether an index is loaded. Unchanged.

**Tests:**
- a ctest (`panvar_source_provenance`) runs the script in a temporary git repository: clean commit A, then a tracked edit
  (dirty), then commit B (clean, new commit); outside git it records `unavailable`/`unknown`;
- a recorded unit commits in this worktree, rebuilds without reconfiguring, builds an index and requires the recorded
  `source_commit` to equal the new HEAD and `source_tree clean`.

## 4. Verification (units)

| unit | checks |
|---|---|
| `build` | all targets |
| `contract` | the CLI and curation contract tests, the e2e unrefined-v2 test, `panvar_source_provenance` |
| `ex.refined_{c4,cyp2d6,acot}` | the ordinary refined commands: every output identical to M15-J (M15-H rule), so predictions, sequences, scores and assessments are unchanged |
| `ex.unrefined_{c4,cyp2d6,acot}` | `--adaptive-repair off`: M15-K binary against M15-J binary. Candidates and scores identical; every other difference must be one listed in §1, reported per file; the v2 report consistent with the scores (the representative's score is the best union score); exit code from the v2 prediction |
| `ex.index_c4` | index build (recorded `source_commit` = HEAD, `source_tree clean`) and load; cold = indexed; payload equal to M15-J's apart from the provenance strings and the header digest |
| `ex.provenance_rebuild` | an empty commit, `cmake --build` without reconfigure, build an index: `source_commit` = the new HEAD |
| `ctest.full` | the full enabled suite on the final revision |

**e2e (in the suite):** an unrefined v2 test over unique, tied (block consensus), sequence-equivalent and phase-ambiguous
(prior-only) fixtures. Each runs the internal legacy contract and v2. Candidates and scores must be identical between them;
the v2 prediction, sequence status and phase must follow the pooled decision; the exit code must be the v2 exit.

The report separates the expected unrefined changes (§1) from any unexpected difference.
