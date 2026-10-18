# S1-C: reporting-integrity fixes for `sequence_first_v1`

**Branch:** `s1-c-reporting-integrity`, from S1-B `6a2d98d`.

**What changed and what didn't.**
- No model, prior, search, default or cap change.
- Legacy outputs are unchanged.
- The S1-B results and its failed comparison records are preserved.

| commit | content |
|---|---|
| `527f54e` | fixes |
| `cead053` | validation tools and frozen manifest |
| `8dad559` | amendment 1: the e2e and closure-control declarations |

**Results:** `tests/results/s1c_reporting_integrity_20261018/README.md`.

## 1. Output-prefix lifecycle

**The five reports are run-owned.** They are `.sequence_locus.tsv`, `.sequence_blocks.tsv`, `.sequence_haplotypes.fa`, `.representation_intervals.tsv` and `.sequence_alternatives.tsv`. Their `atomic_write` staging names (`<file>.writing`) count too: exact suffixes, from `sequence_report_outputs()`.
- **At run start** they are removed with every other run output.
- **On a watchdog stop** they are removed too.
- **On a prepared-locus compatibility refusal** nothing changes. It is raised before the run owns its prefix.

**Publication.**
- **Timing.** Publication happens in the adaptive completion's stop, before the status commit.
- **Order.** Content files first; the locus record, which is the commit record, last.
- **On failure:**
  - every owned file of the prefix is removed;
  - a NO_PREDICTION record with `completion report_publication_failed` is attempted;
  - the status says `sequence_publication FAILED`;
  - the command exits 2.

## 2. Orientation

**The rule.** Each decision-set group's two sequences are spelled once and bound to their SHA-256. Labels, bases, hashes and lengths are canonicalized together:
- `representative_labelling_hap1` spells FASTA hap1, which has `representative_hap1_sha256` and `_bp`;
- `representative_labelling_hap2` spells hap2 likewise;
- an input homologue swap changes nothing.

**Intervals.**
- `hapA` is the group's SHA-256-smaller sequence and `hapB` the other.
- For the representative group, hapA = FASTA hap1, and its labelling is `representative_labelling_hap1`.
- Every interval row carries `hapA_sha256` and `hapB_sha256`.

## 3. Intervals: linear scan

**The rule.** One enclosing range per decision-set group: the first and last block where any explored labelling differs from the first one. Two labellings can differ at a block only if one of them differs from the first there.
- **Why this is the same as before.** With every pair of labellings considered, the S1-B pairwise merge always produced this single range.
- **What is kept.** Intervening context and the whole observed tuples.
- **Cost:** O(members × blocks).

## 4. Scope wording

**Representative.** The canonical member of the numerically inseparable decision set (`representative_scope`), not necessarily the strictly highest floating-point score.

**Dosage.** The report and its FASTA are the caller's decision *before* the dosage stage.
- **Locus record:** `dosage_scope`, `post_dosage_sequence` and `dosage_mode`.
- **FASTA header:** `scope=caller_decision_before_dosage_stage` and `dosage_mode=`.

When dosage repair changes the legacy reported genotype, this FASTA is **not** the post-dosage sequence. Dosage is not redesigned here.

## 5. Schema differences from S1-B

**Added:**
- the status row `sequence_publication` (`complete` / `FAILED`);
- exit code 2 for a failed publication;
- the locus rows `representative_scope`, `post_dosage_sequence` and `dosage_mode`;
- the FASTA header fields `scope` and `dosage_mode`;
- the interval columns `hapA_sha256` and `hapB_sha256`.

**Changed:** `dosage_scope` wording.

**Unchanged:** every other row and column.

## 6. Remaining limitations

- **Projection coverage.** EXPLORED coverage is not completeness. Intervals are conservative and not minimal.
- **Calibration.** There is no calibration and no PASS threshold.
- **Dosage.** The report describes the caller's decision before the dosage stage, and does not follow dosage repair.
- **Unsupported routes.** The unrepaired view and union-only diagnostics are refused.
- **Publication failure.** The NO_PREDICTION record written after a failure is best effort. On a failing filesystem it may itself be absent; the status row and exit 2 still say FAILED.
- **No caller-level coupled fixture.** As in S1-B, a caller-level tiny fixture for coupling across an unchanged context was not constructible. It is covered at the reporting level, and at real scale by CYP2D6.
