# M9-J result: opt-in compound interval annotations integrated under contract v2; with it off nothing changes, with it on only the new table and declared provenance appear; replay agrees with the Python reference 34/34

**Run.** Branch `m9-j-compound-integration`.
- **Commits:** implementation and freeze `ea71e0a`; Amendments 1–4 `cc8d16d`, `9f9582b`, `eb63c71`, `41a7d66`.
- **Results:** `tests/results/m9j_compound_integration_20261020/`.

**Execution: 1,112.7 of 1,800 s,** including 75 s of disclosed pre-freeze development.

**Recorded failures** (each a harness or runner defect of mine, fixed by a declared amendment; no caller or library change
followed any of them):

| amendment | failure | fix |
|---|---|---|
| A1 | the replay harness dropped trailing empty fields, which are empty alleles | keep trailing empty fields |
| A2 | the harness's `stod` rejected subnormal printed masses | use `strtod` |
| A3 | the smoke command builder expected `--out-prefix`, but M1-B used `--out`; an empty `runs/smoke.cyp.off/` is left | accept either flag; CYP2D6 and C4 pairs rerun as `.a3` |
| A4 | some smoke-check status keys were not declared | inspected cache-work split, profile provenance and thread-partition keys declared explicitly; no scientific key excluded |

## Implementation (about 470 lines of library and caller code, plus a 258-line test and replay harness)

| file | lines | content |
|---|---|---|
| `include/panvar/compound_annotation.hpp`, `src/compound_annotation.cpp` | 94 + 272 | model-agnostic helper |
| `src/genotype_command.cpp` | +92 | options and refusals, declaration loading, pooled / provisional / no-prediction hooks |
| `src/genotype_sequence_report.cpp` | +11 | file write before the locus commit record, owned suffixes, status keys |
| headers | 6 | — |

**Profile keys:** `compound_annotation off|contract_v2` (default off) and `compound_interval_declarations <file>`.
- **Refused:** `contract_v2` without `sequence_first_v2`; declarations without the annotation; a malformed declaration
  file (before the prefix).

**Output:** `<prefix>.sequence_compound.tsv`, one row per span. Each row gives:
- boundaries and sources;
- the target scope (exact unordered concatenated sequences, including within-interval phase, not structural);
- status, target resolution, and decision-set projections with digests;
- target SHA-256 and length;
- mass reference and competitor;
- support and error lower and upper bounds;
- group, missing, ambiguous and zero-mass counts;
- original total and normalisation discrepancy;
- the explored-pool limitation string.

**Lifecycle:** the file is an owned sequence-report output. It is removed at run start and on publication failure, and
written before the locus commit record. Provisional reports get `NOT_ASSESSED_PROVISIONAL` rows; no-prediction reports get
no rows.

**Unchanged:** predictions, candidates, scores, priors, search, block rows, phase, GQ, filters, FASTA and dosage.

**Declarations:** `compound_interval_declarations_v1` artifacts, generated from **existing** inventories with the model
fingerprint and ordered block ids.

| model | fingerprint | blocks | intervals |
|---|---|---|---|
| ACOT | `12cbe6be…` | 19 | sv:5..ctx:6 (2,045 equivalences) |
| CYP2D6 | `e41f77fd…` | 19 | sv:3..sv:5 (42 equivalences) |
| C4 | `14851130…` | 11 | none: 0 equivalences, so a genuine control |

- At run time they are validated against the active model and coalesced with the run's representation intervals:
  overlapping spans merge, adjacent spans stay apart.
- There are no locus names in library code and no per-sample search. The persisted index is unchanged.

## Gates

| gate | result |
|---|---|
| `compound_annotation_test` | **21/21**: M1–M8 ported, including the realisable mixed projection and the phase counterexample; homologue exchange; duplicates; decision-set resolved and unresolved; declaration fingerprint and block-id mismatch, non-consecutive and malformed files; overlap coalescing; 4 concurrent invocations identical |
| JP1.1, sequence-report and e2e `sequence_first` gates | pass (publication and failure handling) |
| **replay against Python** (34 saved M9-I cases, saved tables only) | **34/34** within the declared tolerances: status, competitor and counts equal; bounds zero iff zero, otherwise ≤ 1e-12 relative. All 34 targets RESOLVED |

## Smoke pairs (cases frozen before running; off = the source's recorded command, on = the same plus the compound profile)

| case | wall s off / on | identity of pre-existing outputs | annotation |
|---|---|---|---|
| ACOT M9-H L4 (JP1.1 on, as saved) | 170.3 / 161.1 | **identical** | sv:5..ctx:6 ASSESSED, RESOLVED (1 decision-set projection, over both decision-set groups); support / error 0.982 / 0.0178; 86,434 zero-mass groups; **target equals the truth interval** |
| CYP2D6 M1-B cyp_30 | 123.0 / 122.5 | **identical** | sv:3..sv:5 ASSESSED, RESOLVED; error 1.0e-7; **target equals the truth interval** (block labels 118 + 118 bp) |
| C4 M1-B c4_10 (control) | 49.3 / 49.2 | **identical** | declaration VALID, no interval, no representation interval: 0 rows |

**Identity: 76/76 checks.**
- **Compared:** every pre-existing file, with each arm's prefix normalised.
- **Allowed to differ:** in the rounds and branch tables, only non-scientific columns; in the status, dosage and settings
  keys, only `compound_*`, timing, memory and cache-work keys and the profile provenance (A4).
- **New in on runs only:** `sequence_compound.tsv`. Off runs write no compound file.

**Caller against replay of its own saved tables** (full precision against 12-digit printed masses): status, competitor,
resolution and bounds agree, within 1e-9 relative, zero iff zero.

**1 against 4 threads** (ACOT on, 198.4 s at 1 thread): **31/31**. The annotation file is byte-identical, and so is every
scientific output. Only thread-count and work-partition status keys differ.

**Annotation cost:** 0.07 s (ACOT, about 87,000 groups), 0.005 s (CYP2D6) and 0.0004 s (C4).

**Scientific-output differences:** **none.**

## Handoff

**The opt-in compound annotation is implemented, gated and smoke-validated on three loci.** It changes no scientific
output, and its cost is under 0.1 s.

**What it does not do:**
- it implies no PASS threshold or user-facing quality decision;
- its masses remain pool-conditional and uncalibrated.

**Next (as stated):** define and calibrate the production quality target. Wherever a declared interval exists, the
compound unit is the natural evaluation unit.

**Status:** no default promotion.
