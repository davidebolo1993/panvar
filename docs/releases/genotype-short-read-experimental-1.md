# panvar genotype: short-read experimental release 1 (proposed version 0.2.0-exp.1)

**Status: EXPERIMENTAL.** This release genotypes one structural locus at a time from paired short reads, against a
pangenome graph panel.
- It is **not** a validated or high-confidence production caller.
- Nothing in its outputs is calibrated confidence.
- It has **no reliable-call category and no quality threshold.**

Usage is documented in `docs/guides/genotype_short_read.md`, and the method in `docs/algorithms/genotype.md`.

## What is new

**An ordinary command without a profile:**

```
panvar genotype --graph locus.gfa --reads1 S.R1.fq.gz --reads2 S.R2.fq.gz --threads 4 --out S
```

**The accelerated configuration is the default.** The short-read preset carries the configuration the evaluation used:
- the exact computational backends and their cache budgets, which are gated identical on completed runs;
- junction proposals (JP1.1, a candidate-generation setting);
- the early evaluated incumbent;
- the sequence-first v2 report;
- a sampled 6 GiB footprint guard;
- a model-derived refinement-neighbourhood bound, replacing the old locus-specific caps.

**There is no time limit by default.** Locus-specific profiles remain only as an exact record of the evidence
configuration.

**Automatic library estimation** (the default; `--library auto` is the same):
- **Estimated from the reads:** the insert mean and sd, from at least 500 decisive anchored pairs, and the fragment-start
  intensity, from the marker depth fit.
- **Fixed:** the substitution, outlier and gap parameters.
- **On an invalid estimate:** the run refuses (`LIBRARY_REFUSED`), never a fallback.
- **Exact reproduction:** `--library FILE`, a declared library that is validated before anything is written.

**An opt-in dosage check:** `--dosage-check` runs the copy-number audit of the reported alleles.
- It never changes the prediction, and it adds run time.
- Without it, the dosage check is reported as NOT_REQUESTED.
- It is the only public dosage control: no flag or profile key repairs dosage.

**One public reporting contract:** ordinary genotyping reports `sequence_first_v2`, recorded as `output_schema` in
`.settings.tsv`. No flag or profile selects another contract. The unrepaired development view (`--adaptive-repair off`)
still reports `legacy`; whether it should report v2 is an open decision.

**A public panel index:**
- `--build-panel-index FILE` builds the sample-independent preparation (structural model and marker panel) once per
  locus.
- `--panel-index FILE` reuses it for each sample.
- `--graph` is still required, and is still read and hashed. Integrity, versions and identity are verified before any
  output is touched.
- A missing, corrupt or incompatible index is refused, never rebuilt silently.
- The index saves locus preparation (about 1–2.5 s per run on the development loci), not the sample's work.

**Conditional FASTA and extraction references:**
- `.sequence_haplotypes.tsv` gives each predicted haplotype's SHA-256, length, origin and panel-path extraction
  references.
- When both haplotypes are on the panel, no FASTA is written.
- When either is off-panel, `.sequence_haplotypes.fa` holds both.

**The primary assessment, `assessment_v1`,** is written in `.sequence_locus.tsv` and `.sequence_blocks.tsv`, mirrored in
the locus and block tables. Every check concerns the published prediction:

| field | values |
|---|---|
| `sequence_resolution` | sequence ambiguity, as distinct from label ambiguity |
| `genotype_check` | NONE_DETECTED · CONFLICT · NOT_ASSESSED |
| `dosage_check` | AGREES · SUSPECT · DISAGREES · UNAVAILABLE · NOT_REQUESTED · NOT_APPLICABLE; corrections are listed separately, as diagnostics |
| `read_support` | READS · PRIOR_DEPENDENT · NOT_ASSESSED |
| `phase` | locus-level |
| `findings`, `not_assessed` | every finding and every check that could not run |

Each locus value carries coverage counts. A missing check is never a failure, and it never hides a finding.

**Help and inventory:**
- `--help` covers everyday use.
- `--profile` is documented as optional overrides: the literal `key<TAB>value` format, the required `profile<TAB>name`
  row, a minimal working example and the precedence (default < preset < profile < flag). One profile per run.
- `--help-advanced` lists every option by purpose, with values, effective defaults, units and what each can change.
- The guide's option inventory is generated from the parser and checked by the test suite.

**Correctness fixes in this cycle:**
- A declared library is validated before the output prefix is touched.
- The process exit, the sequence record and the run status always state the same exit code (124 / 125 / 3 / 2 / 1 / 0).
- A requested dosage audit that does not complete is reported on stderr.
- **panphorte** (shipped with this integration): approximate folding now shares identical flank fragments instead of
  creating one per haplotype. Identical haplotypes keep identical walks.

## Breaking and default behaviour changes

Relative to the previously published `genotype-mosaic` (`cee58de`). Relative to `main` (`1d87162` plus the panphorte fix),
the whole mosaic-caller line below is new:

| behaviour | before | now |
|---|---|---|
| report contract | legacy (exit 0 call / 1 no-call) | **`sequence_first_v2`**: exit 0 PREDICTION, 3 PROVISIONAL, 1 NO_PREDICTION; `.mosaic_haplotypes.fa` is not written; the unrepaired view still reports legacy |
| a decision that ties only over labellings of one sequence pair | legacy no-call (`distinct_tie`) | a **PREDICTION** (exit 0), with `legacy_no_call_reason distinct_tie` kept |
| configuration | reference routes unless a profile was given | the accelerated preset (**JP1.1 can change predictions**; the early incumbent can turn a budget stop into PROVISIONAL) |
| refinement-neighbourhood bound | 670 (CYP2D6 and ACOT needed profiles) | `model` (the model's own maximum); a numeric cap is still accepted |
| footprint guard | none | 6 GiB, sampled (exit 125, or 3 PROVISIONAL) |
| dosage | off unless `--dosage-mode audit\|repair` | **off unless `--dosage-check`**, an audit only; `--dosage-mode` and `--dosage-audit` removed (refused) |
| reporting contract | `--output-schema legacy\|sequence_first_v1\|sequence_first_v2` | fixed `sequence_first_v2` (`--output-schema` removed, refused; a profile may only restate v2); the unrepaired view still reports legacy |
| removed spellings | accepted | refused by name: `--phase-audit-candidates`, `--repair-origin-cache-bytes`, `--repair-tie-seed-cache-bytes` (their profile keys remain), `--dosage-allocation-test`, `--library-diagnostic` |
| `--profile` given twice | applied in order, the status naming only the last file | refused |
| primary interpretation | the headline filters | `assessment_v1` (the headline fields are kept, unchanged, as **superseded** diagnostics) |
| panphorte with `--min-similarity` < 1 | one fragment node per haplotype | shared fragment nodes: **graphs and node IDs differ** from earlier outputs for approximately folded sites |

**Legacy compatibility retained:**
- `--adaptive-repair off` (the unrepaired development view, which still reports legacy);
- every expert profile key except `output_schema` (which may only restate `sequence_first_v2`), and the locus profiles of
  the evidence, unchanged;
- the `headline_filters_v2` columns and every output file, unchanged;
- the legacy and `sequence_first_v1` contracts and the dosage repair, through the internal API only (regression tests).

**Removed or refused flags (M15-J):** `--output-schema`, `--dosage-mode`, `--dosage-audit`, `--phase-audit-candidates`,
`--repair-origin-cache-bytes`, `--repair-tie-seed-cache-bytes`, `--dosage-allocation-test` and `--library-diagnostic` are
refused by name with what replaces them (exit 2, before the output prefix is touched). `--max-candidates` remains
refused with a pointer to `--mosaic-candidates`. No flag was renamed.

## Evidence

**Reserved ACOT evaluation (M15-B):** 18 reserved families, simulated paired 100 bp reads at 5×, a declared library and
leave-one-out curation.
- **The families are prior-exposed:** their assemblies were used in masked prior scoring.
- **They lie inside the full-panel graph frame:** their paths help define the blocks.
- **So this is not pristine independent validation.**

| measure | result |
|---|---|
| completed predictions | 18 of 18 (14 strict, 4 genotype-equivalent); no no-calls or resource refusals |
| informative blocks | 296 |
| per-block excess over the availability floor: median / q90 / q99 / max | 0 / 8 / 98 / **478** bp |
| mean excess per block | 6.26 bp (family bootstrap 95%: 2.30–11.36) |
| floor attainment | 66.6% |
| catastrophic blocks (≥ 1 kb) | 0. With 18 families, the one-sided 95% upper bound on the family-level rate is about 15–17% |

**Two unflagged errors of about 478 bp:**
- **Where:** block sv:9 in families HG02572 and NA19468, called 902 + 424 bp against a truth of 424 + 424 (excess 478 and
  475 bp, normalized excess 0.36 and 0.35).
- **How they were reported:** both had GENOTYPE_FILTER PASS, no repeated-unit or dosage warning, and support READS
  (HG02572) / PRIOR_ONLY (NA19468).
- **No reported field flagged them.** They were the largest reserved errors.

**C4 and CYP2D6: development evidence only** (M11-D: 12 C4 cells and 6 CYP2D6 donors at 5×, simulated, declared library).
- No independent held-out set exists for them.
- All five development catastrophes were at C4, CYP2D6 or development ACOT.

**The automatic library: development checks only.**
- At the three development cells, its estimates (insert median 348–350 bp, robust sd 40–41.5) gave the same predictions as
  the declared library.
- This is not an accuracy validation: no reserved or independent evaluation used `--library auto`.

**Not established:**
- real-read performance;
- other coverages;
- other loci;
- calibrated confidence.

## Known limitations

**Confidence:**
- GQ, error mass and support are conditional on the explored candidates and uncalibrated.
- There is no reliable-call category.
- PASS and NONE_DETECTED establish nothing beyond the checks that ran.

**Dosage audit:**
- **At the CYP2D6 development cell,** the audit's unit library refuses block sv:9 ("alleles disagree on their flanks or
  anchor trims"). The audit reports FAILED, and dosage is UNAVAILABLE. Other cells may hit the same refusal.
- **It is not run** when the legacy decision is a no-call.
- **Underpowered blocks** are UNAVAILABLE.

**Assessment:**
- **Legacy-tied v2 predictions:** their block checks are NOT_ASSESSED. Computing them would need new scoring and dosage
  work.
- **Phase** is assessed per locus, not per boundary.
- **The superseded headline fields keep known misreadings:**
  - they read NO_CALL beside a legacy-tied v2 prediction;
  - NOT_EVALUATED outranks LOW_QUALITY;
  - DOSAGE_FILTER PASS can describe an unapplied correction.

**Resources:**
- **The footprint guard** is sampled every 100 ms, not an allocation ceiling.
- **There is no default time limit.**
- **Work caps can end a run** as PROVISIONAL or with no prediction.
- **ACOT at 30×** has no reference-route or ordinary-route equivalence evidence; backend equivalence was tested at 5× only.
- **Not implemented:** memory-per-thread budgeting and disk spilling.

**API and interface:**
- One locus and one sample per invocation.
- The caller-only API (`run_genotype` without the command layer) writes no assessment view.
- In-process multi-sample reuse is API-only.
- The index load reads and hashes the graph, within the existing limits (1 GiB file, 3.5 GiB accounted decoding).
- Long reads are not implemented.

**Errors:**
- Only argument, profile and declared-library errors are refused before the output prefix is touched.
- Errors found while running are not transactional.

## Tested platform

All builds and tests ran on **one machine**:
- macOS 26.5.2 (Darwin 25.5.0), Apple silicon (arm64);
- Apple clang 21.0.0, CMake 4.3.3, Python 3.14.7;
- 4 threads, a 6 GiB envelope.

**No other operating system, compiler or architecture has been tested.** Portability is not validated.

**The full enabled regression suite** (64 tests; 2 disabled: `panvar_real_regressions`, `panvar_rebuild_stats`) passes on
the candidate revision.
