# M15-H result: help and documentation describe the implemented interface; one generated option inventory, checked against the parser and the preset; outputs and assessment unchanged

**Run.**
- **Branch:** `m15-h-help`, from `e1b900d`.
- **Commits:** freeze `c4ab40e`; amendments A1 `07f7543` and A2 `f83fd4e` (preservation harness only); result after.
- **Results:** `tests/results/m15h_help_20261022/` (the final help texts are in `runs/ex.help_text/`).
- **Budget:** 1,201.3 of 1,800 s, including 90 s of disclosed development.
- **`panvar` SHA-256:** `3f7e03a41c9d507d8e4600df819140202aa069b8feb717cad67abc96b17b462b`.
- **Not done:** no model, search, filter-rule or default change; no flag added, removed or renamed; no reserved rerun,
  merge, push or tag.

## 1. Results

| unit | result |
|---|---|
| `cli_contract` | pass. Checked: the basic help's order and content; the corrected claims; the removed headline explanation; the six advanced groups and the terms; inventory coverage of every public setting and the profile keys; inventory defaults equal to the parser and the preset; **the guide's inventory equal to the generated one** |
| `ex.help_text` | both texts captured: `--help` 60 lines, `--help-advanced` 322 lines |
| `ex.preserved_c4` → `a1` → `a2` | first run: harness slice bug (compared nothing). A1: 19 of 20 (the exact caches' scheduling-dependent counters were missing from the volatile set; every science file identical). **A2: 20 of 20.** The no-profile C4 command is identical to M15-G's run in every output: settings, sequences, candidate scores, assessment, and the old headline columns |
| `ctest.full` (final binary) | **64 of 64 passed** (1,002.7 s); disabled, reported separately: `panvar_real_regressions`, `panvar_rebuild_stats` |

## 2. The final basic help (`panvar genotype --help`)

```
Usage:
  panvar genotype --graph <locus.gfa> --reads1 <R1.fq[.gz]> --reads2 <R2.fq[.gz]> --threads 4 --out <prefix>
  panvar genotype --graph <locus.gfa> --build-panel-index <locus.pvpl>        (once per locus)
  panvar genotype --graph <locus.gfa> --panel-index <locus.pvpl> --reads1 ... --reads2 ... --out <prefix>

Genotypes one structural locus from paired short reads: a diplotype (two whole-locus haplotypes) over the
graph's panel, projected to one allele pair per block. EXPERIMENTAL. No profile is needed: the short-read
preset is applied, and the library is estimated from the reads (insert size and fragment-start intensity
only; the error parameters stay fixed; an invalid estimate refuses the run, never a fallback).

Everyday options:
  --graph FILE               the locus graph (alias -i)
  --reads1 FILE --reads2 FILE  paired FASTQ (.gz accepted), mates in the same order
  --out PREFIX               output prefix (default genotype; alias -o)
  --reference-path NAME      coordinate reference (default: the single path of the graph's RS:Z sample;
                             refused with the candidates listed when absent or ambiguous)
  --threads N                worker threads (default 1); results do not depend on them
  --dosage-check             optional copy-number audit of the reported alleles. It never changes the
                             prediction and adds run time (not necessarily cheaper than the call); without
                             it dosage_check is NOT_REQUESTED
  --build-panel-index FILE   build the sample-independent panel index for this graph, reference,
                             curation, bubbles and marker settings, then stop; never overwritten
  --panel-index FILE         genotype against a built index: it saves locus preparation, not the rest of
                             the run. --graph is still required and still read and hashed; the index's
                             integrity, versions and identity are verified before the output prefix is
                             touched, and a missing, corrupt or incompatible index is refused (exit 2),
                             never rebuilt. Rebuild it when the graph, reference, curation, bubbles,
                             --min-variant-bp, --kmer-size or --syncmer-s change
  -q, --quiet                no progress messages (warnings are still printed)
  -h, --help                 this help
  --help-advanced            every option by purpose: library file, curation and bubbles; the copying
                             prior and marker model; candidates and bounded refinement; resource limits,
                             caches and backends; reporting; legacy and development controls

Outcomes: PREDICTION (a completed ranking; exit 0), PROVISIONAL (the best fully evaluated candidate after
a budget stop or declared cap: the search is incomplete; exit 3), NO_PREDICTION (exit 1).
Sequences: .sequence_haplotypes.tsv always gives each haplotype's SHA-256 and origin; when both haplotypes
are on the panel, no FASTA is needed (extraction references point into the graph); when either is
off-panel, .sequence_haplotypes.fa holds both.

Primary assessment (assessment_v1), in .sequence_locus.tsv and .sequence_blocks.tsv (mirrored in the locus
and block tables). Every check concerns the published prediction:
  sequence_resolution  RESOLVED, GENOTYPE_RESOLVED_SEQUENCE_SET, SEQUENCE_TIED (sequence ambiguity)
  genotype_check       NONE_DETECTED, CONFLICT or NOT_ASSESSED (+ coverage)
  dosage_check         AGREES, SUSPECT, DISAGREES, UNAVAILABLE, NOT_REQUESTED or NOT_APPLICABLE; corrections
                       are listed apart, as diagnostics (dosage_alternatives)
  read_support         READS, PRIOR_DEPENDENT or NOT_ASSESSED: conditional on the explored candidates;
                       prior-dependent does not mean the prior alone chose the call
  phase                READ_SUPPORTED, PRIOR_ONLY, UNRESOLVED, NOT_ASSESSED or NOT_APPLICABLE, assessed once
                       per locus (not independently at each boundary)
  findings, not_assessed  every concrete finding and every check that could not run, with block ids
Read a locus summary with its coverage and not_assessed list: a missing check is not a failure and never
hides a finding. These assessments are not calibrated reliability, and there is no quality threshold. The
older headline fields (genotype_filter, dosage_filter, phase_filter, call_summary) are kept as superseded
compatibility diagnostics.

Exit status: sequence-first (default) 0 PREDICTION, 3 PROVISIONAL, 1 NO_PREDICTION; legacy schema 0 a
genotype was emitted, 1 no-call; both: 2 argument, input or output error; 124 time limit, 125 footprint
limit (the 6 GiB footprint guard is sampled, not an allocation ceiling). Full guide:
docs/guides/genotype_short_read.md.
```

## 3. `--help-advanced` and the option inventory

`--help-advanced` is generated from `genotype_option_inventory()`. It lists six purpose groups:
1. input and panel overrides;
2. the copying prior and the marker / placement model;
3. candidate generation and bounded refinement;
4. resource limits, work caps, cache budgets and exact backends;
5. reporting and diagnostics;
6. legacy and development compatibility.

It then lists the terms.

**Every entry gives:** its values, the effective short-read default (computed; the built-in value where it differs), its
units, what it changes (inference / computation / reporting / completion), its status, and its aliases and profile key.
For example:
- `--repair-max-neighbours model|N` — default `model (the model-derived bound) (built-in 670)`;
- `footprint_limit_bytes` — default `6 GiB (built-in 0)`;
- `union_scoring_backend` — default `indexed_fine_experimental (built-in reference)`, `computation`.

**The inventory: 86 entries.**

| group | entries |
|---|---|
| everyday | 10 |
| input | 6 |
| model | 6 |
| search | 11 |
| resources | 37 |
| reporting | 10 |
| legacy and development | 6 |

| status | entries |
|---|---|
| everyday | 10 |
| supported | 21 |
| expert profile only | 49 |
| legacy | 2 |
| development | 4 |

The full table is in the guide's §3, between generated markers, and the test fails if it drifts.

**Terms stated explicitly** (in the advanced help, the guide's §3, and in brief in the basic help):
- donor redraw is not a recombination probability;
- the initial mosaic-candidate count does not cap refinement-generated mosaics;
- refinement is bounded, not globally exhaustive;
- `min_variant_bp` constructs sites and does not filter calls;
- library estimation learns the insert size and intensity only;
- support is conditional on the explored candidates;
- prior-dependent does not mean the prior alone chose the call;
- phase is assessed per locus;
- summaries are read with their coverage and `not_assessed`;
- nothing is calibrated, and there is no threshold.

## 4. Documentation

- **Guide:**
  - §2: preset prose without duplicated numbers;
  - **§3: the generated inventory**, columns and terms;
  - §4: the resource kinds, without hand-typed values;
  - §5: reporting (assessment primary, headline superseded);
  - **§10: recommended future removals and aliases** (not implemented).
- **Module page:** option table updated (index flags; the initial-mosaics wording; the `dosage_check NOT_REQUESTED`
  default); the stale `output_schema` default corrected (built-in legacy, preset v2).
- **Code:** the settings table's experimental-key mapping is extracted unchanged (`experimental_setting_value`), so the
  help's built-in defaults and `.settings.tsv` share one source. `.settings.tsv` is byte-identical (preservation unit).

## 5. Recommended future removals and aliases (guide §10; none implemented)

**Remove after the experimental release:**
- the superseded headline columns, or move them behind a flag;
- `--dosage-audit`;
- the refused `--max-candidates` message;
- the `legacy_i2a` store layout;
- the diagnostic flag spellings `--repair-origin-cache-bytes` / `--repair-tie-seed-cache-bytes`.

**Consider:** `--time-limit` / `--memory-limit` for the two run guards (now profile-only), and `--curation` as an alias of
`--panel-curation`.

## 6. Remaining-release checklist

`RELEASE_CHECKLIST.md`.

**Done:** the basic command; preflight and exit agreement; the panel index; `assessment_v1`; help and inventory.

**Open:**
- the headline-column removal timing;
- legacy-tied block checks;
- the caller-only API view;
- the run-guard flags and aliases;
- the removals;
- the known ACOT 30× and CYP2D6 audit limitations;
- the release notes;
- merge, tag and publication (Davide).

**Stopped for review.**
