# M15-E: public panel-index build and load flags (frozen before execution)

**What this is:** the existing persisted prepared-locus index (I3-C format, I3-D panel section 16) exposed through the
basic CLI.
- **Reused unchanged:** the format, writer, reader, integrity checks, compatibility rules and prepared-locus execution
  route.
- **Not done:** no new scoring model, filter change, multi-sample concurrency, reserved rerun, merge, push, tag or cleanup.

**Logistics:**
- **Branch:** `m15-e-panel-index`, from `m15-d1-correctness` `bbc9052`.
- **Results:** `tests/results/m15e_panel_index_20261021/`.
- **Budget:** 2,400 s, including 260 s of disclosed development.
- **Resources:** 4 threads, one process at a time, 6 GiB.

## 1. Interface

| | build | genotype |
|---|---|---|
| command | `panvar genotype --graph G [--reference-path R] [--panel-curation C] [--bubbles B] --build-panel-index locus.pvpl` | `panvar genotype --graph G … --panel-index locus.pvpl --reads1 … --reads2 … --threads 4 --out S` |
| needs | no reads, no library, no output prefix | everything an ordinary run needs; `--graph` is still read and hashed (identity check and extraction references) |
| writes | the index, through the existing atomic publication (`O_EXCL` temporary, fsync, link; never overwrites); a build record on stdout (versions, identity, model fingerprint, size, SHA-256, prepare and write seconds). The one temporary, the staged sorted graph next to the index, is removed. No prefix is cleared and no sample report is written | the ordinary outputs; the status adds `panel_index`, `panel_index_file_bytes`, `panel_index_load_seconds`, `panel_index_panel_rebuild_reason`, besides the existing `prepared_locus_*` rows |
| refused (exit 2) | an existing destination (`exists`); `--panel-index` in the same invocation; sample arguments (`--reads1/2`, `--library`, `--out`/`-o`/`--out-prefix`, `--dosage-check`, `--dosage-mode`, `--dosage-audit`) | before the output prefix is touched: load failures (`io`, `truncated`, `checksum_mismatch`, `bad_magic`, `unsupported_schema`, `incompatible_preparation`, `limit_exceeded`, `reconstruction_mismatch`, …); identity mismatches (graph content, reference, bubbles, `min_variant_bp`, curation, k, s), named. Never a fallback to an in-memory preparation |

**Compatibility:** the schema version and the preparation-algorithm version decide it, plus the input-content identity.
- **Not:** the git revision (provenance is recorded, never gating).
- **Not:** reads, library, threads, `--dosage-check` or any other sample or search option.
- **A panel rebuilt** from the stored model (the I3-D path, for an index without a usable section 16) is reported on stderr
  and in the status, never as a loaded panel. Indexes built here always carry section 16.

**Never serialised:** library estimates, read counts or depth, candidate pools, placements, refinement state and dosage
results. The index saves locus preparation, not the sample's work.

## 2. Tests (frozen; development controls only)

| unit | what | frozen expectations |
|---|---|---|
| `build`, `cli_contract` | every target; the CLI contract | pass; the contract checks build + load and sample-argument refusals, and the help needles |
| `e2e.prepared_locus` | the existing prepared-locus gates plus **`test_panel_index_cli`** on the tiny locus | pass. The new test covers: a build without reads writes only the index; an existing destination refused byte-identical; build + load and four sample arguments refused, touching nothing; cold == indexed (declared library; automatic library; `--dosage-check` with the dosage status equal; a second sample; 2 threads); mismatched curation, k, `min_variant_bp`, reference and graph refused; missing, truncated and corrupted indexes refused; every refusal leaves three sentinels untouched and writes nothing; the index byte-unchanged after every use |
| `ex.build_{c4,cyp2d6,acot}` | build from the development graph (and curation) | exit 0 without reads; the unit's index directory holds only the index and its record; the record states the versions, the graph SHA-256 and the file's size and SHA-256 (verified); a second build is refused (`exists`), byte-identical. **Recorded:** size, prepare and write seconds |
| `ex.index_{c4,cyp2d6,acot}` | ordinary automatic-library run vs the same with `--panel-index` | equal exits; byte-identical science files (sequence tables, equivalence set, block tables, scores, candidates, intervals, alternatives, compound, sensitivity blocks), FASTA and locus table (prefix names aside), library estimate (timing aside), every effective setting and the prediction key; the status records a loaded model **and** panel. **Recorded:** wall time, `seconds_total`, load time and peak footprint |
| `ex.dosage_check_c4` | declared library + `--dosage-check`, cold vs indexed | as above, plus the dosage tables and status equal (timing aside) |
| `ex.two_samples_c4` | a second sample (NA19650) against the same C4 index | cold == indexed as above; the index byte-identical to its build record after every use |
| `ex.refusals_c4` | the C4 index with: a curation; `--kmer-size 27`; another reference path; the CYP2D6 graph; a missing file; half the file; a flipped byte; and build + load, and a build with reads | exit 2; the reason named (`panel curation`, `kmer_size`, `reference path`, `graph content`, `(io)`, `(truncated)`, `(checksum_mismatch)`); four sentinels untouched and nothing written; nothing built for the two refused builds |
| `ctest.full` | the full enabled suite on the final revision | every enabled test passes; disabled tests reported separately |

**Comparisons use narrow allowances only:**
- prefix-derived names (FASTA records, the locus table's file links);
- timing rows of the library estimate and dosage status;
- the status itself, whose provenance and timing are compared only through the named fields above.

**No universal speedup is claimed:** times and sizes are recorded, not checked.

**Failure handling:** a failed check is recorded and never retuned. An amendment is declared before any rerun.

## 3. Deliverable

- the implemented flags, help and guide (§8a);
- the build and load times, index sizes and end-to-end runtimes;
- the remaining limitations.

The batch stops for review.
