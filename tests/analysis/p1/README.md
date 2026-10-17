# P1: consolidated development baseline (performance work, before S1)

**What it is:** a development baseline, not a production release. It is branch `p1-consolidated-baseline`, from `s1-a-sequence-reporting` `5a4285f`; the contract and fixture work there is kept unchanged.

**What it adds:** the reviewed, exact performance work of I1–I4, taken from the research stack `i4-c-seed-reuse` `4cad112`. Every backend stays opt-in.

**Not changed:**
- public defaults: every imported backend defaults to off. On C4 and CYP2D6, the reference route was byte-identical to the production sources;
- the scientific model, priors and parameters;
- S1 reporting, which is not implemented here.

| commit | content |
|---|---|
| `13b2837` | snapshot: `src/`, `include/` and the `tests/*.cpp` regression suites byte-identical to `4cad112`; research-only CMake targets dropped |
| `392fc0b` | removal of the I-series measurement-only mechanisms (below) |
| `1cecf7c` | development configuration, the reviewed comparator chain (unchanged from `4cad112`), validation runner and frozen manifest |
| `093f16c` | validation amendment 1 (declared before further execution; see the results README) |

## 1. Selected-change inventory

**In the baseline.** Each backend is exact, and off unless its key is set.

| component | origin | key (default) |
|---|---|---|
| shared-content placement index, prefilter-decision cache, interior-origin reuse, retained incremental scoring, fine discovery units | I1–I1.4 (`206ba6d`, `aa64ebf`, `8ff4244`, `80e9445`, `aabba41`) | through `union_scoring_backend` |
| indexed union backend | I2-A `f21ff0e` | `union_scoring_backend` (`reference`) |
| fine adaptive window discovery | I2-B `87bc1e6` | `adaptive_window_discovery` (`reference`) |
| capacity-charged exact union store, stage timers | I2-C `e02c00d` | `indexed_union_store_layout` (`exact`), `timing_detail` (0) |
| fused anchor seed pass | I2-D `d665288` | `adaptive_fused_seed` (0) |
| marker-proposal contribution reuse | I2-E `ae390e2` | `marker_proposal_reuse_table_bytes` (0) |
| prepared library terms, with the I2-G tally correction | I2-F `68797bd`, I2-G `15b9148` (tally only) | `combination_prepared_library_terms` (0) |
| fixed-geometry gapped recursion; the scaled reference is kept | I2-H `8058690` (test hooks), `d4cdf70` | `gap_recursion_backend` (`scaled_reference`) |
| interned dosage census, with its census work counters | I2-J `90dc305`, `705504f` | `dosage_census_backend` (`reference`) |
| union-to-adaptive emission handoff: **present, left off** | I2-K `69a8cf1`, `394fe0a` | `adaptive_emission_handoff` (0) |
| prepared-locus reuse, between-sample release, batch resource ledger | I3-A `22ea56c`, I3-A.1 `57787d5`, I3-B `945775c` | C++ API and driver (section 3) |
| persistent prepared-locus index v1 and panel section 16 | I3-C `01dc34b`, I3-D `bf226d8` | C++ API and driver (section 3) |
| factorised fine evidence views | I4-B `62e1d7f` | `indexed_union_evidence_layout` (`reconstructed`) |
| seed-cache accounting fix, typed tie comparator | I4-C `71ada7f`, `ca6db98` | correctness fixes (no key) |

**Kept but parked.** The tie-seed cache predates the I-series. It stays at its default, `repair_tie_seed_cache_bytes 0`.

**Removed from the snapshot (`392fc0b`):**
- the I2-G sampled verification profile (`verification_profile_period`);
- the I2-I verification-request census (`verification_census_period`, `verification_census.{hpp,cpp}`);
- the S0 `add_haplotypes_with_censuses` member.

Each was inert when off.

**Not brought in:**
- the I4-D census and I4-D.1, the I4-E fine-origin aggregates, the I4-F sidecar, the I4-G page release;
- I5 (pruning, dumps, prototypes);
- the `s0-v2` pattern-mass cache;
- every research runner, benchmark and design-measurement tool.

They remain on their branches.

**Already in the production caller, untouched:** the diagnostics that predate the I-series — emission census, union verification census, piece-hit census, union work trace and seeding profile.

**Two release paths: only one is in this baseline.**
- **In:** I3-A.1's between-sample release (`release_between_samples`). It runs only when the driver is called with `CLEANUP=1`.
- **Not in:** the I4-G boundary page-release experiment.

## 2. Development configuration

The configuration is `profiles/genotype_dev_performance_v1.tsv`. It sets every accepted key explicitly, in four sections:

| section | keys |
|---|---|
| performance backend choices | the eight backend keys above, plus `marker_proposal_reuse_table_bytes` |
| cache and workspace budgets | indexed union store 2 GiB, fine store 512 MiB, and the six caches that already exist in production (C4 v6 campaign values) |
| resource capacities and guards | the union and fine consumption caps |
| reporting | `timing_detail 1` |

**It contains no scientific or model setting.** `compose_profile.py` enforces this, and rejects any key outside these sections.

**What it is not:** a public default, or a universally validated configuration. Its cache and capacity values come from C4. They were also exercised on the CYP2D6 and ACOT development cells, but never tuned or validated beyond them.

**Caches and capacities can change whether a run completes, never what it computes.**
- They change what is executed and stored.
- Caps count executed work, so a cap, time limit or footprint limit can stop one configuration and not another.
- A completed result is unchanged.

**Composing one run's profile.** A run records a single profile path and SHA-256, so each comparison composes one file:

```
python3 tests/analysis/p1/compose_profile.py OUT.tsv NAME BASE_PROFILE full \
    "footprint_limit_bytes=4294967296:4 GiB local test envelope" "time_limit_seconds=280:declared timeout"
panvar genotype --preset short-read --profile OUT.tsv --threads 4 --graph G.gfa --reference-path REF \
    --reads1 R1.fq --reads2 R2.fq --library LIB.txt --dosage-mode repair --out PREFIX
```

**Arguments:**
- `BASE_PROFILE` is the analysis's own profile: its scientific settings and its own guards (state-class and neighbour caps).
- An override may only replace a key the base profile already sets, and must declare a reason. Duplicate keys are refused.
- `reference` in place of `full` keeps only the production caches and consumption caps, with every scoring backend at its default. This is the reference arm used in validation.

The 4 GiB footprint limit is the local test envelope. It is an override chosen per comparison, never part of the configuration.

## 3. Prepared-locus index: where each part is available

| functionality | normal CLI (`panvar genotype`) | C++ API (`panvarlib`) | development driver |
|---|---|---|---|
| prepare a locus once, analyse many samples in one process | **no** | `prepare_genotype_locus(args \| options)`, `run_genotype_command_prepared(args, locus)`, `run_genotype_prepared(options, locus)` | `genotype_prepared_locus_driver prepared …` |
| write a persistent index (v1, optionally with panel section 16) | **no** | `prepare_genotype_locus_record(…, persist_panel)`, `encode_prepared_locus`, `write_prepared_locus` (never overwrites; atomic publication) | `… write` / `… write-panel CEILING 0 OUT_PREFIX INDEX ARGS_FILE` |
| load an index (bounded, checksum- and digest-verified) and analyse samples against it | **no** | `read_prepared_locus_record`, `load_prepared_locus(path, PreparedLocusLoadLimits, &report)` | `… loaded CEILING CLEANUP OUT_PREFIX INDEX ARGS_FILE [ARGS_FILE …]` |
| release between samples; batch resource verdict | **no** | `release_between_samples()`, `BatchResourceLedger` | `CLEANUP=1`; `OUT_PREFIX.batch_resource.tsv` |

**The driver** is `tests/genotype_prepared_locus_driver.cpp`. It is built with `PANVAR_BUILD_TESTING` and never installed. Each `ARGS_FILE` holds one `panvar genotype` argument list, one argument per line.

**There is no user-facing CLI integration.**

**Compatibility.** An index or prepared object serves a sample only if the sample's graph content, reference, bubbles, `min_variant_bp`, panel curation (manifest content), k-mer and syncmer settings all match.
- Any mismatch is refused with a typed error **before any output** is written.
- The prepared object holds no read-, library- or scorer-dependent state. Library, scoring options and backends are applied per sample: the e2e gate checks a second declared library and `--library auto` against cold runs.
- Malformed or corrupted files are refused. The SHA-256 covers header and payload, and invariants are checked behind valid checksums.
- Unsupported panel-section versions are rebuilt, and the reason is recorded.

**Format.** v1 with section 16 is unchanged. An index written by the I3-D build (`bf226d8`) loads and analyses correctly on this baseline.

## 4. Reproducing the validation

```
python3 tests/analysis/p1/p1_run.py next      # repeat; one process at a time, budget-checked
```

- The manifest (`MANIFEST`, plus `AMENDMENT_1`) and the pairs are in `p1_run.py`.
- `p1_compare.py` prints the complete strict `perf_compare` problem list, plus the qualified verdict of the reviewed `i4c_compare --profile i4b` chain.
- Results: `tests/results/p1_consolidation_20261017/README.md` (gitignored).

## 5. Limitations

**ACOT 30× does not complete on this baseline under the local caps.**
- Without fine-origin aggregates (I4-B views, the route here), it stops at the declared **500 M physical-origin cap** during the union pass (RESCORE_INCOMPLETE; I4-B, 70.4 s, process peak 4.14 GB).
- The **≈4.37 GB footprint refusal in the repair phase** came from the separate, experimental I4-E fine-origin-aggregate route, which is not in this baseline.
- These are different runs, and the baseline does not complete the union pass under those conditions. ACOT 30× was not rerun here.
- Caps are unchanged. The allocator investigation is closed and not reopened.

**Resource governance:**
- Caps count executed work.
- 4 GiB is the local testing envelope, not a production requirement.
- Footprint limits are detection plus refusal. There is no envelope policy for cluster nodes.

**Executed-work counters vary between runs at more than one thread.** The indexed union's `indexed_union_executed_prefilters` and `…_executed_dps` change between identical runs, because workers race to fill shared caches. Science is unaffected.

**No CLI for the prepared-locus index** (section 3).

**The configuration is derived from C4.** It is not validated beyond the development cells.

## 6. Issues carried to S1 (recorded, not implemented)

1. **Explored projection coverage is not completeness.** The block labellings seen among explored candidates do not prove that a sequence's projection is unique, so they do not give unconditional error bounds. A reported bound must be qualified as conditional on the explored candidates, or come from a complete enumeration.
2. **Representation intervals must include coupled context.** Where alternative labellings of one sequence are coupled across a context block, the interval must include that intervening context. The proposal's maximal runs of adjacent ambiguous blocks are insufficient. An example is CYP2D6 sv:3..sv:5 across ctx:4, where one 59 bp unit can be counted on either side.
3. **Legacy-output identity needs explicit treatment of new settings and provenance fields.** A new `output_schema` key, and any provenance rows, appear in the settings and status tables. "Byte-identical legacy output" must say whether an unset key writes no row, whether a set `legacy` value writes one, and how the comparison treats them.
