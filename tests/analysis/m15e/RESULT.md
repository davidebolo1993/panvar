# M15-E result: public panel-index build and load flags; cold and indexed runs are scientifically identical at C4, CYP2D6 and ACOT; every refusal is named and leaves the prefix untouched

**Run.**
- **Branch and commits:** `m15-e-panel-index`, from `bbc9052`. Freeze `09ec00c`; amendment A1 `4b214a9` (test fixture
  only); result after.
- **Results:** `tests/results/m15e_panel_index_20261021/`.
- **Budget:** 1,033.8 of 2,400 s, including 260 s of disclosed development.
- **Units:** 13 of 14 passed first time. The 14th (`ex.refusals_c4`, a fixture-naming error in the test, not the
  software) passed on its declared A1 rerun.
- **Resources:** 4 threads, one process at a time, peak 2.7 GiB.
- **Not done:** no merge, push, tag, cleanup, reserved rerun or graph regeneration. Existing indexes, graphs, WIP and stash
  are untouched.

**Binaries** (source at `09ec00c`; A1 changed only `tests/analysis/m15e/examples.py`):

| binary | SHA-256 |
|---|---|
| `panvar` | `8d6bfde21e0587158784c6eb9aafc738258dce1510af6da8cf74ad64298cb5f7` |
| `genotype_command_e2e_test` | `2d6189687a56814dc6da68b44b2ae627b0dcd66ecf514339b91eaa7ce83acba9` |

## 1. Interface (implemented)

```
panvar genotype --graph locus.gfa [--reference-path R] [--panel-curation C] [--bubbles B] --build-panel-index locus.pvpl
panvar genotype --graph locus.gfa [same locus inputs] --panel-index locus.pvpl --reads1 S.R1.fq.gz --reads2 S.R2.fq.gz --threads 4 --out S
```

**Build:**
- **Reuses:** the existing I3-C/I3-D writer and its atomic publication; the index carries panel section 16.
- **Needs:** no reads, library or output prefix. An existing file is never overwritten.
- **Refused:** sample arguments and a simultaneous `--panel-index`.
- **Writes:** a build record on stdout (versions, identity, model fingerprint, size, SHA-256, prepare and write seconds).
  Nothing else is left: the staged sorted graph next to the index is removed.

**Load:**
- **Verified first, before the output prefix is touched:** the index's integrity, versions and reconstruction digests;
  then the sample's graph content, reference, bubbles, `min_variant_bp`, curation, k and s against its identity.
- **On failure:** a named refusal (exit 2), never an in-memory fallback.
- **Still read:** `--graph` is still read and hashed.
- **Status:** records `prepared_locus_origin loaded`, `prepared_locus_panel_source` and the `panel_index*` facts. A
  rebuilt panel would be reported, never described as loaded; it did not occur here.

**Documentation:** basic help (usage lines, both flags, what forces a rebuild); guide §8a (build once and genotype many,
what is verified, what is never stored, when to rebuild); module page.

## 2. Results

| unit | result |
|---|---|
| `cli_contract` | pass: build + load refused; eight sample arguments plus `--dosage-check` refused in build mode; a load accepts every sample argument; help needles present |
| `e2e.prepared_locus` (tiny locus; includes the existing I3-A/C/D gates) | pass (27.4 s). **`test_panel_index_cli`:** build without reads writes only the index; an existing destination refused byte-identical; build + load and four sample arguments refused, touching nothing; **cold == indexed** for a declared library, the automatic library, `--dosage-check` (dosage status equal), a second sample and 2 threads; curation, k, `min_variant_bp`, reference and graph mismatches and missing, truncated and corrupted indexes refused, each leaving three sentinels untouched and writing nothing; the index byte-unchanged after every use |
| `ex.build_{c4,cyp2d6,acot}` | 5/5 each: built without reads; only the index and its record; versions, graph SHA-256, size and SHA-256 verified; a second build refused (`exists`), byte-identical |
| `ex.index_{c4,cyp2d6,acot}` (automatic library) | **18/18 each**: equal exits (0); byte-identical science files, FASTA, locus table and library estimate (with the declared name and timing allowances); every effective setting and the prediction key identical; model **and** panel loaded |
| `ex.dosage_check_c4` (declared library, `--dosage-check`) | **21/21**: as above, plus the dosage tables and status identical (timing aside) |
| `ex.two_samples_c4` (NA19650 against the NA19705-graph index) | **19/19**: cold == indexed; the index byte-identical to its build record after both samples |
| `ex.refusals_c4` → `a1.ex.refusals_c4` | first run 8/9 (fixture error, §4); **A1 rerun 9/9**. Curation (`panel curation`), `--kmer-size 27` (`kmer_size`), another reference path (`reference path`), the CYP2D6 graph (`graph content`), a missing file (`io`), half the file (`truncated`) and a flipped byte (`checksum_mismatch`) were all refused with exit 2, four sentinels untouched and nothing written; build + load and a build with reads were refused, building nothing |
| `ctest.full` (final revision) | **63 of 63 passed** (535.2 s); disabled, reported separately: `panvar_real_regressions`, `panvar_rebuild_stats` |

**Sample-independent options never invalidated the index:** reads (two samples), library (`auto` and `FILE`), threads (1,
2, 4) and `--dosage-check`.

## 3. Build, load and runtime (observed; not a promised speedup)

| cell | index size | build (prepare + write) | load | `seconds_total` cold → indexed | wall cold → indexed | peak footprint cold → indexed |
|---|---|---|---|---|---|---|
| C4 NA19705 (auto) | 44.5 MB | 1.79 + 0.20 s | 0.47 s | 12.69 → 11.10 s | 12.7 → 11.6 s | 0.98 → 0.94 GB |
| CYP2D6 HG02769 (auto) | 39.2 MB | 1.99 + 0.18 s | 0.43 s | 22.86 → 21.05 s | 22.9 → 21.6 s | 1.22 → 1.06 GB |
| ACOT HG03209 (auto) | 160.2 MB | 6.18 + 0.77 s | 1.71 s | 34.40 → 30.09 s | 34.5 → 32.0 s | 3.23 → 1.88 GB |
| C4, declared + `--dosage-check` | | | 0.48 s | 11.50 → 9.42 s | 12.3 → 10.6 s | 0.94 → 0.70 GB |
| C4 NA19650 (auto) | | | 0.47 s | 9.92 → 8.44 s | 10.0 → 9.0 s | 0.96 → 0.94 GB |

**How to read this:**
- **Savings:** the index saved 1.0–2.5 s of wall time per run here (about 6–14%; 1.5–4.3 s of `seconds_total`), consistent
  with I3-D. Reading and scoring the
  sample is the rest of the run and is not saved.
- **ACOT's peak footprint fell most** (3.23 → 1.88 GB): the in-memory preparation's peak is avoided.
- **Machine load:** absolute times in this session are about half of M15-D's on the same cells (the full suite also ran in
  535 s against 1,006 s). That shift is outside the batch's control. Cold and indexed were measured back to back.

## 4. Amendment A1

**What happened:** `ex.refusals_c4` wrote its curation manifest as `curation.txt` inside its refusal prefix `curation`,
so the "nothing written" count included the manifest.
- **The refusal itself was correct** in the original run: exit 2, named, sentinels intact.
- **The fix:** the manifest was renamed and the unit rerun once (2.4 s, 9/9). The failed unit stays recorded.

## 5. Remaining limitations

- **Loading still reads and hashes `--graph`** and the other identity inputs. This is not an index-only invocation, and the
  load itself is about 70% SHA-256 (I3-D).
- **Sequential reuse only, one sample per process** from the CLI. In-process multi-sample reuse stays API-only, and there
  is no concurrency.
- **Load limits are the existing ones:** file ≤ 1 GiB and accounted decoding ≤ 3.5 GiB. A much larger locus would refuse
  with `limit_exceeded`.
- **No partial-build recovery or index management** (listing, migration): an incompatible index is rebuilt by the user.
- **The rebuilt-panel route** (an index without a usable section 16) is reported but was not exercised here, because
  indexes built by this command always carry section 16. The I3-D gates in `e2e.prepared_locus` cover it at the API level.

**Stopped for review.**
