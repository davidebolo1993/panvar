# M7-B result: provenance verified; the Locityper comparison stopped at calibration (Ns in the frozen background); Panvar-side rankings only

**Run.** Branch `m7-b-fixed-pool`:
- freeze `e614f13`;
- Amendment 1: λ check at printed precision;
- calibration stop record and Amendment 2: truth-distance parsing.

Results are in `tests/results/m7b_fixed_pool_20261020/`. M7-A's specification and manifest are unchanged.

**Budget.** 66.0 of 900 s, 4 threads, one process at a time. Three recorded failures:
- `background`: a check-precision defect;
- `preproc.c4_lzo_hg00171_5`: **the calibration failure**;
- `truthdist`: a parsing defect.

**Scope.** No downloads, production changes, scoring code, merge, push or cleanup.

## 1. Simulator provenance: all lineages verified byte for byte (`provenance`, 7.6 s)

| lineage | evidence | verdict |
|---|---|---|
| C4 (C4 HG00171 LZO 5×, C4 NA19650 LOO 10×) | Harness `tests/genotype_coverage_c4.sh` at `abb35a0`: `simulate --gfa c4.gfa --reference HG00096#1… --donor <d> --depth 30 --read-length 100 --seed 20260916`. The graph SHA-256 `b13af390…` equals the production `c4.gfa`, and the package's source reads equal the production source reads. The recorded probe SHA (`35bbd627…`) exists nowhere on the machine. The **first declared candidate build** (`s1-d-panel-references/build`) regenerated **both** donors' 30× source reads byte for byte. Each cell's recorded subsample (the `mosaic-11fe681` build) reproduced the cell's reads byte for byte | **PASS** (behavioural; the original binary is unavailable) |
| CYP (CYP HG02769 LZO 5×, CYP NA19650 LOO 5×) | The recorded `simulate.log` commands, with the pilot's own `bin/genotype_probe` (SHA-256 `e9bdf722…`), `cyp2d6.gfa`, `grch38#1#chr22:42031864-42245566`, seed 20260921 and the donor as salt, regenerated cov30 byte for byte. `subsample.5.log` reproduced cov5 byte for byte | **PASS** |

`executables.tsv` records SHA-256 and version for locityper (v1.7.2), strobealign, samtools, jellyfish 2.3.1, the simulator
builds and the scoring build. **The installed Locityper binary's source identity (checkout `7a6f2993`) remains not
established.**

## 2. Background and calibration: STOPPED

**The background construction passed its checks** (`background.a1`):
- GRCh38 ACOT (252,518 bp) and GSTM1 (240,715 bp) spellings match their region lengths, with no non-ACGT bases;
- one contig of 494,233 bp with a 1,000-N spacer;
- a homozygous diploid GFA: exactly 2 identical donor paths plus a reference path;
- four simulations with the verified lineages, each cell's depth and the frozen salts. λ matches (depth/2)/(2·100)/(1 − 0.001)
  at printed precision, and template z-scores lie between −0.73 and 2.07.

**The calibration failed.** `locityper preproc` refused it:

    Runtime error: Cannot count k-mers for sequence with Ns

The N spacer lies **inside** the declared background region (`-b m7bg:1-494233`). The M7-B design wrongly assumed that
Locityper would exclude N windows. Locityper counts k-mers over the whole background region and rejects any N.

**Per the stop rule,** the background was not changed:
- no spacer removal and no region restriction;
- every Locityper `preproc`, `target` and `genotype` unit is **NOT RUN**;
- the three other `preproc` units would fail identically and were not attempted;
- no Locityper ranking, calibration estimate, GC-support analysis, seed check or pool-membership check exists.

**Remedies for review** (none applied; each changes the frozen calibration):
- **(a) Concatenate the two spellings with no spacer:** one artificial ACOT|GSTM1 junction inside the background.
  Locityper's own unique-k-mer window filter would see it as ordinary sequence.
- **(b) Keep the reference but restrict `--bg-region` to one N-free region,** e.g. ACOT, 1–252,518. This halves the
  background.
- **(c) Two separate backgrounds:** not supported, since `preproc` takes one background region.

Whichever is chosen, ~0.25–0.49 Mb may still prove too little for Locityper's per-GC depth fit. The existing gates
(usable windows; GC bins with fewer than 100 windows) would decide that, before any ranking.

## 3. Canonical pools (`pools`; all checks pass; frozen biological pool unchanged)

| cell | frozen candidates → canonical diplotypes / haplotypes | checks |
|---|---|---|
| C4 HG00171 LZO | D 11 → 11 / 11 | orientation (worst forward − reverse containment 0.952); no N |
| CYP HG02769 LZO | D 13 → 13 / 14 | 0.933 |
| C4 NA19650 LOO | D 11 → 11 / 18; D+O 34 → **29** / 41 (sequence-identical constructions merged) | M1-D `reported.1/.2` = spelled reported haplotypes; 0.954 |
| CYP NA19650 LOO | D 13 → 13 / 15; D+X 14 → 14 / 17 | 0.937 |

`mapping.tsv` maps every frozen candidate (labels, tier, sources, block projections, target alleles, unit structure,
`also_truth_selected`) to its canonical sequence pair. The pool FASTAs and the `--priors` files exist for D, D+X and D+O,
but were never handed to Locityper.

## 4. Panvar side (`score-set` on D; gate M2 passes)

**M2.** Every saved expectation was reproduced: 11/11, 13/13, 11/11 and 13/13. The X candidate was scored separately. The
O `score-set` was not run, because it serves only the stopped D+O diagnostic.

Block excess is the target-block edit distance to truth (best homologue assignment) minus the availability floor, and is
evaluation only. Every M6-A row implies the same floor per cell: 0, 0, 36 and 24.

| cell | Panvar F + E top (class, excess) | Panvar total top | total copy number | allocation | correct control |
|---|---|---|---|---|---|
| C4 HG00171 LZO | the caller's 2+2 (excess 0; tied with its target-swap relabelling) | caller 2+2 | true total ranks 1; the best wrong total (2+3) is 216.9 nats behind | best wrong allocation (1+3) is 852.9 behind | **not deteriorated** (rank 1 under both) |
| CYP HG02769 LZO | caller 2+2 (excess 0) | caller 2+2 | the best wrong total is 2,244.4 behind | best wrong allocation is 11,417 behind | **not deteriorated** |
| C4 NA19650 LOO | reported 3+1 (excess 65,444) | reported 3+1 | true total 4 ranks 1 (the best wrong total is 374.6 behind) | **error:** the best true-class 2+2 (caller-saved, excess 24) ranks 5th, 252.5 behind (total: 271.4). The floor-reaching 2+2 (excess 0) ranks 7th, 407.0 behind (total: 425.9) | — |
| CYP NA19650 LOO | **caller-saved 2+3 (excess 17)**, 3.75 ahead of the reported 2+2 | **reported 2+2** (excess 12,159), 3.64 ahead of the 2+3 | under F + E the true total 5 ranks 1 by 3.75; under the total it ranks 2nd by 3.64 | the true class leads all wrong allocations by > 2,900 | — |

**Within true classes,** Panvar orders by sequence only loosely:
- **CYP LOO, F + E, by rank:** 1 (excess 17), 3 (12), 4 (0, also truth-selected), 5 (13), 6 (12), 7 (0).
- **C4 LOO, by rank:** 5 (24), 7 (0), 8 (6).

These reproduce M1-B and M6-A:
- The C4 failure is content-driven allocation misranking: both F + E and total prefer 3+1.
- The CYP failure is a dosage decision of a few nats that the copying prior flips: F + E prefers the true total; the total
  does not.

## 5. What was planned but not run, and why

| item | status |
|---|---|
| Locityper D-only controls, failure cells, D+X / D+O diagnostics | NOT RUN: calibration failed |
| Calibration estimates (insert, error, NB per GC), usable windows, target GC support | NOT RUN |
| Seed stability (s1r, s2) and pool-membership sensitivity | NOT RUN |
| Read-pair counts through recruitment, mapping and scoring | NOT RUN |
| Panvar O `score-set` | not run (serves only the D+O diagnostic; M1-D/M1-E scores exist) |

**Conclusion.**
- **No comparison of the two objectives exists yet.** Nothing here supports, or argues against, the alignment +
  local-depth objective.
- **What is established:**
  - all generator lineages are verified;
  - the background construction and its simulations pass every check except Locityper's no-N requirement;
  - the canonical pools and their provenance mapping are verified;
  - Panvar's rankings of the frozen pools reproduce the saved records.
- **The one decision needed:** which N-free background remedy (§2: a, b or neither) to authorise.
- **Execution budget:** 834 s of the 900 s remain unused. A resumed run would need a new review and a new budget.
