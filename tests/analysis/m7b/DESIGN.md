# M7-B: corrected fixed-pool comparison, Panvar vs Locityper (frozen before any unit)

**Ancestry.** Branch `m7-b-fixed-pool` from `7ba52d4` (M7-A). M7-A's `SPEC.md`, `CANDIDATES.tsv` and `BINDINGS.tsv` are
**unchanged**. Every change below is prospective and implements Sol's corrections.

**Bounds.**
- 900 s counting every attempt; 4 threads; one experiment process at a time; 6 GiB envelope.
- Admission: spent + the unit's declared timeout + 40 s (reserve for compare and reporting) ≤ 900.
- **Stop rule:** any validity or calibration failure stops the affected comparison. Settings are never changed in
  response.
- Not done: no downloads, production changes, new scoring code, merge, push or cleanup.

## 1. Simulator provenance (unit `provenance`; must pass before calibration)

**C4 lineage** (cells C4 NA19650 LOO 10×, C4 HG00171 LZO 5×).
- **Recorded:** `genotype_coverage_c4_held_out_abb35a0/<donor>/source/{simulation.tsv, library.txt}`, harness
  `tests/genotype_coverage_c4.sh` at commit `abb35a0`. Command: `simulate --gfa c4.gfa --reference
  HG00096#1#haplotype1-0000024:31848049-32074186 --donor <d> --depth 30 --read-length 100 --seed 20260916` (no salt).
  Graph SHA-256 `b13af390…` equals the production `c4.gfa`. Recorded probe SHA-256 `35bbd627…`.
- **No binary with the recorded probe hash exists on the machine.** The lineage is verified **behaviourally**:
  - candidate builds, in this declared order: (1) `panvar-worktrees/s1-d-panel-references/build/genotype_probe`,
    (2) `panvar-worktrees/mosaic-11fe681/build/genotype_probe`;
  - the first that regenerates **both** donors' archived source reads byte for byte (the package's `.fq.gz`, which equal
    the production `source/` reads) is the C4 simulator.
- Then each cell's own subsample (recorded in `reads/<d>/cov<k>/subsample.log`) is re-run with its recorded binary and
  must reproduce the cell's reads byte for byte.

**CYP lineage** (cells CYP NA19650 LOO 5×, CYP HG02769 LZO 5×).
- **Recorded:** `cyp2d6_shortread_pilot_v1/data/reads/<d>/simulate.log`. The pilot's own `bin/genotype_probe` (SHA-256
  `e9bdf722…`), `cyp2d6.gfa`, reference `grch38#1#chr22:42031864-42245566`, depth 30, seed 20260921, salt = donor.
- Regeneration of cov30, and of cov5 via `subsample.5.log`, must be byte-identical.

**If a lineage is not reproduced,** its cells stop before any ranking is read (M5 is not waived). All regeneration
writes only to new directories under the M7-B results.

## 2. Background (unit `background`)

- **Sequence, frozen:** one contig `m7bg` = GRCh38 ACOT path (`GRCh38#0#chr14:73416151-73668669`, from
  `genotype_pilot_acot_2x_d24e471/panel.sorted.gfa`) + 1,000 N + GRCh38 GSTM1 path
  (`GRCh38#0#chr1:109561616-109802331`, from `genotype_pilot_gstm1_2x_d24e471/panel.sorted.gfa`).
  - Paths are spelled with minus steps reverse-complemented. Each spelling's length must equal its region length.
  - Non-ACGT bases are allowed only in the spacer.
  - The simulator skips fragments containing N; Locityper excludes N windows.
- **Diploid representation:** a one-segment GFA with three paths, `M7REF#0#m7bg` (reference) and `M7BG#1#m7bg` /
  `M7BG#2#m7bg` (donor), all spelling the same bytes: a homozygous diploid. Ploidy (exactly two donor paths), identical
  spelling and the reference path are checked before simulation.
- **Simulation:**
  - each cell's background uses that cell's verified simulator, at the cell's depth (C4 LOO 10×, the others 5×);
  - read length 100, insert 350 ± 40, 4 sigmas, error 0.002, outlier 0.001;
  - **seed rule:** `--seed <the cell's original seed> --seed-salt m7b-bg-<cell>` (C4 20260916, CYP 20260921);
  - the expected fragment intensity, (depth/2) / (2 · 100) / (1 − 0.001), is checked against `simulation.tsv`.
- **No substitution with random sequence, no enlargement, and no calibration change after rankings are seen.**

## 3. Calibration (unit `calibrate.<locus>` and `preproc.<cell>`)

- **Reference per locus:**
  - C4: `m7bg` + `c4ref`, the spelled `HG00096#1#haplotype1-0000024:31848049-32074186` from `c4.gfa`.
  - CYP: `m7bg` + `cypref`, the spelled `grch38#1#chr22:42031864-42245566` from `cyp2d6.gfa` (the CYP model's GRCh38
    reference choice is preserved).
- `samtools faidx`; `jellyfish count -C -m 25 -s 5M -t 4 --out-counter-len 2 --lower-count 2` (jellyfish 2.3.1, as in
  the October 15 smoke).
- `locityper preproc -i bg_1.fq bg_2.fq -r ref.fa -j ref.jf -b m7bg:1-<len> --seed 20261020 -@ 4 --debug`, run per
  cell.
- **Recorded:**
  - usable windows;
  - per-GC window counts and the fitted NB per GC bin;
  - window and neighbourhood sizes;
  - insert-size and error-profile estimates;
  - the GC support of the **target** windows (from each genotype run's `windows.bed.gz`, weight ≥ 0.001). Target windows
    whose GC falls in a bin with fewer than 100 background windows (Locityper's `min_tail_obs`) are flagged as relying
    on extrapolated tails.
- **Limitation carried throughout:** the reduced reference (≈ 0.49 Mb + locus) cannot establish genome-wide k-mer
  uniqueness.

## 4. Canonical pools (unit `pools`)

- **Spelling:** D and X candidates are spelled with `genotype_probe spell`, in the model frame of each cell's reference.
  O candidates are copied from the M1-D / M1-E FASTAs, which are also in the model frame.
- **Checks:**
  - M1-D's `reported.1/.2` must equal the spelled reported haplotypes;
  - every haplotype's forward 31-mer containment in the locus reference must exceed its reverse-complement containment
    (consistent orientation);
  - no haplotype contains N.
- **Canonicalisation:** by complete sequence (SHA-256) into canonical haplotypes `<cell>.hNN`, and by unordered pair into
  canonical diplotypes. `mapping.tsv` maps every frozen manifest candidate to its canonical sequence pair, keeping all
  labels, tiers, block projections and target-block alleles.
- **Pools:** D, D+X (CYP LOO) and D+O (C4 LOO). Each pool has a FASTA of its canonical haplotypes, a `--priors` file
  (locus, `hA,hB`, log₁₀ prior 0) over its canonical diplotypes, and the Panvar candidate files. **The frozen biological
  pool is unchanged.**
- **Evaluation-only:** target-block edit distance of each D candidate to the truth block (M6-A `blockseq` truth FASTA,
  best homologue assignment) with `genotype_probe edit-distances`. The availability floors are those recorded by
  M6-A / M1-B.

## 5. Executables (recorded in `executables.tsv` by `provenance`)

- SHA-256 and version of: `locityper` (conda v1.7.2; source checkout `7a6f2993`; **binary-to-source identity not
  established**), `strobealign`, `samtools`, jellyfish 2.3.1, both simulator builds, and the scoring build.

## 6. Execution order (frozen) and seeds

| step | units |
|---|---|
| 1 | provenance; background; calibrate.C4 / .CYP; preproc for each cell; pools |
| 2 | Panvar `score-set` on D for each cell (+ X for CYP LOO; + O explicit for C4 LOO, diagnostic) |
| 3 | **Locityper D-only, controls first:** `loc.c4_lzo_hg00171_5.D.s1`, `loc.cyp_lzo_hg02769_5.D.s1` |
| 4 | **Locityper D-only, failure cells:** `loc.c4_loo_na19650_10.D.s1`, `loc.cyp_loo_na19650_5.D.s1` |
| 5 | diagnostics (never support a deployable claim): `loc.cyp_loo_na19650_5.DX.s1`, `loc.c4_loo_na19650_10.DO.s1` |
| 6 | stability: `loc.c4_lzo_hg00171_5.D.s1r` (same seed repeated), `loc.c4_loo_na19650_10.D.s2`, `loc.cyp_loo_na19650_5.D.s2` |
| 7 | compare |

- **Seeds:** s1 = 20261020, s2 = 20261021.
- **Locityper flags:** `target -l <LOCUS> <ref contig>:1-<len> pool.fa`; `genotype --priors priors.tsv --filt-diff 1e9
  --prob-thresh -1e9 --dont-skip -S greedy:i=5k,a=1 -S anneal:i=<pool size>,a=20 --basis none --debug 1 -s <seed> -@ 4`.
- **Per run, verified:**
  - the database's `haplotypes.fa.gz` equals the pool exactly (after Locityper's own deduplication and N filter, with
    nothing discarded);
  - the `res.json.gz` options equal the priors set;
  - every canonical genotype has 20 annealing (stage 2) rows in `sol.csv.br`;
  - read pairs at recruitment, mapping and scoring are taken from the log, and `total_reads` from `res.json.gz`.

## 7. Interpretation (`compare.py`, frozen)

- Rankings and gaps are compared **within each method only.** Converting log₁₀ to nats does not make absolute values
  commensurate.
- **Reported separately:**
  - Panvar fragment + exposure ranking;
  - Panvar total ranking;
  - Locityper `lik_mean` ranking (s1);
  - for each ranking's top D diplotype: unit total, allocation and target-block truth distance;
  - correct-control deterioration (the caller's call, or a sequence-identical relabelling, not ranked first);
  - pool-membership sensitivity (D vs D+X / D+O, per D genotype);
  - seed sensitivity (s1, s1r, s2: means, differences, rank changes).
- **A winning margin smaller than the s1 / s2 (or s1 / s1r) difference of either genotype involved is labelled
  INCONCLUSIVE.**
- A favourable Locityper result motivates attribution only. It does not show that local depth alone caused it, nor
  improved general genotyping.

## Amendment 1 (declared after `background` exited 1; 8.2 s spent)

**What failed.** The fragment-intensity check compared `simulation.tsv`'s λ, which is printed to 6 significant digits, at
a relative tolerance of 1e-6. Spelling, ploidy, template z-scores and simulation exits all passed.

**Fix.** λ is compared at the printed precision (relative 1e-5). The unit is re-run as `background.a1`, into a new
directory, with the same deterministic seeds. Downstream units read `background.a1`. Nothing else changes.

## Calibration stop (recorded after `preproc.c4_lzo_hg00171_5` exited 1; 9.6 s spent)

**What failed.** Locityper `preproc` refused the frozen background with "Cannot count k-mers for sequence with Ns": the
1,000-N spacer lies inside the declared background region. The design wrongly assumed that N windows would be excluded.

**Decision, per the stop rule.** No change to the background construction (spacer removal, restricting `--bg-region` to
one region, or any other) is made without review. **All Locityper preprocessing, `target` and `genotype` units are
stopped and recorded NOT RUN.** The remaining three `preproc` units would fail identically and were not attempted.

**What still runs.** Calibration-independent units on the frozen pools, which change no setting: `pools`, Panvar
`score-set` (D for each cell; X for CYP LOO), `truthdist`, and a compare that reports Locityper as not available. The O
`score-set` is not run: it serves only the D+O Locityper diagnostic.

## Amendment 2 (declared after `truthdist` exited 1; 63.4 s spent)

**What failed.** `m6a.tsv` has two sections with different headers. The block-excess section was not parsed, so no floor
could be implied.

**Fix.** The section whose header carries `block_excess` is read. The units are re-run as `truthdist.a2` and
`compare.a2`.
