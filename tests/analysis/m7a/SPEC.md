# M7-A: real-failure comparison specification (Panvar fragment/exposure objective vs Locityper's alignment + local-depth objective; fixed candidates; not executed)

**Ancestry.** Branch `m7-a-real-failure-spec` from `767c7ec` (M6-G). Commits:
- `f272ede`: M6 closure appended to M6-G; manifest builder;
- Amendment 1: units-table comment line;
- Amendment 2: reachability split.

**Diagnostic execution.** 1.4 of 300 s, all read-only over saved records. The manifest builder ran 3 times; the first
attempt failed and is recorded. Results are in `tests/results/m7a_real_failure_spec_20261020/`.

**Not done:** no scoring, genotyping, download, production change, flag, merge or push.

**Question.** Could an alignment-plus-local-depth objective distinguish the existing real structural alternatives better
than Panvar's current fragment/exposure objective?

**What this answers.** A fixed-pool ranking comparison. **It does not establish truth-free candidate discovery.** The four
cells are inspected development cases, not an independent validation set.

## 1. Development cases and bindings (`BINDINGS.tsv`: 54 files, SHA-256)

| cell | model | reads (library) | saved prediction | known outcome (evaluation only) |
|---|---|---|---|---|
| C4 NA19650 LOO 10× | `genotype_c4_production_11fe681/graphs/NA19650.loo.gfa` | campaign `reads/NA19650/cov10` | M1-B fresh run `call.c4_10` (equivalence set, mosaic scores/candidates/status) | allocation error: reported 3+1, truth 2+2 |
| CYP2D6 NA19650 LOO 5× | `cyp2d6_shortread_pilot_v1/.../loo.NA19650.gfa` (+ curation table) | pilot `reads/NA19650/cov5` | M1-B `call.cyp_5` | dosage error: reported 2+2, truth 2+3 / 3+2 |
| C4 HG00171 LZO 5× (control) | `.../graphs/c4.gfa` | campaign `reads/HG00171/cov5` | M4-B constant arm `run.c4_lzo_hg00171_5.const` | correct 2+2 |
| CYP2D6 HG02769 LZO 5× (control) | `.../graphs/cyp2d6.gfa` (+ curation) | pilot `reads/HG02769/cov5` | M4-B `run.cyp_lzo_hg02769_5.const` | correct 2+2 (sequence-exact) |

The same reads, models and libraries were bound by M6-A; the digests are unchanged.

## 2. Frozen candidate pools (`CANDIDATES.tsv`, SHA-256 prefix `a8cd9516`, 72 rows)

**Tiers** (from `manifest.py`, read-only over saved records):

| tier | contents | C4 LOO | CYP LOO | C4 LZO | CYP LZO |
|---|---|---|---|---|---|
| **D** deployable | the reported call; the caller's own best **saved scored** candidate per target-block unit class (M1-B Amendment 1 rule, truth-free steps); M6-A's frozen truth-free count selection | 11 | 13 | 11 | 13 |
| **X** cross-run | M1-B `cross_30x`: the 30× report of the same donor, i.e. other data | — | 1 | — | — |
| **O** oracle diagnostics | M1-B `sub_truth_a/b`; M1-D truth-block and oracle compositions (8); M1-E patched constructions (a22/e13 × none/U1/T1/U2/T2/all) and its truth block (13) | 23 | — | — | — |

**Columns per candidate:**
- tier, labels and sources (all provenances);
- allele paths, or the explicit FASTA and record names with sequence-hash prefixes;
- target-block alleles;
- per-homologue unit copies (M1-B structure rule, from bytes), unit total and class;
- `caller_scored`: whether the run's own search scored it;
- `target_substitution_of_reported`: whether it differs from the report only at the target block;
- `copying_prior_defined`: yes for panel allele paths; **no** (prior_support none) for O explicit sequences;
- the saved Panvar total, fragment, exposure and prior where one exists.

**Structural coverage of D.**

| cell | unit classes in D |
|---|---|
| C4 LOO | 3+1 (reported) and 1+3; 2+2 (caller-saved and two count-selected catalogue members the caller never scored); 2+1; 3+2; 3+3; 1+1 |
| CYP LOO | 2+2 (reported); 3+2 and 2+3, including the caller-saved 2+3/3+2; 3+1; 2+1; 2+4; 3+3; 3+4; 4+1 |
| controls | the correct 2+2 call plus every other class the caller saved |

**Disclosures.**
- **The CYP LOO count-selected 3+2 / 2+3 pairs coincide with M1-B's truth-closest substitutions.** They stay in D, because a
  truth-free rule reached them, but carry `also_truth_selected = yes`. A favourable result for them must be read with that
  coincidence in mind.
- The CYP LOO X candidate differs from the report at several blocks.
- M1-E's complete true diplotype exists only as FASTA, was never a declared candidate, and is excluded.

## 3. The alternative objective, as implemented in the pinned Locityper

**Source.** `external/locityper` (clean checkout, commit `7a6f2993227c15920df415a19045d5baf98d90e3`, v1.7.2). The
executable is `miniconda3/envs/locityper/bin/locityper`, v1.7.2.

**Binary identity gap.** The executable's source commit is not established. Building from source would need crate
downloads (not permitted). The comparison runs Locityper itself; **no Locityper-like loss is re-implemented.**

| component | implementation (file:line) | defaults and provenance |
|---|---|---|
| **retained paired placements** | Each mate is mapped **independently** with strobealign (`-M = -N = min(25000, 4H) -S 0.5 -f 0.001 -k 15 --eqx`; `command/genotype.rs:962`); samtools keeps `AS ≥ 50` or primary. Per read end, a passable edit distance comes from the error profile's edit-distance distribution (`model/locs.rs:298–315`, `EditDistCache`). For each read pair and contig pair, at most **10** alignments (2 for low-weight reads) within `prob_diff` of the best are kept (`locs.rs:740–797, 876–908`). Mates are reconnected with the insert-size log probability. An **unmapped option** costs 2 · `unmapped_penalty` + the insert penalty (`locs.rs:816–866`) | `unmapped_penalty` = 10⁻¹⁰ (Illumina, `model/mod.rs:55`); `prob_diff` = \|unmapped\| + 1 log₁₀ unit (`genotype.rs:1294`) |
| **alignment / error scoring** | ln P(alignment) = Σ over CIGAR operation counts × ln P(op), from an error profile of operation probabilities plus a beta-binomial edit-distance distribution estimated from background alignments (`bg/err_prof.rs:142–212`). Pair score = ln p(mate 1) + ln p(mate 2) + ln P(insert) (`seq/aln.rs:236`) | estimated per sample at preprocessing |
| **assignment optimisation** | For each genotype, every read pair takes exactly one retained placement or the unmapped option, maximising **L = (1 + s) · Σ_w depth term_w + (1 − s) · Σ_reads ln P(assigned placement)** (`model/assgn.rs:80–81, 235–237`). Greedy, then simulated annealing (`solvers/solve.rs:211–228`). Each attempt re-randomises read-middle "tweaks" and window offsets (`assgn.rs:127–151`). Reported as `lik_mean` / `lik_sd` (log₁₀) over attempts in `res.json.gz` (`solve.rs:732–757`) | s = `lik_skew` = **0.85**, so coefficients 1.85 (depth) and 0.15 (alignment) (`model/mod.rs:121`). The source documents no calibration provenance for 0.85: it is an author default. Stages: greedy (≤ 5,000 genotypes, 1 attempt), then annealing (top 20, 20 attempts) |
| **preliminary filter** | alignment only: genotype score = prior + Σ_reads (best placement on any member haplotype). Genotypes below best − `filt_diff` are dropped (`solve.rs:86–118`) | `filt_diff` = 100 log₁₀ units; `prob_thresh` = −4 log₁₀ |
| **local-depth windows and distributions** | Windows of the background-estimated size along each haplotype, excluding a 200-bp boundary at each end. A read pair counts in the window of its (tweaked) middle. Each window has weight × **ln[P(k \| CN = 1) / (P(k \| CN = 1) + Σ_c P(k \| CN = c))]**, with alternative CN ∈ {0.3, 2, 3, 4, 5} (`math/distr/bayes.rs:27–35`, `model/distr_cache.rs`). P(k \| CN) is a negative binomial for the window's GC bin, ×2 for paired-end, ×CN. Window weight is the product of a linguistic-complexity weight and a locus-specific-k-mer weight; windows below `min_weight` contribute nothing (`model/windows.rs:323–420`) | NB per GC fitted from background windows: LOESS over GC with `frac_windows` 0.5, tail blurring, ploidy 2 (`bg/depth.rs:130–200, 300–340`). Weight calculators (0.5, 4) and (0.2, 4); `min_weight` 0.001; `boundary_size` 200 |
| **calibration requirements** | `preproc` needs a reference FASTA, its jellyfish k-mer counts, and WGS reads (or alignments) over a background region (default `BG_REGIONS`, or `--bg-region`). It estimates the insert distribution, error profile, depth NB per GC bin, and window and neighbourhood size; `--like` copies a similar dataset's calibration. `add` needs the same reference and counts to compute locus-specific k-mer fractions | none of these exist for the cells (§5) |

**Observation-model differences that will remain explicit.**

| aspect | Panvar | Locityper |
|---|---|---|
| fragment treatment | sums every admissible origin of each fragment (marginalised point-process fragment likelihood, exact within its band) | assigns each read pair to one retained placement (maximisation over assignments, bounded alternatives) |
| depth | depth enters through the exposure term | depth enters through a weighted per-window normalised CN posterior |
| errors | substitution error from the declared library, plus an outlier/background mixture | learned operation probabilities, plus an unmapped penalty |
| prior | copying prior | uniform genotype priors (§4) |
| search | its own search produced the pool | sees only the fixed pool |

## 4. Placements and pool-membership control

- **Panvar placements cannot be reused.** Panvar never exports a placement table: it sums origins internally with its own
  banded emission. Locityper needs its own CIGAR-scored, edit-filtered placements. The comparison therefore needs **new
  mapping**, done by Locityper's own pipeline (strobealign on the candidate haplotypes).
- **Genotype set.**
  - A `--priors` file lists exactly the pool's genotypes, at log₁₀ prior 0 (uniform). `generate_genotypes` then uses
    that list as the **complete** genotype set (`genotype.rs:1096–1121`).
  - Pruning is disabled with `--filt-diff 1e9`, `--prob-thresh -1e9` and `--dont-skip`.
  - The scheme `-S greedy:i=5k,a=1 -S anneal:i=<pool size>,a=20` lets every declared genotype reach annealing.
- **Database membership.** The locus database holds exactly the pool's distinct haplotype spellings:
  `locityper add -l <locus> <reference-path region> pool.fa`. Mapping caps and best-relative thresholds depend on the
  database, so membership is tested, not assumed:
  - the **D-only** database is the primary run;
  - **D+X** and **D+O** databases are separate diagnostic runs;
  - a gate reports how much each D genotype's score moves between them;
  - O and X are never ranked together with D.
- **Panvar side.**
  - D and X: `genotype_probe score-set` with allele paths, reproducing every saved expectation.
  - O: `score-set --explicit`, likelihood only, with no prior.
  - Panvar's ranking is reported both by likelihood (F + E) and by total (with the copying prior), because Locityper uses
    no copying prior.

## 5. Calibration inventory and the smallest defensible remedy

| required input | status |
|---|---|
| candidate haplotype sequences | **exists** (models + manifest; explicit FASTA for O) |
| sample reads | **exists** (hash-bound) |
| declared generative library | **exists**: insert 350 ± 40, substitution 0.002, fragment intensity (`library.txt`). The simulator is `genotype_probe simulate`: Poisson templates per haplotype, substitutions only, **no simulated outliers**, uniform coverage |
| reference genome FASTA + jellyfish counts | **missing.** No GRCh38/CHM13 genome is on the machine, and downloads are not permitted |
| background WGS reads for the sample | **missing.** The cells are locus-only simulations |
| Locityper calibration (insert, error profile, depth NB, window size) | **missing.** It must come from `preproc` on background reads |

The October 15 smoke used an artificial uniform-random background chromosome. It was not a matched accuracy benchmark,
and its result is not reused.

**Remedy (smallest defensible; no hand-written calibration, no tuned weights).** Build a **simulation-matched
background** and let Locityper estimate its own calibration from it:
1. **Background sequence:** real human sequence already in Panvar's records and unrelated to C4 and CYP2D6. These are the
   GRCh38 paths of the ACOT pilot model (`GRCh38#0#chr14:73416151-73668669`, ≈ 252 kb) and the GSTM1 pilot model
   (`GRCh38#0#chr1:109561616-109802331`, ≈ 240 kb), spelled from `genotype_pilot_{acot,gstm1}_2x_d24e471/panel.sorted.gfa`.
2. **Background reads:** simulated with **the cells' own simulator** at each cell's depth (5× or 10×), with the cell's
   declared library and the declared seed rule. A one-path GFA wrapper feeds `genotype_probe simulate`.
3. **Locityper reference:** the background contigs plus the cell's model reference path (C4 `HG00096#1…`; the CYP model's
   reference path). Jellyfish counts are computed over that reference.
4. **`locityper preproc`** on the background reads (`--bg-region` set to the background contig) estimates insert, error
   profile and depth itself.

**What the remedy does not fix, stated rather than hidden.**
- No genome-wide paralogs or off-target reads: k-mer uniqueness is computed on the reduced reference.
- The simulator has no GC bias and no indel errors, so Locityper will learn a substitution-only profile and GC-flat depth.
- About 0.49 Mb of background may be too little for Locityper's depth estimation. A gate stops the run if `preproc`
  fails or reports too few usable windows. Random sequence is **not** substituted without a new declaration.

So the comparison calibrates Locityper to the **simulator**, exactly as the cells were generated. It is not a test of
Locityper under real-WGS calibration.

## 6. The executable proposal (M7-B; frozen here, not run)

**Minimal diagnostic implementation needed (no new scoring code):**
- `m7b/pool_fasta.py`: spells each cell's D, X and O haplotypes with the existing `genotype_probe spell` (allele paths)
  or copies the O records, and writes `pool.<tier>.fa` plus the `--priors` file.
- `m7b/bg.py`: writes the one-path background GFA wrappers and runs `genotype_probe simulate`.
- `m7b/compare.py`: parses `res.json.gz` and the score-set tables, and applies the planned comparisons.
- A runner in the established pattern.

**Commands per cell.** `<C>` is the cell, `<D>` its depth, `<L>` its locus.

```text
# Panvar objective (existing tool)
genotype_probe score-set --model-gfa <model> --reference <ref path> --donor <donor> --run <saved run prefix> \
   --reads1 … --reads2 … --library <library.txt> --candidates m7b/<C>.pool.D.tsv --out m7b/<C>/panvar.D
genotype_probe score-set … --explicit m7b/<C>.pool.O.tsv --out m7b/<C>/panvar.O          (C4 LOO only)

# simulation-matched background, then Locityper's own calibration
genotype_probe simulate --gfa m7b/bg.gfa --reference BG --donor BG --depth <D> --insert-mean 350 --insert-sd 40 \
   --error 0.002 --seed <declared> --out m7b/<C>/bg
samtools faidx m7b/<C>/ref.fa;  jellyfish count -m 25 -s 100M -C -o m7b/<C>/ref.jf m7b/<C>/ref.fa
locityper preproc -i m7b/<C>/bg/reads_1.fq m7b/<C>/bg/reads_2.fq -r m7b/<C>/ref.fa -j m7b/<C>/ref.jf \
   --bg-region <BG contig region> -o m7b/<C>/preproc -@ 4

# fixed-pool genotyping, primary database = D only
locityper add -d m7b/<C>/db.D -r m7b/<C>/ref.fa -j m7b/<C>/ref.jf -l <L> <reference-path region> m7b/<C>/pool.D.fa
locityper genotype -i <reads_1> <reads_2> -d m7b/<C>/db.D -p m7b/<C>/preproc -o m7b/<C>/loc.D \
   --priors m7b/<C>/priors.D.tsv --filt-diff 1e9 --prob-thresh -1e9 --dont-skip \
   -S greedy:i=5k,a=1 -S anneal:i=<|D|>,a=20 -s <declared seed> -@ 4
# diagnostics: db.DX (CYP LOO) and db.DO (C4 LOO), same flags, scored separately
```

**Gates (all must pass before any comparison is read):**

| gate | check |
|---|---|
| M1 | every input digest equals `BINDINGS.tsv`; the manifest digest equals `a8cd9516…` |
| M2 | Panvar `score-set` reproduces every saved expectation in the manifest (\|Δ\| ≤ 1e-6) |
| M3 | each Locityper database holds exactly the pool's distinct haplotypes (per-contig SHA-256 equals the spelled manifest), and nothing else |
| M4 | Locityper evaluates exactly the declared genotypes: the `res.json.gz` options set equals the `--priors` set; none is filtered; each has a finite `lik_mean` from 20 annealing attempts |
| M5 | simulator identity: regenerating one cell's archived source reads from its declared donor, seed and library reproduces them byte for byte. This is an evaluation-side check of the simulator only. The campaign build (`panvar-worktrees/mosaic-11fe681/build/genotype_probe`) exists, but M7-A did **not** locate the record of how the archived `source` reads were generated. M7-B must resolve that provenance first. If it cannot, M5 is reported UNRESOLVED, the background is simulated with the campaign build and the declared library, and the gap is disclosed |
| M6 | `preproc` succeeds and reports its usable background windows; failure stops the run |
| M7 | pool-membership sensitivity: each D genotype's `lik_mean` difference between D-only and D+X / D+O databases is reported |
| M8 | reproducibility: one cell repeated with the same seed; differing `lik_mean` values are reported as stochastic search noise |

**Planned comparisons (reported separately; truth labels are evaluation-only):**

| comparison | measure, under each objective (Panvar F + E, Panvar total, Locityper `lik_mean` converted to nats) |
|---|---|
| total copy number | rank and score gap of the best D candidate with the true unit total against the best with a wrong total (CYP LOO: 5 vs 4) |
| allocation between homologues | among D candidates with the true total, rank and gap of the true allocation against wrong allocations (C4 LOO: 2+2 vs 3+1 / 1+3) |
| sequence-content mismatch | block excess over the availability floor of each objective's top D candidate, and whether, within the true class, lower excess ranks higher. Oracle constructions (C4 LOO, O tier) are reported separately as a ceiling, never as competitors |
| correct-control deterioration | LZO cells: rank of the caller's correct call and its gap to each objective's top. Any demotion counts as deterioration |

All differences in observation model, retained alignments, search and calibration stay listed beside every result
(§3, §5).

**Expected cost.**

| step, per cell | estimate |
|---|---|
| Panvar `score-set` (11–14 candidates; M6-A: 6.8 s for 6 at 10×) | 10–25 s |
| background simulation + jellyfish + `add` + `preproc` (smoke: 0.5 + 13 + 2.8 s with 131 haplotypes) | 15–30 s |
| fixed-pool genotyping (≤ 14 genotypes; the smoke's 8,515-genotype greedy stage took 29 s) | 10–40 s |

Total about 6–8 min for four cells plus diagnostics and gates. **Proposed bound: 900 s, 4 threads**, one process at a
time.

**Not claimed in advance.**
- A better fixed-pool ranking under either objective does not establish truth-free candidate discovery, calibration
  under real WGS, or improved genotyping.
- The results will be specific to four inspected development cells.
