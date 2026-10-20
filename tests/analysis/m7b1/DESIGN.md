# M7-B.1: N-free background and resumed fixed-pool comparison (frozen before any unit)

**Ancestry.** Branch `m7-b1-resume` from `8e992d7` (M7-B stopped result). The M7-B package
`tests/results/m7b_fixed_pool_20261020/` (outputs, failures, manifests, records) is **unchanged** and read only. This
continuation writes `tests/results/m7b1_fixed_pool_20261020/`.

**Execution accounting.** M7-B spent **66.0 s** (recorded there). M7-B.1 counts **new execution only**, against
**800 s**.
- Admission: spent + timeout + 40 s reporting reserve ≤ 800.
- 4 threads, one process at a time, 6 GiB.

## The one authorised change (recorded before execution): remedy (a)

**Background.** `m7bg` is the **direct concatenation**, with no spacer, of the same GRCh38 ACOT
(`GRCh38#0#chr14:73416151-73668669`) and GSTM1 (`GRCh38#0#chr1:109561616-109802331`) spellings.
- **Artificial junction** (1-based): between position 252,518 (the last ACOT base) and 252,519 (the first GSTM1 base). It
  is disclosed as a limitation.
- **Checks:** each source spelling's length equals its region; its SHA-256 is recorded; **the complete background is
  asserted to contain only A/C/G/T**.
- **Representation:** the same homozygous diploid (reference path plus two identical donor paths).

Remedy (b) is not used: preprocessing adds 50-kb padding before counting k-mers, which could include the old spacer.

**Rebuilt** (background-dependent only):
- background reads, with the M7-B-verified simulator lineages and unchanged depths, libraries and seed rule
  (`--seed <cell seed> --seed-salt m7b-bg-<cell>`);
- references (`m7bg` + `c4ref` / `cypref`, with the CYP GRCh38 reference choice kept);
- FASTA indexes, jellyfish counts and `preproc` outputs.

No old reads are N-stripped and no old calibration is reused.

## Reused after binding checks (unit `reuse`)

Reused from M7-B: the provenance verdicts, canonical pools and mapping, Panvar `score-set` tables (D for each cell, X for
CYP LOO) and the truth-distance table. `reuse` checks that:
- the M7-B units exited 0;
- the lineages PASS;
- the pool checks passed;
- `score-set` reproduced its expectations;
- every recorded executable's SHA-256 is unchanged.

It also binds the reused files by SHA-256 (`reused_bindings.tsv`). No completed scoring is repeated.

**Unchanged:** the M7-A manifest, D/X/O separation, Locityper score weights (defaults), optimisation settings, flags and
seeds (s1 = 20261020, s2 = 20261021, preproc 20261020).

## Calibration validation (explicit; `calibration_ok` gates every Locityper run of that cell)

| check | requirement |
|---|---|
| `preproc` | exits 0 |
| insert-size NB | n, p finite, n > 0, 0 < p < 1 |
| error profile | all values finite |
| depth NB | 101 GC bins, every n, p finite with n > 0 and 0 < p < 1 |
| depth window and neighbourhood sizes | positive |
| background windows | counted (`depth.csv.br`); total > 0 |
| GC support | at least one GC bin with ≥ 100 background windows (Locityper's `min_tail_obs`) |

- **Reported, not used as a stop rule:** per Locityper run, the target windows (weight ≥ 0.001) whose GC bin has fewer
  than 100 background windows, i.e. that rely on extrapolated GC ranges.
- **Any failed check stops that cell's comparison.** Background size, thresholds, weights and distributions are never
  changed automatically.
- **Limitation carried throughout:** the background is a small, synthetic calibration dataset made from real sequence.
  It provides neither genome-wide k-mer uniqueness nor real-WGS calibration.

## Order (frozen)

| step | units |
|---|---|
| 1 | reuse; background; reference.C4 / .CYP2D6; preproc for each cell |
| 2 | Locityper D-only, controls first: `loc.c4_lzo_hg00171_5.D.s1`, `loc.cyp_lzo_hg02769_5.D.s1` |
| 3 | failure cells: `loc.c4_loo_na19650_10.D.s1`, `loc.cyp_loo_na19650_5.D.s1` |
| 4 | pool membership (diagnostics only): `loc.cyp_loo_na19650_5.DX.s1`, `loc.c4_loo_na19650_10.DO.s1` |
| 5 | seeds: `loc.c4_lzo_hg00171_5.D.s1r`, `loc.c4_loo_na19650_10.D.s2`, `loc.cyp_loo_na19650_5.D.s2` |
| 6 | compare |

## Questions (kept separate; rankings and gaps within each method only)

- **CYP NA19650:** does the Locityper pipeline favour the correct total dosage (5) over the wrong one (4)?
- **C4 NA19650:** does it rank the correct allocation (2+2) above 3+1 / 1+3 despite missing sequence content?
- **Controls:** is the correct call ranked first, after canonical relabelling?

**Interpretation rules.**
- Unstable rankings (a winning margin ≤ the seed difference of the top two genotypes) are INCONCLUSIVE.
- A favourable result motivates attribution only. It is not evidence that local depth caused it, nor of improved
  genotyping.
