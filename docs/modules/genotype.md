# `panvar genotype`

Genotypes a structural locus from paired short reads against the locus pangenome graph. It reports the diplotype
(two whole-locus haplotypes) over the graph's panel, and optionally corrects dosage (unit copy number) errors.

This page states what is frozen at `26401b4`, the C4 short-read freeze. The contracts are in
[`algorithms/genotype_mosaic.md`](../algorithms/genotype_mosaic.md) and
[`algorithms/genotype_dosage.md`](../algorithms/genotype_dosage.md). The planned user-facing command and presets are in
[`algorithms/genotype_cli.md`](../algorithms/genotype_cli.md). It is partly implemented at `c758ae2`:
`--preset short-read`, `--library auto` and the basic flags.

## One adaptive algorithm

There is one production algorithm, **panel-anchored adaptive completion**, selected by a profile.
- It scores every panel diplotype and up to K off-panel mosaics built from the graph's alleles.
- It repairs the incumbent by exact neighbour search, branching exactly at improving-move ties.
- It certifies the winner against every rival beyond floating-point error.
- A representation guard and a combined block quality annotate the call.

Leave-zero-out, leave-one-out and synthetic-mosaic runs are **evaluation scenarios**, not modes. They differ only in
the graph and curation supplied.

The frozen C4 invocation:

    panvar genotype --profile profiles/genotype_c4_v6_campaign_local.tsv --dosage-mode repair \
      --gfa <graph> --reference-path <reference path> --reads1 R1.fq --reads2 R2.fq \
      --library library.txt --threads 8 --out-prefix <prefix> [--panel-curation <manifest>]

**`--library` is required.** It takes one of two values:
- **`--library FILE`:** a declared insert distribution, error, outlier and fragment-start intensity. This is the
  reproducibility path, unchanged.
- **`--library auto`** (`c758ae2`): estimates the insert centre and spread from decisively placed pairs (at most
  20,000 examined, at least 500 required).
  - Every other value comes from `--preset short-read`, including the fragment-start intensity, which is not inferred.
  - A failed validation refuses the run as `LIBRARY_REFUSED` and names `--library FILE`.
  - The sources are written to `<prefix>.library_estimate.tsv`.
  - **It is not the default,** and it changes calls at low coverage: see
    [`algorithms/library_estimate.md`](../algorithms/library_estimate.md) §5.

**`--library-diagnostic`,** off by default, writes an estimate beside a `--library FILE` run and is never used by the
call.

**The basic command:** `--graph G --reads1 R1 --reads2 R2 --out PREFIX --preset short-read --library FILE|auto
--threads N`.
- `--preset short-read` switches on the production algorithm with short-read technology constants. It sets no locus,
  resource or dosage value.
- `--preset long-read` is refused as not implemented.
- Without `--reference-path`, the graph's declared reference sample (`RS:Z`) is used. With none declared, the run
  refuses and lists the candidates. The C4 graph declares none.

## `--dosage-mode off|audit|repair` (default `off`)

| mode | what it does |
|---|---|
| `off` | no dosage work; every output is byte-identical to the caller without the dosage stage |
| `audit` | runs the dosage stage (census, features, classification, repair search, provenance) and writes only `<prefix>.dosage_*`; the call is unchanged |
| `repair` | as audit. A block whose dosage repair is complete and resolved (COMPLETE_DIPLOTYPE_REPAIR) gets the repaired genotype; the caller's outputs are kept as `<prefix>.original.*`, a complete prefix the evaluator scores unchanged |

`repair` changes a genotype only under the four eligibility rules of `genotype_dosage.md` §22.2. Otherwise the block
is marked `DOSAGE_REPAIR_NOT_APPLIED` with the reason.

## Outputs: locus level and block level

- **`<prefix>.locus_genotype.tsv`: the locus-level call,** one row, written in every mode.
  - Run status and no-call reason, call class, locus filter, quality counts, phase support, and candidate and search
    scope.
  - The two haplotypes: IDs, SHA-256 digests, lengths and ordered block paths.
  - The unordered genotype, score and separation (with their basis), and panel-template switch counts and positions.
  - The dosage counts, the whole-locus repair scope, and links to every related file.
  - **By call class:** a genotype-equivalent call reports the unordered genotype only, with no phase and no
    haplotypes. A no-call reports status and reason only.
- **`<prefix>.block_genotypes.tsv` and `<prefix>.genotypes.tsv`: the per-block decomposition.** In repair mode they
  carry trailing dosage and provenance columns.
- **`<prefix>.mosaic_haplotypes.fa`: the two haplotypes of a strict call.** The `hap1`/`hap2` labels assert phase only
  where `phase_support` is `READS_WITHIN_CANDIDATES` or `READS_AUDITED`.
- **`mosaic_scores.tsv`, `mosaic_candidates.tsv`, `equivalence_set.tsv` and the repair trajectory: the search record.**
  They describe what the caller scored. After a dosage repair they still mark the caller's winner, and
  `dosage_search_winner_superseded true` says so.

## Three kinds of result, never mixed

- **Exact genotype.** The reported block alleles and haplotypes: the caller's, or an applied COMPLETE_DIPLOTYPE_REPAIR
  (`dosage_repair_role=RESOLVED`).
- **Dosage content** (`corrected_content`, for `CONDITIONAL_DOSAGE_CONTENT`). The corrected unit count, HERV (I-state)
  count and feature copy vector, authoritative within the declared model.
  - The block keeps the caller's alleles.
  - "Complete" means every admissible candidate shares that content. It does not mean the content is correct.
- **Conditional representative.** One arrangement serialising a set of arrangements the reads cannot distinguish:
  `REPRESENTATIVE_ONLY_NOT_RESOLVED`, only in `dosage_repaired*.fa` and the equivalence table.
  - It is never a call.
  - It never enters the allele columns, the locus digests or the main FASTA.

## Filters

Each filter describes the **reported** genotype. After an applied repair, the caller's values are in `original_*`.

| field | values | meaning |
|---|---|---|
| `filter` (representation) | `PASS`, `LOWQUAL_REPRESENTATION_CONFLICT`, `NOT_EVALUATED` | the marker depth and the exposure term prefer another pair while the fragment term prefers the call, each by ≥ ln 1000 (`representation_conflict_rule_v1`) |
| `genotype_support` | `READS`, `PRIOR_ONLY` | whether the reads alone put every rival of another genotype at this block below the reported diplotype |
| `quality_filter` | `PASS`, `LOW_QUALITY` | `block_quality_rule_v1`: PASS only for representation PASS with READS support |
| `dosage_filter` (repair mode) | `PASS`, `DOSAGE_CONTENT_CORRECTED`, `DOSAGE_TRIGGERED_UNREPAIRED`, `DOSAGE_SUSPECT`, `DOSAGE_UNDERPOWERED`, `DOSAGE_NOT_EVALUATED`, `DOSAGE_REPAIR_NOT_APPLIED` | the dosage stage's verdict on the reported block; information limits are reasons, never merged into `quality_filter` |
| `locus_filter` | `PASS`, `LOWQUAL_REPRESENTATION_CONFLICT`, `NOT_EVALUATED` | fails closed over the blocks |

A repaired run's filters are recomputed on the repaired diplotype. If the caller's own annotation cannot first be
reproduced, the value is `NOT_EVALUATED_AFTER_DOSAGE_REPAIR`: never copied, never PASS.

## C4 timing and memory (M4, 8 threads)

Genotyping: the sealed 484-job campaign, per job.

| scenario | 2× median / max | 30× median / p90 / max | peak footprint at 30× |
|---|---|---|---|
| leave-zero-out | 13.2 / 23.2 s | 72.1 / 96.9 / 112.5 s | 3.45 GiB |
| synthetic mosaic | 26.3 / 35.2 s | 146.4 / 295.4 / 347.3 s | 3.13 GiB |
| leave-one-out | 17.1 / 38.5 s | 117.5 / 231.5 / 306.7 s | 3.58 GiB |

**Dosage stage, repair mode** (649-cell replay at `26401b4`, 4 threads):
- The per-cell mean is 2.7–8.0 s.
- The worst cell took 165.9 s.
- The peak footprint was 2.85 GiB.

## Known limitations

**Scope of validation:**
- **Validated on one locus (C4), with simulated short reads and a declared library.** The simulator writes constant base
  qualities.
- **Resource caps are sized on C4.** The dosage default is `off` until other loci validate it.
- **Leave-one-out accuracy is limited by availability.** The truth is often not spellable from the panel, so read
  excess, not raw distance. The floor is attained in 58/67/75/84% of blocks at 2/5/10/30×.

**Calls:**
- **Two no-calls in the frozen C4 set,** both `tie_branching_incomplete`: `loo.HG01573.30`, and fresh `loo.HG01573.5`.
- **HG00096 is the C4 coordinate reference.** Leave-one-out on it is not a supported scenario.

**Dosage:**
- **Which module carries a HERV (allocation) cannot be seen from short reads.** Conditional representatives are
  non-authoritative.
- **Dosage content can be complete but wrong** when the bootstrap interval excludes the truth (fresh `loo.NA19434.2`).
- **The four applied C4 complete repairs (HG03456 sv:1) have no defined truth.** A synthetic known-truth repair
  passes. The repaired pairs are dosage-supported, not read-supported (`PRIOR_ONLY`).

**Inputs and builds:**
- **Automatic library estimation is a diagnostic only.**
- **Rebuilds are not byte-identical.** A rebuilt `panvar` differs from the frozen binary only in the Mach-O UUID,
  debug map and code signature; the code and data sections are identical.
