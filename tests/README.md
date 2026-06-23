# Tests

Test data, smoke tests, and the GWAS demo/validation drivers for `panvar`. The build registers one fast test
with CTest (`synthetic_smoke.sh`); everything else is run manually. Together the scripts exercise all six
subcommands (`bubble`, `panphorte`, `call`, `describe`, `inspect`, `associate`) — see the coverage table at
the end.

## Data

`real_data/` — small bundled graphs and copy-number ground truth for four example loci:

| file | what it is |
|------|------------|
| `lpa.gfa.gz`, `c4.gfa.gz`, `gstm1.gfa.gz`, `cyp2d6.gfa.gz` | pangenome graphs (GFA) for the four loci |
| `c4.snarls.jsonl`, `gstm1.snarls.jsonl`, `lpa.snarls.jsonl` | precomputed `vg`-style snarls for `bubble --snarls-in` |
| `c4.bed`, `cyp2d6.bed`, `gstm1.bed`, `lpa.repeats.tsv` | copy-number ground truth (gene-count / direct / KIV-2 repeats) |
| `Homo_sapiens.GRCh38.116.gtf.gz` | reference GTF for `--gtf` gene annotation |

`synthetic_data/` — tiny hand-built graphs for the synthetic smoke: `syn.gfa` / `syn_w.gfa` (with their
`.snarls.jsonl`), `expected.tsv` (the asserted records), and `make_synthetic.py` (regenerates them).

Outputs go to an output directory you pass on the command line (e.g. `/tmp/...` or the gitignored
`results/`); there is no checked-in results directory.

## Smoke tests

- `synthetic_smoke.sh` — fast, dependency-free; runs `bubble` + `call` on the synthetic graphs and asserts
  every `call` event type (DEL/INS/INV/DUP, substitution, and the flip case) against exact records.
  Registered with CTest.
- `real_smoke.sh` — fuller integration run of `bubble -> inspect -> describe -> panphorte -> call` on a real
  locus; checks the main handoff files plus the compressed inspect/describe outputs. Not in CTest (slower).

```bash
ctest --test-dir build --output-on-failure                          # runs synthetic_smoke
tests/synthetic_smoke.sh ./build/panvar tests/synthetic_data /tmp/panvar_smoke
tests/real_smoke.sh      ./build/panvar tests/real_data/c4.gfa.gz /tmp/panvar_smoke_c4
```

## GWAS drivers (`tests/gwas/`)

These exercise `describe --samples` and `associate`, which the smoke tests above do not.

- `run_lpa_real.sh` — end-to-end LPA association demo: `bubble -> panphorte -> call -> describe --samples ->
  associate` on the real LPA graph, plus the population-structure-correction demo and the plots. Needs a
  numpy/scipy-capable python and (for plots) `Rscript` + `ggplot2`.
- `make_lpa_phenotype.py` — simulates the LPA/KIV-2 cohort (phenotype, covariates, kinship, and the synthetic
  genome-wide-like panel) consumed by both drivers here.
- `validate_gemma.sh` — checks `panvar associate` (linear and LMM) against GEMMA on identical BIMBAM inputs;
  reports the Pearson r of effect size and of −log10 p, and skips cleanly if GEMMA is not installed. It reuses
  the `panphorte.copies.tsv` from an LPA run (default path under `results/`), so run `run_lpa_real.sh` (or
  `scripts/regen_results.sh lpa`) first.

```bash
# usage: run_lpa_real.sh   <panvar_bin> <out_dir> [python] [Rscript] [--big]
bash tests/gwas/run_lpa_real.sh   ./build/panvar results/real_data/lpa/gwas python3 Rscript

# usage: validate_gemma.sh <panvar_bin> <out_dir> [python] [gemma_bin] [copies.tsv]
bash tests/gwas/validate_gemma.sh ./build/panvar /tmp/gemma_check python3 gemma
```

## Subcommand coverage

| subcommand | synthetic_smoke | real_smoke | run_lpa_real | validate_gemma |
|------------|:---:|:---:|:---:|:---:|
| `bubble`    | ✓ | ✓ | ✓ |   |
| `panphorte` |   | ✓ | ✓ |   |
| `call`      | ✓ | ✓ | ✓ |   |
| `describe`  |   | ✓ | ✓ |   |
| `inspect`   |   | ✓ |   |   |
| `associate` |   |   | ✓ | ✓ |
