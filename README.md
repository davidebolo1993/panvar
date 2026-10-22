# panvar

`panvar` is a modular C++17 toolkit for pangenome-graph SV calling and association testing. It runs an end-to-end pipeline — snarl finding, tandem-repeat normalization, graph-native structural-variant calling, feature description, and association testing: `bubble → panphorte → call → describe → associate`, with an opt-in `rebuild` step (re-induce a locus graph too fragmented to decompose) before `bubble`, an opt-in `refine` step (POA-realign bubble interiors to remove graph-builder alignment artifacts) between `panphorte` and `call`, the `inspect` utility available at any step, and `benchmark` for assessing the quality of the calls.

## Install

```bash
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/panvar --help
```

#### Requirements

A C++17 compiler, CMake ≥ 3.16, `make`, and zlib. Five libraries ship as git submodules (hence the `git submodule update --init --recursive` step — nothing extra to install): `minimap2`, statically linked for INS-subtype realignment; `Eigen`, the header-only linear-algebra library used for linear mixed models; `edlib`, the bit-parallel edit-distance library behind the banded fit alignment; `abPOA` (with its nested `SIMDe`), the partial-order aligner behind `refine`; and `minigraph`, the progressive graph generator behind `rebuild`. Run the enabled regression suite with `ctest --test-dir build --output-on-failure` (several minutes, not just a smoke test). Genotype assessment fixtures require Python 3; plotting tests additionally require NumPy, pandas and matplotlib. See [tests/README.md](tests/README.md) for focused test groups.

## Short-read genotyping (experimental)

Genotype one diploid sample at one locus, using a graph panel and paired reads:

```bash
mkdir -p results
./build/panvar genotype --graph locus.gfa --reference-path REF \
  --reads1 sample.R1.fq.gz --reads2 sample.R2.fq.gz --threads 4 --out results/sample
```

Replace REF with the coordinate-reference path. The command estimates the library automatically and needs no
profile. Add `--dosage-check` for a diagnostic copy-number audit; it never changes the prediction. The optional
`--build-panel-index` / `--panel-index` flags reuse locus preparation across samples.

Read the [genotype manual](docs/modules/genotype.md) for inputs, all flags, output layout, assessment meanings
and evidence limits. Primary reports sit beside the output prefix; supporting reports are in its `.diagnostics/`
directory. A completed prediction is not calibrated reliability. Evidence is limited to simulated short reads,
with reserved evaluation at ACOT and development-only C4/CYP2D6 results. Long reads are not implemented.

## Docker

The image builds `panvar` and bundles a few companion tools from Bioconda for adjacent graph work and for
handling `panvar`'s outputs: `vg` (snarls), `odgi` (check graph sort/flip), `bcftools` (work with the called
VCFs), and `R` + `ggplot2` + `ggrepel` (the plotting scripts):

```bash
docker build -t panvar:latest .
docker run --rm -v "$PWD":/work -w /work panvar:latest panvar --help
```

## Documentation

- [docs/README.md](docs/README.md) — documentation index
- **Modules:** [rebuild](docs/modules/rebuild.md) · [bubble](docs/modules/bubble.md) · [panphorte](docs/modules/panphorte.md) · [refine](docs/modules/refine.md) · [call](docs/modules/call.md) · [describe](docs/modules/describe.md) · [associate](docs/modules/associate.md) · [genotype](docs/modules/genotype.md) · [inspect](docs/modules/inspect.md) · [benchmark](docs/modules/benchmark.md)
- **Algorithms:** [rebuild](docs/algorithms/rebuild.md) · [bubble](docs/algorithms/bubble.md) · [panphorte](docs/algorithms/panphorte.md) · [refine](docs/algorithms/refine.md) · [call](docs/algorithms/call.md) · [describe](docs/algorithms/describe.md) · [associate](docs/algorithms/associate.md) · [genotype](docs/algorithms/genotype.md) · [inspect](docs/algorithms/inspect.md) · [benchmark](docs/algorithms/benchmark.md)
- [docs/walkthrough.md](docs/walkthrough.md) — the full pipeline on the LPA locus, step by step with plots
- [docs/gwas.md](docs/gwas.md) - an example association testing (LPA KIV-2 repeat ~ Lp(a))
- [docs/references.md](docs/references.md) — tools and papers per module

## License

MIT — see [LICENSE](LICENSE).
