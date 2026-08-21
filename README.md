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

A C++17 compiler, CMake ≥ 3.16, `make`, and zlib. Five libraries ship as git submodules (hence the `git submodule update --init --recursive` step — nothing extra to install): `minimap2`, statically linked for INS-subtype realignment; `Eigen`, the header-only linear-algebra library used for linear mixed models; `edlib`, the bit-parallel edit-distance library behind the banded fit alignment; `abPOA` (with its nested `SIMDe`), the partial-order aligner behind `refine`; and `minigraph`, the progressive graph generator behind `rebuild`. A fast smoke test runs with `ctest --test-dir build --output-on-failure`.

## Docker

The image builds `panvar` and bundles a few companion tools from Bioconda for adjacent graph work and for
handling `panvar`'s outputs: `vg` (snarls), `odgi` (check graph sort/flip), `bcftools` (work with the called
VCFs), and `R` + `ggplot2` + `ggrepel` (the plotting scripts):

```bash
docker build -t panvar:latest .
docker run --rm -v "$PWD":/work -w /work panvar:latest panvar --help
```

## Scope and limitations

Worth knowing before reading any number this project reports:

- **`benchmark` reports four levels and they are not interchangeable.** `graph`, `called` and `carrier` are ceilings — they implant the haplotype's true block — while only `genotype` reconstructs the emitted VCF. Every QV figure quoted before the module review pass was the `graph` column.
- **The region VCF is the interpreted output; the allele VCF is the lossless one.** `call --allele-vcf` reconstructs its haplotypes with 0 bp residual on all six reference loci. The region VCF merges records for readability and reconstructs 39–97% of the reference-to-truth distance depending on locus; a merged record hands every carrier the representative's sequence, and a `DUP` is reconstructed by tiling an inferred span.
- **`associate` implements common single-variant association.** Firth and SPA make a rare single-variant test better behaved; there is no burden, collapsing or SKAT-style aggregation, rare binary p-values in the far tail remain about 1.7× nominal, and the LMM is experimental and not numerically validated against a pinned reference.
- **`genotype` (short-read sample genotyping) is not built by default.** It is the one module that has not completed a review pass. `cmake -DPANVAR_ENABLE_EXPERIMENTAL_GENOTYPE=ON` builds it for development.

The complete list, tagged by whether it blocks a release, is [docs/review_followups.md](docs/review_followups.md); numbers this project reported and later corrected are in [docs/reports/module-review-summary.md](docs/reports/module-review-summary.md).

## Documentation

- [docs/README.md](docs/README.md) — documentation index
- **Modules:** [rebuild](docs/modules/rebuild.md) · [bubble](docs/modules/bubble.md) · [panphorte](docs/modules/panphorte.md) · [refine](docs/modules/refine.md) · [call](docs/modules/call.md) · [describe](docs/modules/describe.md) · [associate](docs/modules/associate.md) · [inspect](docs/modules/inspect.md) · [benchmark](docs/modules/benchmark.md)
- **Algorithms:** [rebuild](docs/algorithms/rebuild.md) · [bubble](docs/algorithms/bubble.md) · [panphorte](docs/algorithms/panphorte.md) · [refine](docs/algorithms/refine.md) · [call](docs/algorithms/call.md) · [describe](docs/algorithms/describe.md) · [associate](docs/algorithms/associate.md) · [inspect](docs/algorithms/inspect.md) · [benchmark](docs/algorithms/benchmark.md)
- [docs/walkthrough.md](docs/walkthrough.md) — the full pipeline on the LPA locus, step by step with plots
- [docs/gwas.md](docs/gwas.md) - an example association testing (LPA KIV-2 repeat ~ Lp(a))
- [docs/references.md](docs/references.md) — tools and papers per module

## License

MIT — see [LICENSE](LICENSE).
