# panvar

`panvar` is a modular C++17 toolkit for pangenome-graph SV calling and association testing. It runs an end-to-end pipeline — snarl finding, tandem-repeat normalization, graph-native structural-variant calling, feature description, and association testing: `bubble → panphorte → call → describe → associate`, with the `inspect` and `benchmark` utilities
available at any step.

## Install

```bash
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/panvar --help
```

#### Requirements

A C++17 compiler, CMake ≥ 3.16, `make`, and zlib. Three libraries ship as git submodules (hence the `git submodule update --init --recursive` step — nothing extra to install): `minimap2`, statically linked for INS-subtype realignment; `Eigen`, the header-only linear-algebra library used by the `associate` LMM; and `edlib`, the bit-parallel edit-distance library behind the banded fit alignment. A fast smoke test runs with `ctest --test-dir build --output-on-failure`.

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
- **Modules:** [bubble](docs/modules/bubble.md) · [panphorte](docs/modules/panphorte.md) · [call](docs/modules/call.md) · [describe](docs/modules/describe.md) · [associate](docs/modules/associate.md) · [inspect](docs/modules/inspect.md) · [benchmark](docs/modules/benchmark.md)
- **Algorithms:** [bubble](docs/algorithms/bubble.md) · [panphorte](docs/algorithms/panphorte.md) · [call](docs/algorithms/call.md) · [describe](docs/algorithms/describe.md) · [associate](docs/algorithms/associate.md) · [inspect](docs/algorithms/inspect.md) · [benchmark](docs/algorithms/benchmark.md)
- [docs/walkthrough.md](docs/walkthrough.md) — the full pipeline on the LPA locus, step by step with plots
- [docs/gwas/example.md](docs/gwas/example.md) - an example association testing (LPA KIV-2 repeat ~ Lp(a))
- [docs/references.md](docs/references.md) — tools and papers per module

## License

MIT — see [LICENSE](LICENSE).
