# panvar

`panvar` is a modular C++17 toolkit for pangenome-graph analysis from GFA. It runs an end-to-end pipeline —
snarl finding, tandem-repeat normalization, graph-native structural-variant calling, feature description, and
association testing: `bubble → panphorte → call → describe → associate`, with the `inspect` utility available
at any step.

## Install

```bash
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/panvar --help
```

Requirements: a C++17 compiler, CMake ≥ 3.16, `make`, and zlib. Two libraries ship as git submodules (hence
the `git submodule update --init --recursive` step — nothing extra to install): `minimap2`, statically linked
for INS-subtype realignment, and `Eigen`, the header-only linear-algebra library used by the `associate` LMM.
A fast, dependency-free smoke test runs with `ctest --test-dir build --output-on-failure`.

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
- **Modules (usage):** [bubble](docs/modules/bubble.md) · [panphorte](docs/modules/panphorte.md) · [call](docs/modules/call.md) · [describe](docs/modules/describe.md) · [associate](docs/modules/associate.md) · [inspect](docs/modules/inspect.md)
- **Algorithms (mechanism + worked traces):** [docs/algorithms/](docs/algorithms/)
- **GWAS:** [example](docs/gwas/example.md) — runnable LPA walk-through (concepts explained in context) + GEMMA validation
- [docs/references.md](docs/references.md) — tools (GitHub) and papers per module

## License

MIT — see [LICENSE](LICENSE).
