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

A C++17 compiler, CMake ≥ 3.16, `make`, and zlib. Seven libraries ship as git submodules (hence the `git submodule update --init --recursive` step — nothing extra to install): `minimap2`, statically linked for realigning long alleles and inserted sequence (`call`, `rebuild`); `Eigen`, the header-only linear-algebra library used for linear mixed models; `edlib`, the bit-parallel edit-distance library behind the sequence alignments; `abPOA` (with its nested `SIMDe`), the partial-order aligner behind `refine`; `minigraph`, the progressive graph generator behind `rebuild`; `SSHash` (with its nested PTHash), the k-mer dictionary behind `background --compact`; and `htslib` (with its nested htscodecs), which reads SAM, BAM and CRAM for `recruit`. bzip2, liblzma and libdeflate are used when installed (CRAM codecs, faster BGZF decompression). A fast smoke test runs with `ctest --test-dir build --output-on-failure` (CTest 3.20 or later; with an older CMake, run `ctest --output-on-failure` inside `build/`); the tests also need `python3`.

## Docker

The image holds `panvar` and every tool its scripts and the `dev/` checks call: `vg`, `odgi`, `bcftools` and `htslib`; `minimap2`, `samtools`, `svim-asm`, `truvari` and `gemma`; Python with numpy; and R with the packages of the plots, the viewer and the association validation. The entrypoint is `panvar`; other tools need `--entrypoint`. The scripts are not in the image, so mount the repository and give them the image's binary:

```bash
docker build -t panvar:latest .
docker run --rm -v "$PWD":/work panvar:latest bubble --help
docker run --rm -v "$PWD":/work --entrypoint /usr/local/bin/_entrypoint.sh panvar:latest \
  bash scripts/bubble_examples.sh /usr/local/bin/panvar
```

## Documentation

- [docs/README.md](docs/README.md) — documentation index
- Modules: [rebuild](docs/rebuild.md) · [bubble](docs/bubble.md) · [panphorte](docs/panphorte.md) · [refine](docs/refine.md) · [call](docs/call.md) · [describe](docs/describe.md) · [associate](docs/associate.md) · [inspect](docs/inspect.md) · [benchmark](docs/benchmark.md)

## License

MIT — see [LICENSE](LICENSE).
