# panvar

`panvar` is a modular C++ toolkit for pangenome-graph analysis from GFA.

Modules:

1. `bubble`: site extraction/refinement from precomputed `vg snarls`
2. `inspect`: path FASTA and node/edge traversal matrix for called bubbles
3. `panphorte`: normalize tandem-repeat bubbles into a compact copy-number-explicit GFA
4. `call`: graph-native structural variant calling (DEL/INS/INV/DUP) into a multi-sample VCF
5. `describe`: YYY

## Build

```bash
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Binary:

```bash
./build/panvar
```

## End-to-End Example (C4)

XXX

## Documentation

- [docs/README.md](docs/README.md)
- [docs/modules/bubble.md](docs/modules/bubble.md)
- [docs/modules/inspect.md](docs/modules/inspect.md)
- [docs/modules/panphorte.md](docs/modules/panphorte.md)
- [docs/modules/call.md](docs/modules/call.md)
- [docs/modules/describe.md](docs/modules/describe.md)
- [docs/modules/example.md](docs/modules/example.md)
- [docs/gwas_example.md](docs/gwas_example.md) — worked pangenome k-mer GWAS (LPA KIV-2, multiplicity vs presence/absence)
- [docs/presentation.md](docs/presentation.md) — slide-ready schematics (panphorte + SV types)

## Tests

```bash
ctest --test-dir build --output-on-failure
```

Direct smoke test:

```bash
tests/smoke.sh ./build/panvar tests/real_data/c4.gfa /tmp/panvar_smoke_c4
```

## Docker

The image builds `panvar` and includes related tools from Bioconda:

- `vg`
- `odgi`
- `minimap2`
- `samtools`
- `bcftools`
- `R` + `ggplot2` for helper plotting scripts

Build:

```bash
docker build -t panvar:latest .
```

Run:

```bash
docker run --rm -v "$PWD":/work -w /work panvar:latest panvar --help
```