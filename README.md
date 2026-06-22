# panvar

`panvar` is a modular toolkit for pangenome-graph analysis from GFA. It runs an
end-to-end pipeline — snarl finding, graph normalization, structural-variant calling,
and variant description for association studies.

## Modules

| Command | Role | What it does |
| --- | --- | --- |
| `bubble` | Module 1 | Extract bubble sites from a GFA. |
| `panphorte` | Module 2 | Normalize tandem-repeat bubbles into a compact, copy-number-explicit GFA. |
| `call` | Module 3 | Graph-native structural-variant calling. |
| `describe` | Module 4 | Per-bubble haplotype features for association. |
| `inspect` | Utility | For one or all called bubbles: per-path sequence, allele clustering, and node/edge traversal matrices. |

Typical flow: `bubble` → `panphorte` → `call` → `describe`, with `inspect` available at any point for ad-hoc checks.

## Build

```bash
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Requirements: a C++17 compiler, CMake ≥ 3.16, `make`, and zlib. `minimap2` is bundled as a
submodule and statically linked — no system install needed.

Binary:

```bash
./build/panvar --help
```

## Tests

The build registers a fast, dependency-free smoke test that exercises every `call` event type on
tiny hand-built graphs with exact-record assertions:

```bash
ctest --test-dir build --output-on-failure
```

Run that smoke test directly:

```bash
bash tests/synthetic_smoke.sh ./build/panvar tests/synthetic_data /tmp/panvar_smoke
```

A real-data integration smoke (`tests/real_smoke.sh`) runs the whole `bubble → inspect →
describe → panphorte → call` pipeline on a real locus; it is intentionally not part of `ctest`:

```bash
bash tests/real_smoke.sh ./build/panvar tests/real_data/c4.gfa.gz /tmp/smoke_c4
```

## Documentation

- [docs/README.md](docs/README.md) — documentation index (and the copy-number topology table)
- **Modules (usage):** [bubble](docs/modules/bubble.md) · [panphorte](docs/modules/panphorte.md) · [call](docs/modules/call.md) · [describe](docs/modules/describe.md) · [associate](docs/modules/associate.md) · [inspect](docs/modules/inspect.md)
- **Algorithms (mechanism + worked traces):** [docs/algorithms/](docs/algorithms/)
- **GWAS:** [primer](docs/gwas/primer.md) (concepts) · [example](docs/gwas/example.md) (runnable LPA walk-through + GEMMA validation)
- [docs/references.md](docs/references.md) — tools (GitHub) and papers per module

## Docker

The image builds `panvar` and bundles companion tools from Bioconda for adjacent graph work and for handling `panvar`'s outputs:

- `vg`, `odgi`
- `minimap2`, `samtools`, `bcftools`
- `R` + `ggplot2`

Build:

```bash
docker build -t panvar:latest .
```

Run:

```bash
docker run --rm -v "$PWD":/work -w /work panvar:latest panvar --help
```

## License

MIT — see [LICENSE](LICENSE).
