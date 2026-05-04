# panvar

`panvar` is a modular C++ toolkit for pangenome-graph analysis from GFA.

Current modules:

1. `bubble` (Module 1): site extraction/refinement from precomputed `vg snarls`
2. `allele` (Module 2): allele extraction and clustering from module-1 bubbles
3. `call` (Module 3): SV calling (`INS`/`DEL`/`INV`) from clustered alleles
4. `describe` (Module 4): scaffold command (planned)

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

## Quick Start (C4)

1) Bubble:

```bash
./build/panvar bubble \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/bubble \
  --snarls-in tests/real_data/c4.snarls.jsonl
```

2) Allele:

```bash
./build/panvar allele \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/allele \
  --bubbles-csv-in tests/results/c4/bubble.bubbles.csv
```

3) Call:

```bash
./build/panvar call \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/call \
  --bubbles-csv-in tests/results/c4/bubble.bubbles.csv \
  --clusters-csv-in tests/results/c4/allele.allele_clusters.csv \
  --assignments-csv-in tests/results/c4/allele.allele_assignments.csv \
  --reference-path "$(awk -F '\t' '$1==\"P\" || $1==\"W\" { print $2; exit }' tests/real_data/c4.gfa)"
```

Useful call options:

- `--debug --debug-out-dir <dir>`: write per-cluster FASTA/PAF/dotplot/VCF
- `--minimap-preset asm5|asm10|asm20` (default `asm20`)
- `--dotplot-gtf <gtf.gz>` with repeatable `--dotplot-gene-match <pattern>`

## Documentation

- [docs/README.md](docs/README.md)
- [docs/modules/bubble.md](docs/modules/bubble.md)
- [docs/modules/allele.md](docs/modules/allele.md)
- [docs/modules/call.md](docs/modules/call.md)
- [docs/modules/describe.md](docs/modules/describe.md)

## Tests

Run configured smoke tests:

```bash
ctest --test-dir build --output-on-failure
```

Or run one smoke test directly:

```bash
tests/smoke.sh ./build/panvar tests/real_data/c4.gfa /tmp/panvar_smoke_c4
```

## Docker

The Docker image builds `panvar` and includes related graph/SV tools from Bioconda:

- `vg`
- `odgi`
- `minimap2`
- `samtools`
- `bcftools`

Build:

```bash
docker build -t panvar:latest .
```

Run:

```bash
docker run --rm -v "$PWD":/work -w /work panvar:latest --help
```

Tool binaries from the container (for example `vg`, `odgi`) are available directly in the shell.

## Repository Layout

- `src/`, `include/panvar/`: C++ implementation
- `external/minimap2/`: minimap2 C library source
- `tests/real_data/`: bundled loci and snarl JSONL inputs
- `tests/results/`: local output artifacts
- `docs/modules/`: module docs
