# panvar

`panvar` is a modular C++ toolkit for pangenome-graph analysis from GFA.

Modules:

1. `bubble` (Module 1): site extraction/refinement from precomputed `vg snarls`
2. `inspect` (utility): path FASTA and node traversal matrix for one or all called bubbles
3. `allele` (Module 2): allele extraction and clustering from module-1 bubbles
4. `call` (Module 3): SV calling (`INS`/`DEL`/`INV`) from clustered alleles
5. `describe` (Module 4): per-bubble k-mer feature tables for downstream association

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

1. Bubble:

```bash
./build/panvar bubble \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/bubble \
  --snarls-in tests/real_data/c4.snarls.jsonl \
  --merge-nearby-bp 20
```

2. Inspect one bubble:

```bash
./build/panvar inspect \
  -i tests/real_data/c4.gfa \
  --bubble-prefix-in tests/results/c4/bubble \
  --bubble-id 1 \
  -o tests/results/c4/inspect/bubble_1
```

Inspect all bubbles by omitting `--bubble-id`:

```bash
./build/panvar inspect \
  -i tests/real_data/c4.gfa \
  --bubble-prefix-in tests/results/c4/bubble \
  -o tests/results/c4/inspect/all
```

3. Allele (default walk clustering, node-length weighted):

```bash
./build/panvar allele \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/allele \
  --bubble-prefix-in tests/results/c4/bubble
```

4. Allele (predefined path-group JSON):

```bash
./build/panvar allele \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/allele_json \
  --bubble-prefix-in tests/results/c4/bubble \
  --clusters-json tests/real_data/c4.clusters.json
```

5. Call (standard):

```bash
./build/panvar call \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/call \
  --bubble-prefix-in tests/results/c4/bubble \
  --allele-prefix-in tests/results/c4/allele \
  --reference-path grch38#1#chr6:31891045-32123783
```

6. Call (debug dotplots + INS classification + pangene copy annotations):

```bash
./build/panvar call \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/call_pangene \
  --bubble-prefix-in tests/results/c4/bubble \
  --allele-prefix-in tests/results/c4/allele_json \
  --reference-path grch38#1#chr6:31891045-32123783 \
  --classify-ins \
  --pangene-bed tests/real_data/c4.pangene.bed.gz \
  --pangene-gene-match C4 \
  --pangene-tune-ins \
  --dotplot-gtf tests/real_data/gencode.v49.annotation.gtf.gz \
  --dotplot-gene-match C4 \
  --debug
```

7. Describe:

```bash
./build/panvar describe \
  -i tests/real_data/c4.gfa \
  --bubble-prefix-in tests/results/c4/bubble \
  --bubble-id 1 \
  --out-dir tests/results/c4/describe \
  --kmer-size 31
```

## Useful Options

- `allele`
  - `--bubble-prefix-in <prefix>`: consume bubble outputs using prefix convention
  - `--clusters-json <path>`: use predefined path->cluster labels
  - `--similarity-out-dir <dir>`: per-bubble distance matrices, cluster-separation stats, and heatmap helper
  - `--skip-no-reference-bubbles --reference-path <path>`: drop bubbles not traversed by the reference
- `bubble`
  - `--merge-nearby-bp <N>`: optionally fuse nearby bubble candidates by graph bp distance
- `inspect`
  - `--bubble-prefix-in <prefix> [--bubble-id <N>]`: write bubble path FASTA plus node traversal counts
  - `scripts/plot_node_coverage_heatmap.R`: plot inspect node-count tables as path-by-node coverage heatmaps
- `call`
  - `--bubble-prefix-in <prefix>` + `--allele-prefix-in <prefix>`: simplified module handoff
  - `--debug --debug-out-dir <dir>`: per-cluster FASTA/PAF/dotplot/VCF
  - debug status manifests: `<debug>/debug_summary.tsv`, `bubble_<id>/bubble_status.tsv`, `bubble_<id>/cluster_status.tsv`
  - `--minimap-preset asm5|asm10|asm20` (default `asm20`)
  - `--classify-ins`: INS NOVEL vs DUP-like subtype/CN
  - `--pangene-bed <bed(.gz)>`: event-level gene-copy delta annotation
  - `--pangene-gene-match <expr>` (repeatable): filter pangene genes
  - `--pangene-tune-ins`: set `INS_SUBTYPE=DUP_PANGENE` when pangene copy gain supports duplication
- `describe`
  - `--bubble-prefix-in <prefix>`: consume bubble outputs using prefix convention
  - `--kmer-size <K>`: canonical 2-bit k-mer size (`1..31`, default `31`)
  - `--feature-mode all|minimizer|syncmer`: optionally sample k-mers to reduce feature count
  - `--max-wide-features <N>`: skip very wide matrices above a feature cap
  - `--no-wide-matrix`: keep only compact graph-aware feature map plus sparse JSONL counts

## Documentation

- [docs/README.md](docs/README.md)
- [docs/modules/bubble.md](docs/modules/bubble.md)
- [docs/modules/inspect.md](docs/modules/inspect.md)
- [docs/modules/allele.md](docs/modules/allele.md)
- [docs/modules/call.md](docs/modules/call.md)
- [docs/modules/describe.md](docs/modules/describe.md)
- [docs/modules/example.md](docs/modules/example.md)

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

## Repository Layout

- `src/*_command.cpp`, `include/panvar/*_command.hpp`: CLI wrappers for each subcommand
- `src/cli_utils.cpp`, `include/panvar/cli_utils.hpp`: shared command-line helpers
- `src/graph_utils.cpp`, `include/panvar/graph_utils.hpp`: shared graph/path/sequence helpers
- `src/`, `include/panvar/`: module algorithms and public data structures
- `scripts/plot_distance_heatmap.R`: ggplot2 helper for plotting allele similarity distance matrices
- `scripts/plot_node_coverage_heatmap.R`: ggplot2 helper for plotting inspect path-by-node coverage matrices
- `external/minimap2/`: minimap2 C library source
- `tests/real_data/`: bundled example loci and inputs
- `docs/modules/`: module docs

## Output Paths

If you pass output prefixes/files inside non-existing directories, `panvar` now creates parent directories automatically.
