# panvar

`panvar` is a modular C++ toolkit for pangenome-graph analysis from GFA.

Modules:

1. `bubble` (Module 1): site extraction/refinement from precomputed `vg snarls`
2. `allele` (Module 2): allele extraction and clustering from module-1 bubbles
3. `call` (Module 3): SV calling (`INS`/`DEL`/`INV`) from clustered alleles
4. `describe` (Module 4): per-bubble haplotype feature tables for downstream association

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
  --snarls-in tests/real_data/c4.snarls.jsonl
```

2. Allele (default clustering):

```bash
./build/panvar allele \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/allele \
  --bubbles-csv-in tests/results/c4/bubble.bubbles.csv
```

3. Allele (predefined path-group JSON):

```bash
./build/panvar allele \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/allele_json \
  --bubbles-csv-in tests/results/c4/bubble.bubbles.csv \
  --clusters-json tests/real_data/c4.clusters.json
```

4. Call (standard):

```bash
./build/panvar call \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/call \
  --bubbles-csv-in tests/results/c4/bubble.bubbles.csv \
  --clusters-csv-in tests/results/c4/allele.allele_clusters.csv \
  --assignments-csv-in tests/results/c4/allele.allele_assignments.csv \
  --reference-path grch38#1#chr6:31891045-32123783
```

5. Call (debug dotplots + INS classification + pangene copy annotations):

```bash
./build/panvar call \
  -i tests/real_data/c4.gfa \
  -o tests/results/c4/call_pangene \
  --bubbles-csv-in tests/results/c4/bubble.bubbles.csv \
  --clusters-csv-in tests/results/c4/allele_json.allele_clusters.csv \
  --assignments-csv-in tests/results/c4/allele_json.allele_assignments.csv \
  --reference-path grch38#1#chr6:31891045-32123783 \
  --classify-ins \
  --pangene-bed tests/real_data/c4.pangene.bed.gz \
  --pangene-gene-match C4 \
  --pangene-tune-ins \
  --dotplot-gtf tests/real_data/gencode.v49.annotation.gtf.gz \
  --dotplot-gene-match C4 \
  --debug
```

6. Describe:

```bash
./build/panvar describe \
  --vcf-in tests/results/c4/call_pangene.region.vcf \
  --out-dir tests/results/c4/describe \
  --gtf tests/real_data/gencode.v49.annotation.gtf.gz \
  --gene-match C4 \
  --size-bins 100,1000
```

## Useful Options

- `allele`
  - `--clusters-json <path>`: use predefined path->cluster labels
  - `--similarity-out-dir <dir>`: per-bubble clustering diagnostics
- `call`
  - `--debug --debug-out-dir <dir>`: per-cluster FASTA/PAF/dotplot/VCF
  - `--minimap-preset asm5|asm10|asm20` (default `asm20`)
  - `--classify-ins`: INS NOVEL vs DUP-like subtype/CN
  - `--pangene-bed <bed(.gz)>`: event-level gene-copy delta annotation
  - `--pangene-gene-match <expr>` (repeatable): filter pangene genes
  - `--pangene-tune-ins`: set `INS_SUBTYPE=DUP_PANGENE` when pangene copy gain supports duplication

## Documentation

- [docs/README.md](docs/README.md)
- [docs/modules/bubble.md](docs/modules/bubble.md)
- [docs/modules/allele.md](docs/modules/allele.md)
- [docs/modules/call.md](docs/modules/call.md)
- [docs/modules/describe.md](docs/modules/describe.md)

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

Build:

```bash
docker build -t panvar:latest .
```

Run:

```bash
docker run --rm -v "$PWD":/work -w /work panvar:latest panvar --help
```

## Repository Layout

- `src/`, `include/panvar/`: C++ implementation
- `external/minimap2/`: minimap2 C library source
- `tests/real_data/`: bundled example loci and inputs
- `docs/modules/`: module docs
