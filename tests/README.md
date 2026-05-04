# Tests

## Data

- `real_data/`: bundled loci (`*.gfa`) and precomputed snarls (`*.snarls.jsonl`)
- `results/`: local outputs (kept out of git except `.gitkeep`)

## Smoke Test

Run CTest after build:

```bash
ctest --test-dir build --output-on-failure
```

Run directly:

```bash
tests/smoke.sh ./build/panvar tests/real_data/c4.gfa /tmp/panvar_smoke_c4
```
