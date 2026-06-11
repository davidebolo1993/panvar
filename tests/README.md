# Tests

## Data

- `real_data/`: bundled small example dataset (C4: `c4.gfa`, `c4.snarls.jsonl`)
- `results/`: local outputs (kept out of git except `.gitkeep`)

## Smoke Test

The smoke test runs `bubble -> inspect -> describe -> call` and checks
the main handoff files plus compressed inspect/describe outputs.

Run CTest after build:

```bash
ctest --test-dir build --output-on-failure
```

Run directly:

```bash
tests/smoke.sh ./build/panvar tests/real_data/c4.gfa /tmp/panvar_smoke_c4
```
