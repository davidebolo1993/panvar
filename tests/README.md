# Tests

## Data

- `real_data/`: bundled small example dataset (C4: `c4.gfa.gz`, `c4.snarls.jsonl`)
- `results/`: local outputs (kept out of git except `.gitkeep`)

## Smoke Tests

Two smoke tests:

- `synthetic_smoke.sh` — fast, dependency-free; exercises every `call` event type on tiny
  hand-built graphs with exact-record assertions. This is the one registered with CTest.
- `real_smoke.sh` — fuller integration run of `bubble -> inspect -> describe -> panphorte -> call`
  on a real locus; checks the main handoff files plus compressed inspect/describe outputs. Not part
  of CTest (large/slow); run it manually.

Run CTest after build (runs the synthetic smoke):

```bash
ctest --test-dir build --output-on-failure
```

Run either directly:

```bash
tests/synthetic_smoke.sh ./build/panvar tests/synthetic_data /tmp/panvar_smoke
tests/real_smoke.sh ./build/panvar tests/real_data/c4.gfa.gz /tmp/panvar_smoke_c4
```
