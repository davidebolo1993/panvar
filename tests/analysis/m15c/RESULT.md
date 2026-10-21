# M15-C result: user-facing interface reviewed — one authoritative guide, corrected help and docs, portable complete profiles, every documented example tested

**Run.** Branch `m15-c-interface-review`; code commit `66da796`; documents committed after.
- **Results:** `tests/results/m15c_interface_review_20261021/`.
- **Budget:** 305.6 of 1,800 s, including 150 s of disclosed development (two help rebuilds and a CLI-contract run).
- **Units:** 8 of 8 passed.

**Not changed:**
- no model, search, threshold or default change;
- no reserved family touched;
- no merge, push or tag.

**Deliverables:**
- `docs/guides/genotype_short_read.md`: the authoritative usage guide. It covers:
  - inputs and orientation;
  - preset, profile and precedence;
  - every flag and profile key, with values, defaults, refusals and whether it affects predictions, follows a route or is
    report-only;
  - resource guards, cache budgets and work caps, kept separate;
  - reporting;
  - outputs;
  - a decision guide for filters, ambiguity, audit and support, with the sv:9 example;
  - development-only interfaces;
  - the evidence per locus;
  - the tested commands.
- **Corrected help** (`src/genotype_command.cpp`; strings only):
  - `sequence_first_v2`, its conditional FASTA and its 0/3/1 exits;
  - the dosage default by schema, with v2 refusing repair;
  - PASS is not calibrated reliability (the old "exact genotype =" policy is removed);
  - bounded refinement, not a global optimum;
  - profile semantics;
  - a pointer to the guide.

  `tests/genotype_cli_contract.cpp` now requires these statements and forbids the old claims.
- **Corrected module doc:** `docs/modules/genotype.md` (usage-policy passage, dosage default, profile, repair).
- **Portable complete profiles:** `profiles/genotype_experimental_sr_v1_{c4,cyp2d6,acot}.tsv`. Each is the full tested
  configuration (base + performance + M10-A, with locus-specific neighbour limits), with no absolute paths and no
  archived-layout dependency.
- `RELEASE_CHECKLIST.md`.

**Binary:** the help changes rebuilt `panvar` to sha256 `5197fd975e69a9e7b9f0b4ba51f07f08e9ada4301bf3b945d8e7f046a411d8bc`.
Its sources differ from the release candidate's `47042a5c…` only in help strings.

**Verification** (development controls only):

| example | result |
|---|---|
| standard run with the committed complete profiles (C4 NA19705, CYP2D6 HG02769, ACOT HG03209) | exit 0; **predictions identical to the archived M11-D runs** (rank-1 score, block paths, haplotype SHA-256, status, model fingerprint); **every effective setting equal**, except the archived runs' declared time guard and diagnostics key |
| off-panel FASTA (C4) | `.sequence_haplotypes.fa` with both records, SHA-256 equal to the table; legacy FASTA suppressed |
| panel recovery without FASTA (ACOT, both haplotypes PANEL) | no FASTA; `extract_panel_haplotype.py --check` passes |
| resource-adjusted profile (C4: time 120 s, footprint 4 GiB, emission cache 256 MiB) | exit 0; prediction unchanged; effective values recorded |
| binding time limit (ACOT, 8 s) | exit 124, `NO_PREDICTION` / `TIME_LIMIT`; never reported as a completed prediction |
| 7 invalid combinations (v2 + dosage repair; early_incumbent without v2; repair off with rounds; a duplicate key; an unknown key; an unnamed profile; `--max-candidates`) | all exit 2 with the reason named; **nothing written at the prefix** |
| CLI contract test (extended) | passed |

**Stopped for Davide's interface review.**
