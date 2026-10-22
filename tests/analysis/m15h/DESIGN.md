# M15-H: help organisation and the public option inventory (frozen before the recorded run)

**Logistics:**
- **Branch:** `m15-h-help`, from `m15-g-assessment` `e1b900d`.
- **Results:** `tests/results/m15h_help_20261022/`.
- **Budget:** 1,800 s, including 90 s of disclosed development.

**Not changed:** no model, search, filter-rule or default. No flag is added, removed or renamed.

## Changes

**One option inventory in code** (`genotype_option_inventory()`). Every option has:
- a group (everyday; then 1–6 by purpose: input and panel; copying prior and marker / placement model; candidates and
  bounded refinement; resources; reporting; legacy and development);
- aliases, the profile key, accepted values and units;
- what it can change (inference, computation, reporting, completion);
- a status (everyday, supported, expert profile only, legacy, development).

**The effective short-read default is computed** from the preset and the built-in options:
- public keys through `genotype_setting_value`;
- profile-only keys through `experimental_setting_value`, the settings table's own mapping, extracted unchanged so that
  `.settings.tsv` is byte-identical.

**The basic `--help`** leads with:
- the ordinary no-profile command and automatic library estimation;
- the everyday options (graph, reads, output, reference, threads, `--dosage-check`, the two index flags, quiet and help);
- outcomes (completed, provisional, absent), the FASTA rule, and the primary `assessment_v1` view with its caveats.

The old headline fields are named only as superseded compatibility diagnostics.

**`--help-advanced`** is generated from the inventory, grouped by purpose. It shows the model-derived neighbour policy
and the terms (donor redraw, initial mosaic candidates, bounded refinement, `min_variant_bp`, library scope, conditional
support, prior-dependent, locus-level phase, coverage).

**The guide** (§2–§5 restructured, §10 added):
- §3 carries the generated inventory between markers, and the CLI-contract test fails if it differs from the parser and
  the preset;
- the terms; recommended future removals and aliases (not implemented).

**The module page** has its defaults corrected.

## Verification units

| unit | expectation |
|---|---|
| `cli_contract` | the new help organisation (order and content), the corrected claims, removal of the headline explanation from the basic help, the six advanced groups and the terms, inventory coverage of every public setting and the profile keys, inventory defaults equal to the parser and the preset, and **the guide's inventory equal to the generated one** |
| `ex.help_text` | both help texts captured (the deliverable) |
| `ex.preserved_c4` | the no-profile C4 command: every output identical to M15-G's run of the same command, apart from the run's name and timing and resource rows. This covers settings, sequences, scores, assessment and the old headline columns |
| `ctest.full` | the full enabled suite on the final binary |
