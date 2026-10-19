# M3-F: bounded validation beyond the original joint-search recovery

**Ancestry.** Branch `m3-f-validation` from `f31b9e6` (M3-E.1). There is no algorithm change in this batch: the code is
the M3-E.1 library, caller and probe, built once, and every run uses that one binary.

**Scope.** Development validation on cells from the frozen campaign. This is not independent held-out validation and
not a default-promotion gate.

## 1. Frozen settings (identical within each pair)

**Base.** Profile `genotype_c4_v6_campaign_local.tsv` (sha `74a3723b…`), 4 threads, the cell's own campaign inputs.

| arm | profile | differs by |
|---|---|---|
| `off` | `profiles/support_v2_off.tsv` | v6 + `block_support_rule full_data_likelihood_v2_experimental` |
| `on` | `profiles/support_v2_on.tsv` | the off profile + `repair_joint_search experimental_v1` (exactly one key) |

**Unchanged:** k = 4, both 2,048 bounds, the one-pass schedule, the continuation rules, caps and budgets.

## 2. Deterministic selection (from the frozen campaign `genotype_c4_production_11fe681`)

**Excluded:** the M3-C/D cells (`one_homologue_once_1` at 2× and 5×, LZO HG00096 at 5×).

### Low-coverage mosaic error cells
**Rule.** Plan families are the plan ID without its replicate suffix: `both_homologues`, `dense_switching`,
`one_homologue_once`, `panel_control`, `three_switches`. Families are taken in alphabetical order, and within a family
the replicate ID ascending, then coverage ascending. From each family, take the first cell with
`excess_over_floor_total > 0`, at most one per family, until two cells are chosen.

| family | outcome |
|---|---|
| `both_homologues` | every cell has excess 0: none |
| `dense_switching` | **`dense_switching_1` at 2×** (excess 81) |
| `one_homologue_once` | `_1` is excluded (M3-C/D); `_2` has excess 0 everywhere: none |
| `panel_control` | excess 0 everywhere: none |
| `three_switches` | **`three_switches_1` at 2×** (excess 8) |

Not selected: `dense_switching_2` at 2× (10) and `three_switches_2` at 2× (7), as second replicates.

### Correct higher-coverage twins
**Rule.** The lowest higher coverage of the same plan with an **exact sequence call**: `excess_over_floor_total` 0 and
status COMPLETE, which gives a unique, phased sequence call.

**Unavailable for both selected plans.** Every excess-0 cell of `dense_switching_1` (10× and 30×) and of
`three_switches_1` (5×, 10× and 30×) is CONDITIONAL_GENOTYPE_EQUIVALENT: an unphased genotype with no FASTA.

**Declared replacement rule.** Use the lowest higher coverage with `excess_over_floor_total` 0 whatever the status. These
are labelled **genotype-equivalent twins (not exact sequence calls)** and evaluated on block genotypes and excess only.

- `dense_switching_1` at **10×** (GE; campaign 123.3 s at 8 threads)
- `three_switches_1` at **5×** (GE; campaign 47.5 s)

### Ordinary LZO control
**Rule.** The same coverage as the existing control (5×): the first sample by ID, excluding HG00096, with
`excess_over_floor_total` 0 and status COMPLETE. That is **HG00171 at 5×** (campaign 16.5 s).

### LOO structural-error case
**Rule.** From the existing attribution set (M1-B `AUDIT_SET.tsv`, C4 cases `c4_10` and `c4_30`), the lowest coverage:
**NA19650 at 10×**.
- Campaign: excess 65,506, availability floor 113, COMPLETE, 62.2 s.
- M1-B fresh run: 48.2 s.

### Outside the experiment
CYP2D6 and ACOT exceed the frozen search bounds (M3-D eligibility), so they are not run. A skipped pass is never counted
as cross-locus validation.

## 3. Run order and stop rule

1. **Controls first:** HG00171 5× (off, on), `three_switches_1` 5× (off, on), `dense_switching_1` 10× (off, on).
2. **`control.check`** (truth-free) compares each control's off and on predictions: status, called sequences, block
   genotypes, final score. Its digest freezes the control predictions.
3. **`eval.<control>.<arm>`** (truth), then **`control.gate`** (`accuracy_m3f.py control`).
   - **Stop rule:** if any control's on arm is HARMFUL (higher excess, fewer exact blocks or fewer exact SV blocks),
     MIXED, a refusal or no-call in the on arm only, or a status downgrade (COMPLETE to anything else, or GE to
     incomplete or no-call), the gate exits 1.
   - The runner then refuses every error-cell run. The evidence is preserved and the mechanism reported.
4. **Error cells:** `dense_switching_1` 2× (off, on), `three_switches_1` 2× (off, on), NA19650 LOO 10× (off, on).
5. **`compare`** (truth-free, every cell): freezes all predictions with a digest.
6. **`eval.<error>.<arm>`** (truth), then **`accuracy`**.
   - **Outcomes:** RECOVERED, IMPROVED, HARMFUL, MIXED, CHANGED_EQUAL, UNCHANGED, REFUSAL_OR_NO_CALL; plus a separate STATUS_DOWNGRADE flag.
7. **Diagnostics only** (truth-guided, after the freeze; they never influence predictions or settings), for every error
   cell whose on arm still has positive excess, changed or not:
   - **joint pass not complete** (`joint_search_status` not `complete*`): PASS_SKIPPED_OR_INCOMPLETE, no probe;
   - **mosaic cells:** `two-block-lik` (the frozen rule: k = 4, both 2,048 bounds, repair configuration) from the **off**
     run's incumbent, which is where the on arm's pass started, then `two-block-eval` against the plan's truth;
   - **the LOO case:** the evaluator's pool fields (is the availability floor attained by a scored candidate?). M1-B's
     existing `block-sequences` output (run-independent, hashed in INPUTS.tsv) with M1-B `structure.py` gives sv:7 copy
     allocation for both arms. structure.py's "historical" column is the off arm and its "fresh" column the on arm.
   - **Classes, in precedence** (`attribute_m3f.py`):
     1. PASS_SKIPPED_OR_INCOMPLETE;
     2. EVALUATED_ALTERNATIVE_LOSES: the truth is in the catalogue but ranks below first, or (LOO) a floor-attaining
        candidate was scored and lost;
     3. OUTSIDE_TWO_BLOCK_NEIGHBOURHOOD: the truth differs from the incumbent on both homologues, or at more than two
        blocks of one;
     4. ALTERNATIVE_ABSENT_FROM_SHORTLIST: a truth allele was not kept in the top 4;
     5. UNRESOLVED: anything else, including a truth that ranks first while the call stays wrong, and a LOO floor that no
        scored candidate attains.

   No class claims exhaustive coverage.

## 4. Compared, per cell (resource refusals and no-calls stay in every denominator)

- **Selection:** the selected unordered sequence pair and block projections; the full score and its fragment, exposure
  and prior terms.
- **Search:** the joint decision, move, continuation, ties and skips.
- **Accuracy, kept separate:**
  - sequence excess and distance (evaluator);
  - structural outcomes (SV-block vs context distance, SV blocks exact; copy allocation for the LOO case);
  - genotype-data support (v2 labels) and phase support.
- **Cost:** whole-caller wall time and peak RSS (ledger); joint-pass seconds, retained-census bytes and footprint
  (status).
- **Harmful moves are reported as prominently as recoveries.** A higher score is not by itself a better genotype.

## 5. Runtime projection and budget (1,200 s, counting every attempt)

M3-D measured the current caller at 4 threads at 0.65–1.14× the campaign's 8-thread time. The joint pass took 17–31 s at
2–5×.

| cell | projection, off + on (s) |
|---|---|
| HG00171 5× | 15 + 40 |
| `three_switches_1` 5× | 50 + 90 |
| `dense_switching_1` 10× | 125 + 200 |
| `dense_switching_1` 2× | 30 + 60 |
| `three_switches_1` 2× | 30 + 60 |
| NA19650 10× | 60 + 120 |
| evaluations | 12 × ~10 |
| diagnostics | ~150 |

**Total ≈ 1,100 s.** M3-D/M3-C actuals suggest closer to 900 s: off arms ran at 0.4–1.0× campaign time and evaluations
took about 10 s.

**Declared timeouts (s).**

| unit | off | on |
|---|---|---|
| HG00171 5× | 60 | 120 |
| `three_switches_1` 5× | 100 | 160 |
| `dense_switching_1` 10× | 180 | 240 |
| `dense_switching_1` 2× | 70 | 120 |
| `three_switches_1` 2× | 70 | 120 |
| NA19650 10× | 120 | 180 |

- Each evaluation: 40.
- Each comparison, accuracy or attribution step: 30.
- `two-block-lik`: 120.
- `two-block-eval`: 40.
- structure: 30.

**Rules.**
- A unit is refused when the time spent plus its declared timeout exceeds 1,200 s.
- A refused, timed-out or failed unit is recorded and never substituted or repeated.
- The caller's own 4 GiB footprint guard (v6) stays in force: a refusal is a result, and caps are not raised.

**Manifests** (frozen by `m3f_run.py manifest C|E|D` before the first unit of each):
- C: controls;
- E: error cells, refused unless `control.gate` exited 0;
- D: diagnostics.

**Inputs.** Every input is hashed in INPUTS.tsv and verified before each caller run.
