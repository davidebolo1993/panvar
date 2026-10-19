# M3-A: copying-prior semantics and fixed-candidate audit

## Ancestry

Branch `m3-a-copying-prior` from `6d85cce`.
- `6d85cce` is the M1-E head and descends from the frozen model-correctness baseline `282e64e` (M1-A baseline; S2-L
  combined route).
- `282e64e..6d85cce` changes no library code (`src/` and `include/` are identical). It changes only the diagnostic
  probe and the M1 analysis tools (`score-set`, `--decompose`, `--contrast`, `block-sequences`, `edit-distances`).
- **No M2 code is carried.** The M2 branches descend from `6d85cce`; nothing from them is merged here.

## Scope

- No production prior or default change.
- No genotyping campaign, new sequence construction or long-read work.
- Diagnostic tools only: the probe subcommand `prior-pool`, and `tests/m3a_prior_representation.cpp`, which includes
  `tests/genotype_mosaic_core.cpp` unchanged to reuse its brute-force oracles.

## 1. What the current prior computes (static trace; established by reading code)

**The transition** (`src/ls_hmm.cpp:40-77`, `100-180`). It is Li–Stephens copying: T = (1−r)I + (r/n)·11ᵀ between
consecutive blocks.
- **r is a donor *redraw* probability.** The new donor is drawn uniformly over all n panel templates, the current one
  included.
- Actual donor change is r(n−1)/n; staying on the same donor is (1−r) + r/n.
- An allele change happens only when the redrawn donor carries a different allele at the next block.
- **0.01 is a declared setting** (`switch_probability`, default `0.01` at `src/genotype_command.cpp:2273`, allowed range
  (0, 0.5]). It is not a calibrated biological crossover probability.

**Where r applies.**
- At each of the B−1 boundaries between consecutive blocks, whatever the length of the invariant segment between them.
- There is no transition at the locus ends.
- The start is uniform, 1/n per template.
- Invariant segments must match exactly. They are hard constraints and carry no transition or weight.
- **Consequence:** a boundary's switch opportunity depends on how the graph is cut into blocks, not on distance. A
  neutral boundary (a monomorphic block) composes as 1 − (1−r)², so r = 0.01 becomes 0.0199. The existing counterexample
  (d) in `genotype_mosaic_core.cpp` shows this.

**Emission.** It is hard: a template can copy an allele only if it carries it.
- Block alleles are canonicalised by bytes (`src/genotype_model.cpp:255-285`). Source-ID aliases with identical bytes
  become one allele.
- Templates are **not** deduplicated (`src/genotype_model.cpp:303-343`); there is one row per template name. Panel
  multiplicity is therefore the template count: an allele's weight at the start or after a redraw is
  (carriers / n).
- A duplicated template entry and a genuine second carrier are indistinguishable to the model.

**Summation.**
- Hidden template histories are summed by the forward recursion (`log_haploid_mosaic_prior`).
- The sequence prior additionally sums every allele path that spells the same bytes (`log_haploid_sequence_prior`, a
  forward recursion over (block, offset, template)).
- The unordered pair gets log 2 decided on bytes (`log_unordered_sequence_pair_prior`).

**Caller stages.** Every stage builds its vector as `options.switch_probability` repeated B−1 times, so they all use the
same constant transition vector.

| stage | prior form | where |
|---|---|---|
| proposal beam (grouped; production) | same transition (stay (1−r)+r/n, move r/n). Retains allele-prefix pairs by summed history mass, then best history, scored with **proposal** emissions, not fragment likelihoods | `src/mosaic_beam_grouped.cpp:543-600`, `better_pair` |
| early incumbent and progressive ordering | proposal emissions plus the exact **path** prior | `genotype_command.cpp:3188`, `3606` |
| union ranking (selection) | fragment + exposure + **sequence-pair** prior | `genotype_command.cpp:3443`, `4190` |
| block projection of a sequence | the member allele path with the **highest path prior** | `genotype_command.cpp:4194-4199` |
| coordinate repair (neighbours, ranking) | sequence-pair prior | `src/genotype_repair.cpp:258`, `709`, `957` |
| dosage repair | sequence-pair prior | `src/genotype_dosage_repair.cpp:295`, `403`, `1484` |

## 2. Representation tests (unit `repr`; independent oracles reused)

- **R0.** The existing mosaic-core prior fixtures: brute-force histories, sequence-equivalent aggregation, and
  counterexamples (a)–(d).
- **R1. Neutral subdivision.** A monomorphic block inserted between two blocks preserves sequences and information.
  - The admissible sequence set must be unchanged.
  - Constant model: subdivided equals undivided at r_combined = 1 − (1−r1)(1−r2), within 1e-13, against production and
    the brute oracle, for r pairs (0.01, 0.01), (0.001, 0.05) and (0.2, 0.3).
  - The same check at sequence level: a shared segment written as invariant context vs as a one-allele block.
  - The distance model, with fixed L = 50 bp, must be invariant (1e-13).
- **R2. A non-neutral edit.** Splitting a polymorphic block admits new recombinants: 2 admissible sequences become 4.
  It must not be called neutral.
- **R3. Summed vs best history** at r = 0.01, over 400 random small panels (brute oracle). Production must equal the
  summed oracle. Order flips are counted.
- **R5. Aliases.**
  - A source-ID alias with identical bytes changes nothing (exact).
  - A duplicated template row changes the prior; it is reported as indistinguishable from a genuine carrier.
- **R6. Redraw semantics.** Allele-change probability against the analytic expression; donor-change and stay
  probabilities are reported.

## 3. Fixed-pool prior comparison (units `pool.<cell>`; truth-free pools)

**Pools.** `POOLS.tsv` freezes, for each of seven cells, every sequence candidate a saved caller run scored
(`<run>.mosaic_scores.tsv`, with allele paths from `<run>.mosaic_candidates.tsv`), by SHA-256.
- **Cells:** C4 LZO control, C4 mosaic 2× (M1-A SEARCH error), C4 NA19650 10× sv:7 (M1-B fresh baseline caller), CYP2D6
  label tie, CYP2D6 LOO sv:9 5× (M1-A RANKING error), ACOT 5× and 30×.
- **Scores.** Fragment and exposure terms are read from the run and never recomputed. Only the sequence-pair prior
  changes.
- **Gate.** The current prior, recomputed, must reproduce every saved prior (and total) before anything is reported.

**Prior models** (declared now; no other is evaluated):
- `const_r` for r ∈ {0.001, 0.003, 0.01 (current), 0.03, 0.1}.
- `zero` (r = 0, diagnostic): only sequences that one template spells whole remain finite.
- `dist_matched`: r_b = 1 − exp(−d_b/L).
  - **Coordinates:** d_b is the panel-mean span in bp from the start of block b−1 to the start of block b, i.e. the mean
    over templates of len(allele at b−1) plus len(invariant segment b).
  - **Zero span:** d_b = 0 gives r_b = 0, and such boundaries are counted.
  - **Scale:** L = Σd_b / ((B−1) · −ln 0.99), so the total hazard equals the current model's. This separates
    redistributing switch opportunities from changing prior strength.
  - The model is exactly subdivision-consistent, and L is fixed per locus from panel geometry before any score is read.

**Reported per cell and model:**
- the winner, whether it is the saved call, and the call's rank;
- winner minus call, as total, fragment, exposure and prior;
- best intact pair minus best mosaic (intact means both haplotypes are whole panel templates);
- the winner's probability within the pool, the count within 1 nat, and the gap to second.

**Frozen rule.** A cell is SENSITIVE to a factor when the pool winner differs from the `const_0.01` winner:
- **strength:** the other constants;
- **placement:** `dist_matched`;
- **zero redraw:** `zero`, diagnostic only.

**What this does not show.** Fixed-pool results say nothing about what a fresh search would find. This is not a
calibration or a default selection.

**Evaluator-only rows** (`*.counterfactual.tsv`) go to a separate table:
- the C4 2× truth diplotype, with fragment and exposure from M1-A;
- the M1-B sv:7 alternative chosen by a rule that used the true unit class.

## 4. Flank evidence without double-counting

1. **Spanning molecules (directly informative).** Fragment-side classes come from the **existing** M1-B decompositions:
   block interior, junction-spanning, and outside the block on all haplotypes.
   - Available for C4 sv:7 10× (reported vs saved 2+2) and CYP2D6 sv:9 5× (reported vs saved 2+3).
   - They are descriptive. M1-C later showed the flagged mass mismatch was a reduction artefact.
   - There is none for C4 2×, which is not rescored here.
2. **Donor association through the copying prior.** For each contrast (call vs alternative):
   - Δprior(full path) is split into Δprior(target-marginal) and Δflank association.
   - The target-marginal prior treats every non-target block as monomorphic (transitions kept, emission constraint
     removed), so `log P(full) − log P(target-marginal)` is the flank alleles' association.
   - This is path-level and uses the same reads in nothing: the prior has no read input.
3. **Competing origins outside the locus.** These enter only through the background term (η · P_bg) and through
   fragments placed outside the target on every haplotype (whose delta is reported).
   - Off-locus paralog placement is not modelled, and it is not measured here.

**Proper flank ablation: specified for review, NOT run.**
- **Procedure:**
  - truncate every pool haplotype to the target window plus a declared flank margin;
  - re-score the unchanged fragment model on the **truncated** haplotypes, so exposure is recomputed on the truncated
    lengths and reads that need the removed flank fall to background;
  - compute the prior over the window's blocks only;
  - compare rankings with full and truncated flanks at two declared margins.
- **Why not simply drop reads.** Dropping flank reads while keeping the full-length exposure would create an exposure
  deficit, which is not a valid ablation.
- **Data separation.** No read set may serve both as an empirical prior and as likelihood.
- **Cost:** one score-set per cell per margin, with the reads.

## 5. Units and budget

**Units:**

| unit | declared timeout |
|---|---|
| `repr` | 30 s |
| `pool.<cell>`, C4 and CYP2D6 cells | 60 s each |
| `pool.acot_5`, `pool.acot_30` | 120 s each |
| `summary` | 30 s |

**Limits.** 600 s cumulative, including the 0.91 s M2 closure unit run under this brief. The check is spent time plus the
declared timeout. At most 4 threads, 6 GiB, one process at a time.
