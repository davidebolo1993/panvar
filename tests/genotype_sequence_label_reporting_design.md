# Sequence, block-genotype and representation-ambiguity reporting: prospective design (no caller change yet)

**Status:** a design note, written 2026-10-15 after the CYP2D6/D7 pilot (sealed diagnostic baseline). It is not implemented. Sections 1–5 are version 1. **Section 6 is version 2 (S0 v2, 2026-10-15): the sequence-first output and confidence contract.** Where they differ, section 6 supersedes; nothing in sections 1–5 is withdrawn except where section 6 says so.

**Gate for implementation:** caller behaviour stays unchanged until the adversarial fixtures below exist, their current-behaviour baseline is recorded, and this contract is accepted. Implementation is then a versioned output-contract change with its own baseline; historical freezes keep their recorded implementation.

**Related:**
- review finding C4 (`docs/reports/genotype-shortread-code-review-20261014.md`, cross-block spelling degeneracy);
- review finding C2 (projection-aware support and GQ);
- the pilot's `degeneracy/` inventory and AMENDMENT_5 reporting-layer classification.

## 1. Problem

A haplotype is reported as one allele per block, and its sequence is the concatenation of those alleles. That map need not be injective.
- **In CYP2D6:** a 59 bp tandem unit can be counted in sv:3 or in sv:5, across the 105 bp ctx:4. One sequence then has several block labellings.
- **Current caller behaviour:** when one winning sequence group carries different labellings, the caller refuses to invent labels, which is right, but it turns this into a whole-locus no-call.
- **What is lost:** the sequence was identified (7 of 8 such pilot no-calls spelled the true sequence), and so was every block outside the ambiguous interval.

**Separate from this, and to be preserved:** a genuine tie between distinct sequences is a different situation. In the pilot, mosaic three_2 at 30× tied eight distinct sequence candidates on total score. It must stay unresolved at the sequence level.

## 2. Three reporting units

Each unit is reported separately. None is derived by picking a member of a set.

| unit | identified when | otherwise |
|---|---|---|
| **diploid sequence**: the unordered pair of haplotype sequences (hash, length, FASTA) | every member of the decision set spells the same unordered pair | `SEQUENCE_UNRESOLVED` (distinct sequences tied, or not separated within the declared tolerance) |
| **block genotype**: an unordered allele pair at one block | across every member of the decision set and every admissible projection of each member's sequence pair, exactly one unordered allele pair occurs at that block | `REPRESENTATION_AMBIGUOUS` if the alternatives come only from alternative segmentations of one identified sequence; `GENOTYPE_UNRESOLVED` if they come from distinct sequences |
| **representation interval**: a minimal compound interval of consecutive blocks, i..j, over which segmentation is non-unique | always reported when present | lists every admissible label tuple and the identified interval sequence per haplotype (hash, length). It never chooses a tuple |

**Decision set:** the caller's existing pooled top set, under its current tolerance and tie rules. That rule is unchanged.

**Projection:** an allele tuple whose concatenation equals the member's haplotype sequence. A member can have several projections only inside representation intervals.

**Representation intervals come from the model, not from reads:**
- They are computed at preparation by an exhaustive overhang search over block alleles, as in `tests/analysis/cyp2d6_pilot/degeneracy_inventory.py`, and bound to the model fingerprint.
- The per-candidate projection set is the product, over intervals, of the label tuples consistent with that candidate's interval sequence.
- **Blocks outside every interval have exactly one projection.**

## 3. Rules

1. **No arbitrary label.** A representation-ambiguous block reports the admissible set and no allele pair.
2. **No silent merge.** Block identifiers and the block table keep their meaning. A compound interval is reported in addition, as its own row, not by replacing blocks. A model that would coalesce blocks is a different model, with a different fingerprint, and needs its own validation.
3. **Unambiguous blocks are retained.** A block identified under rule 2.2 is reported even when another block or the sequence is unresolved. Its quality fields say they are conditional on the decision set.
4. **Distinct-sequence ties stay unresolved.** The sequence is `SEQUENCE_UNRESOLVED`. Blocks on which all tied sequences agree, and agree under every projection, are reported. All others are `GENOTYPE_UNRESOLVED`. No member is promoted.
5. **Support and GQ are projection-aware** (review finding C2).
   - **Read support:** must consider every admissible projection of every contributing sequence group.
   - **Block posterior:** for a block identified under 2.2, the posterior uses sequence-group probabilities, each counted once per sequence group, not once per decomposition.
   - **Representation-ambiguous blocks:** a block GQ is not defined. An interval-level sequence quality is reported instead.
   - **Lower-ranked groups:** a lower-ranked group whose projections include a different block genotype contributes error mass to that block.
6. **Phase, dosage and FASTA follow the same distinction.**
   - **FASTA:** written when the sequence is identified, with a header flag listing representation intervals.
   - **Phase:** relations are evaluated between identified units, so an interval counts as one unit.
   - **Dosage eligibility:** evaluated on compound intervals where they apply. Review finding C3, node-representation dependence, is handled separately.
7. **Order independence is required.**
   - **Invariant to:** candidate enumeration order, template (P-line) order, node renumbering and alias order.
   - **Canonical choices:** every set, ranking tie and listed alternative is ordered canonically, by sequence hash and then label tuple, never by discovery order.
   - **Enforcement:** the fixtures below include permuted variants, which must give identical canonicalised outputs.
8. **Versioning.**
   - **New schema version:** the output schema changes to a new version: new columns, and a locus-level `sequence_status` beside `call_class`.
   - **Caller status names:** the existing caller status names keep their meaning for comparison with historical outputs.
   - **Baseline:** a versioned development baseline is recorded before and after, on the fixtures and a small development set.

## 4. Adversarial fixtures (built before any caller change; `tests/analysis/representation_fixtures/`)

Each is a small synthetic locus, a graph plus a panel, with a donor whose reads are simulated.

| fixture | construction | expected under this contract | recorded current behaviour |
|---|---|---|---|
| F1 pure label ambiguity | a 60 bp unit R; bubbles A (∅ / R) and B (∅ / R) around a context block equal to R. Templates spelling the same sequence with labels (R, ∅) and (∅, R) are both in the panel | sequence identified; A and B `REPRESENTATION_AMBIGUOUS`; interval A..B reported with both tuples; flanks identified | `baseline_current.tsv` |
| F2 ambiguity plus an independent block | F1 plus a distant heterozygous bubble E with distinct 80 bp alleles | as F1, and E identified with its own GQ | (as above) |
| F3 phase tie | two heterozygous bubbles farther apart than any fragment; a symmetric panel, so both phasings tie | sequence `SEQUENCE_UNRESOLVED` (distinct sequences); both blocks identified, unphased | (as above) |
| F4 ambiguity plus phase tie | F1 combined with F3 | sequence unresolved; A and B `REPRESENTATION_AMBIGUOUS` under every tied sequence; phase blocks identified | (as above) |
| permuted variants of F1–F4 | P-line order reversed, node ids renumbered, sample names permuted | canonicalised outputs identical to the unpermuted fixture | (as above) |

**Real-data regression cases** (from the sealed pilot, read-only):
- the eight label-ambiguity cells must become sequence-identified with ambiguous sv:3/sv:5;
- mosaic three_2 at 30× must stay a distinct-sequence tie;
- the 27 calls' identified blocks must be unchanged.

## 5. Implementation order (after acceptance)

1. **Preparation:** representation intervals as part of the prepared model, covered by a preflight check, so that `representation_independent` no longer implies injectivity.
2. **Annotations:** projection sets per candidate and per sequence group, plus the C2 support/GQ correction, behind a schema version flag. The fixtures and real-data cases run old against new.
3. **Downstream:** the new outputs (sequence status, interval rows, FASTA flags), then phase and dosage handling of intervals.
4. **Validation:** a separate versioned validation of the changed outputs. This is never mixed with performance patches.

## 6. Version 2: sequence-first output and confidence contract (S0 v2)

**Scope:** an output contract only. No caller, default or scoring change is made in S0. Everything below is conditional on the explored candidate pool unless a completeness certificate says otherwise.

### 6.1 Principles

1. **Predict, then qualify.** Scientific uncertainty keeps the best available prediction and lowers its confidence. It does not by itself produce a no-call. A locus produces no prediction only for a hard failure: invalid input, model refusal, or no admissible candidate.
2. **Labels are not sequence.** Uncertainty about panel labels (which block carries a repeat unit, which alias a donor is) is reported separately and never lowers sequence confidence. Sequence-identical candidates are pooled before any confidence is computed: they are one sequence group, counted once.
3. **Small and material differences are different.** The difference between two diploid sequence predictions is the minimum, over the two homologue pairings, of the per-haplotype distances. It is classified as:
   - `IDENTICAL`;
   - `SMALL`: every difference shorter than a tolerance τ, with no change to copy structure or allocation. **τ is not yet set.** The 50 bp `min_variant_bp` is a graph-construction cutoff, not an approved correctness tolerance (correction C1, 6.8);
   - `MATERIAL`: any SV-scale difference, any copy-number or block-structure change, or any change in which chromosome copy carries a unit.
4. **One headline confidence.** Users read one confidence value and one filter. Diagnostics, alternatives and component scores go to the advanced outputs. Users never combine filters themselves.
5. **Completed and provisional results are never confused.** A search that stopped at a cap or a time limit reports its best prediction so far as `PROVISIONAL`. Its confidence is marked uncalibrated.

### 6.2 Outcome classes (locus `sequence_status`)

| class | meaning | reported prediction | block rows | CYP pilot cells (sealed, 36) |
|---|---|---|---|---|
| `RESOLVED` | one unordered sequence pair in the decision set (labellings pooled) | that pair | every block with a unique projection; one interval row per representation interval with all label tuples and no chosen tuple | 22 current calls, unchanged; **8 current label-ambiguity no-calls** (sv:3..sv:5 interval row) |
| `GENOTYPE_RESOLVED_SEQUENCE_SET` | distinct sequence pairs in the decision set share every block genotype: they differ in phase, or in label-equivalent placement only | the representative pair (6.4) | all blocks, with phase sets; phase across unlinked blocks marked unknown | 5 current genotype-equivalent calls |
| `PARTIALLY_RESOLVED` | the decision set's sequence pairs agree on some blocks and differ at others | the representative pair | agreeing blocks identified; differing blocks `GENOTYPE_UNRESOLVED`, with the alternatives in advanced output | none in the pilot (distinct from the tie below only in how many blocks agree) |
| `SEQUENCE_TIED` | distinct sequence pairs with different block genotypes tie within tolerance | the representative pair, at low confidence | blocks on which every tied pair agrees under every projection; all others `GENOTYPE_UNRESOLVED` | 1 (mosaic three_2 at 30×, 9 tied members; truth among them) |
| `PROVISIONAL` (orthogonal flag) | the search did not complete | best prediction so far, together with the class it would have had | as for that class | none at A5. At A3, all 30× cells would have been `PROVISIONAL`, not no-calls |

`PARTIALLY_RESOLVED` and `SEQUENCE_TIED` are kept apart because users act differently on "most of the locus is known" and "two different structures are equally likely". Both retain a prediction.

### 6.3 Headline confidence: target and calibration

- **Correctness target.** `structure_confidence` is the probability that the reported unordered diploid sequence pair is within `SMALL` of the truth:
  - same copy structure;
  - same chromosome allocation;
  - every difference shorter than τ (undeclared; 6.1, C1).

  A second, advanced-only quantity, `exact_sequence_confidence`, targets byte identity. The headline is the structural one because it is what users act on (star alleles, CNV, hybrids). Exactness is also not attainable under availability limits, such as a truth that cannot be spelled under leave-one-out.
- **Model estimate.** The decision-pool posterior is the sum over sequence groups g of P(g | reads, pool) · 1[d(g, reported) ≤ SMALL]. It is computed from:
  - pooled, label-collapsed sequence-group masses;
  - the existing joint likelihood and prior;
  - the normalisation over the explored pool.

  It is candidate-conditional, and **not a bound in either direction** (C2). Candidates outside the pool can carry mass inside the reported pair's tolerance as well as outside it, so renormalising over the explored pool can overstate or understate the full-space value. It is shown in advanced output as `pool_conditional_posterior`.
- **Calibration.** Map the raw estimate, plus a small set of declared covariates **observable at prediction time** (coverage, completion, score margins, pool and search diagnostics, model features), to `structure_confidence`:
  - LZO, mosaic and leave-one-out are **evaluation strata, not inputs** (C3). A production calibrator never knows which stratum a sample belongs to, so it may not condition on it;
  - fit a monotone map on development cells with truth, pooled over strata in declared proportions;
  - report reliability curves and the calibration error **per stratum**. Leave-one-out is the stratum that governs production claims, because it reflects panel availability;
  - held-out cells are used once, for a final reliability check only, never for fitting.

  Until a calibration exists, the headline field is `NA` with reason `UNCALIBRATED`, and the raw estimate is advanced-only.
- **Filter.** The single `FILTER` value is one of:
  - `PASS`: `structure_confidence` ≥ the declared threshold and the result is complete;
  - `LOW_CONFIDENCE`;
  - `PROVISIONAL`;
  - `NO_PREDICTION` (hard failure only).

  The current headline filters (GQ floors, dosage flags, representation flags) become calibration covariates or advanced diagnostics. They are not separate gates a user must combine.

### 6.4 The representative prediction

**Rule for the first reporting implementation:** keep the current score-based estimator, the caller's pooled decision with its existing tie rules (C4). The confidence is then defined against that estimator.

**Experimental alternative, not to be implemented first:** report the pool member that minimises expected structural loss, the minimum-Bayes-risk (MBR) representative. Over pool sequence groups g′ with posterior w(g′), it is argmin over g of Σ_{g′} w(g′)·loss(g, g′), where:
- loss = 0 for `IDENTICAL`;
- loss = a small value for `SMALL`;
- loss = 1 for `MATERIAL`.

**Why MBR might do better than the highest-scoring pair:**
- When one group dominates, MBR and the top-scoring pair coincide.
- In a tie between structurally clustered alternatives, MBR picks the member closest to the others. That raises the probability of being structurally right without claiming more.
- If it were adopted, the headline confidence would have to be redefined against it.

**Status.** Experimental. Whether MBR beats top-score on real cells is an empirical question. It can be checked post hoc on the sealed pilot and C4 equivalence sets (these outputs list every tied member). That check must be labelled post hoc. The representative is always a pool member, never a synthesised consensus.

### 6.5 Advanced outputs

- **Alternatives table:** the top sequence groups, each with:
  - pool-conditional posterior;
  - distance class and size to the representative;
  - block-genotype projections;
  - representation intervals.
- **Label-equivalence sets** per interval.
- **Per-block conditional GQ** (projection-aware, C2).
- **Phase sets.**
- **Raw components:** fragment, exposure and prior.
- **Completion and caps record.**
- **Status mapping:** the old `call_class` and status names, kept for comparison with history.

### 6.6 What changes for the sealed CYP2D6 pilot (reporting projection only; `tests/results/s0_v2_20261015/cyp_pilot_contract_projection.tsv`)

**Three separate measures (C5).** Every comparison reports them separately, never one figure for all three:
- **prediction yield:** cells with any prediction;
- **high-confidence yield:** cells whose prediction passes the declared confidence threshold;
- **accuracy:** of all predictions, and of high-confidence predictions, each with its own denominator.

The changes below affect prediction yield. High-confidence yield cannot be stated until a calibrated confidence exists.

**Label-ambiguity no-calls.** The 8 label-ambiguity no-calls become `RESOLVED` predictions, reported with an sv:3..sv:5 interval row. Post hoc, 7 of them spell the true sequence. One is wrong at 8 blocks: leave-one-out HG02953 at 30×, where the truth cannot be spelled from that model. This shows why a resolved label must not raise confidence by itself. **No callability gain is claimed** until these cells' calibrated `structure_confidence` exists. A label rule converts a no-call into a prediction. It does not make the prediction right.

**Distinct-sequence tie (1 cell).** It stays tied. It now carries the score-based representative at low confidence, plus the blocks its 9 members agree on.

**Genotype-equivalent calls (5).** Their block genotypes are unchanged. The sequence pair becomes an explicitly chosen representative with its own confidence, no longer an implicit one.

**The 22 strict calls.**
- Blocks and sequences are unchanged.
- Three of the 27 current calls carry 12 kb sv:9 errors (NA19983 ×2, excluded by today's filter; NA19650 at 5×, retained). These are `MATERIAL` errors under the target, so a calibrated confidence must account for them.
- No threshold may be tuned on these three cells. They are development cells, not a test.

### 6.7 Demonstrations and the remaining gate

**Distinctions already covered:**
- resolved sequence with ambiguous labels: F1/F2;
- distinct sequence pairs sharing block genotypes, or phase: F3;
- a combination of both: F4;
- order independence: the permuted variants;
- sum versus best history, non-Markov allele states, neutral subdivision and distant-factor rank: `tests/genotype_mosaic_core.cpp`, `test_conceded_counterexamples`, each with a negative control.

**Missing:** `PARTIALLY_RESOLVED` with a `MATERIAL` alternative, and a `PROVISIONAL` fixture. Add them when the output change is implemented, not before.

**Gate.** Implementation of section 6 is a versioned output-schema change. It follows the order in section 5, plus a calibration study that has its own pre-registered split.

### 6.8 Corrections (Sol's review of S0 v2, 2026-10-15)

Section 6 above has been amended in place. The original statements and their corrections:

- **C1.** "SMALL: shorter than 50 bp, the model's `min_variant_bp`." The 50 bp construction cutoff is not an approved correctness tolerance. τ is a decision still to be made, with its own justification.
- **C2.** "Omitted mass is not bounded, so the raw value is an upper bound." A candidate-conditional posterior is not generally an upper bound. Omitted candidates can add mass inside the tolerance as well as outside it.
- **C3.** Calibration "by scenario" with "leave-one-out exposure" as a covariate. LZO, mosaic and leave-one-out are evaluation strata. That information is not available to a production calibrator. They are used only to report calibration.
- **C4.** MBR as "the rule". Minimum-expected-loss selection stays experimental. The first reporting implementation keeps the current score-based estimator.
- **C5.** Yield and accuracy were not separated. Prediction yield, high-confidence yield and accuracy are now reported separately.
- **C6** (S0 v2 measurements, `tests/results/s0_v2_20261015/README.md`). The distant-evidence "minor-side mass share" is a descriptive statistic of the explored pool. It is not a bound, and it does not permit distant evidence to be discarded.
- **C7** (same README). The proposed CYP2D6 sv:9 depth diagnostic needs expected counts produced by the same read-selection process as its observed counts. The 701 fragments were selected because they discriminate between pool candidates. They cannot be compared with a whole-block exposure term.


## 7. S1-A: the corrected contract for the first reporting implementation (2026-10-17)

**What this section is.** It narrows section 6 to what S1 implements, and corrects assumptions that section 6 left open. Where it differs from section 6, it supersedes. Evidence is in `tests/analysis/s1a/` (fixtures F1–F7 and permutations, current outputs against proposed ones) and `tests/results/s1a_sequence_reporting_20261017/`.

### 7.1 Essential reporting fixes (S1) and later scientific changes (not S1)

| S1: reporting only, no score, prior or search change | later, each with its own study |
|---|---|
| keep a prediction after any completed, valid ranking | calibration (`structure_confidence`) and any PASS threshold |
| pool sequence-identical labellings; label-only ambiguity becomes `REPRESENTATION_AMBIGUOUS` blocks and interval rows, not a locus no-call | a structural tolerance τ, `SMALL`/`MATERIAL` classes, homologue-allocation equivalence |
| distinct-sequence ties keep a canonical representative plus alternatives; agreeing blocks stay identified | MBR or any other representative rule |
| `PROVISIONAL` predictions from a valid, exactly scored candidate when the search stops at a declared cap | projection weighting across alternative block labellings, and projection-aware GQ as a probability |
| conservative projection-aware error-mass bounds for identified blocks (explored pool); the existing GQ kept, labelled as representative-projection GQ | complete projection sets from model-level representation intervals (section 5 step 1) |
| canonical, order-independent outputs | dosage and phase handling of intervals |

### 7.2 Corrections to section 6

- **K1. The posterior is not an upper bound** (restates C2). The candidate-conditional posterior is not an upper (or lower) bound on biological correctness. It is reported only as the advanced `pool_conditional_posterior` of the representative's sequence group, never as confidence. **RESOLVED does not imply high confidence:** F4 is `RESOLVED` with `pool_conditional_posterior` 0.333, because distinct phase alternatives sit just below the decision set.
- **K2. Strata are not inputs** (restates C3). LZO, mosaic and LOO are benchmark strata. Nothing in S1 conditions on them, or names them, at prediction time.
- **K3. No calibrated probability and no new PASS.** S1 introduces no `FILTER=PASS` semantics:
  - `structure_confidence` is written as `NA` with reason `UNCALIBRATED`;
  - the existing headline filter columns keep their current definitions and are computed only on identified blocks;
  - section 6.3's PASS rule is deferred until a calibration study exists.
- **K4. No tolerance classes in S1** (restates C1). S1 compares sequences only as `IDENTICAL` or `DISTINCT` (by digest). There are no `SMALL`/`MATERIAL` classes, and no claim of structural or homologue-allocation equivalence: a 50 bp cutoff defines neither.
- **K5. The representative.** Among the decision set's sequence groups, it is the one with the highest certified score, ties broken by the canonical sequence-digest pair (the group order `decide_pooled_genotype` already uses). It never depends on candidate, template or alias order. MBR stays research.
  - The alternatives table is the decision set plus the top 5 groups by score.
  - Each alternative is compared with the representative only, block by block. No pairwise sequence distance is computed across the pool.
- **K6. Projections (review C2).** Every sequence group carries exactly one weight. For an identified block, S1 reports conservative error-mass bounds over the explored pool:
  - **lower:** groups all of whose projections differ at the block;
  - **upper:** groups with any projection that differs.

  No point probability is invented for a labelling choice, and a group's weight is never duplicated across labellings to form a probability. The existing `conditional_GQ` is kept, labelled `representative_projection_GQ`.
  - **F6 shows why.** The current GQs are 35.2 (sv:1) and 15.8 (sv:2) for a symmetric situation, while both blocks' error mass lies in [1.3e-4, 0.026].
  - **Coverage.** S1 projection sets come from the labellings the caller explored (its sequence groups' members), and are flagged `projection_coverage=EXPLORED`. Complete sets need the model's representation intervals (section 5 step 1). On F1–F7 the explored labellings happened to cover every projection; that is a fixture property, not a guarantee.
- **K7. `PARTIALLY_RESOLVED` is dropped as a separate S1 class.** Section 6 distinguished it from `SEQUENCE_TIED` only by "how many blocks agree", with no rule. S1 reports `SEQUENCE_TIED` whenever distinct decision-set sequence pairs differ at any block genotype. Partial resolution is carried by the block rows: agreeing blocks `IDENTIFIED`, the others `GENOTYPE_UNRESOLVED` with their alternatives (F5).
- **K8. When a result is `PROVISIONAL`.**
  - **Provisional:** only when the search stopped at a declared cap (round cap, cumulative cap) **and** a valid, exactly scored candidate exists — the last committed incumbent. `sequence_status` is then `NOT_DETERMINED`, because no decision set was computed, and no confidence is defined (F7).
  - **No prediction:** stops by the time or footprint watchdog. The process is stopped from another thread with no safe state to report, so S1 writes no prediction there. Incomplete score accumulation never yields a prediction.
- **K9. Hard failures** (invalid input, model refusal, no admissible candidate) remain `NO_PREDICTION`.

### 7.3 Status mapping (old names kept for comparison)

| current (`run_status`, `call_class`, `no_call_reason`) | S1 (`prediction`, `sequence_status`) | fixture |
|---|---|---|
| COMPLETE, strict_call | PREDICTION, RESOLVED | F6, F7full |
| CONDITIONAL_GENOTYPE_EQUIVALENT, genotype_equivalent_call | PREDICTION, GENOTYPE_RESOLVED_SEQUENCE_SET (the representative's FASTA is now written) | F3 |
| RESCORE_INCOMPLETE, no_call, distinct_tie — top set is one sequence group | PREDICTION, RESOLVED, with `REPRESENTATION_AMBIGUOUS` blocks and interval rows | F1, F2, F4 |
| RESCORE_INCOMPLETE, no_call, distinct_tie — distinct sequence groups | PREDICTION, SEQUENCE_TIED, agreeing blocks identified | F5 |
| REPAIR_INCOMPLETE, no_call, round_cap or cumulative_cap | PREDICTION, NOT_DETERMINED, provisional | F7 |
| TIME_LIMIT, RESOURCE_REFUSAL, CANDIDATE_INCOMPLETE, invalid input | NO_PREDICTION | — |

`run_status`, `call_class` and `no_call_reason` keep their current values under the new schema, so history stays comparable.
