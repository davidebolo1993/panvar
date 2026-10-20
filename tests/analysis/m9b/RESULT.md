# M9-B result: the caller's output unit is the reportable block (a bubble leaf), and no SV record is emitted; proposed contract: one conditional, uncalibrated block-allele agreement value per reportable block, with phase separate

**Run.** Branch `m9-b-confidence-contract`.
- **Commits:** freeze `3532125`; Amendment 1 `8c9fa90` (fixture labels); Amendment 2 `2b8554b` (the syntax defect in
  Amendment 1).
- **Results:** `tests/results/m9b_confidence_contract_20261020/`.

**Execution: 5.0 of 300 s,** including 2 s disclosed before the freeze.

**Recorded failures:**
- `fixtures`: J1/J2 competitor labels did not disagree where intended (a fixture defect).
- `fixtures.a1`: my inline comments swallowed closing brackets.
- An empty `runs/evidence/` directory is left from the gate refusal.

**Scope:** read-only. No caller, genotyping, fitting, clustering, threshold, GQ/PASS or default change.

**Appended:** a qualification to M9-A `RESULT.md`. Historical packages are unchanged.

## 1. M9-A qualified

The qualification is appended to `tests/analysis/m9a/RESULT.md`. It states:
- **A and C are different propositions.** In HG01975 LOO:

| run | A | C |
|---|---|---|
| r0.01 | 0.2025 | 0.8378 |
| r0.1 | 0.1551 | 0.8085 |

  This replaces the "understates by 4–5×" wording.
- **C** is joint agreement across the reportable blocks; **D** is agreement at one block.
- **Neither** is a calibrated probability of correctness.

**Fixtures (`fixtures.a2`, 3/3 PASS),** through M9-A's frozen core aggregation, unchanged:

| fixture | D1 | D2 | C (joint) | product of marginals |
|---|---|---|---|---|
| J1 | 0.9 | 0.9 | **0.8** | 0.81 |
| J2 | 0.9 | 0.9 | **0.9** | 0.81 |

- The Fréchet bounds for these marginals are [0.8, 0.9].
- **J3:** C covers only reportable blocks (C = 0.9 while B = 0.6), and D at a reportable block is unaffected by the other
  blocks.

## 2. Output mapping (from code; citations in the agent trace, spot-checked)

**Prediction units**, short-read `panvar genotype`:
- **Locus:** `.locus_genotype.tsv`, one row (status, call class, headline filters).
- **Block:** `.block_genotypes.tsv`, one row per block, context included. `.genotypes.tsv` is one row per **variant
  block only**: the writer skips `kind != Variant` (`genotype_command.cpp:7393`).
- **Whole-locus sequence:** `.mosaic_haplotypes.fa`, two haplotypes, written for a unique call only.
- **Sequence hypotheses:** `.mosaic_scores.tsv` and `.equivalence_set.tsv` (groups); `.mosaic_candidates.tsv` (label
  pairs).
- **Per homologue, dosage stage only:** `.dosage_repair_candidates.tsv` and `.dosage_allocation_equivalence.tsv`
  (`units_a`/`units_b`).

**A reportable block is one bubble leaf, `sv:<index>`.**
- `reportable` is exactly `kind == Variant` (`genotype_command.cpp:7410-7412`); `Variant` is exactly
  `ChainBlockKind::Bubble` (`genotype_locus.cpp:649-653`).
- Bubbles survive `min_variant_bp` (50 here) if **any panel path** carries ≥ 50 bp inside, or there is an inversion
  signal (`bubbles.cpp:964, 1048-1053`).
- **"SV" is therefore a property of the bubble, not of the called alleles.** ACOT sv:5's called and truth alleles are
  20–30 bp.

**No SV record is emitted.**
- The genotype caller writes no VCF, SVTYPE, POS or END (no occurrence in `src/genotype_*.cpp`). Internal coordinates are
  not written to the block tables.
- `panvar call` writes VCF records from **panel-graph walks** against the reference walk, and reads no genotype output.
  `benchmark` states that no genotype is read.
- **So no block-to-record mapping exists.** Matching identifiers across representations (`sv:N`, bubble ids, call
  records) would not establish equivalence.

**What an allele label is:**
- A label is a **distinct sequence within its block**. Identical sequences are merged, so `within_block_duplicates = 0`
  (`genotype_model.cpp:256-286`).
- One label can hold **several graph walks** spelling the same bytes.
- **Example, C4 sv:7** (graph `genotype_c4_production_11fe681/graphs/NA19650.loo.gfa`, reference
  `HG00096#1#haplotype1-0000024:31848049-32074186`, preset short-read, `min_variant_bp` 50, fingerprint
  `14851130f453a071`):
  - label 47 (32,714 bp) is spelled by **2 distinct walks**, which differ in node list but have the same shape (992 steps,
    no node revisits);
  - label 41 is one walk with nodes visited up to 3 times;
  - so a label does not carry the walk, and a structural difference between walks sharing a label is possible but not
    demonstrated here.

**Different block-allele pairs can spell the same sequence.**
- **Example, the CYP2D6 LOO NA19650 model** `e41f77fda080c143`, the same model as the M1-B CYP runs (sealed
  `cyp2d6_shortread_pilot_v1/degeneracy/`, model-only):
  - **42 cross-block equivalences** over sv:3 + ctx:4 + sv:5: different allele tuples spell identical 251-bp sequences;
  - there are no within-block duplicates.
- **Consequence:** a per-block label target can count sequence-identical predictions as disagreement there. The pilot
  reported such cases as the no-call class SEQUENCE_IDENTICAL_BLOCK_LABEL_AMBIGUITY.

**Projection-ambiguous sequence groups:**

| situation | handling |
|---|---|
| a group whose label pairs project to different block genotypes | flagged, representative = highest path prior (`genotype_repair.cpp:1419-1438`) |
| **winner ambiguous, pooled route** | `DistinctTie` and a **no-call**, `distinct_tie`/`divergent_tie`, with no block tables (`genotype_repair.cpp:1452-1467`, `genotype_command.cpp:6592-6599`) |
| non-winner ambiguous groups | enter the block posterior/GQ and `genotype_support` **through their representative only**; the sequence report labels this `legacy_representative_gq … not_projection_aware` |
| non-pooled route | the `projection_ambiguous` column repeats the locus-wide winner flag on every row |
| headline filters | ignore the flag |

In all eight audited ACOT runs and the four M1-B runs, no group was ambiguous.

## 3. Three targets (described, not implemented)

**(a) Exact block-allele genotype agreement.**
- **Defined at:** each reportable block. The reported unordered pair of block-allele labels (= distinct block sequences)
  agrees with the hypothesis.
- **Mass:** existing D (`conditional_genotype_posterior`/`conditional_GQ`), and C jointly.
- **Contract value must be:**
  - projection-aware (M9-A interval semantics), or certified free of MIXED groups;
  - declared at the **minimal compound interval** wherever cross-block degeneracies exist (for example CYP2D6
    sv:3..sv:5), not at a single block.

**(b) Emitted SV-record genotype agreement. Not defined: there is no sound mapping today.** A future mapping would need:
- a deterministic function from block-allele labels to record alleles, bound to one graph, model and export setting;
- explicit many-to-one cases: distinct sequences mapping to one record allele, as in small internal differences;
- explicit one-to-many cases: a label spanning several records;
- aggregation of mass over every label combination that maps to the same record genotype, by summing, with MIXED
  semantics.

**(c) Future structural-representative agreement. Information that must be preserved per haplotype:**
- unit identity class and copy number (per homologue, not pooled);
- allocation of units across homologues;
- unit order;
- orientation of each unit;
- the relevant junctions: flank-to-unit, unit-to-unit, and unit-to-flank at entry and exit.

Where some of this exists today:
- `.dosage_repair_candidates.tsv` / `.dosage_allocation_equivalence.tsv` (`units_a`, `units_b`, `states`) and
  `.dosage_allocation.tsv` (`copies`, `arrangement`, `flank_order`, `shape`) are dosage-stage diagnostics only;
- walk signatures carry order and orientation implicitly (`allele*_sources`);
- no main table reports per-haplotype copy number, order, orientation or junctions.

Small sequence differences may eventually be tolerated, but no distance or clustering is proposed here.

## 4. Illustrative failures (saved evidence only; truth from archived evaluators)

| case (graph / model) | call | conditional support | violates | unassessed |
|---|---|---|---|---|
| **ACOT HG01975 LOO sv:5** (`acot.gfa`, LOO curation, M8-C profiles r0.01 / r0.1, fingerprint `12cbe6be…`) | 89/62 (27/30 bp) against truth-closest 89/113 (27/20 bp); floor 0; excess 10 bp | D = 1.0, GQ 495 / 598, READS | **(a)** exact block-allele agreement | **(c):** whether a 20- vs 30-bp allele difference is structural. (b) is undefined. |
| **C4 NA19650 LOO sv:7** (10× and 30×; `NA19650.loo.gfa`; `14851130…`) | 41/47 = 3+1 units against truth 2+2 (65,452 / 65,449 bp); excess 65,444 | D = 1.0, GQ 149 / 741, READS, PASS | **(a)**, and **(c)** per-haplotype copy number and allocation | — |
| **CYP2D6 NA19650 LOO sv:9, 5×** (`loo.NA19650.gfa`; `e41f77fd…`) | 23/73 = 2+2 against truth 3+2; excess 12,159 | D = 0.986, GQ 18.5, PRIOR_ONLY, quality LOW_QUALITY | **(a)**, and **(c)** total and per-haplotype copy number | (b) undefined |
| **CYP2D6 sv:9, 30×** (same model; contrast) | 22/73 = 3+2, UNIT_PAIR_CORRECT; excess 12 bp | D = 1.0, READS | **(a)** only (small exact-sequence error) | **(c):** not violated at the unit-count level that was assessed; order, orientation and junctions were not assessed |

**ACOT sv:5:**
- The truth-closest pair was **reached**: 37 / 39 evaluated hypotheses carry it.
- As whole hypotheses, the best of them loses by 121.9 / 151.8 nats, of which 119.1 / 132.7 is likelihood.
- These hypotheses differ at other blocks too, so the contrast is not block-isolated, and attributing it to sv:5 is
  unassessed.

**C4 sv:7:**
- The dosage audit reads CONSISTENT (`units=4u/4I`): total content is right, allocation is wrong, and allocation was
  `NOT_TESTED`.
- Earlier work (sv7 self-review; M1-B) showed that allocation is not identifiable from short reads in this case.
- The existing conditional value (D = 1.0, PASS) gives no warning.

**CYP2D6 5×:**
- The existing labels did flag it: PRIOR_ONLY and LOW_QUALITY.
- The dosage audit FAILED, so the dosage stage did not assess total copy number.

**What not to conclude:**
- A small exact-sequence error (sv:5; CYP2D6 30×) is not a demonstrated structural error.
- A high posterior (D = 1.0 in three of these cases) is not demonstrated correctness: all pools are CONDITIONAL.

## 5. Proposed user-facing contract

| element | content |
|---|---|
| **Prediction** | Per reportable block (bubble leaf `sv:N`; at a declared compound interval where label degeneracy exists): the unordered pair of block alleles (label, length, sequence). The two whole-locus haplotype sequences (FASTA) when the call is unique. |
| **Primary genotype quality (one value per output unit)** | Target (a): conditional posterior mass of the reported unordered block-allele pair **within the evaluated candidate pool**, shown as GQ, labelled *conditional, uncalibrated*, with its evidence basis (READS / PRIOR_ONLY) beside it. Not a probability of biological or structural correctness. |
| **Phase (separate)** | `phase_support`, plus phase-given-genotype mass in detailed output. Never folded into genotype quality. |
| **Debug / detailed** | Locus joint agreement C, the exact-sequence mass A, `separation` / `numeric_separation_gap` (named by target), competitors with fragment / exposure / prior gaps, `fragment_contenders`, projection / MIXED counts, pool size and search flags, dosage-stage unit tables. |
| **No new filters** | Existing headline filters are unchanged. The contract defines meaning only. |

Before this is called a contract, the value must be projection-aware or must certify that no MIXED group exists. Today's
GQ assigns non-winner ambiguous groups to their representative.

## 6. Validation still needed, and the smallest next experiment

**Needed before any calibration claim:**
1. An independent check of whether target (a) mass tracks exact block-allele correctness, on samples never used for
   genotyping development.
2. Separately, structural truth (per-haplotype units, allocation, order, orientation, junctions) for target (c).
3. A sound record mapping before (b) can be stated at all.

**Sample status:**

| locus | samples | status |
|---|---|---|
| C4 | all panel families (M8-B inventory) | development only: **0 uninspected families**; nothing held out |
| CYP2D6 | all panel families | development only: **0 uninspected families**; nothing held out |
| ACOT | 13 families | development |
| ACOT | 205 families (training bank) | used for M8-B parameter profiles |
| ACOT | 30 complete-diploid outer-evaluation families | **2 now inspected** (HG03209, HG01975 in M8-C/M8-D) |
| ACOT | the other **28** outer-evaluation families | **never genotyped and never confidence-inspected**; their haplotypes were used once, in M8-B's score-level outer evaluation of r |

**Smallest suitable next experiment** (declared, not run): a descriptive reliability check of target (a) on the 28
held-out ACOT families.
- **Calls:** one LOO call per donor at 5×, legacy defaults (r = 0.01, no promotion), the frozen M8-C read and profile
  recipe.
- **Unit of analysis:** each reportable block whose two truth alleles are available in the calling panel (floor 0), so
  that exact agreement is attainable.
- **Report:** reported D against observed exact agreement, in pre-declared bins with counts and intervals. Separately
  report blocks where agreement is unattainable, and PRIOR_ONLY vs READS.
- **No fitting, thresholds or GQ remapping.**
- **Expected yield:** about 28 × 9 = 252 reportable blocks, mostly at D ≈ 1. The check mainly tests the high-confidence
  region (as sv:5 shows, D = 1 can be wrong). Low-D bins will be sparse.
- **Cost estimate:** about 28 × (55–165 s + 40 s evaluation) ≈ 45–95 min.

C4 and CYP2D6 need genuinely new samples or assemblies outside the inspected panels before either can be calibrated.

**Stopping here** at the contract and mapping checkpoint.
