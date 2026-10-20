# M10-B result: existing model information supports a narrower, named target (structure of model-spellable haplotypes); an off-model sequence parse is the minimum missing component

**Run.** Branch `m10-b-structural-projection`; freeze commit `b674c6f`.
- **Results:** `tests/results/m10b_structural_projection_20261020/`, with `runs/coverage/{model_supply,labelling_invariance}.tsv`
  and `runs/cases/{projections,comparisons}.tsv`.
- **Budget:** 64.7 of 1,200 s, including 45 s of disclosed development. All 10 units ran first time.
- **Not touched:** no reads, no genotyping, no thresholds fitted, no caller change. M10-A's configuration and indexes stay as
  the development checkpoint.

## 1. What the existing model supplies

**Source of each kind of information:**
- **Walks.** Every block allele carries its oriented node walk(s). Over a whole template the walks spell exactly the
  concatenated allele sequences, and blocks share no node.
- **Repeat copies.** A tandem copy re-visits the same nodes, so copies are located and oriented through walk multiplicity.
  The model has no separate unit decomposition: panphorte folding is a different command.
- **Pooled dosage features.** These are unlocated k-mer copy counts per allele. **They are not used as occurrences.**

**Measured per model:**

| measured | C4 (dev / LOO) | CYP2D6 (dev / LOO) | ACOT (dev / LOO) |
|---|---|---|---|
| blocks with a repeated long node; maximum multiplicity | sv:7; 3 | sv:9; 5 | sv:9, sv:13; 4 |
| blocks with a reverse-strand long node | none | sv:1, 3, 5, 9, 13, 15 | none |
| node classes with ≥ 2 nodes of identical sequence | 48 | 68 | 116 / 115 |
| alleles with differing walks for one sequence | 36 | 35 | 2 |
| degeneracy equivalences: raw walks identical | — (0 equivalences) | 0 / 42 | 0 / 2,458 |
| long-node skeletons identical | — | 8 / 42 | 2,458 / 2,458 |
| **parse-set projection** identical (by construction) | — | 42 / 42 (6 resolved, **36 set-valued**) | 2,458 / 2,458 (all resolved) |

**Coverage of the target's structural components:**
- *available* means usable for haplotypes that parse exactly into model alleles;
- *off-model truth* means a truth sequence that does not parse exactly.

| component | C4 | CYP2D6 | ACOT |
|---|---|---|---|
| copy count | available (sv:7 walk multiplicity) | available (sv:9) | available (sv:9, sv:13) |
| per-homologue allocation | available | available | available |
| order | available at long-node resolution | available | available |
| orientation | available, but no inversion is represented in the model | available (reverse long nodes present) | available, none represented |
| junctions | long-node adjacencies only; breakpoints inside runs of short nodes are invisible | same | same |
| unit class | long classes (for example a HERV-length insertion) available; **short class-defining differences (for example C4A/C4B) absent unless declared** | hybrid and paralog classes only where long nodes differ; short sites absent unless declared | short sites absent unless declared |
| ambiguity | 36 alleles have two walks, so some projections are set-valued | **sv:3..sv:5: the same sequence has several structural placements (36 / 42 equivalences set-valued)** | none in the inventory |
| off-model truth | **absent:** NA19650 sv:7 truth does not parse, so UNKNOWN | the saved truths at sv:3..sv:5 parse; none was saved at sv:9 | the saved truths parse |

## 2. The prototype (`project.py`; 17/17 frozen fixtures)

π is defined on a **sequence**. It is the **set** of long-node skeletons over every exact parse into model alleles, and over
every recorded walk of each allele.
- Each step is a **sequence class and its orientation**, never a node id.
- Long means ≥ 50 bp: the configured `min_variant_bp`, used as an **experimental convention**.
- Class-defining short sites are kept only if **explicitly declared**. None is declared for the real cases.

**Properties:**
- Identical sequences have identical projections by construction.
- More than one skeleton is AMBIGUOUS.
- No exact parse is UNKNOWN, and is never equated with anything.
- The diploid comparison is per homologue, so allocation counts.
- Dosage equality is only a flag on a DISJOINT verdict. There is no tolerance and no clustering.

**Fixtures covered:** synonymous split and duplicate node; neutral SNP; a toy class site, undeclared and declared; copy
number; 3+1 vs 2+2; homologue exchange; order; orientation; junction; an allowed mosaic absent from the panel; ambiguous
parse; off-model sequence; dosage-equal but order-different.

## 3. Saved development cases (inspected, not accuracy validation)

| case | comparison | verdict | exact bytes equal | distinction |
|---|---|---|---|---|
| C4 NA19650 allocation failure (sv:7) | reported 3+1 vs saved 2+2 | DISJOINT, `DOSAGE_EQUAL_ALLOCATION_DIFFERS` | no | **preserved**: profiles 3 copies + 1 copy against 2 + 2 (206 long classes) |
| | reported 3+1 vs truth; saved 2+2 vs truth | **UNKNOWN** | no | **unsupported**: the truth does not parse into the LOO model (off-model) |
| | reported 3+1 vs saved 1+2 | DISJOINT | no | preserved (dosage differs) |
| C4 HG00096 control (whole locus) | call vs the donor's panel templates | EQUAL_SET_VALUED | yes | preserved; h2 is set-valued only because an allele has two walks |
| CYP2D6 NA19650 dosage failure (sv:9) | 5× reported 2+2 vs 30× twin 3+2 | DISJOINT | no | **preserved** (3-copy profile on one homologue) |
| | 5× reported vs saved 2+3 | DISJOINT | no | preserved |
| | 30× twin vs saved 2+3 | EQUAL | **no** | small sequence differences **ignored** (below the convention) |
| | vs truth | — | — | **unsupported**: no truth block sequence was saved at sv:9 |
| CYP2D6 NA19650 synonymous (sv:3..sv:5) | 30× call vs truth | EQUAL_SET_VALUED | yes | **labelling difference removed** (block labels wrong by 118 + 118 bp); the structure is **ambiguous**: 3 and 2 placements on the two homologues |
| ACOT L4 (sv:5..ctx:6) | JP on vs truth | EQUAL | yes | synonymous labels removed |
| | JP off (interval distance 14 bp) vs truth | EQUAL | no | the 14 bp **ignored for the declared reason** (short nodes, no declared class site) |
| ACOT L2 | JP on (3 bp) / JP off (5 bp) vs truth | EQUAL / EQUAL | no | the few bases ignored for the declared reason |

**Also observed.** Both ACOT truths parse exactly into the LOO model. Each true haplotype was therefore model-spellable, and
the residual few-base error was a search or ranking outcome, not unavailability.

## 4. Decision: **the existing information supports a narrower, explicitly named target**

**Target that is supported now: "structural genotype of model-spellable haplotypes".** For any haplotype that parses exactly
into model alleles, π gives a labelling-invariant structural projection. It covers:
- copy number, homologue allocation, order and orientation;
- junctions at long-node resolution;
- class-defining short differences only when declared;
- explicit set-valued ambiguity where the model places one sequence in several ways.

**The cases show this already does what the full target needs on the prediction side:**
- C4 3+1 ≠ 2+2 at equal dosage;
- CYP2D6 2+2 ≠ 3+2;
- synonymous labellings (CYP2D6 118 + 118 bp, ACOT 3 + 3 bp) are no longer errors;
- few-base ACOT residuals are not structural.

**Not supported:**
1. **Off-model sequences, the minimum missing component.**
   - **What fails:** C4 NA19650's truth, and any truth or haplotype with novel content, has no exact parse. So π is UNKNOWN
     and the C4 allocation failure **cannot be scored against truth**.
   - **What is needed:** a deterministic parse of an arbitrary sequence onto the model's long-node sequence classes
     (located, oriented occurrences), in which declared small differences do not abort the parse.
   - **What must be approved first:** its detection rule and margin are exactly target-v2 decision item 2 (the
     unit-detection rule), so they need Davide's approval before implementation.
   - **Why the prototype stops here:** it deliberately builds no aligner and does not borrow M1-B's evaluator thresholds.
2. **Canonical placement** where the model spells one sequence in several structural ways (CYP2D6 sv:3..sv:5).
   - **Why it matters:** set-valued projections stay safe, but they weaken decisions. A sequence one neutral SNP away from
     truth there would be COMPATIBLE, not EQUAL.
   - **What is needed:** a canonical placement rule, for example left-normalised tandem expansion. This is secondary.
3. **Short class-defining differences** (C4A/C4B and similar) need a declared rule. That is target-v2 decision item 1.
4. **Untested risk of the 50 bp convention.**
   - **Too lenient:** a structural difference made only of short nodes is invisible.
   - **Too strict:** two parallel long branches differing by a few SNPs would count as different classes.

   Neither occurred in these cases, but the cases do not rule them out.

**Not done:**
- no confidence threshold;
- no calibration evidence;
- no accuracy claim (development cases only);
- no caller integration, new families, clustering, long-read work, merge, push or default change.


## Scope clarification (appended by M10-D, 2026-10-20; the record above is unchanged)

**Davide clarified the intended notion of structural similarity.** It is **sequence similarity**, not equality of an
explicit structural alphabet.
- Small substitutions and indels contribute proportionately to the distance.
- The evaluator is **not** required to identify or classify inversions, rearrangements, junctions or copy-number events.

**What this closes:**
- the structural-alphabet and event-comparator investigation (M10-B, B.1, C);
- M10-C method A was unsuitable for event decomposition, but that does **not** disqualify its exact global edit distance,
  which M10-D adopts as the evaluation metric.

M10-A remains the computational checkpoint.
