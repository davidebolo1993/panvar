# M10-B.1 result: boundary transformation reconciled; most M10-B per-interval claims are not coordinate-supported; no saved comparison is a uniquely resolved agreement

**Run.** Branch `m10-b-structural-projection`; freeze commit `82ef5c1`.
- **Results:** `tests/results/m10b1_projection_hardening_20261020/`, with `runs/audit/{model_audit, block_audit, case_scopes,
  degeneracy_v2}.tsv` and `runs/cases/{projections_v2, comparisons_v2}.tsv`.
- **Budget:** 39.2 of 600 s, including 15 s of disclosed development. 11/11 units passed first time.
- **Unchanged:** M10-B's outputs, its `project.py` and its fixtures (re-gated, 17/17).
- **Not done:** no reads, genotyping, thresholds, caller change or defaults.

## 1. Allele coordinates against source walks

**From the implementation** (`genotype_locus.cpp`; DESIGN_B1 §1):
- **Recorded walks come before normalisation.** A variant walk includes both bubble boundary nodes. Abutting bubbles drop
  their shared node from one side.
- **The divergence trim** moves each variant block's prefix P and suffix S common to all source alleles into the
  neighbouring context blocks. So:
  - variant: `normalized = spell(walk)[P : len − S]`;
  - context: `normalized = S(prev) + spell(walk) + P(next)`.

**Audit, all six models** (extended dump: the report's trims, source kinds, invariant segments):

| model | blocks verified | blocks trimmed | walks verified / total | walks that respell their normalized allele unaided | templates respelled whole | metadata-missing alleles | empty / bypass sources | non-empty invariant segments |
|---|---|---|---|---|---|---|---|---|
| C4 dev / LOO | 11/11, 11/11 | 5, 5 | 736/736, 726/726 | **0**, **0** | 131/131, 129/129 | 0 | 0 | 0 |
| CYP2D6 dev / LOO | 19/19, 19/19 | 9, 9 | 869/869, 860/860 | **0**, **0** | 126/126, 124/124 | 0 | 0 | 0 |
| ACOT dev / LOO | 19/19, 19/19 | 9, 9 | 2,912/2,912, 2,605/2,605 | **0**, **0** | 467/467, 407/407 | 0 | 0 | 0 |

**What follows:**
- **Sol's observation is confirmed and explained.** No individual walk respells its normalized block allele, yet every
  walk satisfies the trim relation, and every whole template respells exactly.

**What can be projected safely (`project_v2.scope`):**

| scope | supported when |
|---|---|
| whole haplotype | always, on these models (relation verified, no clip) |
| compound interval or single normalized block | the relation is verified, **no bytes enter from outside** (an end context block next to a trimmed variant block takes bytes whose nodes are in a neighbour's walk, which the scope does not determine), and **no skeleton node straddles** the clip of an end variant block |
| anything else | UNSUPPORTED; correspondence is not invented |

**The frozen case scopes:**
- C4 sv:7: supported (clips of 197 and 124 bp fall between skeleton nodes);
- C4 whole locus: supported;
- CYP2D6 sv:3..sv:5: supported (clips of 278 and 244 bp);
- CYP2D6 sv:9: **the scope is valid, but an 81 bp node straddles the clip** on three of the six haplotypes;
- ACOT sv:5..ctx:6: **unsupported**, because 35 bytes enter ctx:6 from sv:7.

## 2. Missing metadata (`project_v2.py`)

- A non-empty allele with no walk, or with fewer recorded walks than walk sources, makes its projection INCOMPLETE.
  Comparisons with it are UNSUPPORTED.
- A valid empty record ('empty' or 'bypass', with no walk) contributes an empty walk.
- Fixtures M1–M3 pass. M4 reproduces version 1's defect: absent walks were RESOLVED, and two unrelated sequences compared
  EQUAL.
- **The saved models have no metadata-missing allele and no empty or bypass source.** The defect therefore changed no saved
  verdict.

## 3. The long-node convention under adversarial fixtures (exposed, not tuned)

| fixture | intended distinction | result |
|---|---|---|
| A1: the same unit as one 60 bp node vs ten 6 bp nodes (plus a neutral SNP) | same structure | **not preserved**: DISAGREEMENT |
| A2: a 55 bp insertion made only of 5 bp nodes | structural difference | **not preserved**: UNIQUE_AGREEMENT |
| A3: a neutral SNP inside a long node vs as its own short node | same structure | **not preserved**: DISAGREEMENT |
| A4: shifted block boundaries, same compound sequence | identical | preserved |
| A5: declared multi-base class site vs a neutral short change | only the declared one differs | preserved, **only through the explicit declaration** |
| A5b: a 1 bp site declared by sequence class (G) | G site only | **not preserved**: the neutral C change is captured too. Short sites need positional declaration |

**Class identity** is now the canonical sequence bytes. Gate C1 replaces the hash with a constant: the classes stay
distinct (ids disambiguated) and every verdict is unchanged. **The fixtures pass 34/34.**

**Conclusion:** the skeleton depends on graph node segmentation (A1–A3). That confirms it cannot be the production
structural alphabet.

## 4. Saved-case replay (M10-B's frozen set; version 2 beside version 1)

| category | comparisons |
|---|---|
| **UNIQUE_AGREEMENT** | **0** |
| SET_IDENTICAL | 2: C4 control (call = donor templates, bytes equal); CYP2D6 sv:3..sv:5 synonymous (bytes equal) |
| COMPATIBLE | 0 |
| DISAGREEMENT | 2: C4 NA19650 reported 3+1 vs saved 2+2 (`DOSAGE_EQUAL_ALLOCATION_DIFFERS`); 3+1 vs saved 1+2 |
| UNSUPPORTED | 11 (below) |

**The 11 UNSUPPORTED comparisons:**
- C4 NA19650 truth ×2: the truth is off-model (UNKNOWN). Unchanged from version 1.
- CYP2D6 sv:9 ×3: an 81 bp node straddles the clip. **Changed**: version 1 said DISJOINT, DISJOINT and EQUAL.
- ACOT L4 and L2 ×6: 35 bytes enter the scope from sv:7. **Changed**: version 1 said EQUAL for all six.

**Changed verdicts: 9 of 15**, all through coordinate support, none through missing metadata. Every comparison, unsupported
ones included, stays in the tables.

The degeneracy inventory under version 2:
- CYP2D6: 6 resolved, 36 set-valued;
- ACOT: all 2,458 unsupported (same 35-byte boundary at sv:5..ctx:6).

## 5. Conclusion

**Claims that remain supported:**
- **Labelling invariance.** The projection is a function of the haplotype sequence. Synonymous labellings cannot create a
  structural difference: CYP2D6 sv:3..sv:5 is SET_IDENTICAL with equal bytes.
- **Allocation at C4 sv:7.** 3+1 and 2+2 project differently at equal dosage, among model-spellable hypotheses, on a
  coordinate-verified single-block scope.
- **Whole-haplotype projection** is coordinate-safe on all six models.
- **Missing metadata** now gives INCOMPLETE, never equality.

**Claims withdrawn:**
- the CYP2D6 sv:9 dosage distinction on its block scope;
- every ACOT result at sv:5..ctx:6: synonymous agreement, and few-base residuals being "not structural";
- any reading of SET_IDENTICAL or UNKNOWN as a uniquely correct structural genotype;
- the long-node skeleton as a structural alphabet (A1–A3).

**Boundary metadata needed:**
1. **Per normalized allele, its walk clipped to normalized coordinates:** the walk, plus the offset into its first and last
   node. The trims are now dumped from the report.
2. **For intervals:** either boundaries at representation-independent coordinates (blocks without moved bytes), or the node
   segmentation of the bytes moved across each interval end (the common flank, which the trim makes identical in sequence
   but not necessarily in nodes).
3. **Source-kind records** (walk, empty, bypass) carried with every allele, as dumped here, so that "missing" can be told
   apart from "empty".

**Minimum semantic decisions before an off-model parser** (all for Davide):
1. **A representation-independent structural alphabet:** units and segments defined on sequence, not on graph nodes.
   A1–A3 show why node segmentation cannot define it.
2. **Small-difference semantics:** which differences are ignorable, and the indel length that is structural however it is
   represented (A2).
3. **Positional class-defining sites** in unit coordinates, not sequence classes (A5b).
4. **Scope semantics:** whole haplotype, or intervals with representation-independent ends, and how bytes moved across ends
   are treated.
5. **Ambiguity semantics:** whether identical set-valued projections count as agreement in evaluation, and whether a
   canonical placement rule is required (CYP2D6 sv:3..sv:5).

**Not built:** the off-model parser. No merge, push, storage cleanup, prior change, calibration, clustering or long-read
work.
