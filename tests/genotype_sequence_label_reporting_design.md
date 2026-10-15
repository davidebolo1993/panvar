# Sequence, block-genotype and representation-ambiguity reporting: prospective design (no caller change yet)

**Status:** a design note, written 2026-10-15 after the CYP2D6/D7 pilot (sealed diagnostic baseline). It is not implemented.

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
