# Production quality target, version 2 (M10-A): specification only

**This replaces M9-J.1's `QUALITY_TARGET.md`**, which stays in place as the historical record. Nothing here is implemented,
calibrated or thresholded, and no tolerance value is chosen.

**Priority (Davide).** Correct SV structure first; exact panel identity and every SNP or indel second.

## 1. Why "within tolerance" cannot define the target

Version 1 called two haplotype structures equal when aligned segments were "within the declared tolerance". A pairwise
tolerance is **not transitive**.

- **Counterexample.** Take three copies of one unit. B differs from A at one site; C differs from B at a different site. So
  d(A,B) = 1, d(B,C) = 1 and d(A,C) = 2. With tolerance 1, A≈B and B≈C, but A≉C. "≈" does not partition sequences, so "the
  predicted class" and "the truth class" are undefined.
- **Clustering does not fix it.**
  - Single linkage chains A…Z across any distance.
  - Complete or average linkage depends on which other sequences happen to be present, and on their order.
  - Either way, correctness would depend on the rest of the cohort, not on prediction and truth.
- **A distance threshold is also the wrong kind of rule.** C4A and C4B differ by a handful of bases, yet the difference
  defines the class. A long run of neutral SNPs may define nothing. Whether a difference is class-defining depends on
  **where** it is, not on how many differences there are.

**Correction.** Equivalence is the kernel of a **function** of each haplotype alone: two haplotypes are equivalent iff
π(h1) = π(h2). That relation is reflexive, symmetric and transitive by construction. No comparison between two haplotypes
ever enters the definition. Small differences are ignored by **construction of π** (a declared position set), never by a
distance threshold.

## 2. The structural projection π (deterministic, model-derived, locus-agnostic)

**Input:** the sequence of a haplotype over a structural interval (section 4), and the genotype model.
- **Never the labelling.** Synonymous labellings spell the same sequence and therefore get the same π, as the compound
  target already guarantees for bytes.
- **The same function** applies to explored candidates (at run time) and to truth haplotypes (offline only).

**Output:** an ordered token string, or a **set** of token strings (below).

| token | meaning | derived from (generic; no locus annotation) |
|---|---|---|
| `U(unit, strand, v)` | one occurrence of a model unit, its orientation, and its alleles `v` at that unit's class-defining sites | the model's unit catalogue. Its generic source is the per-sv-block **dosage features** the model already builds for any locus (today they give copy counts per allele, `dosage_feature_copies.tsv`); π extends them to located, oriented occurrences. Presence/absence sv-allele cores that are not repeat units are units too |
| `C(block)` | a context anchor | the model's `context` blocks; they fix order and the frame |
| `J(left, right, w)` | a junction that is **not** one of the model's standard adjacencies: the two flanking unit regions and the breakpoint **window** `w` | the window lies between the consecutive informative sites that distinguish the two flanking units. A window, not a point: deterministic, and blind to differences inside it |
| `X(length class)` | content no unit or context explains (novel sequence of at least SV length) | keeps the projection total; makes the haplotype **unclassifiable** for evaluation |

**Class-defining sites.** Each unit has a fixed set of positions in its own coordinate frame. The set is derived from the
model's panel copies of that unit by **one declared rule** (Davide's choice; section 7).
- π records the alleles at those positions only.
- Substitutions elsewhere, and indels below SV length outside those positions, do not enter π. This is the "declared small
  difference", now a position set.
- A short difference at a class-defining site (for example a C4A/C4B isotype residue) always changes π.

**Representation preserved:** for each haplotype, the U-token count (**copy number**), their sequence (**order**),
`strand` (**orientation**), the J tokens (**meaningful junctions**) and `v` (**class**).

**Determinism and ambiguity.**
- The parse of a sequence into tokens follows a fixed rule. Ties are broken only between parses that give the **same**
  token string.
- If two parses give different token strings, π returns the **set**. Examples: an orientation that cannot be resolved in a
  palindromic arrangement; a breakpoint two unit pairs explain equally; an occurrence on the margin of the unit-detection
  rule.
- Ambiguity is reported, never resolved arbitrarily.

**Locus annotations enrich π; they are not a prerequisite.**
- An annotation may name tokens (`C4A`, `HERV-L`, `CYP2D7-hybrid`) and may add class-defining sites.
- Without one, tokens carry model identifiers (for example `U(sv:7/f0,+,v)`), and any locus the model can build can be
  called and evaluated.

## 3. The genotype and phase quantities

- **Interval genotype:** the unordered pair `{π(h1), π(h2)}` per structural interval. Within-interval **allocation** is
  part of it, so C4 3+1 and 2+2 are different genotypes.
- **Phase between intervals:** which interval genotype components lie on the same homologue. It is assessed separately,
  never part of the genotype, and never able to downgrade it.

## 4. Reporting unit and the locus target

**Structural intervals.**
- **Generated by a fixed rule from the model:** consecutive blocks spanning each sv block with its flanking context, merged
  where units, junctions or declared compound intervals overlap. This is the same model-bound artifact type as
  `compound_interval_declarations_v1`.
- **One row per interval** gives its genotype, support and status (detail).

**Locus event.** E = "every reported interval genotype equals truth's".
- Each interval is compared as an unordered pair, so E does **not** require resolved phase between intervals.

**The primary locus-quality decision** is about E: `locus_quality` ∈ {`CONFIDENT`, `NOT_CONFIDENT`, `NOT_ASSESSED`}, with
one reason code.

**Its input is the joint pool support, never a product of marginals.**
- `support_E_lower`: the pool mass of groups whose **every** recorded labelling projects, **in every interval**, to the
  predicted interval genotype. This is the contract-v2 lower bound applied to the tuple of interval genotypes.
- `support_E_upper`: the matching contract-v2 upper bound.
- Intervals are not assumed independent. Confident intervals do not imply a confident locus.

**Pool support is an input to calibration, not a probability of correctness.**
- `locus_quality` is a calibrated function of truth-free features: `support_E` bounds; missing mass; resolution status;
  search completeness and budget stops; interval count; depth-fit diagnostics; JP1.1 eligibility.
- The function is fitted on approved held-out data (section 6).
- No single support threshold is precommitted. Known errors already sit at extreme support: in M9-C, ACOT sv:5 was wrong in
  all six LOO runs at extreme confidence. If calibration cannot separate such errors from correct calls, that regime stays
  `NOT_CONFIDENT` or `NOT_ASSESSED`, and the report says so.

**The prediction is always kept.** The decision annotates; it never removes, replaces or masks it.

**Phase** gets its own support (the mass of the predicted cis configuration, up to global homologue exchange) and its own
decision, calibrated separately.

## 5. Proposed interface

| output | content |
|---|---|
| `<prefix>.structure.tsv`, locus row | **one prediction** (the interval genotypes as token strings, readable names where an annotation exists); **`locus_quality` with one reason code**; `support_E` labelled uncalibrated until calibrated |
| `<prefix>.structure.tsv`, interval rows | per-interval genotype, support, set-valued or unclassifiable flags |
| `<prefix>.phase.tsv` | one row per adjacent interval pair: predicted cis configuration, `phase_quality`, support |
| existing outputs (unchanged; diagnostics) | sequence prediction and FASTA under the unchanged policy; `sequence_blocks`, `sequence_alternatives`, `sequence_compound`, dosage audit, rounds and ties |

## 6. Evaluation against truth (offline; truth never reaches the caller)

**Per held-out sample and interval:**
1. Project each truth haplotype with the same π.
2. Classify **representability** of the truth diplotype against the **active panel** (after LOO or LZO curation). Each
   haplotype is classified, and the diplotype takes the less representable class:

   | class | definition | how it is determined |
   |---|---|---|
   | R1 `INTACT_PANEL` | π(truth) equals π of an intact active-panel haplotype | direct |
   | R2 `ALLOWED_MOSAIC` | not R1, but equals π of a haplotype in the caller's hypothesis space (a block-allele combination spelled from active-panel alleles) | offline oracle search over the spellable space; truth defines only the target signature; diagnostic only |
   | R3 `UNREPRESENTABLE` | no spellable haplotype has π(truth) | the oracle proves none exists |
   | R4 `UNCLASSIFIABLE` | π(truth) contains `X`, or is set-valued | direct |
   | `UNDETERMINED` | the oracle did not decide within its declared bounds | reported, never dropped |

3. **Search outcome** for R1 and R2, from the run's explored pool:
   - `PREDICTED`: the prediction's signature equals truth's;
   - `IN_POOL_NOT_PREDICTED`: a ranking or likelihood loss;
   - `NOT_EXPLORED`: missed by search, although representable.
4. **Correctness depends only on agreement with truth, never on representability:**
   - `EQUAL`: both signatures are singletons and equal;
   - `COMPATIBLE`: a set-valued side overlaps the other, but equality is not established. Not counted as correct;
     reported;
   - `DISJOINT`: wrong.

   A `CONFIDENT` prediction is an error iff it is not `EQUAL`. Truth outside the panel (R2 or R3) does not make a
   prediction an error by itself. A coarse π can make an off-panel truth `EQUAL` to a prediction, and then the prediction
   is correct.

**Reporting: one full contingency table per locus and coverage.**
- **Rows:** R1, R2, R3, R4, UNDETERMINED.
- **Columns:** `CONFIDENT` × {EQUAL, COMPATIBLE, DISJOINT}; `NOT_CONFIDENT` × the same; `NOT_ASSESSED`; `NO_PREDICTION`.

**Headline numbers:**
- confident-call error rate;
- callability;
- prediction accuracy regardless of the decision.

Each headline number is printed with its numerator and its denominator. Every row stays in the denominators. The counts of
confident calls on R3, R4 and UNDETERMINED cases are printed beside each headline, so unrepresentable or unscorable cases
can never inflate a rate silently.

**Locus level:** the same table for E. Allocation errors (equal multiset of U tokens, different per-homologue split) are
reported as their own sub-count of DISJOINT.

**Calibration protocol, when approved:**
- held-out families frozen beforehand;
- features fixed before fitting;
- the calibration split fits the map;
- a separate test split reports once, with intervals;
- separability of known error classes is reported **before** any threshold is accepted.

## 7. Choices that genuinely need Davide

| # | choice | why it is not a technical default | recommendation (no value chosen) |
|---|---|---|---|
| 1 | **The class-defining-site rule** | it decides which short differences count as structure | paralog-specific sites derived from the panel (sites whose alleles separate panel copies of a unit by lineage or position), enriched by locus annotations. The stricter alternative is every panel-segregating site in the unit |
| 2 | **Unit-detection rule and its margin** (when a partial or divergent copy counts as an occurrence; when content becomes `X`) | it decides copy number for degraded copies; occurrences inside the margin become set-valued | a declared rule on the unit's own sequence, with its margin reported. The value is to be approved, not fitted |
| 3 | **Acceptable confident-call error rate, and the held-out design** (families, LOO vs LZO, coverages) | it determines the calibrated decision | — |

**Already recommended in this specification; Davide may veto them:**
- allocation is part of the interval genotype;
- phase is separate;
- there is one locus decision about E;
- `COMPATIBLE` is never counted as correct.
