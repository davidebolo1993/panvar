# M6-E: panel-grounded sequence-variation model, feasibility before implementation (frozen before any unit)

**Ancestry.** Branch `m6-e-panel-variation` from `2e35609` (M6-D result).

**Status of M6-D, recorded at Sol's review.** M6-D is accepted as a successful diagnostic, not a caller improvement. Both
nuisance priors stay experimental and unintegrated.

**Scope.** Panel-only measurement, then a written model specification.
- **Excluded:** read-likelihood runs, model fitting, production changes, flags, reconstruction, merge, push and cleanup.

## Inputs (panel only; digests frozen in `m6e_run.py`)

**Source.** M1-D's active-panel C4 sv:7 inventory (unit `inventory.a1`, frozen in `27eab4f`), written from the **NA19650
LOO model**:
- `genotype_c4_production_11fe681/graphs/NA19650.loo.gfa`, model fingerprint `14851130f453a071`;
- 117 sv:7 alleles carried by active-panel templates;
- 116 alleles are partitioned at the phase-A anchor from allele_64 into a prefix, 0–2 periods and a tail.

| file | role |
|---|---|
| `inventory.a1/inventory.{tsv,fa,anchors.tsv,log.tsv}` | component table and sequences. The digests must equal `tests/analysis/m1d/INVENTORY_FROZEN.tsv` |
| `seq.panel/seq.sv:7.alleles.fa`, `seq.block_sequences.tsv`, `run.log` | panel allele sequences, carrier templates, and the model binding |

**Verified by `prepare`** (fails if any check fails):
- the frozen digests;
- the LOO-model binding;
- no NA19650 template;
- allele_75 excluded and reported;
- byte-for-byte reconstruction of all 116 included alleles from their components;
- inventory carriers equal the panel block carriers.

**Never opened.** The scripts enforce an allowlist. Excluded:
- M1-D's truth-selected units (`reach.*`, `construct.*`, `cand.dist.*`, `score.*`, `summary`);
- M1-E's oracle patches;
- excluded-donor sequences, reads and evaluation outcomes.

The inventory is panel-only by construction (M1-D froze it before any truth unit). The later constructions were selected
by truth distance and are oracle diagnostics.

**Exclusions and boundaries preserved.**
- allele_75 stays excluded, since its anchor phases disagree.
- Prefixes, periods and tails are analysed **separately**, each against its own reference.
- Periods and tails are related only through one explicit tail-reference vs period-reference alignment, and only where
  that alignment has simple correspondence.

## Units (300 s, 4 threads; one at a time)

| unit | timeout | does |
|---|---|---|
| prepare | 30 s | verify bindings; write the component table (copies, haplotypes, donors, positions, length class) and the copy table; choose references; write the alignment pairs |
| align | 120 s | `unit_correspondence` (new diagnostic, linking the existing edlib library): unbanded global unit-cost alignment of each component to its type's reference, in path mode. Writes SNV, INS and DEL events in reference coordinates |
| analyse | 90 s | correspondence classes, site catalogue, frequencies, co-occurrence (below) |

**Reference rule (frozen, panel-only).** Within each type, the unique component with the most physical copies. Ties go
to more haplotypes, then the type's majority length class, then the smallest SHA-256.

**Pairs:**
- every unique prefix, period and tail against its type's reference;
- the tail reference against the period reference (cross map);
- every period against the most frequent **short** period. This is a reference-sensitivity arm, for correspondence
  classes and site counts only.

**Length class.** Long ≥ 29 kb, otherwise short. The inventory has no component between 26.5 and 32.5 kb.

## Measurements (`analyse.py`, frozen; W = 10 bp, large deletion ≥ 50 bp)

1. **Components.** Distinct sequences with length, position (prefix / period1 / period2 / tail), length class, and
   carriers counted separately by **physical copy** (component occurrences × carrier templates), **haplotype** (templates)
   and **donor**.
2. **Correspondence**, per reference position:
   - **SIMPLE:** every component has a base, and no indel or large-deletion edge lies within W.
   - **CONDITIONAL:** inside a ≥ 50-bp deletion in some components, and simple in the components that have it.
   - **UNCERTAIN:** everything else.

   Also reported: merged intervals, components with end gaps, and SNVs clustered within W.
3. **Sites.** Each SNV position, with alleles counted by copy, haplotype and donor. Also: multiallelic sites, and the
   minor-allele spectra by each counting unit.
4. **Co-occurrence.**
   - **(a) Within copies:** r² across physical copies at informative SIMPLE sites (minor allele ≥ 3 copies), overall and
     by distance. The null permutes each site independently across copies (20 permutations, seed 20261019). Also:
     perfectly linked groups and distinct copy patterns.
   - **(b) Within haplotypes:** sites must have simple correspondence in both period and tail coordinates and be mapped.
     Copy pairs in the same haplotype are compared against all pairs from different donors with the same type
     composition: Hamming distance, median percentile (0.5 = no excess similarity) and the fraction identical. Length-class
     concordance likewise.
   - **(c) Unique-sequence structure:** pairwise Hamming at SIMPLE sites, and single-linkage clusters at 0, 5, 20 and 50
     differences.
5. **Multiplicity vs duplicate representation.** Alleles with more than one carrier (same donor vs different donors);
   components recurring across alleles; the HG00096 templates.

## Deliverables written after the measurements (no execution)

- **One candidate-independent generative model.** It separates structure and copy number, shared sequence families and
  residual variation. It states which latent choices are shared and which are copy-specific, the normalised distribution
  and its induced copy-number effects, and its invariance to base annotation. It distinguishes multiplicity from
  duplicate representation, and states how it interacts with the existing copying prior. No parameter comes from F2.
- **A donor-held-out, panel-sequence predictive check** (future work, not run). A held-out donor's copies are all
  excluded from fitting; donors used in development are excluded from the check.
- **A state-count and runtime estimate for C4** from existing measured costs. It identifies the dependencies created by
  alternative fragment placements and shared sequence states. Any approximation is named, with how its error would be
  assessed. Exact site-local factorisation is not assumed.
- **A decision:** what the panel supports, the model, whether one bounded implementation experiment is practical, and the
  smallest decisive next test (or parking).

Unresolved correspondence and resource limits are recorded, not expanded.
