# M10-B.1: hardening and qualification of the structural-projection prototype (frozen before execution)

**Context.** M10-B stays diagnostic. Its long-node skeleton is **not** the production structural target. M10-A remains the
computational checkpoint.
- **Budget:** 600 s, including 15 s of disclosed development.
- **Scope:** no reads, genotyping, families, fitted thresholds, caller integration or default changes.
- **Branch:** `m10-b-structural-projection`. **Results:** `tests/results/m10b1_projection_hardening_20261020/`.

## 1. The boundary transformation, from the implementation (`src/genotype_locus.cpp`)

1. **Recorded walks are the pre-normalisation walks.** `walk_block` records a template's walk across a chain block. For a
   variant (bubble) block, the walk **includes the bubble's source and sink nodes**. A context (backbone) block takes only
   what lies strictly between the boundaries.
2. **Abutting bubbles.** When two bubbles share a boundary node, there is no backbone between them. That node is dropped
   from one side of the walk (`trim_node`), so consecutive walks never share a node.
3. **Divergence trim** (`trim_variant_block_to_divergence`).
   - Every allele of a variant block is cut by the prefix P and suffix S common to **all its source alleles** (curated ones
     included). The cut bytes are appended to the previous context block and prepended to the next.
   - A block carrying a bypass allele is not trimmed. It is flagged `representation_dependent` if it had common bytes.
   - The amounts are the report's `trimmed_prefix_bp` / `trimmed_suffix_bp`.
4. **The relation, therefore:**
   - variant block: `normalized = spell(walk)[P : len − S]`;
   - context block: `normalized = S(prev) + spell(walk) + P(next)`;
   - whole template: unchanged, because bytes only move between blocks.

   Individual walks therefore **do not** respell their normalized alleles when P or S > 0 (Sol's observation), while whole
   haplotypes do.

**Extended dump** (`m10b_walk_dump`, additive): `blocks.tsv` (the report's trims and flags), `sources.tsv` (per allele:
walk / empty / bypass / other source counts) and `segments.tsv` (invariant segments). The five M10-B files must stay
byte-identical.

**Audit** (`m10b1.py audit`, all six models). Per block, the relation is verified for every allele and walk. Also reported:
- raw per-block respelling;
- whole-template respelling;
- source kinds;
- non-empty invariant segments;
- metadata-missing alleles;
- case-scope support;
- degeneracy projections under version 2.

**Projection scopes** (`project_v2.scope`):

| scope | projectable when |
|---|---|
| whole haplotype | every block's relation is verified (no clip) |
| compound interval or single block | relation verified **and** no bytes enter from outside (an end context block next to a trimmed variant block takes bytes whose nodes are in a neighbour's walk, which the scope does not determine) **and** no skeleton node straddles the clip of an end variant block |

Otherwise the scope is **UNSUPPORTED**. Correspondence is never assumed.

## 2. `project_v2.py` (version 1 unchanged)

- **Missing metadata gives INCOMPLETE.** Recorded walks fewer than the walk sources, or no walk and no valid empty record,
  make any projection using the allele INCOMPLETE. A valid empty record ('empty' / 'bypass', no walk) contributes an empty
  walk.
- **Class identity** is the canonical sequence **bytes**. The display ids come from a registry that disambiguates hash
  collisions. The hash function is injectable.
- **Comparison categories:** UNIQUE_AGREEMENT, SET_IDENTICAL, COMPATIBLE, DISAGREEMENT, UNSUPPORTED.

## 3. Fixtures (`test_project_v2.py`; version 1's `test_project.py` must still pass)

- **Ported:** the M10-B fixtures.
- **Missing metadata:** M1 (fully missing → INCOMPLETE; two unrelated → UNSUPPORTED), M2 (valid empty), M3 (partially
  missing alternative walks), and M4 (version 1's defect reproduced and documented).
- **Adversarial, reported as preserved or LIMITATION** and never tuned:
  - A1: the same unit as one long node against short nodes;
  - A2: a 55 bp insertion made of short nodes;
  - A3: a neutral SNP inside a long node against a separate short node;
  - A4: shifted boundaries;
  - A5: a declared multi-base class site against a neutral short change;
  - A5b: a 1-bp site declared by sequence class.
- **Forced collision:** C1, a constant hash, under which every verdict must be unchanged.
- **Scope:** S0–S4 (whole haplotype over trims; clip on a node boundary; clip cutting a long node; bytes entering a
  context scope; an unreproducible trim).

## 4. Units (manifest P)

`build`, `gate.fixtures_v1`, `gate.fixtures_v2`, `dump.<6>`, `audit`, `cases`.

`cases` replays M10-B's **frozen** case set with version 2, beside M10-B's version-1 verdicts, and records each change. Not
counted as establishing one structural genotype: SET_IDENTICAL, COMPATIBLE and UNSUPPORTED. UNSUPPORTED comparisons stay in
the tables, and so in any later evaluation denominator.
