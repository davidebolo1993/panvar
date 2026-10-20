# M10-B: bounded offline structural-projection prototype (frozen before the real cases)

**Purpose.** Can the information the existing model holds define the structural genotype of the quality target (version
2), before anything is implemented or calibrated? This is a diagnostic, not a caller.

- **Branch:** `m10-b-structural-projection`, from `m10-a-consolidation` `469c7a8`.
- **Results:** `tests/results/m10b_structural_projection_20261020/`.
- **Scope:** no reads, no genotyping, no new families, no thresholds fitted. Saved models and saved sequences only.
- **Kept:** M10-A's configuration and indexes stay as the development checkpoint.

## 1. What the model supplies (established during pre-freeze inspection; disclosed)

| model information | what it is | located / oriented? |
|---|---|---|
| block allele sequences | sequence-unique per block (the model canonicalises them); shared prefixes and suffixes are moved into context blocks (`boundary_bases_moved_to_context`), which is why synonymous labellings exist | sequence only |
| allele walks (`source_ids` "walk:") | the oriented graph-node path(s) of each allele. Over a whole template, the walks spell exactly the concatenated allele sequences, and consecutive blocks share no node | **yes**: every step is a located node visit with a strand |
| repeat units | not a separate decomposition in the model. A tandem copy re-visits the **same** graph nodes: C4 sv:7 alleles visit about 26 kb of nodes twice. Panphorte folding is a separate graph-rewriting command, not part of the genotype model | located copies through walk multiplicity |
| dosage census | per-block route (repeated unit / pair / not sensitive) from node multiplicity and repeated bases | per block, not per occurrence |
| pooled dosage features | groups of varying canonical k-mers that share one copy-count column per allele | **no**: unlocated counts; **not** equated with oriented occurrences |

**Two inspected obstacles for a walk-based projection:**
1. **Distinct nodes with identical sequence.** In CYP2D6, node `1829+` and node `1825-` spell the same 12 bp.
2. **One sequence, different structural placements.** CYP2D6 degeneracy `cyp:1`: one labelling places a tandem expansion
   as a 59 bp node in sv:3, the other as a 58 bp node in sv:5. The total sequence is identical.

**Measured invariance across the dev models' degeneracy equivalences:**

| model | raw walks identical | long-node skeletons identical |
|---|---|---|
| ACOT | 0 / 2,458 | 2,458 / 2,458 |
| CYP2D6 | 0 / 42 | 19 / 42 |

## 2. The projection (`project.py`, frozen)

π(H) is defined on a haplotype **sequence** H over blocks f..l, never on a labelling.
1. **Parses.** Every tuple of model alleles over f..l that spells H exactly.
2. **Skeleton** of a parse. The chosen alleles' walks, every recorded walk of each allele, reduced to their long-node steps,
   each step written as **(sequence class, orientation)**.
   - The sequence class is the sha256 of the node sequence's canonical strand. Node ids never appear.
   - **Long** means at least `ANCHOR_BP = 50`. This is the configured `min_variant_bp`, reused as an **explicitly
     experimental diagnostic convention**, not as approved production semantics.
   - **Separately and explicitly**, nodes whose sequence class is declared class-defining are kept whatever their length.
     The declaration is **empty** for every real case: no locus annotation is assumed. So a panel SNP is never structural
     by default.
3. **π(H)** is the **set** of skeletons over all parses and walk choices.
   - Identical sequences have identical projections **by construction**.
   - A set larger than one is **AMBIGUOUS**: the model places the same sequence in more than one structural way.
   - No exact parse gives **UNKNOWN**. Exceeding the computational guard (100,000 parses or skeletons) gives
     **OVERFLOW**. Both are reported, never truncated.

**Diploid comparison** (unordered per homologue, so allocation counts):

| verdict | condition |
|---|---|
| UNKNOWN | any side UNKNOWN or OVERFLOW; two UNKNOWNs are never equated |
| EQUAL | equal as multisets; EQUAL_SET_VALUED when a projection is set-valued (identical sets) |
| COMPATIBLE | some pairing intersects on both homologues |
| DISJOINT | otherwise. `DOSAGE_EQUAL_ALLOCATION_DIFFERS` flags an equal multiset of skeleton classes summed over the two homologues. Dosage alone never gives EQUAL |

No pairwise tolerance and no clustering anywhere.

**Descriptive copy profile** per skeleton: steps, repeated classes, maximum multiplicity, and multiplicity signature. It is
for reporting only.

## 3. Fixtures (`test_project.py`, 17 checks, a toy model in the dump format; the toy classes are fixture assumptions)

| fixture | expectation |
|---|---|
| F1 synonymous block split, plus a duplicate-sequence node | identical, RESOLVED |
| F2 neutral SNP | EQUAL (below the convention, undeclared) |
| F3 toy class-defining site, undeclared | EQUAL |
| F3 toy class-defining site, declared | DISJOINT |
| F4 copy number 2 vs 3 | DISJOINT |
| F5 3+1 vs 2+2 | DISJOINT with the dosage-equal flag |
| F5 homologue exchange | EQUAL |
| F6 order; orientation | DISJOINT each |
| F7 junction J1 vs J2 | DISJOINT |
| F8 allowed mosaic absent from the panel | RESOLVED, EQUAL to itself |
| F9 one sequence, two structural parses | AMBIGUOUS, labelling-invariant. Identical sets are EQUAL_SET_VALUED; intersecting sets are COMPATIBLE |
| F10 off-model sequence | UNKNOWN; UNKNOWN against UNKNOWN gives UNKNOWN |
| F11 equal class counts, different order | not EQUAL |

## 4. Units (manifest P) and the frozen case set

**Units:** `build` (`m10b_walk_dump`), `gate.fixtures`, `dump.{c4dev, c4loo, cypdev, cyploo, acotdev, acotloo}`,
`coverage`, `cases`.
- **Dump checks.** Each dump's fingerprint must equal the frozen value, or, for the two LOO models without one, the saved
  run's.
- **`coverage`** tabulates, per model: walk multiplicity and orientation, duplicate-sequence node classes, and alleles
  with differing walks. For every dev-model degeneracy equivalence it reports raw-walk, skeleton and parse-set invariance,
  and the RESOLVED/AMBIGUOUS split.

**Frozen case set** (inspected development cases, **not** accuracy validation):

| case | model | interval | haplotypes compared |
|---|---|---|---|
| C4 NA19650 allocation failure (M1-B `c4_10`) | `c4loo` | sv:7 | reported 3+1, saved 2+2, saved 1+2, truth (M1-B block sequences) |
| C4 HG00096 correct control (M10-A `D.c4_5`) | `c4dev` | whole locus | call vs the donor's intact panel templates |
| CYP2D6 NA19650 dosage failure (M1-B `cyp_5`) | `cyploo` | sv:9 | 5× reported 2+2, 30× twin 3+2, saved 2+3. **No truth block sequence was saved at sv:9** |
| CYP2D6 NA19650 synonymous labelling (M1-B `cyp_30`) | `cyploo` | sv:3..sv:5 | call vs truth (M9-I units) |
| ACOT L4 synonymous labelling (M9-H) | `acotloo` | sv:5..ctx:6 | JP on (interval exact) and JP off vs truth |
| ACOT L2, a few bases (M9-H) | `acotloo` | sv:5..ctx:6 | JP on (distance 3) and JP off (5) vs truth |

**Reported per distinction:** preserved, ignored for a declared reason, ambiguous, or unsupported.

## 5. Budget

**1,200 s,** including 45 s of disclosed pre-freeze development. Admission: spent + declared timeout ≤ 1,200. Failed units
stay recorded. No amendment relaxes a rule.
