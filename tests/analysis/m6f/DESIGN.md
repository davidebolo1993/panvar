# M6-F: bounded predictive test of within-unit sequence copying (frozen before any unit)

**Ancestry.** Branch `m6-f-copying-predictive` from `2e7b6e6` (M6-E result), with the M6-E qualifications appended to its
RESULT.md.

**Scope.** Only a panel-sequence predictive comparison.
- **Excluded:** read likelihood, reads, truth-selected constructions, production changes, flags, reconstruction, merge,
  push and cleanup.
- **Not modelled:** cross-copy coupling, a structural prior, segment escape, indels and read placement.

**What is measured.** The predictive probability of a held-out component's bases, **conditional on** its supplied role,
length class and the declared alignment mask. It is not full-sequence probability, and not genotype accuracy.

## Inputs (digests in `m6f_run.py`; read only through allowlists)

- **M1-D panel inventory** (NA19650-LOO model): `inventory.{tsv,fa}` and the panel block carriers. The same digests as in
  M6-E.
- **Pedigree (source-backed):** IGSR `20130606_g1k_3202_samples_ped_population.txt`. URL, retrieval time and SHA-256
  (`4e164b71…`) are in `inputs/SOURCE.tsv`. It documents the HGSVC trios SH032 (HG00512, HG00513 → HG00514), PR05
  (HG00731, HG00732 → HG00733) and Y117 (NA19239, NA19238 → NA19240).
  - Families are connected components of FamilyID co-membership plus parent links. For example, NA19238 is also NA18913's
    mother, under FamilyID Y028.
  - Donors absent from the pedigree (NA21487, NA24385, chm13v2, grch38) are their own families.
  - **Relationships are never inferred from sequence identity.**
- **`DEVELOPMENT_DONORS.tsv`** (committed here, with sources): 54 donors previously inspected in development. These are
  donors named in a memory note or a `tests/analysis` design or result, plus the template donors of all ten campaign
  11fe681 mosaic plans, resolved from their `stage.log` manifests.

## Holdout (`prepare`; chosen before any scoring)

- **Training-only families:**
  - every family containing a development donor;
  - the reference genomes chm13v2 and grch38;
  - HG00096, the coordinate-reference donor (anchor source and reference path).
- **Evaluation:** the first **8** eligible families ordered by sha256("M6-F:" + family id).
- **Coordinate references.** One per stratum (prefix; period long and short; tail long and short): the unique component
  with the largest positive weight among copies of training-only families. Ties go to more haplotypes, then the smallest
  SHA-256. Each component is aligned only to its stratum's reference, with the M6-E `unit_correspondence` tool (edlib,
  global, unit cost, path mode). An evaluation copy identical to a reference is genuine multiplicity and is kept.
- **Weighting policy (family-aware, founder; used everywhere, including inner validation targets):** a copy has weight 0
  if its donor's father and mother are both panel donors (the trio children HG00514, HG00733 and NA19240), and weight 1
  otherwise. Identical sequences in unrelated donors keep their full multiplicity.
- **Per fold** (held-out family F, stratum s):
  - **Library:** the unique positive-weight components of s from families other than F, with weights summed. All of F's
    copies are removed from every library, mask, site inventory and parameter selection.
  - **Training mask:** positions with simple correspondence in every library component (no indel or large-deletion edge
    within 10 bp, base present).
  - **Scored positions:** training mask ∩ the held-out copy's own simple positions. Excluded positions are reported by
    reason (training mask; held-out alignment). Training-variable and training-monomorphic positions are **both** scored,
    so novel substitutions count.
  - A copy is **unassessable** if there is no training library or no scored position.

## Models (`m6f_core.py`)

**Shared residual emission:** e(y | b) = 1 − ε if y = b, else ε/3.

| model | probability of the observed bases |
|---|---|
| C0, independent sites | Π_x Σ_k w_k e(y_x \| b_k(x)) |
| C1, whole-component mixture | Σ_k w_k Π_x e(y_x \| b_k(x)), summed, never the nearest component |
| C2, Li–Stephens copying | at distance d between consecutive scored positions, the donor is redrawn from w with probability 1 − exp(−ρd), otherwise kept |

- **Exact compression.** At training-monomorphic positions the emission is the same for every hidden donor, so it factors
  out. The redraw kernels compose. Only training-variable positions enter the forward pass, at their reference distances.
- **Grids (frozen):** ε ∈ {1e-4, 3e-4, 1e-3, 3e-3, 1e-2}; ρ ∈ {1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3} per bp.
- **Parameter selection (training families only):** within each fold and stratum, leave-one-family-out over the training
  families. Masks and libraries are rebuilt from the inner training set. The validation score is weighted by the founder
  policy. Each model takes the grid maximum; ties within 1e-9 go to the smaller ε, then the smaller ρ.

## Gates (`gates.py`; after `align`, before `evaluate`; all must pass)

| gate | check |
|---|---|
| G1 | normalisation: probabilities over all 4⁵ sequences of a 3-donor fixture (with monomorphic positions) sum to 1 (≤ 1e-12), every model and grid point |
| G2 | compressed forward = exhaustive hidden-path enumeration (3⁵ paths; ≤ 1e-10); compressed = uncompressed forward on a 60-position fixture (≤ 1e-9) |
| G3 | ρ = 0 reproduces C1 (≤ 1e-12) |
| G4 | redraw at every scored site reproduces C0 (≤ 1e-12) |
| G5 | held-out isolation: for every evaluation family × stratum fold, the fit (library, mask, inner validation totals, selected parameters) is identical when that family's copies are replaced by scrambled synthetic components with arbitrary weights |
| G6 | aliases: splitting a donor into two identical halves preserves every model's probability (fixture ≤ 1e-12; one real training-only library ≤ 1e-9) |
| G7 | in every `score_all` call, C0, C1 and C2 see one observation digest |

## Evaluation and decision (`evaluate.py`, frozen)

**Per copy:** C0, C1 and C2 log scores at the selected parameters, the paired differences, per-base scores, coverage,
exclusions, novel substitutions at training-monomorphic positions, and the observation digest. Results are shown per
family × stratum, per family, per stratum and pooled, with runtime.

**Flags,** per comparison (C1 − C0, C2 − C0, C2 − C1), pooled per family:
- **meaningful:** mean paired difference ≥ 2.0 nats per evaluated copy;
- **consistent:** positive in ≥ 75 % of evaluated families, **and** the largest family's share of the summed positive
  differences ≤ 0.5.

**Decision:**
- **C2_IMPROVES_CONSISTENTLY** if C2 − C1 is meaningful and consistent;
- otherwise **C1_SUFFICES** if C1 − C0 is meaningful and consistent (prefer the simpler model);
- otherwise **PARK**.

Strata get the same flags, descriptively. A positive result validates a sequence-prior ingredient only.

## Bounds

- 600 s counting every attempt; 4 threads.
- **Declared timeouts:** prepare 30 s, align 120 s, gates 120 s, evaluate 240 s.
- Failures and incomplete units stay recorded.
- Before the freeze, the scoring core was dry-checked once on fabricated fixtures in the scratchpad (normalisation and
  timing; no panel data). This was not a unit.
