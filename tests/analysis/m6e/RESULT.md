# M6-E result: the panel supports a lineage-copying unit model, not independent per-copy variation; one panel-only test is decisive before any read work

**Run.** Branch `m6-e-panel-variation`: freeze `916235a`, Amendment 1 (a guard defect in the final print), Amendment 2
(state-count and arrangement measurements). Results in `tests/results/m6e_panel_variation_20261019/`.

**Budget.** 6.2 of 300 s, 4 threads, peak RSS 0.26 GiB. One recorded failure: `analyse` exited 1 on its own read-back
guard and was rerun as `analyse.a1`.

**Scope.** No reads, likelihood runs, fitting, production change, flag, reconstruction, merge or push.

**Status of M6-D, recorded at Sol's review.** M6-D is accepted as a diagnostic, not a caller improvement. Both nuisance
priors stay experimental and unintegrated.

## 0. Inputs and bindings (`prepare`: all checks pass)

- **Source:** M1-D's panel-only C4 sv:7 inventory from the **NA19650 LOO model** (`NA19650.loo.gfa`, fingerprint
  `14851130f453a071`).
- **Digests** equal `INVENTORY_FROZEN.tsv` (frozen in `27eab4f`, before any truth unit).
- **No NA19650 template.** allele_75 is excluded and reported (its anchor phases disagree).
- **All 116 included alleles reconstruct byte for byte** from their components, and the carriers match the panel block.
- **Panel:** 128 templates of included alleles (129 with allele_75), 66 donors.
- **Never opened** (enforced by an allowlist): M1-D's truth-selected units (`reach`, `construct`, `cand.dist`, `score`) and
  M1-E's oracle patches. The inventory predates them, and they are oracle diagnostics.
- **Typed separately:** prefixes, periods and tails are aligned to their own references. Periods and tails are related only
  through one explicit tail-reference vs period-reference alignment.

## 1. What a sequence-variation prior must represent

**Components** (copy order: prefix, period1, period2, tail; carriers counted by physical copy, haplotype and donor):

| type | physical copies | unique sequences | appear in one copy only | max copies | in > 1 allele | length classes by copy (haplotypes; donors) |
|---|---|---|---|---|---|---|
| prefix (78 bp) | 128 | 10 | 3 | 62 | 6 | — |
| period | 127 | 106 | 88 | 4 | 13 | long 117 (111; 65), short 10 (10; 10) |
| tail | 128 | 104 | 85 | 4 | 10 | long 53 (53; 38), short 75 (75; 48) |

- Long and short differ by the ~6.4-kb presence/absence segment: 32.7 kb against 26.3 kb.
- Positions and carriers per component are in `runs/prepare/components.tsv` and `copies.tsv`.

**Structure by haplotype.**
- Copies per haplotype: 1 copy in 17 haplotypes, 2 in 95, 3 in 16.
- Length-class arrangement, by haplotype (donors): LS 63 (43), LL 32 (27), L 14 (13), LSS 7 (7), LLL 4 (4), S 3 (3),
  LSL 3 (3), LLS 2 (2).
- Diploid totals by donor: 2 copies in 1 donor, 3 in 14, 4 in 35, 5 in 10, 6 in 2.

**Correspondence.** W = 10 bp. Unit-cost global alignment to the frozen references:

| reference | SIMPLE | CONDITIONAL (absent in some) | UNCERTAIN | non-simple intervals |
|---|---|---|---|---|
| period_u13 (long, 32,738 bp) | **74.4 %** | 2.5 % | 23.1 % | 27 uncertain (7,551 bp) plus 5 conditional. Largest 6,113–10,066 (3,953 bp; 49 components), 4,082–6,066, 10,220–11,277 |
| tail_u9 (short, 26,268 bp) | **89.5 %** | 0 | 10.5 % | 20 (2,750 bp). Largest 3,402–5,699 (2,297 bp; 45 components) |
| period vs a short reference (sensitivity) | 88.7 % | 0 | 11.3 % | site count 264, against 270 with the long reference |
| tail reference ↔ period reference | 74.3 % of period positions mapped with simple correspondence | | | distance 6,833, of which 6,600 bp is deleted |
| prefix_u6 (78 bp) | 69 % | 0 | 31 % | one 24-bp indel region (2 components) |

- **Where simple site correspondence fails:** the long/short segment (~3.4–11.3 kb in period coordinates). Unit-cost
  alignment breaks that presence/absence into fragments surrounded by indel-dense flanks, so it is **not one simple event**
  here. There are also about 30 short indel intervals.
- **Within non-simple intervals:** 2–36 distinct local event signatures among the unique components. This is an upper
  bound, because co-optimal alignments can split one sequence into several signatures.
- **End gaps:** none at unit boundaries. The anchor-defined boundaries align cleanly.
- **Clustered sites:** 53 of 270 SIMPLE period sites and 54 of 239 tail sites have another SNV within 10 bp, which suggests
  complex events.

**Sites.**
- Periods: 270 SIMPLE + 4 CONDITIONAL + 76 UNCERTAIN SNV positions (3 multiallelic).
- Tails: 239 + 0 + 31.

Minor-allele spectrum of the non-uncertain sites:

| type and count | 1 | 2 | 3–5 | 6–20 | > 20 |
|---|---|---|---|---|---|
| period, by copy | 74 | 37 | 66 | 55 | 42 |
| period, by haplotype | 75 | 36 | 68 | 54 | 41 |
| period, by donor | 75 | 36 | 68 | 57 | 38 |
| tail, by copy | 83 | 20 | 27 | 45 | 64 |
| tail, by donor | 83 | 21 | 27 | 60 | 48 |

Copy and haplotype counts nearly coincide because most haplotypes carry one period. Donor counts diverge at common sites,
because both haplotypes of a donor share alleles.

**Co-occurrence.**
- **(a) Strong local linkage within copies.**
  - Periods: 161 informative sites. r² ≥ 0.8 for 3.1 % of pairs overall, **15.2 % within 1 kb**, 2.5 % at 1–5 kb and 0.7 %
    beyond 5 kb. Under the independent-site null the fraction is **0.0000** (20 permutations, max 0).
  - Perfectly linked groups: 108 for 161 sites; the largest has **24 sites**.
  - Distinct copy patterns: 89 over 127 copies.
  - Tails: r² ≥ 0.8 for 5.5 % of pairs within 1 kb against a null of 0.
- **Local diversity is low.** Among 106 unique periods, a 400-bp window has a **median of 3 distinct local patterns** (90th
  percentile 9, max 17). An 800-bp window has a median of 6 (90th percentile 17, max 25). Tails: 4 / 9 / 23 and
  7 / 16 / 34.
- **Unique-sequence structure.** Pairwise Hamming distance at SIMPLE sites ranges from 0 to 114 (median 34).
  Single-linkage clusters: 103 at ≤ 0 differences, 33 at ≤ 5, 8 at ≤ 20 (7 for tails) and 2 at ≤ 50.
- **(b) Copies within a haplotype are neither identical nor independent.** These comparisons use 201 informative sites with
  simple correspondence in both coordinate frames.
  - **Period–tail pairs** (127): mean Hamming 81.3 against 85.7 for between-donor pairs; median percentile **0.32**; none
    identical.
  - **Period–period pairs** (16): 27.1 against 29.8; median percentile **0.20**; none identical (0.8 % of between-donor
    pairs are identical).
  - **Length-class concordance** of period–period pairs in the same haplotype is **0.375**, against 0.855 between donors.
    The two periods of 3-copy haplotypes are usually of different classes, so the structure is arranged, not independently
    drawn.
  - These are descriptive statistics with small n for period–period pairs; no significance claim is made.

**Multiplicity vs duplicate representation.**
- 10 alleles have more than one carrier, and 29 unique components occur in more than one allele.
- Several shared alleles join donors that, **from public pedigree knowledge rather than panel metadata**, are HGSVC trio
  members:
  - YRI NA19238/NA19239 → NA19240: allele_4, allele_10;
  - PUR HG00731/HG00732 → HG00733: allele_45 and part of allele_79;
  - CHS HG00513 → HG00514: allele_73.
- These are **transmitted copies of one chromosome lineage**, not independent recurrences.
- Other sharing is between unrelated donors and is empirical multiplicity: for example allele_79 in NA18534 and NA21487,
  and allele_28, allele_47, allele_65 and allele_77.
- The panel carries no pedigree, so it cannot separate these two cases by itself.

**Is independent per-copy variation an adequate approximation? No.**
- Sites are strongly linked locally: r² ≥ 0.8 in 15 % of pairs within 1 kb, against a null of 0.
- The local sequence space is a handful of patterns, not 2^sites.
- Copies in a haplotype are mildly more similar than random, and their length classes are arranged.
- The copies are not identical either: no within-haplotype pair is identical at the informative sites.

## 2. Proposed model: a panel-grounded lineage-copying unit model (specification only; nothing fitted)

For each haplotype h of a diplotype, the model has three layers.

**L1. Structure and copy number (shared by all copies of h).**
- Latent choices: the copy count n_h ∈ {1, 2, 3} and the arrangement a_h ∈ {L, S}^{n_h}. The last copy takes the tail
  role; earlier copies take period roles.
- Prior: π(n, a) ∝ c(n, a) + α. Here c counts distinct **chromosome lineages** carrying the arrangement, with transmitted
  duplicates counted once; α = 1/2 over the 14 declared patterns with n ≤ 3, fixed in advance.
- The two haplotypes are independent given the panel (an assumption; the held-out check tests it against diploid totals).
- This **replaces** the uniform catalogue weights, which induced 2/9, 3/9 and 4/9 on totals 3, 4 and 5.
- **Induced copy-number prior:** P(n_h) is about 0.13 / 0.74 / 0.13 for n = 1 / 2 / 3. The diploid total is its
  convolution.

**L2. Shared sequence families or ancestral content (copy-specific paths with one shared coupling).**
- Each copy j carries a lineage path Z_j(x) over its role's unit coordinate. The path is a Li–Stephens copying process over
  the library of that role's panel components.
- The start is drawn from lineage weights w_k: distinct-lineage frequencies, with duplicates counted once.
- The path switches between consecutive positions with probability 1 − exp(−ρ·d), choosing the new component by w.
- Inside the long/short segment, the library is restricted to components of the copy's class from L1.
- **Shared across copies of h:** one coupling κ. With probability κ, copy j + 1 starts in copy j's lineage cluster
  (single-linkage at ≤ 20 differences, panel-defined); otherwise it starts from w. κ = 0 is the independent-copy special
  case.
- **Copy-specific:** the path itself and the residuals.

**L3. Residual variation (copy-specific).**
- At SIMPLE and CONDITIONAL positions the base equals the donor component's base with probability 1 − ε_s; otherwise it is
  one of the other three bases, uniformly. ε_s takes one of two declared rates: ε_poly at panel-polymorphic positions and
  ε_mono elsewhere.
- Inside non-simple intervals (**named approximation: segment copy**), the copy takes its donor's whole segment exactly with
  probability 1 − η. Otherwise it takes an explicit **novel-segment escape** state, which is reported, never spelled.

**Normalised distribution.**

  P(G, X) = Π_h π(n_h, a_h) · Π_j P_LS(X_hj | role_j, class_j, start coupling to copy j − 1; w, ρ, κ, ε, η)

- Each P_LS is a normalised HMM over that copy's declared sequence space: SIMPLE and CONDITIONAL alleles × segment choices ×
  escape.
- π is normalised over the declared arrangements. The joint therefore sums to one.
- The support is the panel coordinate system. Indels outside the segment catalogue, and sequence at uncertain positions
  beyond the escape, are **not** represented.
- **Invariance to base annotation.** The probability is a function of the sequence state X and the panel library only. A
  candidate is just a point X, so changing its base annotation cannot change the mass of the same biological state.
  M6-D's gate HB holds by construction, and M6-D's common-content control is the degenerate case of one donor with no
  switching.

**Induced copy-number effects (stated, to be reported in every dosage decision).**
1. **The structural prior π** is the explicit copy-number prior.
2. **The sequence layer is normalised per copy.** An extra copy whose content the reads do not identify costs nothing in
   expectation.
3. **An extra copy whose content the reads do identify pays its description length** (the KL term of M6-D's identity):
   - about −log w(lineage), which is 3.4–4.8 nats for one unique component and about log 8 ≈ 2 nats if only the ≤ 20-difference
     cluster is identified;
   - plus a cost per switch from ρ;
   - plus log(3/ε) per residual difference.

   This is the principled generalisation of M6-D's per-copy charge for rare content. It is a real dosage pressure that is
   larger for rare or novel content, so the log m = E_q log L − KL decomposition must accompany every dosage call.

**Empirical multiplicity vs duplicate representation.**
- Lineage weights w and arrangement counts c count each **distinct chromosome lineage** once.
- Identical sequences in unrelated donors keep their multiplicity, since that is population frequency. Repeated identical
  copies within a haplotype are real copy-level multiplicity.
- Transmitted duplicates (parent–child) count once. This needs a declared **pedigree table**; the panel metadata does not
  provide one.

**Interaction with the existing copying prior (no double counting).**
- The existing Li–Stephens copying prior over panel templates factorises at the block boundary:
  P_copy(path) = P_copy(outside sv:7) · P_copy(sv:7 allele | entry/exit templates).
- The proposal **removes the second factor** for sv:7 and substitutes L1–L3 there, so panel content inside sv:7 is counted
  once, by the unit model.
- The copying prior may enter only as **coupling**: the entry template's own components receive the copying prior's
  no-switch probability in the first copy's start distribution, and the exit transition conditions on the last copy's
  lineage.
- For a first experiment, the flanks are fixed and the copying prior's sv:7 term is set to zero.

**No parameter comes from F2.** w, c, ρ, κ, ε and η are estimated or tested only in the panel-sequence check below.

## 3. Independent validation and computational feasibility

**Donor-held-out, panel-sequence predictive check** (future; not run; no reads).
- **Folds:** hold out one **donor family** at a time. A held-out donor is removed together with its known pedigree
  relatives, otherwise transmitted alleles leak into fitting. All of that family's copies are removed from w, c, ρ, κ, ε, η
  and the cluster definitions.
- **Scored:** the held-out log predictive probability of (i) each haplotype's (n, a); (ii) each copy's alleles at SIMPLE
  and CONDITIONAL positions, plus segment choices, with escape counted as a separate failure category; (iii) diploid-total
  calibration against the observed donor distribution (2: 1, 3: 14, 4: 35, 5: 10, 6: 2).
- **Comparators, fitted on identical folds:**
  - C0: independent per-site frequencies, i.e. the "independent per-copy variation" approximation;
  - C1: nearest whole component (ρ = 0);
  - C2: lineage copying with κ = 0;
  - C3: lineage copying with κ fitted.
- **Decision rules, declared in that batch before execution:** C2 or C3 must beat C0 per copy, with a consistent sign across
  folds. κ is retained only if C3 beats C2 consistently.
- **Separation from development.** Donors already inspected in development are never scored, though they may stay in the
  fitting library: NA19650 (not in this panel), HG03248, HG00171, HG00096 (reference donor and anchor source), HG01573,
  NA19347 (allele_75), and the template donors of the development mosaic cells (three_switches_1, dense_switching_1,
  ts1_2), to be resolved from campaign metadata in that batch.

**State count and runtime for C4 with reads** (from measured panel structure and existing measured costs).

*Exact inference is infeasible.*
- Even restricted to whole panel components with no switches, a 2+2 diplotype has about (106 × 104)² / 2 ≈ 6 × 10⁷
  configurations.
- At the measured ≈ 0.86 ms per oriented candidate at 2× (sv:7 reachability cost model), that is about 15 h at 2×, and more
  at 10–30×.
- Mosaic paths make the space unbounded.

*Dependencies that prevent exact site-local factorisation:*
1. **Alternative placements:** a fragment from a repeated unit places on any copy in either haplotype at homologous
   positions. Its likelihood is a sum over copies, coupling all copies' local contents.
2. **Fragment span** of about 350–500 bp couples the sites in a window: a median of 2–3 SIMPLE sites per 400 bp, max 16–18.
3. **Lineage switching** couples windows along the unit, and **κ** couples copies within a haplotype.
4. **Junction fragments** (prefix–period–tail) couple adjacent copies. They are the only allocation evidence (M1-E, M6-C).

*Named approximation: windowed joint-copy forward pass.*
- Windows of about 800 bp (at least one fragment span) with overlap. The joint state in a window is the tuple of each copy's
  local pattern, restricted to panel-observed local patterns plus one residual escape per copy.
- Each copy's switch transitions factorise. Fragment emissions are computed once per distinct local sequence (≤ 34 per
  window per role) and combined per state as a mixture over copy placements.
- Joint states per window:

  | copies | median window | 90th percentile window | worst window |
  |---|---|---|---|
  | 4 | 6⁴ ≈ 1.3 × 10³ | 17⁴ ≈ 8.4 × 10⁴ | 34⁴ ≈ 1.3 × 10⁶ |
  | 5 | 7.8 × 10³ | 1.4 × 10⁶ | 4.5 × 10⁷ |
  | 6 | 4.7 × 10⁴ | 2.4 × 10⁷ | 1.5 × 10⁹ |

- About 41 windows cover the period coordinate, with roughly 160 fragments per window at 30× with 4 copies.
- Estimated cost: seconds to minutes for 4 copies in typical windows. Worst windows at 5–6 copies, and the segment-copy
  states of the long/short region (up to 36 signatures per copy), need beam pruning with an omitted-mass bound (I5-C style),
  which is a second named approximation.

*How the approximation error would be assessed:*
- against exact enumeration on reduced instances: M6-C-style scaffolds seeded with real local patterns, and one real C4
  window with the library restricted to ≤ 4 lineages, where exact enumeration is feasible;
- with escape posterior mass reported per window, which flags content the panel cannot spell;
- with the beam's omitted-mass bound reported per window.

*Engineering dependency:* the current engine scores whole diplotypes. Per-fragment, per-copy local emissions with confined
placement sets do not exist yet, and this is the main implementation cost.

## 4. Decision

- **What the panel supports.**
  - A structural layer with clear empirical frequencies: 1 / 2 / 3 copies in 17 / 95 / 16 haplotypes; arrangements
    dominated by LS and LL; the periods of a 3-copy haplotype usually of different classes.
  - Strong local linkage and low local diversity, which favour lineage copying over per-site independence.
  - Modest within-haplotype similarity: copies neither identical nor independent.
  - Simple correspondence for 74–90 % of unit positions. The long/short segment and about 30 short intervals need segment
    treatment.
  - Duplicate representation through trio relatives, which needs a pedigree table.
- **The proposed model** is L1–L3 above. Its assumptions:
  - the two haplotypes are independent given the panel;
  - the segment-copy approximation holds in non-simple intervals;
  - escape is reported, never spelled;
  - the two-rate residuals are adequate;
  - the flanks are fixed in a first experiment.
- **Practical for one bounded implementation experiment?** **Yes for a panel-only predictive experiment:** about 106–116
  components and about 270 sites per role, where a Li–Stephens forward pass is trivial. **Not yet for read-likelihood
  integration:** that needs a windowed per-copy emission decomposition that does not exist, and the worst windows at 5–6
  copies reach 10⁷–10⁹ states.
- **Smallest decisive next test: the donor-family-held-out panel-sequence predictive check** (C0–C3 above). It uses no
  reads and no caller.
  - **If lineage copying does not beat independent per-site frequencies on held-out copies, park this direction.**
  - If it does, the next test is an exactness check of the windowed approximation on one restricted real window, against
    exact enumeration, before any read-likelihood implementation.
- **Not claimed.** Nothing here addresses allocation: reads still inform allocation only through junction adjacency. Nothing
  here shows the model would recover C4 NA19650: its true periods and tails are absent from the panel (M1-D: 1–23 bp from
  the nearest components), which the residual/escape layer must absorb.

**Unresolved, recorded rather than expanded:**
- Unit-cost alignment fragments the long/short presence/absence. An affine-gap or segment-aware alignment would represent it
  as one event; the bundled minimap2 has no built binary, and building one was out of scope.
- Signature counts in non-simple intervals are upper bounds.
- There is no pedigree metadata in the panel.
- One reference per type: the sensitivity arm changed the SNV site count by 6 of 270.

## Qualifications appended at M6-F (2026-10-20; Sol's review; the record above is not rewritten)

M6-E's measurements are accepted, with four qualifications.

1. **Within-copy linkage supports testing a linked model, but does not show κ is needed.** Within-copy site linkage supports
   testing a linked sequence model. It does **not** establish that cross-copy coupling κ is necessary.
2. **The copy-count frequencies are descriptive.** They are not a fitted, pedigree-adjusted structural prior.
3. **Two pieces remain unspecified for production:** the escape distribution, and the coupling to the existing block-level
   copying prior.
4. **The exhaustive-enumeration runtime is not a lower bound.** Exhaustive enumeration is expensive, but its estimated
   runtime is not a lower bound for other inference methods.
