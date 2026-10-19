# M4-A: distance-aware copying-prior contract (proposal; not implemented in any caller)

**Status.** Design. The contract is exercised only by the prototype in `tests/m4a_distance_prior.cpp`, together with the
read-only probe `prior-geometry`. Production keeps the constant prior, and nothing here changes a default.

**Builds on M3-A, which is not repeated.** M3-A established four things:
- the current donor-redraw semantics;
- the dependence of constant r on neutral subdivision;
- a distance-matched fixed-pool comparison, using panel-mean spans;
- the sequence-level vs labelled-path distinction.

None of that is rerun here, and none of it is fresh-search validation.

## 1. The transition

**Form.** Between consecutive blocks b−1 and b of the canonical chain:

  **T_b = (1 − r_b) I + r_b U**,  r_b = 1 − exp(−d_b / λ)

- **d_b ≥ 0** is the declared distance between the two copying states (section 2).
- **λ > 0** is the copying length in bp. Equivalently, κ = 1/λ is the redraw hazard per bp (section 3).
- **U** is the declared redraw measure.
  - Version 1: uniform over the n curated template rows. Identical rows count separately; the multiplicity is explicit
    and visible (G8).
  - U includes the current donor. A redraw therefore changes donor with probability r_b(n−1)/n, stays put with
    probability 1 − r_b + r_b/n, and changes the allele only when the redrawn donor carries a different allele.

**Unchanged.** The start distribution (uniform over rows), emissions (a row emits exactly the allele it carries) and the
sequence-level aggregation, which sums every allele path that spells the bytes (G7). The existing recursions
(`log_haploid_mosaic_prior`, `log_haploid_sequence_prior`, `log_unordered_sequence_pair_prior`) already take one r per
boundary, so only the construction of that vector changes.

**Exact composition (the reason for the exponential form).** U is idempotent and 1ᵀ-stochastic (U·U = U, I·U = U), so

  T(d₁) T(d₂) = e^{−(d₁+d₂)/λ} I + (1 − e^{−(d₁+d₂)/λ}) U = T(d₁ + d₂).

- Monomorphic (all-emitting) blocks therefore compose exactly whenever the distances add.
- No other per-boundary rule composes. M3-A R1, and G3 here, show that constant r changes the prior by 0.69 nats on
  one neutral insertion.

**Why a donor redraw is not a biological crossover probability.**
1. **It is a copying-model event.** In Li & Stephens' product of approximate conditionals, a "switch" is a change of
   the template the new haplotype copies from. It approximates, through a genealogy, the effect of historical
   recombination accumulated over many generations, not one meiosis. In their model the switch rate is
   ρ·d/n with ρ = 4N_e·c, so it scales with effective population size and **inversely with the number of templates**:
   larger panels contain closer relatives and longer shared segments.
2. **The definition includes the current donor.** A redraw to the same donor, or to one carrying identical alleles, is
   invisible. The observable change rate is not r.
3. **Redraws absorb more than crossovers.** Gene conversion, recurrent mutation, panel incompleteness (an unsampled
   lineage looks like a mosaic of sampled ones), structural rearrangement at RCCX/CYP2D6, and the switches of simulated
   mosaics all appear as redraws.
4. **Physical scale.** Per-generation crossover probability is about 10⁻⁸ per bp (roughly 1 cM/Mb). It is non-uniform,
   with hotspots, and sex- and population-specific. A copying model fitted to a panel is a different, effective
   quantity.

So r_b, λ and κ are **copying-model parameters**. Calibrating them biologically is not established here (section 3).

## 2. Geometry: where a copying state sits, and d_b

**Recommendation: reference-coordinate core midpoints, trimmed against the frame.**

- **Axis.** Cumulative position along the reference path in canonical chain order: invariant segment 0, then the
  reference allele of block 0, then segment 1, and so on. The probe checks that this tiles the chain and equals the
  reference path's bytes (all 7 buildable models: pass).
- **Copying-state position.** The midpoint of block b's **polymorphic core**: the reference allele's interval with the
  maximal prefix and suffix shared by every allele of the frame removed. A monomorphic block has an empty core and
  contributes only distance.
- **d_b** = core_mid(b) − core_mid(b−1). It is always ≥ 0 (G1); zero only for two empty-core states at one point.
- **What the position means.** A block's allele is atomic: one donor copies the whole block, and redraws inside a block
  cannot be represented. The copying state is the donor of the block's informative content. Its centre is the core
  midpoint, which does not depend on how much shared sequence a representation happens to put inside the block.

**Why not the alternatives** (measured on C4, CYP2D6 and ACOT; RESULT.md):

| geometry | problem |
|---|---|
| invariant gap only | Every inter-block invariant segment is empty (C4 10/10 boundaries, CYP2D6 18/18, ACOT 18/18). The alternating chain puts all sequence in blocks, so this geometry is degenerate. |
| reference start-to-start (and M3-A's panel-mean span, which is also start-to-start) | Assigns a whole block's length to the boundary after it. After a zero-reference-span insertion block (C4 sv:5) the following 48 kb context block gets distance 0 bp (start) or 110 bp (panel-mean), so it can essentially never redraw. The preceding boundary carries all the length. The asymmetry is an artefact of the anchor. |
| untrimmed block midpoints | Moves when shared bytes are written into a block rather than outside it (G4 control). Numerically within 0.5% of core midpoints on the real models, but not representation-invariant. |
| M3-A panel-mean span | Depends on panel composition. Curating 2 templates out of the same graph moves it by up to 227 bp (ts1), while reference geometries move by 0. It also depends on where shared bytes are written (G4) and on donor removal (G9). |

**Required treatments:**

| case | treatment |
|---|---|
| Context blocks and invariant intervals | Treated identically, as reference sequence. A context block is polymorphic (sub-threshold variation), so it gets a core. An invariant segment, or a monomorphic block, has no core and adds length only (G3, G4). |
| Insertions (zero reference span) | The core is an empty point at the insertion site. Distances to neighbours come from the neighbours' half-cores. That is positive unless the neighbour is also an empty point at the same position. |
| Deletions relative to the reference | The reference allele carries the long sequence, so the block's reference span is the long allele. |
| Zero distance | r = 0, T = I. Only panel-observed combinations survive (G2). This happens only for two empty-core states at one point, never for merely adjacent blocks (G2: adjacent single-base cores give d = 1). No floor is added, because a floor would break exact composition. |
| Reverse-oriented paths | Geometry is defined on the reference axis in canonical chain order. A chain recorded as reverse relative to the external frame mirrors the positions, and distances are unchanged (G5). Template source orientation plays no part: positions come from the reference path, not from donors. |
| Equivalent graph representations | The adapter already divergence-trims variant blocks and refuses representation-dependent blocks. Core trimming removes any remaining shared bytes, so node splits and where shared sequence is written do not move positions (G4). |
| Unavailable or ambiguous reference coordinates | If the reference path does not tile the chain, is Uncovered at a block, or visits the chain inconsistently, the distance-aware mode **refuses** with a typed reason. It never falls back to constant r or to panel-mean. |
| Leaving out panel donors | **Same-graph curation must not move geometry.** The frame is fixed before curation, and cores are trimmed against the **frame's** alleles, including curation-inactive ones. That requires the adapter to emit the core (section 5). Measured here: the reference geometries are unchanged by curation (ts1), but cores trimmed against **callable** alleles would move (G9 control). **A rebuilt leave-one-out graph is a different model.** Its divergence trim can move a block boundary: C4 LOO NA19650 moves sv:3 from 166 to 120 bp, which changes two core-midpoint distances by 23 bp each, with the total unchanged. That is a frame change of the model itself, recorded by the geometry digest, not hidden. |
| Independence | d depends only on the graph frame and the reference path: never on reads, truth or the selected candidate. |

**Unresolved geometry choices** (declared, not decided here):
- **G-U1. Structural blocks carry the reference allele's length.** C4 sv:7 counts the reference's 59 kb RCCX allele
  (panel mean 61 kb). A different reference haplotype (more or fewer modules) would change d around sv:7. This matches
  standard practice, since imputation HMMs use reference-genome positions, but the reference dependence at large
  structural blocks is real. Alternatives would be a declared structural-block length or a donor-averaged length, and
  the latter reintroduces panel dependence.
- **G-U2. Frame cores need adapter support.** The prototype probe trims against callable alleles.

## 3. Scale: three separate things

**The scale is one length λ (equivalently κ = 1/λ per bp).** It is never recomputed after subdivision. Because
distances add, the total hazard Σ d_b/λ is unchanged by neutral subdivision (G10).

| role | definition | status |
|---|---|---|
| **(a) Hazard-matched diagnostic control** | λ_HM = Σ_b d_b / ((B−1)·(−ln 0.99)) per model: the same total hazard as the legacy constant model on that representation, redistributed by distance | **Diagnostic only.** It depends on block count by construction (G10 control: one neutral insertion changes it by 497 bp on the fixture). It separates the effect of where the hazard goes from how much there is. Measured: C4 1.87 Mb, CYP2D6 1.13 Mb, ACOT 1.22 Mb. These differ because legacy constant r implies a different hazard per bp at each locus. Freeze λ_HM per model fingerprint for a comparison, and never re-derive it after a representation change. |
| **(b) Physically interpretable copying model** | one declared κ (redraws per bp of reference), the same across loci and representations, recorded with n | **Value not chosen here.** Choosing it is a parameter selection and needs a **later, independent development study** on panel data, not on the evaluation cells and never on the known 2× example. The natural truth-free design is Li & Stephens' own: maximise the leave-one-template-out copying likelihood of panel haplotypes over κ, on development loci. Declaring κ (not ρ) keeps panel size changes from silently rescaling the hazard; under the LS reading κ = ρ/n, and n is recorded. |
| **(c) Biological calibration** | relating κ to recombination maps, N_e and panel size | **Unestablished.** No claim is made. For orientation only (not a calibration), nominal LS values, ρ = 4N_e·c with N_e ≈ 10⁴ and c ≈ 10⁻⁸ per bp, at n ≈ 130 give a copying length of order 10⁵ bp. The MHC class III region is not average, and panel ancestry and structural variation make any such number an effective quantity. |

## 4. Mathematical consistency (oracle gates; RESULT.md)

**Method.** `tests/m4a_distance_prior.cpp` checks the prototype against the brute-force history oracle from
`tests/genotype_mosaic_core.cpp`. λ and U are frozen in every comparison.

**Result.** All contract gates pass: G1 validity, G2 zero distance, G3 neutral insertion at every split point, G4
shared-affix placement, G5 orientation, G6 polymorphic split classified non-neutral (4 new admissible recombinants), G7
equivalent spellings, G8 donor multiplicity and order, G9 donor removal with frame cores, and G10 scale. Controls confirm
the dependences the contract removes: constant r, untrimmed midpoints, panel-mean span, callable-allele cores and the
hazard-matched length. M3-A's test passes unchanged.

## 5. Integration surface (summary; detail in INTEGRATION.md)

**One shared construction.** `build_copying_transitions(model, geometry, spec) -> {r vector, provenance digest}`.
- The geometry is a new adapter output: per-block reference core interval against the frame, plus the reference path
  name and its digest.
- Every stage that today builds `std::vector<double>(B−1, options.switch_probability)` receives the one vector
  instead: proposal beam, union ranking, path priors, refinement and tie seeds, joint search, representation guard,
  dosage audit and repair, and the probes.

**Provenance.** Mode, λ (or the control's λ_HM with its model fingerprint), the U declaration, the geometry digest and
the vector digest go into settings and status. They replace the bare `switch=` hex in cache keys.

**Safety.** Mixing constant and distance-aware vectors across stages is detected by digest equality and refused. A
prepared-locus v1 index carries no geometry, so the distance mode refuses it rather than recomputing from callable
alleles.

## References (primary sources; standard literature, cited from knowledge and not re-fetched in this batch)

- Li N, Stephens M (2003). Modeling linkage disequilibrium and identifying recombination hotspots using single-nucleotide
  polymorphism data. *Genetics* 165:2213–2233. The copying model: switches at rate ρd/n, with a uniform choice of new
  template.
- Scheet P, Stephens M (2006). A fast and flexible statistical model for large-scale population genotype data. *Am J Hum
  Genet* 78:629–644. Distance-dependent jump probabilities between hidden states.
- Marchini J, Howie B, Myers S, McVean G, Donnelly P (2007). A new multipoint method for genome-wide association studies
  by imputation of genotypes. *Nat Genet* 39:906–913. LS-type transitions from genetic distance, N_e and panel size.
- Howie BN, Donnelly P, Marchini J (2009). A flexible and accurate genotype imputation method for the next generation of
  genome-wide association studies. *PLoS Genet* 5:e1000529.
- Lunter G (2019). Haplotype matching in large cohorts using the Li and Stephens model. *Bioinformatics* 35:798–806.
- Ebler J, et al. (2022). Pangenome-based genome inference allows efficient and accurate genotyping across a wide spectrum
  of variant classes. *Nat Genet* 54:518–525. PanGenie: an LS-style HMM over graph bubbles with distance-dependent
  transitions.
- Myers S, Bottolo L, Freeman C, McVean G, Donnelly P (2005). A fine-scale map of recombination rates and hotspots across
  the human genome. *Science* 310:321–324.
- Kong A, et al. (2010). Fine-scale recombination rate differences between sexes, populations and individuals. *Nature*
  467:1099–1103.
- Halldorsson BV, et al. (2019). Characterizing mutagenic effects of recombination through a sequence-level genetic map.
  *Science* 363:eaau1043. Crossovers and non-crossover gene conversion.
