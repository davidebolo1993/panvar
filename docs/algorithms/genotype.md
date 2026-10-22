# Module `genotype` — algorithm

For commands, options, files and assessment values, use the [genotype manual](../modules/genotype.md).
Scientific references are in [references.md](../references.md). This page describes the implemented short-read
model, not a proposed long-read caller or a claim of exhaustive inference.

The caller proposes diplotypes cheaply from markers, then ranks evaluated sequences by paired-read likelihood
plus a copying prior. The search is bounded. Marker proposal scores are not added to the final objective.

## How it works

### 1. Build the blocks

The reference defines the locus orientation. Variant sites come from the graph's snarl decomposition or a
supplied bubbles CSV. A site is retained when some path spans at least `--min-variant-bp` bases inside it;
this threshold constructs the model, rather than filtering final predictions.

Variant blocks (`sv:*`) alternate with context blocks (`ctx:*`). Shared sequence is trimmed from variant
alleles into neighbouring context, and identical allele sequences collapse within a block. The adapter verifies
that block paths re-spell the original haplotypes. Partial paths, inconsistent spelling/representation, or
insufficient blocks/panel haplotypes can refuse the model rather than being silently approximated.

### 2. Define the active panel and prepare reusable data

Every graph path participates in the frame, divergence trimming and representation checks.
`--panel-curation` excludes named templates from the callable panel, not from this framing graph.
An allele survives only if an active template carries it; shared alleles retain active provenance only.
No subsequent search stage can select an allele absent from the active model.

The optional panel index persists this sample-independent model and marker panel, with versions, fingerprints
and input identities. It does not retain read evidence or sample search state. Loading verifies compatibility;
it does not change the inference model.

### 3. Count markers, fit depth and estimate the library

Marker candidates are canonical closed syncmers: the minimum s-mer lies at the first or last position of a
k-mer. A marker is variable in one block and has known per-allele multiplicity `n_m(a)` and fixed multiplicity
`f_m` elsewhere. Syncmers fixed throughout the model supply invariant anchors. Both read mates are counted.

Depth and background are fitted by a Poisson model over anchor multiplicities. With fitted rate `r` per
haploid copy and background `bg`, the proposal emission for block pair `(x,y)` is:

```text
S_b(x,y) = Σ_m [c_m ln μ_m(x,y) − μ_m(x,y) − c_m ln μ0_m + μ0_m]
μ_m(x,y) = bg + r [2 f_m + n_m(x) + n_m(y)]
μ0_m     = bg + 2 r f_m
```

Under the default automatic library, decisive concordant anchored pairs estimate insert mean and standard
deviation after robust validation. The marker-fit rate supplies fragment-start intensity, with read-window
and observed invariant-site mismatch correction:

```text
λ = r / [2 (L̄ − k + 1) (1 − ε_observed)^k]
```

The observed mismatch correction here does not replace the preset's substitution-error parameter in the
read-emission model. That parameter, the outlier weight, insert support width and gap parameters remain fixed.
Too few decisive pairs, contamination, multiple modes or an invalid fit cause automatic-library refusal.
A declared library bypasses parameter estimation, not the marker proposal construction.

### 4. Compute the copying prior

The default Li–Stephens model has `n` active donors and per-boundary redraw probability `ρ`:

- initial donor probability `1/n`;
- stay on a given donor with probability `1 − ρ + ρ/n`;
- change to each other donor with probability `ρ/n`;
- emit a block allele only when the donor carries it.

This applies at every block boundary, including context boundaries. The forward recursion sums compatible
donor histories. Sequence-level prior mass also sums alternative block labellings spelling identical DNA:

```text
ℓ(s) = log Σ_{p spells s} P_LS(p)
log P({s1,s2}) = ℓ(s1) + ℓ(s2) + 1[s1 ≠ s2] log 2
```

`ρ` is a model regularizer, not a biological recombination-rate estimate. The optional experimental
distance mode uses reference-core geometry and `ρ_b = 1 − exp(−d_b / copying_length_bp)`; it is not calibrated
and is not the default. Neither mode creates novel block alleles.

### 5. Construct the initial candidates

Every unordered pair of active panel haplotypes forms the baseline. Additional mosaics are constructed by
a grouped beam over blocks. The proposal combines marker emissions with copying-prior prefix mass:

```text
proposal(p1,p2) = Σ_b S_b + log Σ_i u_p1(i) + log Σ_i u_p2(i) + 1[p1 ≠ p2] log 2
```

`u_p(i)` is the compatible prefix mass ending on donor `i`. `--proposal-state-limit` bounds surviving partial
allele-prefix pairs at a layer; `--mosaic-candidates` bounds completed non-panel proposals admitted to the
initial scoring union. They limit different stages. Work/memory refusals are explicit, not silent approximations.

### 6. Score the read likelihood

For a pair `{a,b}`, the likelihood has an expected-count penalty (exposure) and a fragment-placement term:

```text
log L(a,b) = −(1−η) λ (E_a + E_b)
             + Σ_f ln[(1−η) λ (M_fa + M_fb) + η P_bg(f)]
P_bg(f)   = 0.25^(L1+L2)
```

`M_fh` sums the admitted placements of read pair `f` on haplotype `h`. Each placement contributes an orientation
prior `1/2`, discrete-normal insert probability and the two mate emissions. Verification uses the declared
affine-gap emission model with at most `--mismatch-band` edits per mate. Insert support is truncated to
`mean ± support_sigmas × sd`. `η` is the outlier weight and `λ` the per-copy fragment-start intensity.

Exposure averages valid fragment starts over the insert distribution and observed mate-length floors:

```text
E_h = Σ_F (n_F/N) Σ_i P_lib(i) (L_h − i + 1)⁺
```

An extra copy increases placement mass and exposure; the score balances both. Reads are paired molecules,
not two independent votes. The calculation is exact within its declared placement/emission model and caps,
not a claim that every biological read origin or variation is modeled.

The accelerated preset uses indexed/fine and blocked scoring, prepared library terms, shared placements,
bounded caches and incremental block-replacement scoring. Reused terms must match cold/reference scoring;
mandatory-store or work-cap failures remain explicit. These optimizations differ from proposal heuristics,
which can change which candidates are evaluated. The index saves locus preparation, not per-sample scoring.

### 7. Refine, branch and decide over the evaluated pool

Final ranking uses `log L + log P`. Identical unordered sequence pairs are grouped. Numerical bounds determine
which groups cannot be separated; this is the decision set, not a fixed biological quality threshold.

Refinement starts from the best-ranked group containing a panel pair. Each round generates a bounded local
allele-pair list for every block, substitutes both homologue assignments and scores new neighbours. It moves
to a certified improving neighbour across this combined neighbourhood, not simply to the first weak block.
Round, neighbour, work and tie-branch caps can prevent convergence.

Ordinary lists use block-internal syncmer emissions. When all pair emissions tie and more pairs exist than
the local limit, the preset's junction-aware rule (JP1.1) uses all canonical 31-mers in allele windows with
up to 30 incumbent bases on either side. Its Poisson proposal accounts for copies elsewhere in the locus;
both assignments are considered. Residual ties use content-based ordering, with scores compared to the tie
class head. The cutoff tie is reported. This chooses proposals only: it adds no term to the objective.

Improving ties can branch breadth-first within the declared branch bound. All evaluated union, refinement
and branch candidates are pooled for the v2 report. A completed pool can yield one sequence pair, multiple
pairs with matching block genotypes, or pairs with unresolved block genotypes. All are completed predictions
when the reporting decision completes; sequence ties are not automatically a missing prediction.

The representative is chosen deterministically within the inseparable set by sequence SHA-256 ordering.
It need not have the largest last floating-point bit. Search completion is conditional on the candidates
actually explored, not certification of a global optimum.

The early checkpoint retains a fully scored candidate if a later budget/cap stop prevents completion.
It is reported as `PROVISIONAL`, without a completed decision set or confidence. With refinement off/zero
rounds, the existing union alone supplies the same v2 reporting contract; no additional scoring is performed.

### 8. Annotate evidence and representation

The representation guard compares the reported block against its marker-best alternative when the marker
margin is at least `ln 1000`. A conflict is recorded when markers and exposure prefer that alternative but
fragment placement prefers the call, with the declared minimum margin. It is a component-disagreement check,
not proof that a particular genotype is wrong. It does not run on the unrefined route.

Automatic-library sensitivity reranks the evaluated pool and reruns refinement at 0.9λ and 1.1λ. A changed
unordered block genotype is library-sensitive; incomplete checks remain unavailable.

For the current default support rule, a fragment contender is an explored rival whose fragment-placement
term is not certifiably below the winner's. Exposure and the copying prior are excluded from that comparison.
A block is read-supported when no such contender has a different unordered genotype. Phase support is
computed once per locus, not separately for every boundary. Both depend on the explored pool and representative
labellings; absence of a useful rival can make support vacuous.

Projection-aware error-mass bounds summarize evaluated sequence groups with definite or possible differing
block genotypes. The older representative-based GQ remains a diagnostic. Neither is a calibrated probability
of error, and no automatic reliable-call threshold follows from it.

### 9. Audit dosage without changing the prediction

`--dosage-check` enables an independent count audit of reported copy features. A census identifies sensitive
blocks; large length differences alone are not a universal dosage classifier. Features group block k-mers
with shared copy patterns, compared to single-copy reference k-mers:

```text
ĉ_f = 2 (T_f / |f|) / (T_R / |R|)
```

Read-pair bootstrap resamples give quantiles `q05, q25, q50, q75, q95`. With `m` evaluated features,
the default classification uses:

```text
robust_sd  = (q75 − q25) / 1.349
z          = |reported_count − q50| / robust_sd
z_critical = Φ⁻¹(1 − 0.001 / (2m))
trigger    = z ≥ z_critical and reported_count is not admissible
suspect    = z ≥ 2.5
underpower = robust_sd > 1 / [z_critical + Φ⁻¹(0.90)]
```

Admissible counts are integers in `[ceil(q05), floor(q95)]` intersected with the feature's feasible range.
If that set is empty, the implementation retains the nearest feasible integer to the interval. At block
level a trigger takes precedence over suspicion, which takes precedence over underpower, then consistency.
A trigger is reported as `DISAGREES` whether or not a correction is found; a block classified underpowered
is `UNAVAILABLE`. Failed unit-library construction is also unavailable. These are declared rules, not
calibrated reliability. Count features require at least 500 k-mers under the default audit settings.
When the robust standard deviation is zero, `z` is zero for an exact median match and infinite otherwise.

Diagnostic alternatives can enumerate panel pairs or repeated-unit content/allocations under their own
bounds. Complete count content does not imply unique homologue allocation. Alternative sequences and legacy
file names containing `repaired` are audit diagnostics only. No command-line dosage mode applies a repair;
the primary prediction, its scores and haplotypes remain unchanged. Internal repair machinery is retained
for regression coverage, not offered as a supported prediction-changing mode.

### 10. Publish sequences and independent assessments

The sequence report separates one DNA pair from ambiguity in block labelling, sequence sets and phase.
Panel membership uses exact byte identity in the locus frame. Two panel-represented haplotypes use extraction
references instead of FASTA; if either is off-panel, the FASTA contains both representative sequences.

Only the four primary sequence report names remain beside the output prefix. Other reports are under its
diagnostics directory. The primary `assessment_v1` reports sequence resolution, genotype conflicts, dosage,
read support, phase, findings and missing checks. A fact from a compatibility block table is used only when
verified against the published alleles. Findings never disappear behind missing checks; locus coverage is
reported explicitly. Phase is locus-level. The old headline filters remain superseded diagnostics.

No new scoring is performed by this assessment layer. These annotations do not alter selection and do not
claim calibrated confidence. See the [manual's assessment tables](../modules/genotype.md#assessments-independent-checks-not-one-passfail)
for every outcome and its conditions.

## Worked trace

Consider two variant blocks separated by context. The panel carries `(a,c)` and `(b,d)`; the sample's reads
suggest `(a,c)` and `(a,d)`. This is an illustrative mechanism, not a real-data performance claim.

1. Build a block chain in the reference orientation and verify that its alleles re-spell every full path.
2. Apply curation to callable templates without changing that frame; optionally reuse compatible preparation.
3. Count markers, fit depth and estimate/validate the automatic library, or use declared parameters.
4. Sum copying histories; the mosaic `(a,d)` pays for a donor change under the default prior.
5. Evaluate panel diplotypes and admit a bounded beam shortlist that may include `{ac,ad}`.
6. Score paired-read placement and exposure for each candidate, adding its sequence-level copying prior once.
7. Search bounded single-block neighbours; if supported, move from `{ac,bd}` to `{ac,ad}`. Pool all evaluated candidates.
8. Assess representation, library sensitivity and fragment contenders. These are distinct from sequence uniqueness.
9. If requested, audit implied copy features without substituting any dosage alternative into the prediction.
10. Write the completed prediction with its independent assessments. Here only the second variant block is
    heterozygous, so resolved phase is `NOT_APPLICABLE`; sequence uniqueness alone does not certify accuracy.
