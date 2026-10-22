# Module `genotype` - algorithm

Mechanism for the `genotype` module. For usage, outputs and status values see [modules/genotype.md](../modules/genotype.md); for every tunable constant see `panvar genotype --help-advanced`; references in [references.md](../references.md).

`genotype` calls a diplotype, two whole-locus haplotypes, as a pair of paths through the locus graph. The locus is cut into blocks: variant sites alternating with invariant context. A haplotype takes one allele per block and may switch between panel haplotypes at block boundaries.

Marker counts and a Li–Stephens model propose candidates cheaply. Each evaluated candidate is then scored by a paired-read likelihood, computed without approximation under the declared read model up to floating-point bounds, plus the Li–Stephens prior. Marker evidence only selects candidates and is never added to the score, so no read is counted twice. The candidate set, the repair search and tie branching are all bounded: the call is the best of what was evaluated, not a proof over every possible diplotype.

## How it works

### 1. Build the blocks

The graph's paths are validated, and variant sites are taken from its top-level snarls (or from `--bubbles`). A site is kept when some path spans at least `--min-variant-bp` bases inside it. Nested sites reduce to their innermost bubbles.

The locus becomes a chain `ctx:0, sv:1, ctx:2, …` of context and variant blocks. Each variant block is trimmed to where its alleles diverge, and the shared bases move into the neighbouring context. Identical allele sequences in a block collapse to one canonical allele.

The model is refused (`MODEL_INCOMPLETE`) in these cases:
- a path covers only part of the locus;
- re-spelling a path from its block alleles does not reproduce it;
- a block depends on how the graph happens to represent it;
- fewer than two blocks or two panel haplotypes remain.

### 2. Define the active panel

The active panel is the set of paths not named by `--panel-curation`. Block structure (frames, the divergence trim, representation checks) uses every path, so curation never changes the blocks. A block's source allele stays in the model only if a retained path carries it, so an allele carried only by curated paths is reachable by no later stage. A sequence-identical allele carried by a retained path keeps only its retained sources. A block left with no active allele refuses the run.

### 3. Count syncmer markers and fit the depth

A k-mer is a marker candidate when its canonical form is a closed syncmer, meaning its smallest s-mer sits at the first or last position. Within a block, a syncmer that every allele carries equally often is fixed; otherwise it is variable. A syncmer variable in exactly one block is that block's marker, with per-allele copy number `n_m(a)` and fixed multiplicity `f_m` elsewhere. A syncmer fixed everywhere is an anchor, with haploid multiplicity `m_k`. Every closed syncmer of every read mate is then counted.

Anchor counts fit the depth: `c_k ~ Poisson(bg + r·2·m_k)`, where `r` is the expected count per haploid copy and `bg` a background rate, by maximum likelihood. A failed fit stops the run (`CANDIDATE_INCOMPLETE`). The marker evidence for allele pair `(x, y)` at block `b` is a Poisson log-ratio against "no variable copies":

```text
S_b(x,y) = Σ_m [ c_m·ln μ_m(x,y) − μ_m(x,y) ] − [ c_m·ln μ0_m − μ0_m ]
μ_m(x,y) = bg + r·(2·f_m + n_m(x) + n_m(y))        μ0_m = bg + 2·r·f_m
```

**Automatic library.** Under `--library auto` the library is estimated at this point:
- **Inserts:** read pairs that anchor 31-mers place identically on every panel haplotype give the insert distribution (median and robust spread, then the mean and sd of the inliers).
- **Intensity:** the depth fit gives the fragment-start intensity, `λ = r / (2·(L̄ − k + 1)·(1 − ε)^k)`. Here `L̄` is the mean mate length, `k` the marker k-mer size and `ε` the observed mismatch rate.
- **Refusal:** an estimate that fails validation refuses the run (`LIBRARY_REFUSED`).

### 4. Compute the Li–Stephens prior

A haplotype is a path of block alleles copied from the `n` active panel haplotypes:
- **Start:** `1/n` on each haplotype carrying the first allele.
- **Transitions between blocks:** stay with probability `1 − ρ + ρ/n`, move to a given other haplotype with probability `ρ/n`, where `ρ` is `--switch-probability`. This applies at every boundary, context boundaries included.
- **Emission:** 1 when the copied haplotype carries the path's allele, and 0 otherwise.

The forward recursion sums the probability mass `P_LS(p)` of an allele path `p` over its copying histories. Ranking works with log probabilities. The log sequence prior `ℓ(s)` sums that mass over every allele path spelling the same bytes `s`, and the log prior of an unordered pair adds `log 2` when the two sequences differ, for the two orders in which they can be drawn:

```text
ℓ(s)            = log Σ_{p spells s} P_LS(p)
log P({s1, s2}) = ℓ(s1) + ℓ(s2) + 1[s1 ≠ s2]·log 2
```

A path through an allele that no active haplotype carries has probability mass 0 (log prior `−∞`), which is why such an allele can never be proposed.

### 5. Construct candidates

Candidates come from two sources:
- **Panel baseline:** every unordered pair of active panel haplotypes, always evaluated.
- **Mosaics:** a beam over blocks keeps the best `--proposal-state-limit` allele-prefix pairs. A pair is scored by `Σ_b S_b + log Σ_i u_p1(i) + log Σ_i u_p2(i) + 1[p1 ≠ p2]·log 2`. Here `u_p(i)` is the Li–Stephens probability mass of the paths spelling prefix `p` that end on haplotype `i`.

The first `--mosaic-candidates` complete diplotypes that are not panel pairs join the baseline to form the evaluated union. A beam that exceeds its work or memory budget stops the run (`CANDIDATE_INCOMPLETE`).

### 6. Score the paired-read likelihood

Each evaluated candidate `{a, b}` is scored by

```text
log L(a,b) = −(1−η)·λ·(E_a + E_b) + Σ_f ln[ (1−η)·λ·(M_fa + M_fb) + η·P_bg(f) ]
```

- **`λ`:** the fragment-start intensity per base and haplotype copy.
- **`η`:** the outlier weight, the fraction of fragments from outside the model.
- **`P_bg(f) = 0.25^(L1+L2)`:** the background probability of fragment `f`'s two mates.
- **`M_fh`:** the placement mass of fragment `f` on haplotype `h`, summed over placements. Each placement contributes `½ · P_lib(insert) · e(mate1) · e(mate2)`:
  - `½` is the orientation prior;
  - `P_lib` is a discrete normal insert distribution on `μ ± support_sigmas·σ`;
  - `e(·)` is a mate's emission under a banded affine-gap model. A match emits `1−ε`, a mismatch `ε/3`, and the gap model is fixed. At most `--mismatch-band` edits (substitutions plus gap bases) are allowed per mate.
- **`E_h`:** the exposure, the expected number of fragment starts on `h`. Here `L_h` is the haplotype length, `F` a fragment's length floor (its longer mate) and `N` the number of fragments:

  ```text
  E_h = Σ_F (n_F/N) · Σ_i P_lib(i)·(L_h − i + 1)⁺
  ```

  The exposure term is how depth enters the score: a haplotype that predicts sequence the reads never cover pays `(1−η)·λ` per expected fragment start. A homozygote counts two copies.

Placements are seeded by pigeonhole pieces of each mate, then verified by dynamic programming.

### 7. Rank, repair and branch at ties

**Ranking.** The score is `log L + log P({s1, s2})`. Candidates spelling the same unordered pair of sequences are ranked once. Two scores tie when their difference is within `δ_c = max((F + 16)·ε_double·(Σ|terms_c| + |prior_c|), 0.5e-9)`. The best candidate and every rival within the bounds form the top set.

**Repair.** The search starts from the best-ranked group containing a panel pair (the anchor). Each round does the following:
- for every block, takes the `--repair-local-pairs` best allele pairs by marker evidence;
- substitutes each into the incumbent, in both assignments;
- scores every new neighbour;
- moves only on a certified improvement.

It stops when nothing improves, or at `--repair-rounds`. A neighbourhood larger than `--repair-max-neighbours` ends in a no-call before any work.

**Tie branching.** An improving round whose top set excludes the incumbent branches on every tied neighbour, breadth-first, up to `--repair-tie-max-branches`. Beyond that, or if a branch fails to converge, the run no-calls (`tie_branching_incomplete`).

**Decision.** Every candidate evaluated by the union, the repair and the branches is ranked together:

| top set | call |
|---|---|
| one group | strict call (`COMPLETE`) |
| several groups, the same unordered genotype at every block | genotype-equivalent call |
| groups that disagree at a block | no-call (`distinct_tie`) |

**Support and quality.** A rival whose fragment term is not certifiably below the winner's is a fragment contender. `phase_support` is computed once per locus: `UNRESOLVED` when the decision is not unique, read-supported with no contender, and `PRIOR_ONLY` otherwise. A block's `genotype_support` is `PRIOR_ONLY` when a contender differs there. The block's conditional posterior and GQ are:

```text
log_err_b = logsumexp{ s_c : genotype_c(b) ≠ called(b) } − logsumexp{ all evaluated s_c }
posterior = 1 − exp(log_err_b)        GQ = −10·log10(exp(log_err_b))
```

They are conditional on the evaluated candidates and are not calibrated against truth.

### 8. Annotate representation and library sensitivity

These annotations never change the caller's decision.
- **Representation check.** At each block, the marker-best pair `T` is substituted into the call when its marker margin `m` is at least `ln 1000`, then scored. The block is `LOWQUAL_REPRESENTATION_CONFLICT` when `min(m, F_call − F_comp, X_comp − X_call) ≥ ln 1000`. Here `F` is the fragment term and `X` the exposure term. The conflict is that the markers and the exposure prefer the competitor while the fragment term prefers the call. The threshold is declared and conservative, not a calibrated Bayes factor.
- **Library sensitivity** (automatic library only). The evaluated pool is reranked at 0.9λ and 1.1λ, and the search is rerun at each. A block whose unordered genotype changes is library-sensitive; a view that cannot complete leaves the audit NOT_EVALUATED.
- **Block quality.** `quality_filter` is PASS only for a representation PASS with read support and no library sensitivity.

### 9. Audit and repair dosage

With `--dosage-mode audit` or `repair`, a called locus is tested for copy-number contradictions.

**Census.** With molecule span `S = ⌈μ + 4σ⌉`, a block is dosage-sensitive when two of its active alleles differ in length by at least `S` while the shorter adds at most 50 bp of nodes absent from the longer.

**Estimate.** For each tested feature `f` (a group of at least 500 block k-mers sharing one copy pattern), read hits against a single-copy reference set `R` give a copy-number estimate:

```text
ĉ_f = 2·(T_f/|f|) / (T_R/|R|)
```

Bootstrap resamples of the read pairs give quantiles `q05…q95` and a robust sd `(q75 − q25)/1.349`. With `m` tested features at the locus, the tests are:

```text
z = |call − q50| / robust_sd          z_c = Φ⁻¹(1 − 0.001/(2m))
trigger ⇔ z ≥ z_c and the call is outside the admissible copy numbers [⌈q05⌉, ⌊q95⌋]
suspect ⇔ z ≥ 2.5        underpowered ⇔ robust_sd > 1/(z_c + Φ⁻¹(0.90))
```

**Repair.** A triggered block searches dosage-admissible alternatives:
- **Pair route:** panel allele pairs.
- **Repeated-unit route:** alleles built from unit multisets.

Each alternative is spliced into the call and scored by the likelihood of stage 6. The likelihood alone chooses the content: the top set is every genotype within 1e-9 nats of the best.
- **Pair route:** more than one genotype in the top set leaves the block unrepaired. A tie between the two orders of one genotype is not a genotype tie.
- **Repeated-unit route:** top-set genotypes with different unit contents leave the block unrepaired. Otherwise every arrangement of that content is ranked by its panel-structure prior. The arrangements at the best prior within 1e-9 nats of the best likelihood among them form the allocation-equivalence set. When the set holds more than one distinct unordered genotype, the reported arrangement (the smallest canonical SHA-256 pair) only serialises it: the arrangement carries `serialization_representative=true` and is not a resolved diplotype.

Three coverage fields are recorded separately:
- `dosage_content_scope`: whether the corrected copy-number content is established.
- `diplotype_scope`: whether the candidate enumeration was complete.
- `allocation_equivalence_scope`: whether the allocation-equivalence set was enumerated completely.

Coverage is not uniqueness: a set enumerated completely can still have several members. The outcome is classified as follows:
- **Complete diplotype repair** (`COMPLETE_DIPLOTYPE_REPAIR`): the content is complete, the enumeration is complete, and on the repeated-unit route the allocation-equivalence set has a single member, so the repaired diplotype is uniquely resolved within the model and search. It is the only class whose arrangement has the role `RESOLVED`.
- **Conditional dosage content** (`CONDITIONAL_DOSAGE_CONTENT`): the content is complete but the diplotype is not resolved, either because the enumeration was conditional or because the allocation-equivalence set has several members (`dosage_reason=unresolved_allocation_equivalence_set`; the coverage scopes may all be `COMPLETE`). The arrangement is written only to the dosage diagnostics, labelled `REPRESENTATIVE_ONLY_NOT_RESOLVED`.
- **Triggered but unrepaired:** no admissible correction was found, the top set held a tie that the route cannot resolve, or the content itself is conditional.

In repair mode, only a complete, uniquely resolved diplotype repair can replace the block, under the eligibility rules in [modules/genotype.md](../modules/genotype.md#dosage-modes). Conditional content never does. As a second guard, integration refuses any result marked `serialization_representative=true`, whatever its class, with the reason `unresolved_allocation_equivalence_set`. After an applied repair, the stage-8 annotations are recomputed on the reported genotype. Each must first reproduce the caller's own result, else it is `NOT_EVALUATED_AFTER_DOSAGE_REPAIR`.

### 10. Derive the outputs and headline filters

**Outputs.** The block tables report the decided pair per block. A strict call also gets its two haplotypes, spelled as concatenated block alleles. The locus summary records:
- the haplotypes' SHA-256 digests;
- the separation from the best other genotype;
- for each haplotype, the fewest switches between panel haplotypes needed to spell it.

Sequences produced by the dosage stage are compared with the active panel by global unit-cost edit distance.

**Headline filters.** Three filters are derived from the component fields: GENOTYPE_FILTER per block, DOSAGE_FILTER per block, and PHASE_FILTER per adjacent pair of heterozygous blocks. Each unit takes the most severe applicable status, ranked `NO_CALL > NOT_EVALUATED > LOW_QUALITY > PASS > NOT_APPLICABLE`, and keeps every reason. The locus takes the worst applicable unit. A no-call is NO_CALL everywhere, and a run with no blocks gets a single `locus` unit. The per-status meanings are in [modules/genotype.md](../modules/genotype.md#reading-the-call).

## Worked trace

The locus has five blocks, with variant blocks `sv:1` (alleles `a`, `b`) and `sv:3` (alleles `c`, `d`). The four panel haplotypes are `(a,c)`, `(b,d)`, `(a,c)` and `(b,d)`; the sample carries `(a,c)` and `(a,d)`. Numbers are illustrative.

1. **Blocks and panel.** Five blocks. The two `(a,c)` paths spell the same bytes.
2. **Markers.** An `a`-only marker seen 20 times at `r = 11.28`, `bg = 0.01` gives `S(a,a) = 131.9`, `S(a,b) = 129.3` and `S(b,b) = 0`.
3. **Prior.** `ℓ(ac) = −0.70` and `ℓ(ad) = −5.30`. The mosaic `(a,d)` needs one switch, so `{ac,ad}` trails `{ac,bd}` by 4.6 nats in the prior.
4. **Candidates.** The panel pairs, plus the mosaic `{ac,ad}` from the beam.
5. **Likelihood.** If `b` is 1,700 bp longer than `a`, `{ac,bd}` predicts about 85 nats more unseen fragment starts at `λ = 0.05`, which outweighs the prior deficit.
6. **Decision.** The repair from anchor `{ac,bd}` reaches `{ac,ad}`, and the top set is one group: a strict call. Phase across `sv:1` is not informative, because `sv:1` is homozygous.

```text
sv:1  a/a  HOM  GENOTYPE PASS
sv:3  c/d  HET  GENOTYPE PASS
call_summary=GENOTYPE_PASS;DOSAGE_NOT_EVALUATED;PHASE_NOT_APPLICABLE
```
