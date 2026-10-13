# Module `genotype` - algorithm

Mechanism for the `genotype` module. For usage/flags see [modules/genotype.md](../modules/genotype.md); references in [references.md](../references.md).

`genotype` calls a diplotype, two whole-locus haplotypes, as a pair of paths through the locus graph. The locus is cut into blocks: variant bubbles alternating with invariant context. A haplotype is one allele per block, and it may switch between panel haplotypes at block boundaries. The method has two halves. A fast, approximate proposal uses syncmer marker counts and a Li–Stephens model to find promising candidates. An exact ranking then scores every candidate with a paired-read likelihood plus the Li–Stephens prior. The marker evidence only chooses which candidates to score; it is never added to the final score, so no read evidence is counted twice. An exact repair search around the best panel pair, and a decision over every candidate scored, protect the result from proposal misses.

## How it works

The stages run in this order. Under the default settings, stages 11–13 are skipped (`--dosage-mode off`).

### 1. Read and validate the inputs

The graph is read and its paths validated: unique names, every node present, an oriented link for every step. Variant sites come from `--bubbles` or, by default, from the graph's top-level snarls (cactus decomposition). A site is kept when at least one path spends `--min-variant-bp` (50) bp or more inside it, or when it carries an inversion signal. Nested sites are reduced to their innermost bubbles.

The locus becomes a chain of blocks, `ctx:0, sv:1, ctx:2, …`, alternating invariant context with variant blocks. B blocks give B−1 boundaries. Each variant block is trimmed to the span where its alleles diverge, and the shared bases move into the neighbouring context. Every graph path is a panel template, unless `--panel-curation` names it.

The model fails closed, with status `MODEL_INCOMPLETE`, if any of the following holds:
- a template covers only part of the locus;
- re-spelling a template from its block alleles does not reproduce its path's bytes;
- a block's alleles depend on how the graph happens to represent them;
- fewer than two blocks or two templates remain.

Identical allele sequences within a block are collapsed to one canonical allele.

A declared library (`--library FILE`) is validated:
- `μ > 0` and `σ > 0`;
- `0 < substitution_error < 0.75`;
- `0 ≤ outlier_weight < 1`;
- `λ > 0`.

The gap model is fixed: `gap_open` 1e-4, `gap_extend` 0.5, at most 3 bp net displacement per mate.

### 2. Select syncmer markers and count them

A k-mer (`--kmer-size` 31) is a marker candidate if it is a closed syncmer. Its canonical form has its smallest s-mer (`--syncmer-s` 15) at the first or the last position. About 2/(k−s+1) ≈ 12% of k-mers qualify. Syncmers are collected within each allele, so none spans a block boundary.

Within each block, a syncmer is classified as follows:
- **fixed** if every allele carries it equally often;
- **variable** otherwise. A syncmer variable in exactly one block is that block's marker, with per-allele copy number `n_m(a)` and fixed multiplicity `f_m` from the other blocks.
- **dropped** if it is variable in two or more blocks.
- **anchor** if it is fixed everywhere, with haploid multiplicity `m_k`.

Every closed syncmer of every read mate is then counted against this universe.

### 3. Fit marker depth and estimate the library

The anchors give the sequencing depth. Each anchor count `c_k` is modelled as Poisson with mean `bg + r·2·m_k`:
- `r` is the expected count per haploid copy;
- `bg` is a background rate, floored at 1e-6.

`(bg, r)` is fitted by maximum likelihood with a damped Newton method (at most 200 iterations), which needs at least two distinct anchor multiplicities. If the fit fails, the run stops with `CANDIDATE_INCOMPLETE`.

The marker emission of an allele pair `(x, y)` at block `b` is a Poisson log-ratio against "no variable copies":

```text
S_b(x,y) = Σ_m [ c_m·ln μ_m(x,y) − μ_m(x,y) ] − [ c_m·ln μ0_m − μ0_m ]
μ_m(x,y) = bg + r·(2·f_m + n_m(x) + n_m(y))        μ0_m = bg + 2·r·f_m
```

For example, take `bg = 0.01`, `r = 11.28`, `f = 0`, and one marker carried once by allele `a` and observed `c = 20` times:
- `S(b,b) = 0`;
- `S(a,b) = 129.30`;
- `S(a,a) = 131.88`.

So `(a,a)` is preferred by 2.6 nats.

Under `--library auto`, the library is estimated at this point.
- **Inserts.** Take read pairs placed by anchor 31-mers (single-copy, same strand, in every template). Up to 20,000 pairs are sampled by a salted hash of the read name, so the sample does not depend on file order. A pair is decisive when it is inward (FR) and implies the same insert on every template.
  - The robust centre is the median `m`, with robust sd `s = 1.4826·MAD`.
  - `insert_mean_bp` and `insert_sd_bp` are the mean and sd of the pairs within `6s` of `m`.
- **Intensity.** The fragment-start intensity comes from the depth fit:

  ```text
  λ = r / (2·(L̄ − k + 1)·(1 − ε)^k)
  ```

  - `L̄` is the mean of the two mates' median lengths.
  - `ε` is the mismatch rate at invariant sites over the first 2,000 accepted pairs.

  Example: `r = 11.28`, `L̄ = 150`, `k = 31` and `ε = 0.002` give `λ = 11.28 / 225.6 = 0.050`.
- **Refusal.** The run refuses with `LIBRARY_REFUSED` on:
  - fewer than 500 decisive pairs;
  - a bimodal insert distribution (an Otsu split with groups ≥ 4 sd apart and the smaller ≥ 5%);
  - more than 5% insert outliers;
  - more than 5% non-FR pairs;
  - an unconverged fit or a non-positive `λ`.
- **Everything else** (`support_sigmas` 4, `substitution_error` 0.002, `outlier_weight` 0.001) comes from the preset.

### 4. Propose candidate diplotypes

The **panel baseline** is every unordered pair of templates, deduplicated by allele path. It is always scored.

**Off-panel mosaics** come from a beam over blocks. Its unit is a pair of allele prefixes. Each prefix `p` carries `u_p(i)`, the probability mass of template paths that spell `p` and end on template `i`, under the transition of stage 5:

```text
u'(i) = (1 − ρ)·u_p(i) + (ρ/T)·Σ_j u_p(j)      for each template i carrying the next allele
mass(p1,p2) = Σ_b S_b + ln Σ u_p1 + ln Σ u_p2 + ln 2·[p1 ≠ p2]
```

`T` is the number of templates and `ρ` the switch probability. The best `--proposal-state-limit` (256) prefix pairs survive each block. The run stops `CANDIDATE_INCOMPLETE` if the beam exceeds its work or memory budget.

The first `--mosaic-candidates` (128) complete diplotypes that are not panel pairs join the panel baseline. Together they form the **union** to be scored.

### 5. Compute the Li–Stephens prior

A haplotype is a path of block alleles. Its prior follows a Li–Stephens copying model over the `n` panel templates:
- **Start:** `1/n` for each template carrying the first block's allele.
- **Transition between blocks:** stay on the same template with probability `1 − ρ + ρ/n`; move to any other with probability `ρ/n`.
- **Emission:** 1 if the template carries the path's allele at that block, 0 otherwise.

`ρ` is `--switch-probability` (0.01). It applies at every boundary, context boundaries included, whatever the distance.

The forward recursion sums over all copying histories. Ranking uses the **sequence prior** `h(s)`, which sums this over every allele path that spells the same bytes. An unordered pair of distinct sequences counts both orders:

```text
prior(s1, s2) = h(s1) + h(s2) + ln 2·[s1 ≠ s2]
```

For example, take 4 templates and 4 boundaries:
- A haplotype spelled exactly by two templates and needing no switch has `ln h ≈ −0.703`.
- One that needs one switch has `ln h ≈ −5.303`. Here the switch costs 4.6 nats.

### 6. Place read pairs and score each candidate

Each candidate `{a, b}` is scored by an exact paired-read log-likelihood:

```text
log L(a,b) = −(1−η)·λ·(E_a + E_b) + Σ_f ln[ (1−η)·λ·(M_fa + M_fb) + η·P_bg(f) ]
```

The terms are:
- **`λ`:** the fragment-start intensity per base and haplotype copy.
- **`η`:** `outlier_weight`, the fraction of fragments from outside the model.
- **`P_bg(f) = 0.25^(L1+L2)`:** the background probability of fragment `f`'s two mate sequences.
- **`M_fh`:** the placement mass of fragment `f` on haplotype `h`, summed over every placement. One placement contributes `½ · P_lib(insert) · e(mate1) · e(mate2)`:
  - `½` is the mate-orientation prior.
  - `P_lib` is a discrete normal insert distribution on `[μ − 4σ, μ + 4σ]` (`support_sigmas` 4), normalised over that support. The insert must be at least the longer mate's aligned span.
  - `e(·)` is each mate's emission under a banded affine-gap pair-HMM, normalised over reads of that length:
    - a match emits `1 − ε` and a mismatch `ε/3` (`ε` = `substitution_error`);
    - a gap opens with probability 1e-4 and extends with probability 0.5;
    - a placement may use at most `--mismatch-band` (6) edits, substitutions plus gap bases, per mate.
- **`E_h`:** the exposure, the expected number of fragment starts on haplotype `h`:

  ```text
  E_h = Σ_F (n_F/N) · Σ_i P_lib(i)·(L_h − i + 1)⁺
  ```

  - `L_h` is the haplotype length.
  - The outer sum averages over the observed distribution of fragment length floors `F`, each the longer mate's length.

  A homozygote counts two copies.

Placements are found by splitting each mate into `band + 1` pieces and indexing them. By the pigeonhole principle, a placement with at most `band` edits leaves at least one piece exact. Each candidate start is then verified by dynamic programming.

For example, take `μ = 400`, `σ = 50`, `λ = 0.05` and `η = 0.001`:
- **Insert.** An insert of 450 bp has `P_lib = 4.84e-3`.
- **One placement.** With two error-free 150 bp mates (`ln e ≈ −0.33` each), `ln M = ln ½ − 5.33 − 0.65 = −6.68`.
- **Fragment contributions:**
  - placeable on one haplotype copy: `ln(0.999·0.05) − 6.68 = −9.67`;
  - placeable on both copies: `−8.98`;
  - placeable on neither: `ln 0.001 + 300·ln 0.25 = −422.8`.
- **Exposure.** A haplotype of 10,000 bp has `E = 9,601`. Every extra 1,700 bp of haplotype costs `0.04995 · 1,700 ≈ 84.9` nats of exposure. Reads must pay for that length.

### 7. Rank candidates and form equivalence sets

The ranking score is `log L + prior`. Candidates spelling the same unordered pair of sequences are grouped and ranked once.

Two scores are treated as tied when floating-point error cannot separate them. Each candidate `c` carries a bound

```text
δ_c = max( (F + 16)·ε_double·(Σ|terms_c| + |prior_c|), 0.5e-9 )
```

- `F` is the number of fragments.
- `ε_double ≈ 2.2e-16`.
- `Σ|terms_c|` is the sum of the absolute values of the candidate's per-fragment log terms.

A rival is inseparable from the best when `(s_best − δ_best) − (s_c + δ_c) ≤ 0`. The best candidate and its inseparable rivals form the **top set**. With 20,000 fragments and `Σ|terms| ≈ 241,000`, `δ ≈ 1.1e-6` nats.

### 8. Repair the incumbent and branch at ties

The union is ranked first. Repair then starts from an **anchor**: the best-ranked group containing a panel pair, represented by its member with the highest path prior. Without a panel anchor the run stops with `REPAIR_INCOMPLETE` (`no_panel_anchor`).

**Neighbourhood bound.** Before any work, the size of one repair neighbourhood is checked:

```text
Σ_b 2·min(32, A_b(A_b+1)/2) ≤ --repair-max-neighbours (670)
```

`A_b` is the number of alleles in block `b`. Otherwise the run stops with `neighbourhood_bound`. For example, blocks with 1, 3, 12, 2 and 40 alleles give `2 + 12 + 64 + 6 + 64 = 148`.

**Each repair round:**
1. For every block, take its `--repair-local-pairs` (32) best unordered allele pairs by marker emission.
2. Substitute each pair into the incumbent, in both assignments to the two haplotypes.
3. Score every unvisited neighbour exactly, as in stages 5 and 6. The incumbent must reproduce its own score within 1e-6.
4. Move only if the certified winner is a different candidate scoring above the incumbent.

The search converges when no neighbour improves. A move in the last of the `--repair-rounds` (32) rounds stops it as `round_cap`, which is not converged.

**Tie branching.** When an improving round ends in a top set that excludes the incumbent, every tied neighbour becomes a branch (`--repair-tie-branching`). Branches run breadth-first in canonical order, each with its own repair search, and split again at their own ties. If more than `--repair-tie-max-branches` (16) starts would be queued, or a branch ends unconverged, the run no-calls with `REPAIR_INCOMPLETE` (`tie_branching_incomplete`).

**Pooled decision.** Every candidate scored by the union, the repair and the branches enters one pool. A candidate scored twice must give bit-identical scores. The pool is ranked as in stage 7:

| top set | status | call class |
|---|---|---|
| one group | `COMPLETE` | strict call: phased haplotypes |
| several groups, the same unordered genotype at every block | `CONDITIONAL_GENOTYPE_EQUIVALENT` | genotype-equivalent call: no phase |
| groups that disagree at some block | `RESCORE_INCOMPLETE` (`distinct_tie`) | no-call |

**Read and phase support.** Both come from the same decision. A rival is a **fragment contender** unless its fragment term (the sum over fragments) is certifiably below the winner's.
- **Phase.** `phase_support` is `READS_WITHIN_CANDIDATES` with no contender, `PRIOR_ONLY` otherwise, and `UNRESOLVED` for a non-strict call.
- **Read support.** A block's `genotype_support` is `PRIOR_ONLY` when some contender has a different genotype there, and `READS` otherwise.

**Conditional posterior and GQ.** Each block's posterior is conditional on the pool:

```text
log_err_b = logsumexp{ s_c : genotype_c(b) ≠ called(b) } − logsumexp{ all s_c }
posterior = 1 − exp(log_err_b)        GQ = −10·log10(exp(log_err_b))
```

For example, scores −1000 (called), −1002 and −1005, where both rivals differ at `b`, give a posterior of 0.876 and a GQ of 9.05. GQ is reported as uncalibrated, because the pool is not every possible diplotype.

### 9. Annotate representation, read support and library sensitivity

These annotations never change the call.

**Representation check** (`--representation-guard`). This asks whether the call fits the fragments while the markers and the exposure point elsewhere. For each block:
- **Competitor.** `T` is the marker-best pair. If its marker margin `m = S_b(T) − S_b(called)` is at least `ln 1000 = 6.908` nats, `T` is substituted into the call (the better assignment) and scored exactly.
- **Statistic.** With `F` the fragment term and `X` the exposure term, the block is flagged `LOWQUAL_REPRESENTATION_CONFLICT` when all three differences are at least `ln 1000`:

  ```text
  statistic = min( m, F_call − F_comp, X_comp − X_call ) ≥ ln 1000
  ```

  The conflict is that the markers and the exposure prefer the competitor while the fragment term prefers the call. For example, `min(12, 10, 10) = 10` is flagged.
- **Not evaluated.** A competitor that cannot be scored leaves the block `NOT_EVALUATED`.

`ln 1000` is a declared, conservative threshold, not a calibrated Bayes factor.

**Library sensitivity** (under `--library auto` only):
- **Two views.** The scored pool is reranked at 0.9λ and at 1.1λ. The whole search is also rerun at each intensity.
- **Classification.** A block whose unordered genotype changes under either view is `LIKELIHOOD_RANK` (the rerank), `SEARCH_PATH` (the rerun) or `BOTH`. The locus is `SENSITIVE` if any block is, and `NOT_EVALUATED` if a view cannot complete.

**Block quality.** `quality_filter` is PASS only for a representation PASS with `READS` support that is not library-sensitive. Otherwise it is LOW_QUALITY, with every reason that applies.

### 10. Spell and write the haplotypes

**Spelling.** A haplotype's sequence is the concatenation of its block alleles. For a strict call, the representative member of the winning group is spelled into `hap1` and `hap2`. The FASTA header records `phase_support` and `representation=` (`read_supported_within_candidates` or `prior_selected`). A genotype-equivalent call writes no FASTA, and its block tables carry the unordered pair.

**Block tables.**
- `block_genotypes.tsv` has every block, context included.
- `genotypes.tsv` has the variant blocks, with each allele's panel source alleles.

**Template switches.** For each haplotype, the locus table counts the fewest switches between panel templates needed to spell it, found greedily from the left: keep the set of templates that carry every allele so far, and count a switch where that set becomes empty.

### 11. Audit dosage

The audit asks one question: do the reads contradict the number of copies the call implies? It runs on a call (strict or genotype-equivalent) with `--dosage-mode audit` or `repair`. A no-call is `NOT_RUN`.

**Census.**
- **Molecule span:** `S = ⌈μ + 4σ⌉`.
- **Allele length:** `L_a` is the allele's length, counting each node visit.
- **Dosage-only pair:** an ordered allele pair `(i, j)` where `L_j ≥ L_i + S`, and at most 50 bp of `i`'s distinct nodes are absent from `j`.
- **Sensitive block:** one with a dosage-only pair. Its route is `repeated_unit` if a node is visited at least twice and its repeated bases exceed `S`; otherwise it is `pair`. Other blocks are `NOT_SENSITIVE`.

For example, `μ = 400` and `σ = 50` give `S = 600`. Alleles of 1,000 and 2,700 bp sharing all of the shorter one's nodes form a dosage-only pair.

**Features and reference.**
- **Features.** A sensitive block's canonical 31-mers are grouped by their copy count in each allele; a group of at least 500 k-mers that vary between alleles is a feature. A feature is tested if at least 500 of its k-mers remain after removing those found outside the block in either reported homologue.
- **Reference.** The single-copy reference `R` is the set of k-mers from non-sensitive blocks that occur exactly once in each homologue and nowhere in the block.

**Estimate.** Read-pair hits give each feature `f` a copy-number estimate:

```text
ĉ_f = 2·(T_f/|f|) / (T_R/|R|)
```

`T_f` and `T_R` are hit counts, and `|f|` and `|R|` are k-mer counts. For example, `|R| = 1,000` with `T_R = 30,000`, and `|f| = 600` with `T_f = 36,000`, give `ĉ = 4`.

**Bootstrap.** Two hundred bootstrap resamples of the read pairs, with a fixed seed per block, give the quantiles `q05…q95` and a robust sd `(q75 − q25)/1.349`. The admissible copy numbers are the integers in `[⌈q05⌉, ⌊q95⌋]` that the panel's alleles can produce. The call's copy number is `copies(a1) + copies(a2)`.

**Classification.** With `m` tested features over the locus:

```text
z_c = Φ⁻¹(1 − 0.001/(2m))           z = |call − q50| / robust_sd
trigger      ⇔  z ≥ z_c  and  call not admissible
suspect      ⇔  z ≥ 2.5
underpowered ⇔  robust_sd > 1/(z_c + Φ⁻¹(0.90))
```

The block status is the first that applies: `AUTO_REPAIR` (any trigger), then `SUSPECT`, then `UNDERPOWERED`, then `CONSISTENT`.

For example, `m = 4` gives `z_c = 3.662` and an underpowered threshold of 0.202. With `q25/q50/q75 = 3.8/4.0/4.1` the robust sd is 0.222. A call of 2 copies then has `z = 9.0` and is not admissible (`{4}`), so the block triggers.

### 12. Repair dosage: complete diplotype or conditional content

Each `AUTO_REPAIR` block searches for dosage-admissible alternatives. Each alternative is spliced into the reported homologues and scored with the exact likelihood of stage 6. The reported diplotype must first reproduce its own terms within 1e-6.
- **Pair route.** Every panel allele pair whose copy numbers are admissible for every feature, scored in both orientations (at most 4,096 pairs).
- **Repeated-unit route.** New alleles are built from multisets of the repeated units. A negative-binomial screen ranks them under construction caps of 64 multisets, 256 haploids and 4,096 diploids.

**Choice.** The likelihood alone chooses; the prior is recorded but does not rank. The top set is every candidate within 1e-9 nats of the best:
- More than one genotype in it is a tie (`top_genotype_tie`), and nothing is repaired.
- On the repeated-unit route, arrangements of the same content that the reads cannot separate form an allocation-equivalence set. Its representative is the pair with the smallest canonical identifiers, labelled `REPRESENTATIVE_ONLY_NOT_RESOLVED`.

**Scope.** Two separate facts are recorded. `dosage_content_scope` is COMPLETE when every admissible candidate shares one copy-number content. `diplotype_scope` is COMPLETE when the search space was not truncated.

| content | diplotype | class | status |
|---|---|---|---|
| COMPLETE | COMPLETE | `COMPLETE_DIPLOTYPE_REPAIR` | `REPAIRED` |
| COMPLETE | CONDITIONAL | `CONDITIONAL_DOSAGE_CONTENT` | `CONTENT_REPAIRED` |
| CONDITIONAL | any | `TRIGGERED_UNREPAIRED` (`dosage_content_conditional`) | `TRIGGERED_UNREPAIRED` |

**Application.** In repair mode a complete repair replaces the block's genotype only when all of the following hold; otherwise the failing reason is recorded:
- the run is a strict call (`not_a_strict_call`);
- it is the locus's only complete repair (`several_complete_repairs_not_jointly_scored`);
- the reported call and haplotypes are reproduced;
- both repaired sequences are panel alleles (`repaired_allele_not_a_panel_allele`).

Conditional content never changes the reported alleles.

### 13. Recompute the annotations after an applied repair

When a repair is applied, the caller's outputs are first copied to `<prefix>.original.*`. The annotations of stage 9 are then recomputed on the reported genotype:
- **Representation check:** first rerun on the original call, where it must reproduce the caller's result exactly, then on the repaired diplotype.
- **Read support:** reproduced from the search record, then re-evaluated against the caller's candidates and the other repair alternatives.
- **Block quality.**

If the caller's own annotation cannot be reproduced, the value is `NOT_EVALUATED_AFTER_DOSAGE_REPAIR`; it is never copied and never PASS. The repaired block's posterior and GQ become NA. The FASTA is re-spelled with `dosage_repaired_blocks=` and `dosage_repair_role=RESOLVED`. The status file keeps each superseded value as `original_<key>`.

### 14. Write locus, provenance and status outputs

`locus_genotype.tsv` gathers the call into one row. It has the status and call class, the haplotypes' SHA-256 digests and block paths, the unordered genotype, the score terms of the winner (or of the repaired diplotype), the separation from the best other genotype, the template switches, the dosage summary and links to every file.

The dosage stage's repaired sequences are compared with the panel in `dosage_provenance.tsv`. It reports the nearest panel haplotype and the nearest panel block allele, by global unit-cost edit distance (all tied nearest sequences, ordered by SHA-256), with the alignment's CIGAR and the distance divided by the alignment columns.

The run's status is written to `mosaic_status.tsv` and every effective setting to `settings.tsv`. A no-call writes no genotype tables. The exit status is 0 for a strict or genotype-equivalent call and 1 for a no-call.

### 15. Derive the headline filters

Each of the three filters gives every unit the most severe applicable status: `NO_CALL > NOT_EVALUATED > LOW_QUALITY > PASS > NOT_APPLICABLE`. Its reasons are kept. The locus takes the worst applicable unit, with reasons written `REASON@ids`, and is NOT_APPLICABLE only if every unit is.

**GENOTYPE_FILTER, per block** (all reasons that apply):

| condition | status, reason |
|---|---|
| no call | NO_CALL `RUN_NO_CALL:<reason>` |
| representation not evaluated | NOT_EVALUATED `REPRESENTATION_NOT_EVALUATED` |
| library audit incomplete | NOT_EVALUATED `LIBRARY_SENSITIVITY_NOT_EVALUATED` |
| representation conflict | LOW_QUALITY `REPRESENTATION_CONFLICT` |
| library-sensitive block | LOW_QUALITY `LIBRARY_SENSITIVE` |
| dosage ran, and the block is content-corrected, triggered-unrepaired, or has a complete repair not applied | LOW_QUALITY `EXACT_DIPLOTYPE_UNRESOLVED` |
| none of the above | PASS |

**DOSAGE_FILTER, per block** (the first row that applies):

| condition | status, reason |
|---|---|
| no call | NO_CALL `RUN_NO_CALL` |
| `--dosage-mode off` | NOT_EVALUATED `DOSAGE_NOT_EVALUATED:dosage_mode_off` |
| the stage did not complete | NOT_EVALUATED `DOSAGE_NOT_EVALUATED:audit_<status>` |
| `CONSISTENT` | PASS `DOSAGE_CONSISTENT` |
| complete repair applied | PASS `COMPLETE_REPAIR_APPLIED` |
| complete repair, audit mode | PASS `COMPLETE_REPAIR_AVAILABLE_NOT_APPLIED` |
| complete repair not applied | LOW_QUALITY `COMPLETE_REPAIR_NOT_APPLIED:<reason>` |
| content corrected | PASS `CONTENT_CORRECTED_EXACT_DIPLOTYPE_UNRESOLVED` |
| `SUSPECT` / `TRIGGERED_UNREPAIRED` | LOW_QUALITY `DOSAGE_SUSPECT` / `TRIGGERED_UNREPAIRED` |
| `UNDERPOWERED` | NOT_EVALUATED `DOSAGE_UNDERPOWERED` |
| not dosage-sensitive | NOT_APPLICABLE `NOT_DOSAGE_SENSITIVE` |

**PHASE_FILTER, per adjacent pair of table rows `L>R`.** A pair in which either block is homozygous is NOT_APPLICABLE (`NO_PHASE_INFORMATIVE_BOUNDARY`). Otherwise the status follows the call's `phase_support`:

| `phase_support` | status, reason |
|---|---|
| genotype-equivalent call or `UNRESOLVED` | NO_CALL `PHASE_UNRESOLVED` |
| read-supported | PASS |
| `PRIOR_ONLY` | LOW_QUALITY `MODEL_ONLY_PHASE` |

A boundary next to an applied repair is NOT_EVALUATED (`PHASE_NOT_ASSESSED_AFTER_DOSAGE_REPAIR`).

`call_summary` is `GENOTYPE_<s>;DOSAGE_<s>;PHASE_<s>`. The per-status unit counts beside each filter show how many blocks or boundaries are usable when the locus value is not PASS. `genotype_support`, `quality_filter` and GQ do not enter the headline.

## Worked trace

The steps below follow the fifteen above. The locus has five blocks. Variant blocks `sv:1` (alleles `a`, 1,000 bp; `b`, 2,700 bp) and `sv:3` (alleles `c`, `d`) alternate with invariant context. The panel has four templates, `A = (a,c)`, `B = (b,d)`, `C = (a,c)` and `D = (b,d)`. The sample carries `(a,c)` and `(a,d)`; the second haplotype is in no panel template. The library is declared, `μ = 400`, `σ = 50`, `λ = 0.05`, and the dosage stage is off. Numbers are illustrative.

```text
blocks     ctx:0   sv:1      ctx:2   sv:3     ctx:4
template A   ·     a          ·      c          ·
template B   ·     b          ·      d          ·
template C   ·     a          ·      c          ·
template D   ·     b          ·      d          ·
sample       ·     a / a      ·      c / d      ·
```

1. **Inputs.** Five blocks and four boundaries. `A` and `C` spell identical bytes, so their allele paths collapse into one panel haplotype, `(a,c)`.
2. **Markers.** A syncmer present only in allele `a` is a `sv:1` marker with `n(a) = 1`, `n(b) = 0`. Syncmers in every allele of every block are anchors.
3. **Depth.** The anchors fit `r = 11.28`, `bg = 0.01`. If the `a`-marker is seen 20 times, `S(a,a) = 131.88` against `S(a,b) = 129.30` and `S(b,b) = 0`.
4. **Candidates.** The panel pairs are `{ac,ac}`, `{ac,bd}` and `{bd,bd}`. The beam adds the off-panel mosaic `{ac,ad}`: `a` from `A`, then a switch to `B` for `d`.
5. **Prior.** `ln h(ac) = −0.703`, `ln h(bd) = −0.703`, `ln h(ad) = −5.303`. So `prior{ac,bd} = −0.713` and `prior{ac,ad} = −5.313`: the mosaic pays 4.6 nats for its switch.
6. **Likelihood.** `{ac,bd}` predicts reads across 1,700 extra bases of allele `b`, and none are seen. Its exposure term alone is 84.9 nats worse than `{ac,ad}`'s, and reads that fit `a` but not `b` add more. `{ac,ad}` wins on the likelihood by far more than its prior deficit.
7. **Ranking.** The union's best group is `{ac,ad}`, separated beyond every error bound.
8. **Repair.** The anchor is the best panel group, `{ac,bd}`. Its neighbour that sets `sv:1` of the second haplotype to `a` is `{ac,ad}`, which improves the score. The next round finds no improvement. The pool's top set is one group, so the run is `COMPLETE`: a strict call. Suppose the only other groups score −1002 and −1005 against the winner's −1000, both differ at `sv:3`, and both have fragment terms certifiably below the winner's. Then no rival is a fragment contender: `phase_support = READS_WITHIN_CANDIDATES`, `genotype_support = READS`, and `sv:3` has posterior 0.876 and GQ 9.05.
9. **Annotation.** At each block the marker-best pair is the called pair, so the representation check passes. Read support is `READS` at both variant blocks, so `quality_filter` is PASS.
10. **Spelling.** `hap1 = ctx:0 + a + ctx:2 + c + ctx:4`, `hap2 = ctx:0 + a + ctx:2 + d + ctx:4`. `hap1` needs 0 template switches and `hap2` needs 1.
11.–13. **Dosage.** Skipped (`--dosage-mode off`). With `--dosage-mode audit`, `sv:1` would be dosage-sensitive (`S = 600`), and a call of `b` on one copy would be flagged by the copy-number test.
14. **Outputs.** The locus row, the two block tables, the FASTA, the status and the settings are written. Exit status 0.
15. **Headline.** Every block passes GENOTYPE. DOSAGE is NOT_EVALUATED (`dosage_mode_off`). Every adjacent pair of rows includes a homozygous block, context blocks included, so PHASE is NOT_APPLICABLE.

```text
block  allele1 allele2 zygosity genotype_support genotype_filter dosage_filter
sv:1   a       a       HOM      READS            PASS            NOT_EVALUATED
sv:3   c       d       HET      READS            PASS            NOT_EVALUATED

call_class=strict_call  phase_support=READS_WITHIN_CANDIDATES
call_summary=GENOTYPE_PASS;DOSAGE_NOT_EVALUATED;PHASE_NOT_APPLICABLE
```
