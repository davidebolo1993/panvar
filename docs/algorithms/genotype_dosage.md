# Dosage audit and repair: normative specification (approved with amendments, 2026-09-28)

This document is the contract for the C++ port of the dosage audit and repair add-on to `panvar genotype`. It freezes
the design Sol accepted on 2026-09-28, with Sol's ten amendments of the same day (§17). The Python prototype is the
executable reference: `tests/results/genotype_c4_production_11fe681_dosage_repair/`.
- **The sealed reference outputs** are the policy v2 tables `contiguity/class_baseline_{campaign,fresh}.*` and the
  repairs in `phase1/` and `policyC_test/`.
- **Normative words.** "Must" and "never" are requirements. A value in `code` is a frozen constant, and changing it is
  a new design, not a tuning step.

The add-on is truth-free. It reads the prepared model, the call the caller reports, the reads and the library, and
nothing else. It is a **short-read** specification. Nothing in it applies to a long-read model; a run whose library
does not declare the short-read keys of §2 is refused (§11).

**Exact tie.** Two scores are tied when they differ by at most `1e-9` nats. This is the production
`CertificationOptions::minimum_tie_tolerance_nats`; the prototype's top sets use the same tolerance.

## 1. Scope and states

For every block of a genotyped locus, the audit answers one question: do the reads contradict the copy numbers the
reported genotype implies?

When they contradict it beyond a family-wise error rate, the block is repaired inside the dosage the reads allow.
Otherwise the reported genotype is kept and annotated.

Repair never happens for any other reason. In particular:
- a correct-looking likelihood gain, a weak dosage deviation, low read support, representation conflict or an
  unresolved allocation never trigger it.

Blocks that the census (§2) does not mark as copy-number-sensitive are never examined: `dosage_status NOT_SENSITIVE`.

## 2. Census: copy-number-sensitive blocks (graph and library, once per run)

**Molecule span.** S = ⌈`insert_mean_bp` + 4 × `insert_sd_bp`⌉, from the run's short-read library.
- C4's libraries declare 350 and 40, so S = 510, the prototype's value.
- S is recorded in the status (`dosage_molecule_span_bp`) and in every dosage block row.

For each block, over its panel alleles as the prepared model holds them (each allele's node path; node weight = node
length in bp):

- **len(a):** the allele's length, i.e. the sum of node weights with multiplicity.
- **novel(a, b):** the bp of distinct nodes present in `a` and absent from `b`.
- **Dosage-only pair (a, b):** len(b) − len(a) ≥ S, and novel(a, b) ≤ `50` bp.
  - The shorter allele differs from the longer one essentially by missing sequence, which is what copy number
    measures.
  - `50` bp is the minimum structural variant size. It is declared separately from S.
- **Copy-number-sensitive:** the block has at least one dosage-only pair.
- **Repeated-unit route** (§8.2): the block is copy-number-sensitive, some allele visits a node at least twice, and
  that allele's repeated sequence exceeds S bp. The repeated sequence is the bp of nodes visited ≥ 2 times in one
  allele. Otherwise a sensitive block uses the **pair route** (§8.1).
- **Census-invariant:** every other block. These blocks provide the single-copy reference (§4).

The census uses no reads, calls or truth. C4 yields sv:1 (pair route) and sv:7 (repeated-unit route) in all 66
models.

## 3. Pooled dosage features (graph, then per call)

- **k-mers.** Canonical `31`-mers of every panel allele of the block, counted with multiplicity: A[a, x] = number of
  occurrences of k-mer x in allele a.
- **Varying k-mers:** the column A[·, x] is not constant across alleles.
- **Groups:** varying k-mers with an identical column. A group with ≥ `500` k-mers is a **pooled feature** with copy
  vector s = that column. Features are ordered by the first occurrence of any of their k-mers, scanning alleles in
  ascending allele ID and positions left to right.
- **Per-call exclusion.** For the reported genotype, remove from every feature the k-mers occurring anywhere in
  either reported homologue outside this block's interval. A feature left with ≥ `500` k-mers is **tested**.
- **No tested feature.** A sensitive block with no tested feature is `NOT_EVALUATED` / `no_tested_feature`.
- **m** for a sample = the number of tested features summed over all copy-number-sensitive blocks of the locus.
  - m counts only blocks that reach the estimate; blocks that fail closed in §3–§5 contribute nothing.
  - If m = 0, no z_critical is computed: `z_critical NA`, and every sensitive block keeps its fail-closed status.
- **Call copy number** of feature j: c_j = s_j[allele1] + s_j[allele2], using the two reported block alleles.
- **Feasible range:** [2·min s_j, 2·max s_j].

These features are final:
- The contiguity split and DOSAGE_MODEL_ALIASED were rejected (§13).
- No other feature family is part of this design.

## 4. Single-copy reference (per block, per call)

The reference R is the set of canonical k-mers that satisfy all of the following:
- they occur in both reported homologues within census-invariant block intervals;
- they occur exactly once in each full reported homologue;
- they are absent from every panel allele of the block being audited.

If R is empty, the block is `NOT_EVALUATED` / `no_reference`.

## 5. Estimate and bootstrap

- **Hits.** Every canonical k-mer of both mates of every read pair, with multiplicity. Each is assigned to R (index
  0) or to the tested feature holding it. A read pair contributes a count vector v = (v_0, v_1, …).
- **Retained pairs.** Pairs whose vector is all zero are dropped. The retained pairs, in input order, are the
  bootstrap population P of size n.
- **Point estimate.** T_j = Σ v_j over P. If T_0 = 0 (this includes n = 0), the block is `NOT_EVALUATED` /
  `no_reference_hits`. Otherwise CN_j = 2 · (T_j / |F_j|) / (T_0 / |R|), evaluated in this order in IEEE double.
- **Bootstrap.**
  - `200` replicates; each draws n indices uniformly with replacement and recomputes CN_j by the same formula.
  - The generator is the Python `random.Random(20260925)` stream (MT19937, `init_by_array` seeding, `randrange` by
    `getrandbits` rejection). It is re-seeded for every block and drawn sequentially: replicate r, then index i.
  - This is frozen so the port reproduces the prototype bit for bit. Draw order is independent of thread count.
- **Zero reference in a replicate.** If any replicate has zero reference hits, the block is `NOT_EVALUATED` /
  `bootstrap_zero_reference`.
  - No replicate value is computed from a zero denominator, and no NaN or infinity reaches a quantile.
  - The prototype has no rule for this case: it would have raised ZeroDivisionError, and no sealed run raised it.
    Failing closed is therefore a new, declared rule, not a change to replayed behaviour.
- **Quantiles.** q(p) = sorted replicates[⌊p·(200 − 1)⌋] for p ∈ {0.05, 0.25, 0.50, 0.75, 0.95}.
- **Robust sd** = (q(0.75) − q(0.25)) / 1.349. **Robust z** = |c_j − q(0.50)| / robust sd. If robust sd = 0, then
  z = 0 when c_j = q(0.50), else +∞. That infinity is an explicit rule, never the result of a division.
- **Admissible set** (the declared 90% interval):
  - the integers in [⌈q(0.05)⌉, ⌊q(0.95)⌋] that lie in the feasible range;
  - if there are none, the single feasible integer nearest to the interval (distance to the nearer end, then the
    smaller integer).
  - c_j is **outside** when c_j is not in this set.

## 6. Classification (policy v2, frozen)

- **Per sample:** α = `0.001` family-wise; z_critical(m) = Φ⁻¹(1 − α / (2m)).
- **Per feature:**
  - trigger ⇔ robust z ≥ z_critical(m) and c_j is outside;
  - suspect ⇔ robust z ≥ `2.5`;
  - underpowered ⇔ robust sd > 1 / (z_critical(m) + Φ⁻¹(`0.90`)). A one-copy error would reach z_critical with
    probability below 0.90.
- **Per block**, by precedence:
  1. any feature triggers: `AUTO_REPAIR`, which then ends as `REPAIRED` or `TRIGGERED_UNREPAIRED` (§8);
  2. else any feature suspect: `SUSPECT`;
  3. else any feature underpowered: `UNDERPOWERED`;
  4. else `CONSISTENT`.

**What the statuses mean:**
- SUSPECT is a graded diagnostic. It never triggers repair and is not a prediction of error.
- UNDERPOWERED means the reads cannot verify dosage at this coverage. It is an information limit, not a warning
  about the genotype.

## 7. What is recorded per feature

Record α, m, z_critical, feature index, k-mers after exclusion, copy vector values, c_j, CN_j, the five quantiles,
robust sd, robust z, the admissible set, whether c_j is outside, and the trigger, suspect and underpowered decisions.

## 8. Repair (AUTO_REPAIR blocks only)

**Dosage-admissible candidate:** a diploid block genotype whose copy number of every tested feature lies in that
feature's admissible set.

**Exact score.** Every candidate is spliced into the reported homologues at the block's intervals.
- **Pair route:** both orientations (one if the two block sequences are equal).
- **Repeated-unit route:** one orientation, as the prototype scored it: the lexicographically smaller unit tuple goes
  in the first reported homologue (§18). It is scored by the production likelihood: fragment term +
exposure term, with the run's own scoring settings. The **likelihood** is that sum. The prior is recorded but never
ranks.

**Winner guard.** Before any candidate is scored, the reported diplotype alone must reproduce its persisted total,
fragment term and exposure term within 1e-6 in the same scoring context. Otherwise the repair fails closed:
`TRIGGERED_UNREPAIRED` / `winner_not_reproduced`.

**Resource refusal.** If scoring hits a cumulative scorer cap or the run's watchdog, the block ends
`TRIGGERED_UNREPAIRED` / `resource_refusal`, and the reported genotype is kept.

**Top-genotype rule** (both routes). Numeric IDs and ranks may order enumeration, but never choose between distinct
unordered block genotypes.
1. **Merge.** Scored candidates are grouped by their unordered pair of spelled block sequences. A group's score is
   its best orientation.
2. **Top set.** The top set is every group within the tie tolerance of the maximal likelihood.
3. **Distinct genotypes tied.** If the top set holds more than one group, the reported genotype is kept:
   `TRIGGERED_UNREPAIRED` / `top_genotype_tie`, with every tied candidate reported.
   - **Exception:** on the repeated-unit route, top-set groups that share one content multiset are resolved by the
     allocation rule (§8.2). That rule never resolves a tie between different content multisets.
4. **Orientation ties.** A tie between the two orientations of one unordered block genotype is not a genotype tie.
   Orientation AB (the first block sequence in the first reported homologue) is taken and `orientation_tie` is
   recorded. Orientation carries no block-genotype information, and `phase_support` is unchanged.

In the sealed prototype repairs, every top-set tie (39 blocks) lay within one content multiset, so this rule changes
none of them.

**Repair scope.** Every repair records `repair_scope` and `repair_scope_reason`.
- **COMPLETE:** every dosage-admissible candidate of the declared candidate model was constructed and exactly scored.
- **CONDITIONAL** otherwise. The reasons, in canonical order, comma-joined:
  - `pool_truncated`: some internal- or terminal-capable unit fell outside the 36-unit pool;
  - `cap_multisets`, `cap_haploid`, `cap_diploid`: a construction cap stopped expansion;
  - `proposal_ranked`: admissible multisets beyond those constructed exist, ordered only by the proposal
    likelihood, which is not an upper bound on the exact score.
- A proposal rank or a screening weight is never treated as proof that an omitted candidate cannot win.
- The shadow implementation reproduces the prototype, conditional repairs included. The integration decision reports
  how many repairs are conditional.

### 8.1 Pair route

- **Candidates:** every unordered pair {a, b} (a = b allowed) of the block's panel alleles with s_j[a] + s_j[b]
  admissible for every tested feature.
- **Cap:** `4096` pairs. Above it: `TRIGGERED_UNREPAIRED` / `cap_pairs`.
- **None:** `TRIGGERED_UNREPAIRED` / `no_admissible_candidate`.
- In both cases the reported genotype is kept.
- **Choice:** the top-genotype rule. When the pairs are scored completely, `repair_scope COMPLETE`.
- **Allocation:** `allocation_status NOT_APPLICABLE`. The block alleles are whole panel alleles.

### 8.2 Repeated-unit route

**Unit library (graph only).** Over the block's panel allele walks:
- **Entry node:** the first node visited twice, scanning walks in order, in the allele with the most repeated
  visits.
- **Flanks:** the node-level longest common prefix (stopping before the entry node) and the longest common suffix.
- **Units:** the walk between the flanks, cut at every visit of the entry node, canonicalised by complete sequence.
- **Unit roles:**
  - internal-capable: its last node has a panel edge into the entry node;
  - terminal-capable: its last node has a panel edge into the right flank;
  - anchor-compatible: it leaves unchanged the bases the model trims into the left anchor.
- **Unit state:**
  - optional run = a maximal run of nodes present in some unit copies and absent from others;
  - state `I` if the unit carries an optional run of ≥ `1000` bp, else `O`.
- **KMAX** = the maximum number of units in one panel allele. It is 3 in all 66 C4 models, equal to the
  prototype's constant.
- The block sequence of a unit tuple is left flank + units + right flank, trimmed by the model's anchor overhang.

**Feature constraint:**
- **u_j(unit)** = the median, over the feature's tested k-mers, of each k-mer's count in the unit's sequence.
- **M range:** every M in [2, 2·KMAX] for which each feature's admissible set meets
  [⌊2·min_T u_j + (M−2)·min_I u_j⌋, ⌈2·max_T u_j + (M−2)·max_I u_j⌉]. Here T are the terminal-capable units, I the
  internal-capable units, and the I term is 0 when M = 2.
- **Admissible multiset** (2 terminal + M − 2 internal units): round(Σ u_j) is admissible for every feature.

**Screening (count model).** This matches the prototype, `phase1/proposal_general.py`.
- **Profile classes:** unit k-mers with a non-constant count across units. Excluded: k-mers that can span a unit
  junction (left flank | unit, unit | entry, unit | right flank), and k-mers in the reported homologues outside the
  block.
- **Exposure U[S, unit]** = the union of read-start windows over the class's occurrences in the unit.
- **y[S]** = read pairs with ≥ 1 k-mer of class S.
- **β and φ** come from the reference k-mers of the **context blocks adjacent to the audited block in the chain**
  (one each side). For C4 sv:7 these are ctx:6 and ctx:8, which the prototype names explicitly. The rule reproduces
  them without naming them.
- **Error term:** e = k × substitution error / 3, with the substitution error taken from the library.
- **Mean and likelihood:** λ[S](n) = β Σ n_u U[S, u] + e β max_u U[S, u] M. The composite NegBinomial log-likelihood
  (Poisson when not overdispersed) is taken over classes.
- **Screening:** Poisson EM over unit weights (≤ `3000` iterations, stop at max change < 1e-10, start M0 = min M
  range).
- **Pool:** the top `36` internal-capable and top `36` terminal-capable units by weight (ties by unit ID).
- **Candidates:** every admissible multiset over the pool with M in the M range and M − 2 ≤ 2(KMAX − 1). They are
  ranked by composite log-likelihood descending, then multiset ascending.

**Construction.** Multisets are taken in rank order and expanded into arrangements. Caps: `64` multisets, `256`
haploid sequences, `4096` diploid structures.
- Construction stops at the first cap.
- A cap reached at the first multiset, or at an exact proposal-likelihood tie across the boundary, refuses:
  `TRIGGERED_UNREPAIRED` / `cap_refused`.
- Arrangements are both terminal assignments × every split of the internal units between homologues × every order
  within a homologue. Each homologue holds ≤ KMAX units, and each homologue's first unit is anchor-compatible.
  Structures with identical spelled sequences are merged.
- These caps and the pool make the scope `CONDITIONAL` (§8) unless the reasons listed there are all absent.

**Choice.**
1. **Content:** the top-genotype rule over all scored structures.
   - A top set within one content multiset fixes the content.
   - A top set spanning several multisets ends in `top_genotype_tie`.
2. **Allocation, the split between homologues:** among the scored structures of that multiset, maximise the **panel
   structure log-prior** log((f(σ_a) + 0.5)/N) + log((f(σ_b) + 0.5)/N).
   - σ is a homologue's unit-state string (e.g. `IO`), f(σ) the number of panel alleles with that state string, and
     N the number of panel alleles.
   - Ties go to likelihood, then structure rank, then orientation.
   - Reason: the transposition test shows that 350 bp-insert fragments cannot see the split. The likelihood and the
     original call's split are therefore uninformative about it. Of the rules compared, panel frequency was best on
     both datasets.
   - This is the only place a declared non-likelihood rule chooses between distinct block genotypes, and it applies
     only within one content multiset.
3. The reported block genotype is that structure, oriented by the orientation rule above.

### 8.3 Allocation observability test (annotation only; shadow-only)

It applies to the reported genotype of every repeated-unit block (repaired or not) whose two homologues hold ≥ 3
units in total. Every other block is `NOT_APPLICABLE`.

- **Switch point.** Take a common window: the first 2·`600` bp window of copy 0, centred at offsets
  `250`–`20000` in steps of `25`, that occurs exactly once in every copy. None: `INCONCLUSIVE`.
- **Transpositions.** For every alternative allocation A′ of the same copies (the construction rules of §8.2, both
  homologue orders), build T(A → A′):
  - each copy keeps its body after the switch point;
  - each position keeps the start region that followed the same left neighbour in A.
  - So T has A's junction windows everywhere, with a different allocation.
- **Scoring.** T and A are scored exactly. Take their per-fragment contribution vectors on a 1e-9-nat grid.
  - `UNRESOLVED`: some T has a vector identical to A's and an exposure within 1e-9. No fragment can distinguish the
    allocations.
  - `RESOLVED`: every T differs.

**Meaning.** ALLOCATION_UNRESOLVED is an identifiability annotation.
- It is never evidence that the genotype is wrong.
- It never triggers repair and never changes a call.

**Cost.** The test stays in shadow and off by default until Sol reviews its measured cost: incremental wall time, peak
memory and the number of blocks tested. Its switch is `--dosage-allocation-test`.

## 9. Repaired-sequence provenance (chosen repaired haplotypes)

This covers the **chosen** repaired haplotypes of every REPAIRED block. Losing candidates appear only in the candidate
table. Comparisons use only the panel supplied to the caller (in LOO validation, the LOO panel).

**Emitted for each chosen repaired haplotype:**
- the repaired full-haplotype FASTA, and the original full-haplotype FASTA (the reported homologue it replaces);
- **all equally nearest model-panel sequences, at two levels:**
  - full haplotypes: the full repaired haplotype against every spelled panel template haplotype;
  - block alleles: the repaired block allele against every panel allele of the block.
- **Distance:** raw global edit distance (unit costs; edlib, the existing submodule).
- **Normalised distance:** raw / global alignment length, where the length counts the `=`, `X`, `I` and `D` columns.
- **Deterministic extended CIGAR** (`=`, `X`, `I`, `D`) against each reported nearest sequence.
  - The alignment is edlib's global path, and the gate replays it.
  - Among equal-cost paths, edlib's traceback is deterministic for fixed inputs.
- **Ties:** every panel sequence at the minimal distance is reported. They are ordered by full SHA-256 of the
  sequence, then by panel ID, so the output is invariant to panel ordering.
- **Unit provenance:** each unit in order, with its unit ID, its sequence SHA-256 and every panel allele and position
  that carries it (repeated-unit route); or the two panel allele IDs (pair route).

**Generated-sequence ID.** The complete lowercase SHA-256 of the sequence.
- A display ID (the first 12 hex digits) may be shown alongside it.
- If two distinct generated sequences share a display ID, the shadow audit fails closed: `dosage_audit_status FAILED`,
  note `display_id_collision`.

## 10. Genotype quality after repair

- **Unrepaired blocks:** the production `conditional_GQ` keeps its production meaning. It is conditional on the
  scored candidate set and uncalibrated: `gq_status UNCHANGED_CONDITIONAL`.
- **REPAIRED blocks:** the ordinary GQ is **not** reused. The repaired genotype was chosen inside a different,
  dosage-constrained candidate set, and a repaired non-panel haplotype may have no spellable atomic prior.
  - `original_GQ` and `original_gq_status`: the production value and its conditional status, retained. Never
    overwritten.
  - `repaired_GQ NA`; `gq_status UNAVAILABLE_REPAIRED`.
  - `conditional_repair_gap` = likelihood(chosen) − likelihood(best admissible candidate with a different unordered
    block genotype), in nats; NA if there is none. It is conditional on the repair's candidate set and its scope.
  - `repair_separation` = likelihood(chosen) − likelihood(the original call), in nats. It can be negative: the
    likelihood alone is not a reliable judge of dosage, which is why repair is dosage-constrained.

## 11. Output schema (shadow; separate fields, never collapsed)

Existing outputs are unchanged. The audit writes new files next to them, and each row carries the fields below.
`quality_filter` is not recomputed or merged.

**`<out>.dosage_blocks.tsv`**, one row per block:

| field | values |
|---|---|
| block_id, census_route | pair / repeated_unit / not_sensitive |
| molecule_span_bp | S (§2) |
| dosage_status | CONSISTENT, REPAIRED, TRIGGERED_UNREPAIRED, SUSPECT, UNDERPOWERED, NOT_SENSITIVE, NOT_EVALUATED |
| dosage_reason | see the reasons list below |
| allocation_status | RESOLVED, UNRESOLVED, INCONCLUSIVE, NOT_APPLICABLE, NOT_RUN (test switched off) |
| read_support_status | the call's `genotype_support`, copied |
| representation_status | the call's representation `filter`, copied |
| alpha, m, z_critical, features_tested, max_robust_z | z_critical NA when m = 0 |
| original_allele1, original_allele2 | the reported block alleles |
| repaired_allele1, repaired_allele2 | NA unless REPAIRED; generated alleles by full SHA-256 |
| repair_route, candidate_source | pair / repeated_unit; panel_pair / unit_multiset |
| repair_scope, repair_scope_reason | COMPLETE / CONDITIONAL; reasons per §8 |
| candidates_constructed, candidates_scored, top_set_size, orientation_tie, cap_reached | |
| original_fragment_term, original_exposure_term, original_likelihood, original_prior, original_total | |
| repaired_fragment_term, repaired_exposure_term, repaired_likelihood, repaired_prior, repaired_prior_spellable | |
| repair_separation, conditional_repair_gap | §10 |
| original_GQ, original_gq_status, repaired_GQ, gq_status | §10 |

**Reasons** (`dosage_reason`): `none`, `no_tested_feature`, `no_reference`, `no_reference_hits`,
`bootstrap_zero_reference`, `no_admissible_candidate`, `cap_pairs`, `cap_refused`, `winner_not_reproduced`,
`top_genotype_tie`, `resource_refusal`, `root_evidence_noninvariant` (§12.1), `unsupported_library`.

**Other files:**
- **`<out>.dosage_features.tsv`**, one row per tested feature: the §7 fields.
- **`<out>.dosage_repair_candidates.tsv`**, one row per scored candidate orientation:
  - the route, label, rank, the unit tuples or allele pair, the state strings;
  - the proposal log-likelihood (repeated-unit route), fragment term, exposure term, likelihood, prior;
  - the panel-structure log-prior;
  - `in_top_set`, and whether it was chosen for content and for allocation.
- **`<out>.dosage_allocation.tsv`**, one row per transposition: A, A′, the fragment-vector hashes (bitwise and 1e-9
  grid), the exposure delta, and whether it is identical.
- **Provenance:**
  - `<out>.dosage_repaired.fa` and `<out>.dosage_original.fa`;
  - `<out>.dosage_provenance.tsv`: §9, one row per (chosen repaired haplotype, level, nearest panel sequence), with the
    CIGAR.

**Status keys:**
- `dosage_audit` (off/shadow), `dosage_audit_status` (COMPLETE / REFUSED_LIBRARY), `dosage_policy`
  (`dosage_policy_v2`), `dosage_molecule_span_bp`;
- `dosage_census_sensitive_blocks`, `dosage_m`, `dosage_z_critical`;
- per-status and per-reason counts, and `dosage_repairs_conditional`;
- `dosage_allocation_test` with its blocks tested;
- `seconds_dosage_{census,features,estimate,repair,allocation,provenance}`.

A library without `insert_mean_bp` and `insert_sd_bp` (e.g. a future long-read model) gives
`dosage_audit_status REFUSED_LIBRARY`, and every block is `NOT_EVALUATED` / `unsupported_library`.

## 12. Architecture and configuration

- **Library code.**
  - `include/panvar/genotype_dosage.hpp` with `src/genotype_dosage.cpp`: census, features, reference, estimate,
    classification.
  - `src/genotype_dosage_repair.cpp`: candidates, exact rescoring through the existing scoring API
    (`score_diplotypes`), allocation, provenance.
  - All of it is callable from `panvar genotype`.
- **`panvar genotype --dosage-audit off|shadow`, default off.**
  - Off: the run is byte-identical to today's.
  - Shadow: the audit runs after the call and writes only the §11 files and status keys; every existing output stays
    byte-identical.
  - Repaired calls are **not** integrated into the reported genotype in this phase. Integration is a later,
    separately gated step, never a second required executable.
- **Replay entry (validation only).** `genotype_probe dosage-replay` rebuilds a persisted run's call, guarded by the
  winner-reproduction check, and calls the same library on the persisted reads. It does not re-genotype.
- **Threads.**
  - Hit counting splits read pairs into contiguous ranges, reduced in input order (integer sums, so exact).
  - Bootstrap indices are drawn sequentially; replicates are evaluated in parallel.
  - Scoring uses the existing deterministic scorer.
  - Output must be byte-identical at 1, 2, 4 and 8 threads.
- **Truth isolation.** The library API has no truth parameter, and the replay harness writes a SHA-256 manifest of
  every output before any evaluation script is allowed to run. The evaluation script refuses a directory whose
  manifest does not verify, and it records the manifest digest.

### 12.1 Tie branches without a completed call (targeted gate, `loo.HG01573.30`)

When tie branching is active (v2+ profiles), the procedure is:

1. **Root representatives.** Every tied neighbour at the root improving-move tie, plus the root incumbent.
2. **Evidence per representative.** §3–§6 are computed for each representative, in place of a reported call:
   exclusion and reference from its homologues. This gives its feature sets, reference, estimates and admissible
   sets.
3. **Invariance check.**
   - **Identical across representatives:** the tested features (their k-mer sets), the reference set, the hit
     totals, the quantiles and the admissible sets. Then the constraint is that common admissible system.
   - **Otherwise:** `root_evidence_noninvariant` is reported, and the branch search runs unchanged.
4. **Constraint.** Each branch start, and each neighbour a branch would move to, is kept only if its block genotype
   at every sensitive block is dosage-admissible under the common system.
5. **Cap.** The existing cap of `16` branches is unchanged. The gate records whether the constrained tree completes
   within it.
6. If it does not complete, the fail-closed no-call is kept. A branch is never selected merely to obtain a call.

This mode is a separate switch (`--dosage-tie-constraint`, default off) and is reported only as a gate result until
reviewed.

## 13. Known limitations and rejected alternatives

- **HG01596 sv:1 (fresh reads; a 2 kb segment on one homologue).**
  - In the donor's LOO panel, the 2 kb segment and the neighbouring 4.4 kb region share one copy-count column. The
    only allele that separates them is the donor's own.
  - The pooled feature therefore reads 1.73 at 5×, a non-integer only 2.1 robust sd from the call: silent to the
    dosage check.
  - This is a known limitation. It needs richer evidence (local assembly or long reads).
- **Contiguity-split features, rejected.** They caught HG01596.5, but 13 of 918 campaign and 13 of 301 fresh correct
  calls crossed z_critical(m): small components measure local sequence identity as well as copies. As a generic
  SUSPECT source they were also too nonspecific (166 of 918 and 75 of 301 correct calls).
- **DOSAGE_MODEL_ALIASED, rejected.**
  - Each pooled feature's component matrix has exact rank 1 by construction.
  - The annotation therefore marked 962 of 966 campaign and 326 of 328 fresh block calls, every correct call among
    them.
- **Allocation between homologues is unobservable** from 350 bp-insert reads. It is labelled, not repaired.
- **Also rejected earlier:**
  - likelihood arbitration of repairs (Policy C): its reverts were right 3 times and wrong 4;
  - fixed numeric thresholds, replaced by the family-wise policy;
  - likelihood-only sequence polishing: it overfits.

## 14. Replay, identity and fixture gates (acceptance for the shadow port)

These use existing artifacts only: no genotyping and no campaign rerun.

1. **Classification.**
   - **Coverage:** every cell of the campaign (483 called; the no-call HG01573.30 has no call to audit) and the 164
     fresh cells.
   - **Integer and set quantities, exact:** per-feature k-mer counts, reference size, admissible sets and every
     block status.
   - **Real quantities:** equal to the sealed Python tables at their printed precision (4 decimals). Any decision
     within that precision of its boundary is listed; the smallest observed margins are 0.059 (z_critical), 0.009
     (2.5) and 0.0003 (power).
   - **Molecule span:** derived as 510 from every C4 library, with an identical census.
2. **Counts, content errors:**
   - campaign: 31 AUTO, 5 SUSPECT, 3 UNDERPOWERED, 0 silent;
   - fresh: 19 AUTO, 1 SUSPECT, 2 UNDERPOWERED, 1 silent (HG01596.5);
   - **zero** correct calls AUTO among 918 campaign and 301 fresh block calls.
3. **Repairs.** For every AUTO block, the prototype's results must be reproduced:
   - the admissible pairs or M range, the multisets, the pool, the retained structures, the caps;
   - the exact scores, bitwise, since they come from the same scorer;
   - the top sets, the content choice, the allocation choice and the refusals;
   - `repair_scope` for each repair, with the count of CONDITIONAL repairs reported.
   - Numpy-computed screening quantities (EM weights, proposal log-likelihoods) must agree within 1e-9 relative, and
     every rank and set derived from them must be identical.
   - Any divergence caused by an amendment is listed by cell. For the top-genotype rule, none is expected (§8).
4. **Allocation test** on the self-review subset (`sv7_selfreview/allocation_test.tsv`): identical statuses and
   identical-transposition counts. It also reports the incremental wall time, peak memory and number of blocks
   tested over the full replay.
5. **Byte identity.**
   - `--dosage-audit off` and the flag absent are byte-identical to the 11fe681 outputs on the end-to-end fixtures.
   - In shadow mode every pre-existing file is byte-identical, and the only additions are the §11 files.
6. **Threads.** Every §11 file is byte-identical at 1, 2, 4 and 8 threads.
7. **Provenance.**
   - A fixture with two equally nearest panel sequences reports both, at both levels.
   - Permuting the panel order leaves every provenance row identical.
   - The LOO replay's nearest-sequence search visits exactly the LOO model's panel: the set of compared panel IDs
     equals the model's panel, recorded as a digest. The withheld donor is never among them.
   - The CIGAR re-derives both sequences and the raw distance.
   - Display-ID collision detection is exercised on a forced collision.
8. **Truth isolation.** A test that the audit's sources and replay harness read no evaluation or truth path, and
   that the evaluation step refuses unsealed inputs.
9. **Undefined-estimator fixtures:**
   - `no_tested_feature`, `no_reference`, `no_reference_hits`, `bootstrap_zero_reference`;
   - m = 0 (z_critical NA);
   - robust sd = 0 with the call at and off the median;
   - `unsupported_library`.
10. **Tie fixtures:**
    - distinct genotypes tied (→ `top_genotype_tie`, all reported);
    - duplicate spelled genotypes merged;
    - orientation-only tie;
    - a same-multiset allocation tie, resolved by the prior;
    - a cross-multiset tie, refused.
11. **HG01573.30 targeted gate** (§12.1), reported as:
    - the root representatives and whether their evidence is invariant;
    - completes / does not complete within 16, with the branches required and the admissible starts.
    - Its companion gate permutes the root representatives and requires identical results.

## 15. Quality reporting (prepared, not the final scientific table)

The evaluation is truth-joined after sealing.

**Correction outcome.** For every block call, four facts are reported separately:
1. whether the original call was erroneous;
2. whether the dosage trigger detected it;
3. whether the final call (repaired or kept) is correct;
4. the change in distance to truth after repair (improved / unchanged / worsened, with the bp).

A triggered or repaired original error counts as a successful correction only when fact 3 holds.

**Per individual annotation:**
- `dosage_status` values;
- `allocation_status`;
- `read_support_status` (PRIOR_ONLY);
- `representation_status`;
- the production `quality_filter`.

**Useful combinations:**
- dosage (REPAIRED ∪ TRIGGERED_UNREPAIRED ∪ SUSPECT);
- production filters;
- their union.

**Reported per annotation and combination:**
- error-detection precision and recall for incorrect calls (content errors, and all errors ≥ 1 kb);
- correct-call retention;
- stratified by coverage and by block;
- repaired, warned (SUSPECT), underpowered and unresolved counts, each separately;
- COMPLETE and CONDITIONAL repairs separately.

UNDERPOWERED, PRIOR_ONLY and ALLOCATION_UNRESOLVED describe information limits. They are **never** counted as false
positives when the genotype is correct. Two views are reported:
1. error detection, using only the annotations that predict error (REPAIRED, TRIGGERED_UNREPAIRED, SUSPECT,
   representation conflict);
2. high-confidence retention: the correct calls carrying none of the information-limit or error annotations.

## 16. Review record

**2026-09-28: Sol approved** the specification and its five design choices:
- the extra statuses;
- GQ unavailable on repair;
- the β reference generalised as the adjacent context blocks;
- the expanded allocation test;
- §12.1 as a gate only.

The ten amendments are listed in §17.

## 17. Amendments (Sol, 2026-09-28)

1. **Exact top ties** (§8, top-genotype rule): IDs never choose between distinct tied genotypes, and a
   cross-multiset tie refuses.
2. **Candidate completeness:** `repair_scope` COMPLETE / CONDITIONAL, with reasons (§8).
3. **HG01573 root-tie invariance:** evidence is computed for every root representative, and the constraint applies
   only if it is invariant (§12.1).
4. **Undefined estimator cases** fail closed (§3, §4, §5).
5. **Molecule span** derived from the library; short-read only (§2).
6. **Provenance** only for the chosen haplotypes; normalised by alignment length; full SHA-256 IDs (§9).
7. **Original GQ** retained on repair (§10).
8. **Evaluation:** the correction outcome is reported as four separate facts (§15).
9. **Reasons list** extended (§11).
10. **Allocation test** shadow-only, with measured cost (§8.3).

## 18. Implementation record (shadow port, 2026-09-28/29; for review)

**Code.** Branch `genotype-dosage`, worktree `panvar-worktrees/dosage-dev`.
- **Library:** `genotype_dosage.{hpp,cpp}` (census, features, estimate, policy v2), `genotype_dosage_repair.cpp`
  (routes, choice, allocation test), `genotype_dosage_provenance.cpp`, `genotype_dosage_shadow.cpp` (the one entry
  point).
- **`panvar genotype --dosage-audit shadow [--dosage-allocation-test]`** calls the shadow after every existing output
  is committed.
- **Validation entry:** `genotype_probe dosage-replay` calls the same function on a persisted run.
- **Fixtures:** `tests/genotype_dosage.cpp` (ctest `panvar_genotype_dosage`).
- **Gates:** `tests/genotype_dosage_shadow_identity.sh` and `tests/genotype_dosage_thread_identity.sh`.

**Findings that bind the contract:**
1. **Orientation on the repeated-unit route.** The prototype scored structures in one orientation, and §8 now says
   so. The pair route scores both.
2. **Exact allocation ties.** In 35 of 50 LOO unit repairs (campaign 18 of 31, fresh 17 of 19), several distinct
   genotypes of the chosen content tie exactly on both panel prior and likelihood at the allocation pick. §8.2's last
   key, structure rank, then chooses among them.
   - Every such block reports `allocation_tie` = the number of tied genotypes, and flags each tied candidate.
   - Post hoc, the tied alternatives differ from truth by 6–140 bp, and none changes solved versus unsolved (< 1 kb).
   - The prototype broke the tie by label string instead. That differs from structure rank in one block
     (`loo.NA20847.30`: prototype 108 bp, structure rank 26 bp).
   - **Open for review:** keep structure rank; choose by sequence (smallest SHA-256 pair), which makes the choice
     independent of IDs; or report such blocks as allocation-unresolved representatives.
3. **Numerics.** The dosage sources are compiled without floating-point contraction (`-ffp-contract=off`), so the
   proposal likelihoods equal numpy's bit for bit. Two multisets with mathematically equal likelihoods otherwise
   swapped ranks in `loo.HG04217.2`. Only the normal quantile contracts, because CPython builds it with contraction.
4. **Outputs.**
   - Status keys go to `<prefix>.dosage_status.tsv`, never into `mosaic_status.tsv`, so shadow mode leaves every
     existing file byte-identical.
   - `<prefix>.dosage_feature_copies.tsv` gives each tested feature's per-allele copies and admissible set.
   - The audit runs on the equivalence set's rank-1 representative. A run without one (a no-call) gets
     `dosage_audit_status NOT_RUN`.
5. **Identity normalisation.** The same rule as the tie-branch gate: status keys for wall time, memory, worker gauges
   and scheduling-dependent cache counters are excluded; repair rounds are compared on columns 1–28; the FASTA header
   path is replaced.
6. **HG01573.30 (§12.1).** Decided from the frozen run without changing the branch search.
   - The root evidence is identical across its three representatives and all six of their orderings.
   - Every tie is a phase tie, and all 18 starts the early cap stop required are dosage-admissible.
   - The constraint therefore removes nothing: at least 18 branches are still needed against the cap of 16, and the
     no-call stands.

## 19. Second review record (Sol, 2026-09-28): allocation semantics, scoring lifecycle, candidate space

**Decisions.**
- Allocation ties use the 1e-9-nat tolerance, and every tied allocation is reported.
- The allocation observability test stays off by default.
- CONDITIONAL repeated-unit repairs never replace production calls, and repaired genotypes are not integrated
  (`--dosage-audit shadow` is the only enabled mode).
- HG01573.30 stays a fail-closed no-call.
- Nearest-panel provenance is accepted.

**Allocation semantics (normative; replaces the last key of §8.2).**
1. **Equivalence set.** Within the chosen content multiset, take the maximal panel-structure prior (exact). Every
   arrangement within 1e-9 nats of the best likelihood among those is one **allocation-equivalence set**. It is
   reported in full, one row per distinct unordered spelled genotype, in `<prefix>.dosage_allocation_equivalence.tsv`.
2. **Representative.** The member with the smallest canonical SHA-256 pair is the **serialization representative**
   only. It is labelled `REPRESENTATIVE_ONLY_NOT_RESOLVED` whenever the set has more than one member, and it is never
   described as the called diplotype.
   - Repaired FASTA headers carry `repair_scope`, `allocation_status` and, for such a set,
     `allocation=REPRESENTATIVE_ONLY` with a link to the equivalence set.
3. **allocation_status.**
   - `UNRESOLVED` when the observability test proves it, or when the set has more than one member (basis
     `tied_allocation_equivalence_set`).
   - `NOT_TESTED` when an eligible repeated-unit genotype (≥ 3 units) was not tested.
   - `NOT_APPLICABLE` only for the pair route or fewer than three units.
4. **Declared corrections, not replay failures.** The prototype's choice in this step depended on label strings and on
   pandas' float parser. Two cells showed it:
   - `loo.NA20847.30`: an exact tie broken by label string;
   - fresh `loo.NA19384.5`: likelihoods 1.5e-11 nats apart, merged by the parser.
   - Across all 54 sealed repairs the approved rule leaves the content multiset unchanged and always contains the
     prototype's pick.
   - The representative differs from that pick in 18 repairs, a serialization choice only.
5. **Scope.** `repair_scope` (COMPLETE for the pair route, CONDITIONAL otherwise) is carried in every table and in
   every FASTA header.

**Scoring lifecycle (a lifecycle change only; scores bit-identical).**
- **One context per repair.** Each repair is scored in one `PreparedScoringContext`:
  - the reported homologues are placed once with retained censuses;
  - the reproduction guard is scored first;
  - each distinct repaired haplotype is placed once by the exact block delta and reused by every diplotype holding
    it;
  - columns are released once no pending candidate needs them, and released slots are reused.
- **Caps.** The profile's repair caps (`repair_max_*`) are cumulative over the repair. A cap refuses the repair with
  no partial output.
- **Origin cache.** The stage-6 origin path is on with retention budget 0: nothing is retained, every window is
  recomputed exactly, and the path's merge is kept. With the profile's 768 MiB budget, NA19434.30 peaked at 3.74 GiB
  for no hits.
- **Construction.** Compact fixed-width multiset records, a sparse EM over nonzero entries in the dense index order,
  and parallel proposal likelihoods into per-candidate slots. All three leave every value bit-identical. The earlier
  one-shot path stays selectable (`--repair-scoring cold`) as the reference.

**Candidate space (exact counts, nothing scored).**
- **Pool truncation.** The 36-unit pools omit 71–75 of 107–111 internal-capable and 69–70 of 105–106
  terminal-capable units in every repair.
- **Full-library sizes.** Admissible multisets number 1.0e5–1.5e10, and the theoretical diploid structures
  2.0e5–3.3e11. Only 123–1,128 structures are constructed.
- **Proposal ranking is not a bound.** The chosen content's proposal rank is 1 in only 11 of 50 repairs, and up to 61.
- **Dosage content is determined by the constraint alone** in 48 of 50 repairs: every admissible candidate of the
  full library has the same rounded feature-copy vector and the same (units, I-state units) count.
- **Junction-signature classes are exact.** Classes of equal content and equal unit-junction multiset agree in
  likelihood within 6e-11 nats, but they collapse the constructed structures by only 1%.

## 20. Integration contract: `panvar genotype --dosage-mode off|audit|repair` (written 2026-09-28; repair not enabled)

**Scope.** One genotyping command; no second executable.
- The dosage stage runs inside `panvar genotype` after the call is decided.
- It uses the same scoring context settings, the profile's repair caps and the same model.
- `--dosage-audit shadow`, the validation mode of §12, becomes `--dosage-mode audit` with identical outputs.
- **Default `off`.**
- **`repair` is specified here and not enabled.** The command refuses it until reviewed.

### 20.1 Modes

| mode | dosage work | reported genotype | new files |
|---|---|---|---|
| `off` | none | the caller's | none; every output byte-identical to the build without the stage |
| `audit` | census, features, classification, repair search, provenance | the caller's, unchanged | `<prefix>.dosage_*` (§11, §19) |
| `repair` | as audit | changed only by COMPLETE_DIPLOTYPE_REPAIR blocks (20.2) | as audit, plus `<prefix>.original_*` (20.5) |

### 20.2 What each repair class does in `repair` mode

| repair class | block genotype written | haplotype FASTA | block row `dosage_status` |
|---|---|---|---|
| (not AUTO_REPAIR) | the caller's | the caller's | CONSISTENT / SUSPECT / UNDERPOWERED / NOT_SENSITIVE / NOT_EVALUATED |
| COMPLETE_DIPLOTYPE_REPAIR | the repaired pair replaces the block's alleles (pair route: panel allele IDs) | the repaired full haplotypes, headers `role=RESOLVED` | REPAIRED |
| CONDITIONAL_DOSAGE_CONTENT | **the caller's**, kept, with the corrected content reported beside it (20.3) | the caller's; the representative only in `<prefix>.dosage_repaired.fa`, headers `role=REPRESENTATIVE_ONLY_NOT_RESOLVED` | CONTENT_REPAIRED |
| TRIGGERED_UNREPAIRED | the caller's | the caller's | TRIGGERED_UNREPAIRED (with reason) |

A conditional representative is never written as the block's `allele1`/`allele2`, never enters the main haplotype
FASTA, and never raises a call class. The exact allele sequences and allocation stay non-authoritative.

### 20.3 Block tables (`genotypes.tsv`, `block_genotypes.tsv`)

The existing columns keep their meaning. New trailing columns, so positional readers of the old columns are
unaffected:
- `dosage_status`, `dosage_reason`, `repair_class`;
- `dosage_content_scope`, `diplotype_scope`, `allocation_equivalence_scope`, `serialization_representative`;
- `original_allele1`, `original_allele2`: always the caller's, even after a COMPLETE_DIPLOTYPE_REPAIR;
- `original_content`, `corrected_content`, `admissible_content_vectors`;
- `allocation_status`;
- `original_GQ`, `gq_status`: UNCHANGED_CONDITIONAL, or UNAVAILABLE_REPAIRED with `conditional_GQ` NA;
- `repair_separation`, `conditional_repair_gap` (diagnostics).

**For a CONDITIONAL_DOSAGE_CONTENT block:**
- `allele1`/`allele2` are the caller's;
- `corrected_content` is authoritative within the declared model (units, I-state units, feature copy vector);
- `conditional_GQ` stays the caller's value, labelled as describing the original genotype;
- the representative's identity is in `dosage_blocks.tsv` and in the equivalence set only.

### 20.4 FILTER, quality and call class

- **`filter` (representation) and `quality_filter` (read support)** are unchanged and computed as today. Dosage never
  rewrites them.
- **A new per-block `dosage_filter`:**
  - `PASS`: CONSISTENT, NOT_SENSITIVE, or REPAIRED;
  - `DOSAGE_CONTENT_CORRECTED`: CONTENT_REPAIRED, i.e. the original alleles contradict the complete dosage content;
  - `DOSAGE_TRIGGERED_UNREPAIRED`;
  - `DOSAGE_SUSPECT`, `DOSAGE_UNDERPOWERED`, `DOSAGE_NOT_EVALUATED`;
  - information limits and warnings are listed as reasons, never merged into LOW_QUALITY.
- **The run's call class** (`strict_call`, `genotype_equivalent_call`, `no_call`) is unchanged by audit. In repair
  mode only COMPLETE_DIPLOTYPE_REPAIR changes a genotype.
  - A strict call whose repaired block came from a pair-route complete repair stays `strict_call`, with
    `dosage_repaired_blocks` listed.
  - A call with any CONDITIONAL_DOSAGE_CONTENT block keeps its class. The status records
    `dosage_content_corrected_blocks`, and that block's `dosage_filter` is DOSAGE_CONTENT_CORRECTED.

### 20.5 Status file and recoverability

- **`mosaic_status.tsv` in repair mode** gains:
  - `dosage_mode`, `dosage_policy`, `dosage_audit_status`;
  - per-class counts;
  - `dosage_repaired_blocks`, `dosage_content_corrected_blocks`;
  - the molecule span and m.
- **In audit mode these stay in `<prefix>.dosage_status.tsv` only**, so audit remains byte-identical in every existing
  file.
- **Original call and sequences stay recoverable:**
  - `<prefix>.original_genotypes.tsv`, `<prefix>.original_block_genotypes.tsv` and `<prefix>.original_haplotypes.fa`
    are written whenever repair mode changes any block;
  - the equivalence set and `<prefix>.dosage_original.fa` hold the pre-repair homologues;
  - every repaired sequence is identified by its full SHA-256.

### 20.6 Gates required before `repair` is enabled

- **Identity:** `off` and `audit` byte-identical as today. `repair` identical to `audit` in every block that is not
  COMPLETE_DIPLOTYPE_REPAIR.
- **Round trip:** the original files reproduce the audit-mode outputs exactly.
- **Threads:** 1/2/4/8-thread identity of every file.
- **Downstream:** a schema gate over the new trailing columns.
- **Evaluation:** the frozen evaluator scores `original_*` and the repaired outputs separately.
