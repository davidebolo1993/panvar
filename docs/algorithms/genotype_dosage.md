# Dosage audit and repair: normative specification (draft for review, 2026-09-28)

This document is the contract for the C++ port of the dosage audit and repair add-on to `panvar genotype`. It freezes
the design Sol accepted on 2026-09-28. The Python prototype is the executable reference:
`tests/results/genotype_c4_production_11fe681_dosage_repair/`.
- **The sealed reference outputs** are the policy v2 tables `contiguity/class_baseline_{campaign,fresh}.*` and the
  repairs in `phase1/` and `policyC_test/`.
- **Normative words.** "Must" and "never" are requirements. A value in `code` is a frozen constant, and changing it is
  a new design, not a tuning step.

The add-on is truth-free. It reads the prepared model, the call the caller reports, the reads and the library, and
nothing else.

## 1. Scope and states

For every block of a genotyped locus, the audit answers one question: do the reads contradict the copy numbers the
reported genotype implies?

When they contradict it beyond a family-wise error rate, the block is repaired inside the dosage the reads allow.
Otherwise the reported genotype is kept and annotated.

Repair never happens for any other reason. In particular:
- a correct-looking likelihood gain, a weak dosage deviation, low read support, representation conflict or an
  unresolved allocation never trigger it.

Blocks that the census (§2) does not mark as copy-number-sensitive are never examined: `dosage_status NOT_SENSITIVE`.

## 2. Census: copy-number-sensitive blocks (graph only, once per model)

For each block, over its panel alleles as the prepared model holds them (each allele's node path; node weight = node
length in bp):

- **len(a):** the allele's length, i.e. the sum of node weights with multiplicity.
- **novel(a, b):** the bp of distinct nodes present in `a` and absent from `b`.
- **Dosage-only pair (a, b):** len(b) − len(a) ≥ `510` bp, and novel(a, b) ≤ `50` bp. The shorter allele differs
  from the longer one essentially by missing sequence, which is what copy number measures.
- **Copy-number-sensitive:** the block has at least one dosage-only pair.
- **Repeated-unit route** (§8.2): the block is copy-number-sensitive, some allele visits a node at least twice, and
  that allele's repeated sequence exceeds `510` bp. The repeated sequence is the bp of nodes visited ≥ 2 times in one
  allele. Otherwise a sensitive block uses the **pair route** (§8.1).
- **Census-invariant:** every other block. These blocks provide the single-copy reference (§4).

The census uses no reads, calls or truth. C4 yields sv:1 (pair route) and sv:7 (repeated-unit route) in all 66
models. `510` is the maximum fragment span of the library (insert 350 + 4 × 40). `50` is the minimum variant size.

## 3. Pooled dosage features (graph, then per call)

- **k-mers.** Canonical `31`-mers of every panel allele of the block, counted with multiplicity: A[a, x] = number of
  occurrences of k-mer x in allele a.
- **Varying k-mers:** the column A[·, x] is not constant across alleles.
- **Groups:** varying k-mers with an identical column. A group with ≥ `500` k-mers is a **pooled feature** with copy
  vector s = that column. Features are ordered by the first occurrence of any of their k-mers, scanning alleles in
  ascending allele ID and positions left to right.
- **Per-call exclusion.** For the reported genotype, remove from every feature the k-mers occurring anywhere in
  either reported homologue outside this block's interval. A feature left with ≥ `500` k-mers is **tested**.
- **m** for a sample = the number of tested features summed over all copy-number-sensitive blocks of the locus.
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

If R is empty, or the reads hit R zero times, the block is `NOT_EVALUATED` (reason `no_reference`). That is a fail
closed state, never CONSISTENT.

## 5. Estimate and bootstrap

- **Hits.** Every canonical k-mer of both mates of every read pair, with multiplicity. Each is assigned to R (index
  0) or to the tested feature holding it. A read pair contributes a count vector v = (v_0, v_1, …).
- **Retained pairs.** Pairs whose vector is all zero are dropped. The retained pairs, in input order, are the
  bootstrap population P of size n.
- **Point estimate.** T_j = Σ v_j over P. Then CN_j = 2 · (T_j / |F_j|) / (T_0 / |R|), evaluated in this order in
  IEEE double.
- **Bootstrap.**
  - `200` replicates; each draws n indices uniformly with replacement and recomputes CN_j by the same formula.
  - The generator is the Python `random.Random(20260925)` stream (MT19937, `init_by_array` seeding, `randrange` by
    `getrandbits` rejection). It is re-seeded for every block and drawn sequentially: replicate r, then index i.
  - This is frozen so the port reproduces the prototype bit for bit. Draw order is independent of thread count.
- **Quantiles.** q(p) = sorted replicates[⌊p·(200 − 1)⌋] for p ∈ {0.05, 0.25, 0.50, 0.75, 0.95}.
- **Robust sd** = (q(0.75) − q(0.25)) / 1.349. **Robust z** = |c_j − q(0.50)| / robust sd. If robust sd = 0, then
  z = 0 when c_j = q(0.50), else +∞.
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

**Exact score.** Every candidate is spliced into the reported homologues at the block's intervals, in both
orientations (one if the two block sequences are equal). It is scored by the production likelihood: fragment term +
exposure term, with the run's own scoring settings. The **likelihood** is that sum. The prior is recorded but never
ranks.

**Winner guard.** Before any candidate is scored, the reported diplotype alone must reproduce its persisted total,
fragment term and exposure term within 1e-6 in the same scoring context. Otherwise the repair fails closed:
`TRIGGERED_UNREPAIRED`, reason `winner_not_reproduced`.

### 8.1 Pair route

- **Candidates:** every unordered pair {a, b} (a = b allowed) of the block's panel alleles with s_j[a] + s_j[b]
  admissible for every tested feature. Cap: `4096` pairs; above it, `TRIGGERED_UNREPAIRED`, reason `cap_pairs`. None:
  `TRIGGERED_UNREPAIRED`, reason `no_admissible_candidate`. In both cases the reported genotype is kept.
- **Choice:** maximal likelihood. Ties go to (smaller allele ID pair, then orientation AB before BA). An exact tie at
  the top is recorded (`repair_top_ties`).
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

**Screening (count model).** This matches the prototype, which is itself documented in `phase1/proposal_general.py`.
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

**Construction.** Multisets are taken in rank order and expanded into arrangements (below). Caps: `64` multisets,
`256` haploid sequences, `4096` diploid structures.
- Construction stops at the first cap.
- A cap reached at the first multiset, or at an exact likelihood tie across the boundary, refuses:
  `TRIGGERED_UNREPAIRED`, reason `cap_refused`.
- Arrangements are both terminal assignments × every split of the internal units between homologues × every order
  within a homologue. Each homologue holds ≤ KMAX units, and each homologue's first unit is anchor-compatible.
  Structures with identical spelled sequences are merged.

**Choice.**
1. **Content:** the multiset of the maximal-likelihood exactly scored structure (ties: structure rank, then
   orientation).
2. **Allocation, the split between homologues:** among the scored structures of that multiset, maximise the **panel
   structure log-prior** log((f(σ_a) + 0.5)/N) + log((f(σ_b) + 0.5)/N).
   - σ is a homologue's unit-state string (e.g. `IO`), f(σ) the number of panel alleles with that state string, and
     N the number of panel alleles.
   - Ties go to likelihood, then structure rank, then orientation.
   - Reason: the transposition test shows that 350 bp-insert fragments cannot see the split. The likelihood and the
     original call's split are therefore uninformative about it. Of the rules compared, panel frequency was best on
     both datasets.
3. The reported block genotype is that structure in its best-scoring orientation.

### 8.3 Allocation observability test (annotation only)

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

## 9. Repaired-sequence provenance

This covers every haplotype the repair constructs that differs from both reported homologues. Comparisons use only
the panel supplied to the caller (in LOO validation, the LOO panel).

- the repaired and the original haplotype sequences (FASTA; the original is the reported homologue it replaces);
- **Nearest panel sequences at two levels:**
  - the full haplotype against every spelled panel haplotype;
  - the repaired block allele against every panel allele of the block.
- **Distance:** global edit distance (unit costs; edlib, the existing submodule). Raw, and normalised by the
  repaired sequence's length.
- **Ties:** every panel sequence at the minimal distance is reported. They are ordered by sequence SHA-256, then by
  panel ID, so the output is invariant to panel ordering.
- the extended CIGAR (`=`, `X`, `I`, `D`) against each reported nearest sequence;
- **Unit provenance:** each unit in order, with its unit ID, its sequence digest and every panel allele and position
  that carries it (repeated-unit route); or the two panel allele IDs (pair route).

## 10. Genotype quality after repair

- **Unrepaired blocks:** the production `conditional_GQ` keeps its production meaning. It is conditional on the
  scored candidate set and uncalibrated: `gq_status UNCHANGED_CONDITIONAL`.
- **REPAIRED blocks:** the ordinary GQ is **not** reused. The repaired genotype was chosen inside a different,
  dosage-constrained candidate set, and a repaired non-panel haplotype may have no spellable atomic prior. So:
  - `gq_status UNAVAILABLE_REPAIRED` and `repaired_GQ NA`;
  - `conditional_repair_gap` = likelihood(chosen) − likelihood(best admissible candidate with a different unordered
    block genotype), in nats; NA if there is none;
  - `repair_separation` = likelihood(chosen) − likelihood(the original call), in nats. It can be negative: the
    likelihood alone is not a reliable judge of dosage, which is why repair is dosage-constrained.

## 11. Output schema (shadow; separate fields, never collapsed)

Existing outputs are unchanged. The audit writes new files next to them, and each row carries the fields below.
`quality_filter` is not recomputed or merged.

**`<out>.dosage_blocks.tsv`**, one row per block:

| field | values |
|---|---|
| block_id, census_route | pair / repeated_unit / not_sensitive |
| dosage_status | CONSISTENT, REPAIRED, TRIGGERED_UNREPAIRED, SUSPECT, UNDERPOWERED, NOT_SENSITIVE, NOT_EVALUATED |
| dosage_reason | none, no_admissible_candidate, cap_pairs, cap_refused, winner_not_reproduced, no_reference |
| allocation_status | RESOLVED, UNRESOLVED, INCONCLUSIVE, NOT_APPLICABLE |
| read_support_status | the call's `genotype_support`, copied |
| representation_status | the call's representation `filter`, copied |
| alpha, m, z_critical, features_tested, max_robust_z | |
| original_allele1, original_allele2 | the reported block alleles |
| repaired_allele1, repaired_allele2 | NA unless REPAIRED; generated alleles get content IDs `r<sha256 prefix>` |
| repair_route, candidate_source | pair / repeated_unit; panel_pair / unit_multiset |
| candidates_constructed, candidates_scored, repair_top_ties, cap_reached | |
| original_fragment_term, original_exposure_term, original_likelihood, original_prior, original_total | |
| repaired_fragment_term, repaired_exposure_term, repaired_likelihood, repaired_prior, repaired_prior_spellable | |
| repair_separation, conditional_repair_gap | §10 |
| gq_status, reported_GQ | UNCHANGED_CONDITIONAL with the production value, or UNAVAILABLE_REPAIRED with NA |

**`<out>.dosage_features.tsv`**, one row per tested feature: the §7 fields.

**`<out>.dosage_repair_candidates.tsv`**, one row per scored candidate orientation:
- the route, label, rank, the unit tuples or allele pair, the state strings;
- the proposal log-likelihood (repeated-unit route), fragment term, exposure term, likelihood, prior;
- the panel-structure log-prior;
- whether it was chosen for content and for allocation.

**`<out>.dosage_allocation.tsv`**, one row per transposition: A, A′, the fragment-vector hashes (bitwise and 1e-9
grid), the exposure delta, and whether it is identical.

**Provenance files:**
- `<out>.dosage_repaired.fa` and `<out>.dosage_original.fa`;
- `<out>.dosage_provenance.tsv`: §9, one row per (repaired haplotype, level, nearest panel sequence), with the CIGAR.

**Status keys:** `dosage_audit` (off/shadow), `dosage_policy` (`dosage_policy_v2`), `dosage_census_sensitive_blocks`,
`dosage_m`, `dosage_z_critical`, per-status counts, and `seconds_dosage_{census,features,estimate,repair,allocation}`.

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
  - Integration into the reported genotype (replacing the call) is a later, separately gated step, never a second
    required executable.
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

When tie branching is active (v2+ profiles):
- **Evidence without a call.** Dosage evidence is computed from the **root tie's incumbent**: the diplotype at which
  the repair stopped, in place of a reported call. This uses §3–§6 unchanged, including exclusion and reference from
  the incumbent's homologues.
- **Constraint on branches.** Each branch start, and each neighbour a branch would move to, is kept only if its block
  genotype at every sensitive block is dosage-admissible (§8, every tested feature).
- **Cap.** The existing cap of `16` branches is unchanged. The gate records whether the admissible branch tree
  completes within it.
- If it does not complete, the fail-closed no-call is kept. A branch is never selected merely to obtain a call.
- This mode is a separate switch (`--dosage-tie-constraint`, default off) and is reported only as the gate result
  until reviewed.

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

## 14. Replay and identity gates (acceptance for the port)

These use existing artifacts only: no genotyping and no campaign rerun.

1. **Classification.**
   - **Coverage:** every cell of the campaign (483 called; the no-call HG01573.30 has no call to audit) and the 164
     fresh cells.
   - **Integer and set quantities, exact:** per-feature k-mer counts, reference size, admissible sets and every
     block status.
   - **Real quantities:** equal to the sealed Python tables at their printed precision (4 decimals). Any decision
     within that precision of its boundary is listed; the smallest observed margins are 0.059 (z_critical), 0.009
     (2.5) and 0.0003 (power).
2. **Counts, content errors:**
   - campaign: 31 AUTO, 5 SUSPECT, 3 UNDERPOWERED, 0 silent;
   - fresh: 19 AUTO, 1 SUSPECT, 2 UNDERPOWERED, 1 silent (HG01596.5);
   - **zero** correct calls AUTO among 918 campaign and 301 fresh block calls.
3. **Repairs.** For every AUTO block, the prototype's results must be reproduced:
   - the admissible pairs or M range, the multisets, the pool, the retained structures, the caps;
   - the exact scores, bitwise, since they come from the same scorer;
   - the content choice, the allocation choice and the refusals.
   - Numpy-computed screening quantities (EM weights, proposal log-likelihoods) must agree within 1e-9 relative, and
     every rank and set derived from them must be identical.
4. **Allocation test** on the self-review subset (`sv7_selfreview/allocation_test.tsv`): identical statuses and
   identical-transposition counts.
5. **Byte identity.**
   - `--dosage-audit off` and the flag absent are byte-identical to the 11fe681 outputs on the end-to-end fixtures.
   - In shadow mode every pre-existing file is byte-identical, and the only additions are the §11 files.
6. **Threads.** Every §11 file is byte-identical at 1, 2, 4 and 8 threads.
7. **Provenance.**
   - A fixture with two equally nearest panel sequences reports both.
   - Permuting the panel order leaves every provenance row identical.
   - The LOO replay never references a withheld allele.
8. **Truth isolation.** A test that the audit's sources and replay harness read no evaluation or truth path, and
   that the evaluation step refuses unsealed inputs.
9. **HG01573.30 targeted gate** (§12.1), reported as: completes / does not complete within 16, the branches
   required, and the admissible starts.

## 15. Quality reporting (prepared, not the final scientific table)

The evaluation is truth-joined after sealing.

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
- repaired, warned (SUSPECT), underpowered and unresolved counts, each separately.

UNDERPOWERED, PRIOR_ONLY and ALLOCATION_UNRESOLVED describe information limits. They are **never** counted as false
positives when the genotype is correct. Two views are reported:
1. error detection, using only the annotations that predict error (REPAIRED, TRIGGERED_UNREPAIRED, SUSPECT,
   representation conflict);
2. high-confidence retention: the correct calls carrying none of the information-limit or error annotations.

## 16. Decisions requested at review

1. **Two extra dosage_status values beyond the four listed:**
   - `TRIGGERED_UNREPAIRED`: a trigger with no admissible candidate, a cap or a failed guard (the prototype's
     NA24385.5 sv:1 case);
   - `NOT_EVALUATED`: no reference.
   - With `NOT_SENSITIVE` for census-invariant blocks. None of the three can be expressed as one of the four listed
     statuses without losing meaning.
2. **Tie-breaks.** §8.1 and §8.2 fix numeric tie-breaks where the prototype sorted string labels. They matter only at
   an exact top tie; the replay lists any such case.
3. **β reference.** §8.2 generalises it as the adjacent context blocks, with no named block IDs. For C4 it reproduces
   the prototype exactly.
4. **Allocation test scope.** It covers every reported repeated-unit genotype with ≥ 3 units, not only repaired
   ones. That is extra exact scoring, measured before any adoption.
5. **§12.1 is a gate only.** The tie constraint changes calls, so it stays a separate, default-off switch until
   reviewed.
