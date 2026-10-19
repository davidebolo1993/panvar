# M3-E: separate genotype-data support from phase support (experimental, default off)

**Ancestry.** Branch `m3-e-support-semantics` from `49e211c` (M3-D).

**Scope.** There is no broader joint-search validation, campaign, cap increase, shortlist change or default promotion.
No change to candidate selection, likelihoods, priors, GQ, representation checks, dosage assessment or genotype-filter
membership.

## 1. Evidence contract (frozen)

Three quantities per evaluated candidate:

| quantity | definition |
|---|---|
| **fragment term** F | placement and emission contribution of the fragments (`fragment_term`) |
| **full data likelihood** L = F + E | fragment plus exposure, excluding the copying prior (`log_likelihood`) |
| **final ranking score** S = L + P | full data likelihood plus the existing sequence-level prior P; this alone ranks and selects |

**Candidate unit.** A sequence-level group: the members spelling one unordered sequence pair. A group's genotype at a
block is the **set** of its members' unordered allele pairs there. A group whose members disagree at a block is
projection-ambiguous there, and the ambiguity is **never resolved by choosing a representative**.

**Numerical uncertainty.** Each group carries its certification bound u (the scorer's `numeric_error_bound`; the prior
is exact).
- X is certifiably above Y when L_X − u_X > L_Y + u_Y.
- Rounded output tables are enough for attribution, not for certifying a near tie.

### Experimental block-genotype support (`full_data_likelihood_v1_experimental`)

Per block b, with the called group w:

1. **Projection.** If w is projection-ambiguous at b, the label is **AMBIGUOUS_PROJECTION** and nothing further is
   claimed.
2. **Groups.** Let g be w's unordered genotype at b.
   - **Same** groups: every member carries g at b; w is one.
   - **Disagreeing** groups: no member carries g at b.
   - **Mixed** groups: members both carrying g and not carrying it.
   - Groups that carry g but differ elsewhere, or only in phase, are **same** at b and never weaken b's support.
3. **Comparison.** s is the same group with the largest L; d is the disagreeing group with the largest L.
   - No disagreeing group: **NO_EVALUATED_ALTERNATIVE**. This is vacuous: support only within what was evaluated.
   - L_s certifiably above L_d, and every mixed group certifiably below L_s: **DATA_SUPPORTED**.
   - L_d certifiably above L_s: **PRIOR_DEPENDENT**. The data prefer another genotype at b, and the copying prior selects
     the call.
   - Otherwise, if some mixed group is not certifiably below L_s and d is not above: **AMBIGUOUS_COMPETITOR**.
   - Otherwise: **DATA_TIED**.
4. **Reported per block.**
   - The label and the margin L_s − L_d.
   - Whether it is certified.
   - The strongest disagreeing group, d, as call-minus-competitor differences in F, E, L and P.
   - The counts of disagreeing and mixed groups.

**Scope of the label.** It is support conditional on the evaluated candidate set: not a calibrated probability of
correctness, not a confidence level, and not a phasing claim.

### Phase (kept separate; the label is unchanged)

`phase_support` keeps its existing rule and labels: fragment contenders against the winner; `READS_WITHIN_CANDIDATES`,
`PRIOR_ONLY`, `UNRESOLVED`. The experiment only **classifies** the existing fragment contenders, reported in the status:
- **genotype-changing:** every member's unordered genotype differs from the call somewhere;
- **genotype-preserving:** every member carries the call's genotype at every block, so it differs in phase or
  arrangement only;
- **mixed.**

Exposure-based structural support never becomes a molecular-phasing claim. Equal aggregate likelihoods are not read as
identical per-fragment evidence or as biological unidentifiability.

### Output labels and quality (frozen)

**Option.** Profile key `block_support_rule`: `legacy` (the default) or `full_data_likelihood_v1_experimental`.
`settings.tsv` records it only when a profile sets it.

**Unchanged in both modes:** `genotype_support` (READS / PRIOR_ONLY, legacy fragment-only), `phase_support`,
`quality_filter`, `quality_reason`, GQ, posterior, `numerically_separated` and every filter.

**Appended to `.block_genotypes.tsv` when experimental:**
- `genotype_support_v1x`
- `quality_filter_v1x`, `quality_reason_v1x`
- `data_support_margin_v1x`, `data_support_certified_v1x`
- `competitor_group_v1x`
- `competitor_delta_fragment_v1x`, `competitor_delta_exposure_v1x`, `competitor_delta_likelihood_v1x`,
  `competitor_delta_prior_v1x`
- `disagreeing_groups_v1x`, `mixed_groups_v1x`

**`quality_filter_v1x`.** The legacy combination with only the genotype-support term replaced.

| `genotype_support_v1x` | genotype-support term in `quality_reason_v1x` | effect |
|---|---|---|
| DATA_SUPPORTED | none | may PASS |
| NO_EVALUATED_ALTERNATIVE | none | may PASS |
| DATA_TIED | `data_tied` | LOW_QUALITY |
| PRIOR_DEPENDENT | `prior_dependent` | LOW_QUALITY |
| AMBIGUOUS_PROJECTION | `ambiguous_projection` | LOW_QUALITY |
| AMBIGUOUS_COMPETITOR | `ambiguous_competitor` | LOW_QUALITY |

The representation and library terms are unchanged. Nothing is called "prior only", "high confidence" or calibrated.

**Status rows, when experimental:**
- `block_support_rule`
- the label counts
- `phase_competitors_genotype_changing`, `phase_competitors_genotype_preserving`, `phase_competitors_mixed`

Block tables written by the consensus or plain fallbacks carry `NOT_EVALUATED` in the experimental columns.

**Implementation.** One library function, `block_data_support_v1` (`include/panvar/genotype_block_support.hpp`), is
called by both the pooled (`decide_pooled_genotype` inputs) and the unpooled summary path. Legacy code is untouched.

## 2. Gates (`tests/genotype_block_support.cpp`, synthetic groups through the shared function)

1. A rival wins on the fragment term but loses on L because of exposure: DATA_SUPPORTED, where legacy would say
   PRIOR_ONLY.
2. A rival wins on L, but the prior selects the call: PRIOR_DEPENDENT.
3. Same block genotypes with competing phase: DATA_SUPPORTED at every block, with the phase competitor classed
   genotype-preserving.
4. A rival disagreeing at one block only: only that block is affected.
5. Equivalent sequences (one group), and a projection-ambiguous winner (AMBIGUOUS_PROJECTION) or competitor
   (AMBIGUOUS_COMPETITOR).
6. Numerical ties and bounds: DATA_TIED inside the bounds, certified outside.
7. Group order invariance (permuted input) and thread invariance (pure function; the caller's 1- and 4-thread runs are
   compared in the integration check).

## 3. Evidence audit, then minimum integration

1. **`audit.saved`, before any rerun.** An independent Python implementation of section 1 runs on the saved M3-D tables
   (`run.<cell>.{off,on}/call.mosaic_scores.tsv` + `mosaic_candidates.tsv`). Per block it reports legacy support,
   experimental support, phase, and the strongest disagreeing competitor with its F, E, L and P differences. Near ties
   are flagged as uncertifiable from rounded tables.
2. **Integration runs.** All 3 cells, joint off and joint on, with `block_support_rule full_data_likelihood_v1_experimental`.
   - **Requirement.** Every pre-existing output (selected sequences, scores, candidates, projections, block columns,
     rounds, status apart from the new rows) must equal the matching M3-D run's, with the volatile set of M3-D Amendment
     1.
   - **Cross-check.** The experimental columns must agree with `audit.saved`, except near ties flagged there.
3. **Budget.** 600 s, counting every attempt; 4 threads.
