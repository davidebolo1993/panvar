# M3-E.1: experimental block-support contract v2 (`full_data_likelihood_v2_experimental`)

**Ancestry.** Branch `m3-e1-certification` from `9cd255c` (M3-E).

**v1 is kept unchanged** (`full_data_likelihood_v1_experimental`, `block_data_support_v1`). The historical M3-E results
and Sol's reproducer refer to it. v2 is a separate function (`block_data_support_v2`), a separate profile value, and a
separate column suffix (`*_v2x`). The default stays `legacy`.

**Unchanged.** Selection, scores, priors, projections, phase, GQ, legacy columns and every existing filter's membership.
v2 is still support conditional on the evaluated candidates, not calibrated, and not a phasing claim.

## Reproduced v1 defects (Sol, `/private/tmp/panvar-m3e-review.uDIeES/repro.cpp`, run unchanged as `repro.v1`)

1. **Order dependence.** v1 picks the largest **point** likelihood and checks only that group's bound, so equal points
   with different bounds certify in one order and not the other: DATA_SUPPORTED vs DATA_TIED.
2. **Mixed competitor only.** v1 reports NO_EVALUATED_ALTERNATIVE with an empty quality reason, which would PASS.
3. **Mixed competitor with the highest likelihood.** v1 reports PRIOR_DEPENDENT although the strongest evidence is
   projection-ambiguous.

All three are kept as regression fixtures.

## v2 rule (frozen before implementation)

**Intervals.** Each group k has the interval [L_k − u_k, L_k + u_k], where L = fragment + exposure and u is the
certification bound. With the called group w and its genotype g at block b, the groups split into **S** (every member
carries g), **D** (no member carries g) and **M** (mixed). Over each class:

| envelope | definition |
|---|---|
| S_lo | the largest L − u over S (S always contains w) |
| S_hi | the largest L + u over S |
| D_lo, D_hi | the same over D |
| M_lo, M_hi | the same over M |

All are maxima over sets, so every label and certification is independent of input order.

**Labels, in this precedence:**
1. **AMBIGUOUS_PROJECTION:** w's members disagree at b. No representative is chosen.
2. **NO_EVALUATED_ALTERNATIVE:** D and M are both empty. It is **non-affirmative**.
3. **DATA_SUPPORTED:** S_lo > D_hi (if D is non-empty) **and** S_lo > M_hi (if M is non-empty). A same-genotype lower
   bound sits above every relevant competing upper bound.
4. **PRIOR_DEPENDENT:** D is non-empty, D_lo > S_hi, and (M is empty or D_lo > M_hi). Some definite disagreeing group is
   certifiably above every same-genotype group and every mixed group. Prior dependence is never asserted when a mixed
   group could be as strong as the best definite competitor.
5. **AMBIGUOUS_COMPETITOR:** M is non-empty and M_hi ≥ S_lo, so an unresolved mixed competitor could change the
   conclusion. This covers mixed groups alone, and mixed groups above both definite classes.
6. **DATA_TIED:** otherwise. D overlaps S, and any mixed group is safely below S_lo.

**Certification.** `certified` is true for DATA_SUPPORTED and PRIOR_DEPENDENT only.

**Margins** (order-independent):
- `support_gap_v2x` = S_lo − max(D_hi, M_hi), the envelope gap; it is positive exactly when the competitors are
  certifiably below;
- `point_margin_v2x` = max L over S − max L over D ∪ M.

**Competitor identity** (reporting only; it never affects a label). The strongest **definite** competitor is the D group
with the largest L + u. Exact ties are broken by the stable, content-based rule (largest L, smaller u, larger P, then the
lexicographically smallest genotype rows), which is also independent of input order. Its call-minus-competitor F, E, L
and P differences are reported. The mixed-competitor count and `mixed_upper_minus_same_lower_v2x` (M_hi − S_lo) are
reported.

**Experimental quality** (`quality_filter_v2x`): PASS, LOW_QUALITY or **UNASSESSED**.

| term | quality reason | effect |
|---|---|---|
| representation conflict / not evaluated | as legacy | LOW_QUALITY |
| library sensitivity | as legacy | LOW_QUALITY |
| DATA_TIED | `data_tied` | LOW_QUALITY |
| PRIOR_DEPENDENT | `prior_dependent` | LOW_QUALITY |
| AMBIGUOUS_PROJECTION | `ambiguous_projection` | LOW_QUALITY |
| AMBIGUOUS_COMPETITOR | `ambiguous_competitor` | LOW_QUALITY |
| NO_EVALUATED_ALTERNATIVE | `no_evaluated_alternative` | **UNASSESSED** when it is the only non-PASS term; LOW_QUALITY if any other term fails |
| DATA_SUPPORTED | none | PASS if every other term passes |

An absent comparison is never PASS. It adds no scientific no-call and changes no genotype-filter membership: the
prediction stands.

## Gates (`tests/genotype_block_support.cpp`, v2 section; the v1 gates stay as they are)

- **The reproducer's three cases** under v2:
  - **tied, different bounds:** the same label in both orders; it is DATA_TIED, because the larger-bound competitor's
    upper bound reaches above S_lo;
  - **mixed only:** AMBIGUOUS_COMPETITOR with reason `ambiguous_competitor`;
  - **mixed highest:** AMBIGUOUS_COMPETITOR, not PRIOR_DEPENDENT.
- **Heterogeneous bounds** on both sides, equal point likelihoods, and near ties just inside and just outside the
  envelopes.
- **Exhaustive permutations** of every 3–4-group fixture: labels, certification, gaps and the reported competitor are
  identical.
- **Mixed groups** alone, above both definite classes, between them, and safely below (DATA_SUPPORTED or
  PRIOR_DEPENDENT kept).
- **Absent comparison:** NO_EVALUATED_ALTERNATIVE gives UNASSESSED, never PASS. With a representation conflict it is
  LOW_QUALITY.
- **The v1 gates and the v1 reproducer** behave as before (historical behaviour preserved).
- **Caller (focused):** the synthetic e2e with v2 is annotation-only against the legacy run, with 1 and 4 threads
  identical, and the v1 annotation columns are unaffected.

## Evidence (no genotyping reruns)

1. **`audit.v2`.** The independent Python audit, updated separately to v2, on the saved M3-D and M3-E tables.
2. **`audit.cpp_v2`.** The C++ `block_data_support_v2` on the same saved (rounded) tables, through a probe subcommand.
   It must agree with the Python on every label.
3. **Report:** which annotations change from v1 to v2, and every decision whose envelope gap is within the rounded
   tables' precision (1e-4). Those require full-precision evidence and are not certified from rounded tables.

**Budget:** 300 s, counting every attempt.
