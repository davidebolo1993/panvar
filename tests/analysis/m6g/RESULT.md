# M6-G result: the exact C1 prototype behaves correctly but showed no structural benefit (the missing-variation premise failed), and the real-locus route needs a new approximate engine: PARK

**Run.** Branch `m6-g-c1-structural`: freeze `e9d8002`, Amendment 1 (runner input paths, before any unit ran). Results in
`tests/results/m6g_c1_structural_20261020/`.

**Budget.** 87.3 of 600 s, 4 threads. All 6 units exited 0 on the first attempt.

**Scope.** Synthetic fixtures built from two active-panel windows. No real-sample run, caller integration, flag,
reconstruction, sweep, merge or push.

**M6-F status, recorded at Sol's review.** M6-F's decision is accepted. C1 is carried forward as a research ingredient; C2
switching and κ are parked.

## Read model and library (as frozen)

**Layers.**
- **Structure:** arrangement G over the catalogue 1+1, 1+2, 2+1, 1+3, 2+2, 3+1, with uniform structural weights (explicit,
  identical for both methods). These induce total-copy prior masses of 1/6, 2/6 and 3/6 for totals 2, 3 and 4.
- **Component pattern:** library entry k with weight w_k.
- **Residual variation:** C1 at ε_res = 1e-4. It is distinct from the **sequencing error 0.002** of the existing
  paired-fragment likelihood (placements summed, exposure, background 0.001).
- **One shared sequence state per copy**, used for the whole read set.

**Prior.** C1 conditioned exactly on each method's admitted set and renormalised.
- **RES:** 8 distinct sequences per copy (2 windows × ancestral or novel allele at 2 admitted sites).
- **PANEL:** 2 per copy (the windows only).

**Library** (frozen rule; no reads existed when it was chosen):
- period_u13 and period_u26, weights 0.6 / 0.4 (founder weights 3 and 2 among M6-F training-only families).
- 200-bp window at period_u13 coordinate 2,600; 600-bp window at 2,200.
- In each window the two components differ at **one** position, near the window end (offset 187 of 200; 587 of 600).
- Admitted sites at L/3 and 2L/3. Novel alleles T,A (200 bp) and T,T (600 bp).

## Checks (`gates`: all pass)

| gate | result |
|---|---|
| GA prior normalisation (per copy; joint over every state) | 1; 4.4e-16 |
| GB prototype = independent exhaustive marginal (1+1, 1+2; RES and PANEL) | 0 difference (all four) |
| GC equivalent labels (alias split) | per-sequence prior unchanged (3.3e-16), 8 states |
| GD base annotation (re-annotated enumeration) | 0 difference |
| GE negative control (per-read independent marginalisation) | differs by −5.40 (detected) |
| GF reconstructed per-fragment terms = engine totals | 5.5e-12 |
| GG KL identity | 0 |
| GH constant likelihood over 512 states | returns the constant (no state-count penalty) |

## Outcomes (`summary`)

| fixture (unit) | method | allocation posterior | total-copy posterior (2 / 3 / 4) | frozen verdict | negative-control top |
|---|---|---|---|---|---|
| **A** missing variation, truth 1+2 (600 bp) | RES | 1+2 0.50, 2+1 0.50 | 0 / **1.000** / 5.6e-9 | **TOTAL_CORRECT** | 1+2 |
| | PANEL | 1+2 0.4997, 2+1 0.4997 | 0 / **0.9994** / 6.4e-4 | **PREMISE_NOT_MET** | 1+2 |
| **B** unbalanced control, truth 3+1 (200 bp) | RES | **3+1 1.000** | 0 / 0 / 1 | **CORRECT** | 3+1 |
| | PANEL | **3+1 1.000** | 0 / 0 / 1 | CORRECT | 3+1 |
| **C** uninformative, truth 2+2 (600 bp) | RES | 1+3 / 2+2 / 3+1 **0.3333 each** (exact tie) | 0 / 1.8e-5 / 0.99998 | **RETAINED** | 1+3 (tie) |
| | PANEL | same | same | RETAINED | 1+3 (tie) |

- **A's allocation (1+2 vs 2+1)** is a homologue-exchange tie at 600 bp. The total copy number is resolved separately.
- **B:** 3+1 beats 2+2 by 3,066 nats, through spanning adjacency (d E_q log L = 3,066.05; d KL = 0.00002). This is
  genuine read evidence for an unbalanced allocation, with no preference for balance.
- **C:** the total-4 allocations tie to 0.000000 in log marginal, so uncertainty is retained. Total 4 beats total 3 by
  10.5 nats (d E_q log L 11.0, d KL 0.51).

## Attribution: A (truth 1+2) against each alternative; d log m = d E_q[log L] − d KL

| versus | RES | PANEL | per-read negative control (RES) |
|---|---|---|---|
| 1+1 (total 2) | +4,555.3 = 4,576.3 − **21.0** | +4,532.0 = 4,532.6 − 0.5 | +4,527.2 |
| 2+1 (total 3) | 0 | 0 | 0 |
| 1+3 (total 4) | **+19.0** = 18.9 − (−0.17) | **+7.8** = 7.3 − (−0.5) | +7.25 |
| 2+2 (total 4) | +21.6 = 21.4 − (−0.23) | +7.8 | +7.25 |
| 3+1 (total 4) | +19.0 | +7.8 | +7.25 |

**Cumulative prior cost of the shared novel difference.** The truth's 6 novel site-copies cost **−61.85 nats** under RES:
log prior −63.79 against −1.94 for the same backbones without them. That is −20.62 per fully novel copy. It is reported,
not cancelled.
- Total 2 carries 4 novel site-copies, so the cost shows as **d KL = +21.0 against the truth**.
- Read evidence of 4,576 nats overrides it. Most of that comes from copy–copy junction fragments that a 1+1 arrangement
  cannot place.

**Against total 4, KL is essentially equal.**
- The posterior keeps the extra copy **ancestral** (posterior-mean novel site-copies 6.0, not 8). It pays 18.9 nats in
  likelihood rather than 20.6 in prior.
- So the residual prior **caps** what one extra copy could absorb at about 20.6 nats. Here it did not change the decision.

**RES fits the reads far better** (log m −4,229.2 against −4,620.6 for the truth structure). It **widens** the correct
total's margin over total 4 from 7.8 to 19.0 nats. No support for a wrong copy number was manufactured: no wrong total
exceeds 6.4e-4 under either method.

**The negative control shows why per-read marginalisation is invalid here.** It loses RES's benefit entirely: a margin of
7.25 under both methods, a RES log score of −4,638 against PANEL's −4,641, and a score about 409 nats below the correct
RES marginal. The shared sequence state is what lets the novel alleles be explained once for all reads.

## Computational report

- **States:** 8 distinct sequences per copy (RES), 2 (PANEL). Per structure 64 / 512 / **4,096** distinct states for
  2 / 3 / 4 copies, all enumerated with no truncation.
- **Per fixture:** 13,376 candidate diplotypes and 1,168 unique haplotypes in one engine call. 10.2–25.7 s, peak RSS
  ≤ 0.16 GiB.

**What grows on real C4.**
- **Per-copy states:** (components) × (admitted residual alleles). Here 2 × 2². On C4: about 100 components per role and
  class, with residual variation needing admission across about 32 kb (≥ 270 polymorphic positions per M6-E, plus novel
  positions anywhere). Per-copy states become astronomically many even when restricted.
- **Joint states:** per-copy states raised to the copy count (4–6 copies per diplotype).
- **Each likelihood evaluation** scales with fragments × haplotype length (≈ 32 kb per copy).
- **Placements:** alternative placements across paralogous copies couple every copy's content, and junction fragments
  couple adjacent copies. **Exact site-local factorisation is not available.**

**The only apparent route to real loci is a major new approximate inference engine** (windowed joint-copy forward with
pruning, or sampling). This is stated, not implemented.

## Decision (declared mapping): **PARK**

- **The prototype is correct as an inference object.** All gates pass, the cumulative novel-variation cost was reported
  and did **not** manufacture a wrong copy number, the unbalanced control was recovered, and uncertainty was retained where
  the observation models tie.
- **But the decisive test did not occur.** The frozen missing-variation fixture did **not** make fixed panel sequences
  favour the wrong structure (PANEL P(total 3) = 0.9994). The two methods made **identical** structural decisions on all
  three fixtures. The benefit seen was a wider margin and a better fit, not a corrected decision.
  - Redesigning the fixture until the premise holds would be post-hoc fixture repair, which the brief excludes. So it was
    not done.
- **The real-locus route needs a major new approximate inference engine.**
- **By the frozen mapping, PARK:** the premise was not met and no structural decision improved. The integration
  condition is also unmet.

**If it is revisited:**
- The prerequisite is an independently motivated fixture, or real-locus evidence, in which panel-only content
  demonstrably misleads structure under this read model.
- Separately, any implementation needs an approximate engine whose error against exact enumeration can be bounded on
  instances like these.

**Not claimed.** Nothing here bears on C4 NA19650, genotyping accuracy, or the full C1 model on real sequences.

## M6 closure (appended at M7-A, 2026-10-20; Sol's review; the record above is not rewritten)

1. **C1 remains a supported sequence-prior research ingredient** (M6-F).
2. **Parked:** C2 switching, κ, and full nuisance-marginalised read inference.
3. **M6-G's recovery question is inconclusive,** not negative. The PANEL premise was unmet, so the fixtures could not test
   recovery.
4. **Content-driven structural misranking is already demonstrated** by the existing C4 NA19650 evidence (M1-B to M1-E,
   M6-A/B).
5. **The ~20.6-nat residual cost is not a universal cap on evidence.** It is the prior cost of two novel site-copies at
   ε_res = 1e-4 in that fixture. Stronger data can overcome a prior penalty, so the sentence above calling it a cap on
   what one extra copy could absorb applies only to that fixture's evidence.

No further toy redesign or parameter sweep follows from M6.
