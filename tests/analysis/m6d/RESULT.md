# M6-D result: candidate centring of the nuisance prior caused both F2 failures; the common-content control is not dosage-neutral

**Run.** Branch `m6-d-prior-audit`: freeze `758eca3`, Amendment 1 `dc8ddc4` (numerics only). Results in
`tests/results/m6d_prior_audit_20261019/`.

**Budget.** 53.4 of 300 s, 4 threads, peak RSS ≤ 0.027 GiB. One recorded failure:
- `gates` exited 1. HA failed on naive summation of 1,024 masses: `long double` is 64-bit on arm64.
- HG failed because the posterior was formed through the rounding of log m.
- Amendment 1 changed the numerics only; tolerances were unchanged. `gates.a1` passed.

**Scope.** Synthetic only. No real loci, caller, flags, tuning, factorisation, merge or push. The M6-C qualifications are
appended to `tests/analysis/m6c/RESULT.md`.

## Gates (`gates.a1`): all pass

| gate | result |
|---|---|
| HA normalisation (9 allocations × 2 annotations × 2 priors) | max error 0 |
| HB common-content: base annotation only | Δ = 0. The original prior shifts by 11.78 (contrast, as expected) |
| HC ref/alt renaming with transformed probabilities | Δ = 0 (both priors) |
| HD synonymous labels (16 labels → 4 diplotypes; split ×2) | marginal Δ = 0; mass 1; KL Δ 3.2e-13 |
| HE constant likelihood, 1,024 states | returns the constant exactly; no −log(count) penalty |
| HF independent exhaustive oracle (1+1 base ref, 2+1 base alt) | Δ = 0 (both priors) |
| HG log m = E_q log L − KL | residual ≤ 9.1e-13 |

**Summary checks** (four runs, 9 scaffolds × 2 priors):
- C1, recomputation from the state tables: ≤ 4.6e-12.
- **C2, the original prior reproduces M6-C's recorded marginals: ≤ 4.8e-9.**
- C3, the identity: ≤ 1.8e-12.
- C0, prior mass: 0.

## Outcomes

Uniform weights over the nine allocations are unchanged. They induce total-copy prior masses of **2/9, 3/9 and 4/9**
(totals 3, 4 and 5).

| run | prior | allocation posterior | total 3 / 4 / 5 | allocation class | dosage class |
|---|---|---|---|---|---|
| F2a seed 23 (truth 2+2, all alt; 600 bp) | original | 1+3 0.489, 3+1 0.489, **2+2 2.9e-11** | 1e-4 / **0.978** / 0.022 | SPLIT (wrong allocation) | correct |
| | common-content | 1+3 0.3245, **2+2 0.3245**, 3+1 0.3245 | 0.026 / **0.974** / 4e-5 | SPLIT by the frozen rule. Log spread within total 4 is 1.2e-5, so not tied within 1e-6; this is retained allocation uncertainty | correct |
| F2b seed 11 (truth 2+2, all alt; 200 bp) | original | 2+3 0.500, 3+2 0.500, **2+2 1.8e-8** | 0 / 1.8e-8 / **1.000** | SPLIT (wrong dosage) | **wrong** |
| | common-content | **2+2 0.99999** | 0 / **0.99999** / 8.8e-6 | **CORRECT** | **correct** |
| F4 true 2+2 seed 23 | both, identical | **2+2 0.9997** | 0 / 0.9997 / 2.8e-4 | CORRECT | correct |
| F4 true 3+1 seed 23 | both, identical | **3+1 0.9967** | 0 / 0.9967 / 0.0033 | CORRECT | correct |

## Exact attribution: truth − alternative, d log m = d E_q[log L] − d KL

| run | versus | original: d log m = dE − dKL | common: d log m = dE − dKL | centring effect |
|---|---|---|---|---|
| F2a | 1+3 (original winner) | −23.555 = −0.0006 − 23.555 | −0.000003 = 0.00005 − 0.00005 | +23.555 |
| F2a | 2+3 (total 5) | −19.078 = 4.374 − 23.452 | +10.360 = 4.432 − (−5.928) | +29.438 |
| F2a | 1+2 (total 3) | −14.461 = 9.196 − 23.657 | +3.205 = 9.196 − 5.991 | +17.667 |
| F2b | 2+3 (original winner) | −17.114 = 6.339 − 23.453 | **+12.330 = 6.339 − (−5.991)** | +29.444 |
| F2b | 3+1 (same total) | +3,564.2 | +3,587.8 | +23.6 |
| F4 (both) | every alternative | identical under both priors | | 0 |

Per scaffold:
- **The truth's own KL is unchanged:** 23.97 under both priors. Its base is ref, so the two priors coincide for it.
- **The alternatives' KL rose:** the base-alt scaffolds paid 0.31–0.51 under the original prior and 17.97, 23.97 and
  29.96 under the control (totals 3, 4, 5). Their content had been free because the prior was centred on it.
- **Under both priors the posterior settles on all-alt content** (posterior-mean alt count of 8 or 10). The E_q log L
  differences are therefore the same under both priors in F2b, and nearly the same in F2a.

## Decision

**1. How much of each F2 failure comes from centring the prior on different candidate sequences?**

*F2a: all of the allocation error.*
- Against the winning 1+3/3+1, the truth's deficit of 23.555 is entirely KL. The E_q log L difference is −0.0006.
- Without centring, the three total-4 allocations are within 1.2e-5 nats of each other: the unobservable allocation is
  left uncertain.
- Dosage was correct under both priors (0.978 and 0.974).

*F2b: 23.45 of the 17.11-nat deficit came from the prior; the likelihood term favoured the truth throughout.*
- Under the original prior: 6.34 − 23.45.
- Under the control: 6.34 + 5.99.

**2. Does removing that dependence alter the wrong-dosage result?** Yes. F2b moves from total 5 at 1.000 to total 4 at
0.99999, with 2+2 at 0.99999. Caveat:
- **The control is not dosage-neutral.** It charges each copy for rare content: about 6 nats per extra all-alt copy
  (2·log 0.05), against 0.10 per all-ref copy.
- Of F2b's 12.33-nat margin, 5.99 is this charge. The other 6.34 is the posterior-mean likelihood difference, which here
  equals the all-alt likelihood difference because both posteriors concentrate there.
- In F2a the same charge raises total 3 from 1e-4 to 0.026. The truth still leads total 3 by 3.2 nats.

**3. Are both allocation controls preserved?** Yes, identically. This is **mechanical**, as declared: every F4 base is
ref, so the two priors are the same array. It is not an independent test of the control.

**4. What remains unsupported even though the control improves the results?**
- **The control is not a prior to adopt.** Its 0.05 is declared, not grounded, and its per-copy content charge creates a
  copy-number preference that depends on allele rarity. Any panel-grounded prior will make some such dosage commitment,
  and that commitment needs its own study.
- **One seed per fixture, two sites, synthetic scaffolds.** Missing content at undeclared sites (indels, unknown
  positions) is not addressed.
- **Allocation evidence still comes only from adjacency** (F4, and F2b's exclusion of 3+1 and 1+4). F2a's allocation is
  correctly left unresolved, not resolved.
- **Same-total allocations are not exact observation-model equivalents.** Their log marginals differ by 1.2e-5 (F2a) and
  3e-4 (F2b 2+3 vs 3+2). Nothing here assumes equivalence from equal pooled content.
- **Still open:** candidate-independent site construction on real loci, exact site-local factorisation, and the induced
  total-copy prior (2/9, 3/9, 4/9).
- **Nothing here bears on real C4.**

**Verdict.** The improvement is real and attributable: candidate centring explains both F2 failures, and the allocation
controls are unchanged. By the declared rule, this **justifies investigating a panel-grounded generative prior**, with
its copy-number implications as a first-class question. **This control is not adopted.**
