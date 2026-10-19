# M6-G: bounded bridge from C1 to read-based structural inference (frozen before any unit)

**Ancestry.** Branch `m6-g-c1-structural` from `ddd88ed` (M6-F result).

**M6-F status, recorded at Sol's review.** M6-F's decision (C1_SUFFICES) is accepted. C1 is carried forward as a research
ingredient. C2 switching and κ are parked. There is no further panel-predictive tuning.

**Question.** Can a C1-based sequence prior accommodate missing variation in read-based structural inference without
manufacturing support for the wrong copy number?

**Scope.** One exact, restricted feasibility prototype on synthetic fixtures built from active-panel sequence windows.
- **Excluded:** real-sample runs, caller integration, flags, reconstruction, parameter sweeps, merge, push and cleanup.
- This is **not** the full C1 model for real C4.

## Read model (`tests/c1_structural.cpp`)

| layer | definition |
|---|---|
| structural arrangement G | homologue A = F_A + n_A copies + G_A; homologue B = F_B + n_B copies + G_B. Flanks are synthetic and unique (1,200 bp each, seeded). Catalogue: 1+1, 1+2, 2+1, 1+3, 2+2, 3+1. **Uniform structural weights** (explicit, identical for every compared method), which induce total-copy prior masses of 1/6, 2/6 and 3/6 for totals 2, 3 and 4 |
| panel component supplying a pattern | library entry k with weight w_k |
| residual sample-sequence variation | C1 emission e_res(X_x \| b_k(x)) = 1 − ε_res on a match, ε_res/3 otherwise, with **ε_res = 1e-4** (M6-F's C1 selection; it sat at its grid floor, a known limitation). Admitted set S (below); prior P(X \| S) = P_C1(X) / Σ_{X'∈S} P_C1(X'), **exact and renormalised**, with P_C1(X) = Σ_k w_k Π_x e_res(X_x \| b_k(x)) summed over every library entry |
| sequencing errors and read placements | the existing paired-fragment likelihood (`score_diplotypes`): placements summed, **sequencing substitution error 0.002**, exposure, outlier/background 0.001; insert 350 ± 35 bp; 100-bp reads; diploid depth 20 |

- **One shared sample sequence.** Each copy's sequence X_j is **one** state shared by the whole read set. Copies are
  independent given G (κ parked), and the joint prior is the product over copies.
- **Exact marginal.** log P(reads | G) = log Σ_X P(X | G) P(reads | G, X), over every distinct admitted state. It is never
  replaced by the nearest component, the best edited sequence or a per-read sum.

**Two methods, same reads, same structural weights:**
- **RES:** S = {library window} × {ancestral or declared novel allele at each of 2 admitted sites}, giving 8 distinct
  sequences per copy. That is ≤ 8⁴ = 4,096 distinct states per structure.
- **PANEL:** S = {library window}, i.e. fixed panel sequences marginalised over components; 2 per copy.

**Negative control.** Per-read independent marginalisation, never used for inference:

  NC(G) = Σ_f log E_X[exp(term_f(X))] + log E_X[exp(−expected signal(X))]

The per-fragment terms are rebuilt from the engine's own `log_mass`, `exposure` and `log_background_probability`. Gate GF
certifies that they sum to the engine's totals.

**Decomposition.** For every structure and method: log m = E_q[log L] − KL(q ‖ P(X | G)), with compensated, centred
evaluation (as in M6-D). Contrasts are truth − alternative. d E_q[log L] is read evidence; d KL is the prior/complexity
side. A best-state log likelihood is reported as a labelled diagnostic only.

## Frozen library (`prepare`; outcome-independent; no reads exist when it runs)

- **Components:** the 2 period/long components with the largest founder weight among M6-F's training-only families. Ties
  go to more haplotypes, then the smallest SHA-256. Coordinates are M6-F's period/long reference (period_u13).
- **Weights:** founder weights renormalised over the two (the restriction is explicit).
- **Windows,** for L = 200 and L = 600: the first start s ∈ {0, 50, 100, …} where:
  - neither component has an indel or uncertain correspondence within [s − 10, s + L + 10);
  - the components differ at ≥ 1 position in the window;
  - the admitted sites s + ⌊L/3⌋ and s + ⌊2L/3⌋ carry the same base in both components, with no component difference
    within 10 bp.
- **Novel allele** at each admitted site: the next base in ACGT, absent from both components.

## Fixtures (frozen; one seed each)

| fixture | unit | truth | copy contents (entry 0 = heavier component) | role |
|---|---|---|---|---|
| A_missing_variation | 600 bp | 1+2 | A: [0 + novel at both sites] / B: [0 + novel at both, 1 + novel at both] | missing variation; **the same novel differences in all 3 copies** (6 novel site-copies, cumulative prior cost reported, not cancelled) |
| B_unbalanced_allocation | 200 bp | 3+1 | A: [0, 1, 0] / B: [1] | informative, genuinely unbalanced allocation control |
| C_uninformative_allocation | 600 bp | 2+2 | A: [0, 0] / B: [0, 0] | allocation-uninformative (no fragment spans a copy): uncertainty should be retained |

## Checks (`gates`, on fixture B's world and reads; all must pass before any fixture run)

| gate | check |
|---|---|
| GA | per-copy priors (RES, PANEL) and the joint prior over every enumerated state sum to 1 (≤ 1e-12) |
| GB | prototype marginal = independent exhaustive marginal (each state scored alone, prior recomputed in-loop) for 1+1 and 1+2, RES and PANEL (≤ 1e-6) |
| GC | equivalent labels: splitting library entry 0 into two aliases of half weight leaves every per-sequence prior unchanged (≤ 1e-12), with 8 distinct states |
| GD | base annotation: enumerating the same states from a different base annotation gives the same marginal (≤ 1e-9) |
| GE | negative control: per-read independent marginalisation differs from the correct marginal (> 1e-3) |
| GF | reconstructed per-fragment terms + exposure = engine totals for all candidates (≤ 1e-6) |
| GG | KL identity (≤ 1e-9) |
| GH | constant likelihood over 512 states returns the constant (no state-count penalty) |

## Criteria (`summary.py`, frozen)

| fixture | criterion |
|---|---|
| A | **PANEL premise met** if PANEL's P(total 3) < 0.5. **RES:** TOTAL_CORRECT if P(total 3) ≥ 0.95; MANUFACTURED_WRONG_TOTAL if another total has ≥ 0.95; otherwise TOTAL_UNRESOLVED. Allocation is reported separately |
| B | CORRECT if P(3+1) ≥ 0.95; BALANCED_PREFERENCE if 2+2 is the top allocation; otherwise OTHER |
| C | RETAINED if the largest conditional allocation posterior within total 4 is < 0.95; otherwise OVERCONFIDENT |

## Decision mapping (declared now)

- **Bounded integration route justified:** only if RES is TOTAL_CORRECT in A with the PANEL premise met, B is CORRECT, C
  is RETAINED, **and** a route to real loci exists without a major new approximate inference engine. M6-E already
  estimated that exact enumeration is infeasible for C4, with up to 10⁶–10⁹ windowed states. Unless this batch shows
  otherwise, the third condition is expected to fail, which will be stated rather than implemented.
- **Model needs a specific revision:** if RES manufactures or fails to resolve copy number in A because of the cumulative
  residual prior cost, or another named defect appears.
- **Park:** if RES offers nothing over PANEL (the premise is not met and nothing improves), or RES fails B or C.

## Bounds

- 600 s counting every attempt; 4 threads.
- **Declared timeouts:** prepare 30 s, gates 120 s, each fixture 120 s, summary 30 s.
- Failures are retained.
