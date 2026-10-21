# M13-A: variation-tolerant structural scoring — stopped before substantive execution (computational obstacle)

**Decision: INCONCLUSIVE.** The defensible prototype does not fit the tested contrasts, so the question is unanswered:
- **No real-case comparison was run,** and no truth comparison was opened.
- **Under the brief's stop rule,** this report gives the specific obstacle and a concrete implementation estimate.
- **No other heuristic was substituted.**

**Run.** Branch `m13-a-variation-tolerant`; result commit after.
- **Development evidence:** `tests/results/m13a_variation_tolerant_20261021/development/`.

**Execution time: about 85 s of the 2,400 s,** all of it development:
- builds, about 34 s;
- two gate runs, 4.2 s and 4.4 s;
- dry admission over the six contrasts, 38 s;
- candidate preparation.

**Not done:**
- no protocol or runner freeze (substantive execution never started);
- no caller integration, FASTA, validation family, simulation campaign, default, merge, push or cleanup.

## 1. Carried forward (unchanged)

- **M12's hard dosage-first rule stays parked.**
- **Rates, exposure arithmetic and placement calculations are consistent** (M12-A).
- **Supplying true content reversed the copy-number preference** in M12-A's examined cases.
- **Forcing the pooled feature state did not solve within-state selection,** and it worsened K1 and K5 (M12-B/C).
- **A pooled k-mer feature vector is not automatically biological copy number.**
- **Wrong predictions do not establish short-read non-identifiability.**

## 2. What was built (diagnostic only, tested)

**`tests/vt_marginal.hpp`: the truncated, group-factorised marginal.**
- **Content:** each homologue's block = its spelled block with substitutions at **admitted** sites. Each site has one state,
  shared by every read.
- **Site admission** (truth-free; the same rule for every candidate): ≥ 2 fragments whose best gapless placement shows the
  same non-reference base at a block position.
- **Prior:** the M6-F C1 whole-component mixture over panel block alleles. It is a function of the sequence alone,
  identical across candidates, with ε_res = 1e-4 and active-panel template-count weights. It is substitution-only, so only
  equal-length components contribute.
- **Likelihood:** the existing paired-fragment likelihood (`score_diplotypes`). The residual rate is separate from the
  sequencing error.
- **Exactness argument:** given the component pair, the prior factorises over positions, and a fragment's term depends only
  on sites its placements touch. So the sum over admitted states is a product of exact per-group sums, with groups defined
  by shared fragments.
- **Interaction check:** a joint MAP-state rescoring measures any coupling the grouping misses.
- **Reporting:** as a truncated sum (a lower bound), never as an exact full marginal.

**`genotype_probe vt-score`:** the real-case driver.
- the run's model, reads, library and scoring settings;
- the ordinary fixed-content score with the run's copying prior (the expected total must reproduce);
- m_panel (fixed content under the C1 prior), m_res (variation-tolerant), its decomposition and site/group counts;
- `--dry 1` runs admission only.

**`tests/m13a_vt_gates.cpp`: 9 of 9 checks pass** (`development/gates/gates.log`).

| gate | result |
|---|---|
| fixed content | equals an independent scorer call; m_res = m_panel |
| constant likelihood | equals the analytic admitted prior mass (−2.397 ≤ 0); no state-count bonus |
| alias split, homologue exchange | invariant |
| shared state vs per-read redraw | differ (−3,108.7 against −3,180.8) |
| **exactness** | **engine = brute-force enumeration over every joint state** (M6-G's GB construction): 5 sites, coupled across tandem copies |
| informative unbalanced 3+1 with one variant absent from the panel | the marginal prefers 3+1 |
| uninformative 2+2 / 1+3 | exact tie |

**Synthetic observation** (gate world, not evidence about real loci):
- **Fixed panel content scores 3+1, 2+2 and 1+3 identically,** because the units are longer than every insert.
- **The variation-tolerant marginal separates 3+1 from 2+2 by 61 nats** through the variant's flanking context.
- **1+3 hit the state cap:** one variant was admitted at all 4 copies plus noise, 9 coupled sites.

## 3. The obstacle (dry admission on the frozen candidates; truth-free)

The frozen contrasts are K4, K5, K1 and controls R1–R3. Each includes the incumbent, the M12 selection where different, and
the best-available comparator (controls: the run's own best recorded alternative at the block), in both whole-haplotype
assignments (`development/candidates/`).

| contrast | admitted sites per candidate | largest fragment-coupled group | joint states to enumerate | status |
|---|---|---|---|---|
| K4 | 111–133 | 21–26 sites | 2.1M – 68M | intractable |
| K5 | incumbent 8, M12 13; **best-available 83** | 4 / 4 / **22** | 20 / 27 / **6.3M** | the key comparator is intractable |
| K1 | 79–175 | 12–62 | 5,638 – 4.6e18 | intractable (best-available exceeds the 300-variant cap) |
| R1 | 162–170 | 34–36 | 1.7e10 – 6.9e10 | intractable |
| R2 | 216–217 | 24–26 | 18M – 68M | intractable |
| R3 | 99–100 | 14 | about 17k | intractable under the caps |

**Why the groups are so large:**
- **Copy coupling:** reads place across near-identical tandem copies, so a read-supported difference is admitted at every
  copy where the read fits equally well, and those sites share fragments. Which copy carries the variant **is** the
  allocation question, so it cannot be factored away.
- **Neighbour coupling:** neighbouring sites within one insert share fragments as well.

**Even correct controls admit 99–217 sites.** The admitted set therefore absorbs more than missing variation, including
likely mis-placed paralogous reads. A scaled engine would have to model that, or its scores would be dominated by it.

**Exact enumeration within groups is infeasible on every tested contrast,** and on all but two candidates (K5's incumbent and
M12 selection). Scoring each joint state by re-running the whole-haplotype scorer (about 0.25 s per variant haplotype)
cannot reach these sizes.

## 4. Implementation estimate for a version that could answer the question

**1. A local fragment-term engine.** It would recompute each fragment's placement mass under site states from its recorded
origins, instead of rescoring whole haplotypes.
- **Gapless mates:** the emission changes by an exact factor per site.
- **Gapped mates:** re-run the existing banded pair-HMM on the ≤ 110-bp window.
- **Placements entering the band under a variant state:** collect them with a band-widened base scan (band + 1). This must
  be validated.
- **Cost:** microseconds per fragment-state, against 0.25 s per haplotype now.
- **Size:** about 400–600 lines.

**2. Approximate inference over the site–fragment factor graph.** Fragments are higher-order factors over the sites they
touch, typically ≤ 10–16, so each factor table is enumerable. Two options:
- **loopy belief propagation**, with the Bethe free energy as an **approximate, unbounded** marginal estimate;
- **annealed importance sampling**, with an unbiased marginal estimate and a reported variance. This one is preferable for
  ranking because its error is quantified.

Either must be named as an approximation, never as an exact marginal. Size: about 400–600 lines.

**3. Validation against the exact oracle that already exists** (`vt_marginal.hpp` + brute force) on synthetic worlds with
coupled groups of up to about 12 sites, and on K5's two tractable candidates as a real-data cross-check.
- **Bias rule:** if the estimate's discrepancy or variance could reverse a contrast, that contrast is classified UNRESOLVED.

**4. A design decision needed first:** how to treat the 100–200 sites admitted even under correct content. Options are a
stricter but still truth-free admission rule, or an explicit paralog/mis-placement component. This is a modelling decision,
not a tuning step, and it should be made before any real contrast is opened.

**Expected cost:**
- one implementation batch (engine + gates; under 900 s of execution);
- one evaluation batch (6 contrasts at about 30–60 s each including base scans; under 600 s).

**The question remains open.**

## Closure (2026-10-21): qualifications appended; the historical record is unchanged

**M13-A is closed as INCONCLUSIVE.** No real-case scoring and no accuracy claim.

**Qualifications:**
1. **Factorisation is verified only on the exhaustive fixtures** (engine = brute force). It is **not** established in
   general by the joint-MAP interaction check, which tests one joint state.
2. **Candidate-specific admission creates different truncated state spaces** for different candidates. Rankings of these
   lower bounds do not certify rankings of the full marginals.
3. **High admission counts in correct controls (99–217 sites) do not establish contamination or incorrect placement.**
4. **A widened placement scan (band + 1)** would need a completeness argument covering combined substitutions, where
   several sites change together.
5. **Monte Carlo variance alone would not certify** that an approximate engine explored all important modes.

**The proposed approximate engine is not implemented.** The diagnostic code (`tests/vt_marginal.hpp`,
`tests/m13a_vt_gates.cpp`, `genotype_probe vt-score`) is preserved.
