# M9-A result: SV-block genotypes are far more certain than the exact sequence pair, and phase is a separate, prior-only quantity; recommend the reportable-block genotype as the primary confidence target

**Run.** Branch `m9-a-uncertainty-audit`, freeze `d48e2c5`, with no amendments. Results are in
`tests/results/m9a_uncertainty_audit_20261020/`.

**Execution: 142.3 of 300 s.**
- Before the freeze (disclosed): 98.5 s, mostly one 91-s dry run before the audit was optimised.
- Gates: 0.1 s.
- Audit: 43.7 s, with 0.72 GiB peak.

**Scope:** read-only. No placement, genotyping, model, clustering, threshold, filter, PASS rule or GQ calibration.
Historical packages are unchanged. M8-D received an appended qualification only (§1).

**M8 is closed:** no prior or default promotion, no further fitting, no broader campaign.

## 1. M8-D qualification (reproduced)

| run | locus genotype separation | mass carrying the reported unphased genotype | rank-2 sequence pair (the 0.04-nat gap) carries the same unphased genotype |
|---|---|---|---|
| HG01975 LOO r0.01 | 0.70213715071 (Sol: 0.7021371507) | 0.405311 (Sol: ≈ 0.405311) | **yes**: likelihood gap 0, prior gap 0.0443 |
| HG01975 LOO r0.1 | 0.90184718165 (Sol: 0.9018471817) | 0.313111 (Sol: ≈ 0.313111) | **yes**: likelihood gap 0, prior gap 0.0429 |

**The 0.04-nat near tie is a phase tie.** The qualification is appended to `tests/analysis/m8d/RESULT.md`.
- **Unchanged:** the search conclusion. Both objectives prefer the closer final candidate.
- **Not inferred:** global optimality or better prior calibration.

## 2. Pools

| run | hypotheses (= label pairs) | multi-label | projection-ambiguous | top set | scope | het blocks |
|---|---|---|---|---|---|---|
| HG03209 LZO r0.01 / r0.1 | 84,533 / 84,533 | 0 | 0 | 1 | CONDITIONAL, search not exhaustive | 13 |
| HG03209 LOO r0.01 / r0.1 | 83,719 / 83,719 | 0 | 0 | 1 | CONDITIONAL | 13 |
| HG01975 LZO r0.01 / r0.1 | 84,524 / 84,525 | 0 | 0 | 1 | CONDITIONAL | 14 |
| HG01975 LOO r0.01 / r0.1 | 103,206 / 104,039 | 0 | 0 | 1 | CONDITIONAL, tie-branched | 13 |

**What holds in every run:**
- The reported state is the caller's top hypothesis, unrepaired (dosage changed nothing).
- Every sequence hypothesis has exactly one label pair, so in these data:
  - no MIXED hypothesis occurs;
  - the mass intervals are degenerate (lower = upper);
  - the mixed-projection and same-sequence-under-different-labels paths are exercised only by the gates.

**Pools differ between runs.** Mass differences across runs are not calibration evidence.

## 3. Gates (all 8 PASS) and reproductions

**Gates:**

| gate | fixture |
|---|---|
| F1 | whole-haplotype swap counted once |
| F2 | same block genotypes with different phase: B = 1; A and E split, with the phase competitor |
| F3 | genuinely different genotypes: competitor, total gap 2, likelihood gap 1.5 |
| F4 | identical sequence under different labels: one hypothesis, mass once; MIXED only at the disagreeing block |
| F5a | a label pair in two hypotheses is rejected |
| F5 | mixed projection: interval [0.507, 0.814], and the competitor excludes the MIXED hypothesis |
| F6 | absent alternatives give NO_EVALUATED_ALTERNATIVE |
| F7 | normalisation, and D equals the caller's block posterior formula |

**Reproductions** (`reproductions.tsv`, all eight runs):

| quantity | result |
|---|---|
| normalisation | 1 ± 6e-12 |
| `candidate_set_probability` | within 5.4e-9 |
| `carries_called_genotype` mass | exact |
| `fragment_contenders` | exact (84 … 1,304) |
| `phase_support` | exact |
| per-block `genotype_support` | 152/152 labels |
| per-block posterior | within 7.7e-9 |
| locus separation and `numeric_separation_gap` | within 1e-7 |

**Explained differences:**
- **Rounding.** `mosaic_scores.tsv` prints 12 significant digits (about 1e-7 at |score| ≈ 8e4), against 17 in the locus
  table.
- **Definitions.** The caller's separation and block posterior assign each hypothesis its representative genotype, and
  count MIXED hypotheses as "other". The audit uses intervals and genuine competitors instead. With no MIXED hypotheses
  these coincide here.

## 4. Targets per run (mass within the run's evaluated pool; gaps in nats to the best genuinely disagreeing hypothesis)

| run | A exact sequence pair | B unphased, all 19 blocks | C reportable-block (9 SV) agreement | E phase given B |
|---|---|---|---|---|
| HG03209 LZO r0.01 | 0.99962 (gap 8.52) | 0.99962 (8.52) | 1.000000 (16.3) | 1.000000 (19.0) |
| HG03209 LZO r0.1 | 0.98585 (4.52) | 0.98585 (4.52) | 0.999982 (11.5) | 0.999999 (14.1) |
| HG03209 LOO r0.01 | 0.99327 (5.06) | 0.99327 (5.06) | 1.000000 (16.9) | 1.000000 (19.8) |
| HG03209 LOO r0.1 | 0.91054 (2.63) | 0.91054 (2.63) | 0.999978 (11.4) | 1.000000 (15.0) |
| HG01975 LZO r0.01 | 0.97938 (4.16) | 0.97938 (4.16) | 1.000000 (18.4) | 1.000000 (18.4) |
| HG01975 LZO r0.1 | 0.80004 (1.73) | 0.80004 (1.73) | 0.999994 (13.6) | 0.999999 (13.6) |
| **HG01975 LOO r0.01** | **0.2025** (0.044) | **0.4053** (0.702) | **0.8378** (1.19) | **0.4996** (0.044) |
| **HG01975 LOO r0.1** | **0.1551** (0.043) | **0.3131** (0.902) | **0.8085** (1.10) | **0.4955** (0.043) |

**D, individual blocks** (`blocks.tsv`). Blocks below 0.99:

| run | below 0.99 | weakest |
|---|---|---|
| HG03209 LZO r0.1 | ctx:10 0.989 | — |
| HG03209 LOO r0.1 | ctx:0 0.932, ctx:10 0.990 | — |
| HG01975 LZO r0.01 / r0.1 | ctx:0 0.980 / 0.807 | — |
| HG01975 LOO r0.01 | ctx:0 0.889, ctx:2 0.965, ctx:4 0.721, ctx:6 0.901, **sv:9 0.838**, ctx:12 0.939 | SV: sv:9 |
| HG01975 LOO r0.1 | ctx:0 0.818, ctx:2 0.971, ctx:4 0.675, ctx:6 0.878, **sv:9 0.819**, ctx:12 0.947 | SV: sv:9 |

Every other SV block is ≥ 0.994 in every run. Every block below 0.99 has PRIOR_ONLY genotype support.

**E, phase** (`phase_pairs.tsv`, adjacent het pairs, given B):
- **LZO and HG03209 LOO:** every pair ≥ 0.99.
- **HG01975 LOO:** ctx:8–sv:9–ctx:10 at 0.97–0.98; **ctx:10–ctx:12–sv:13 at 0.52 / 0.53**, a coin flip in one junction.

## 5. What the audit shows

1. **The targets separate cleanly and differ by large factors where it matters.**
   - In six runs, A ≈ B (one B-carrying hypothesis holds the mass) and C ≈ 1.
   - In HG01975 LOO, the exact sequence pair holds 0.20 / 0.16. Within it:
     - phase alone takes about half (E | B ≈ 0.50, a single junction);
     - context blocks (ctx:0/2/4/6/12) take most of the rest of B;
     - the 9 SV blocks jointly keep 0.84 / 0.81, almost all of the shortfall at sv:9.
   - A whole-locus sequence-pair probability therefore understates SV genotype certainty by about 4–5× there.
2. **B is pulled down by non-reportable context blocks.** For example, HG01975 LZO r0.1 has B = 0.80 from ctx:0 = 0.81,
   while C = 0.999994.
3. **Phase is prior-only in all eight runs.**
   - The best phase competitor always has a likelihood gap of exactly 0.
   - `phase_support = PRIOR_ONLY` everywhere.
   - E is a statement about the copying prior, not the reads, and must not be merged into genotype confidence.
4. **Mass and likelihood support disagree often.**
   - The best disagreeing competitor frequently has the *better* likelihood. Examples: HG03209 LZO r0.1 B, −6.89 nats;
     HG01975 LOO sv:9, −0.19; sv:13, −3.46.
   - The prior then decides the reported state.
   - Conditional mass (prior-weighted) and likelihood support (`genotype_support` READS / PRIOR_ONLY, likelihood gaps)
     answer different questions.
5. **Mass is conditional on the pool, not a correctness probability.**
   - In HG01975 LOO, sv:5 has D = 1.000000 and READS support. Its gap to the best *evaluated* disagreeing pair is 113 /
     136 nats.
   - Yet M8-D's evaluation found 10 bp of excess there: called 89/62 against a truth-closest 89/113, with floor 0.
   - All pools are CONDITIONAL with non-exhaustive search.
   - This is an illustration of conditionality, not a calibration measurement.

## 6. Recommendation (production-facing; no implementation, thresholds, filters, PASS rules or GQ calibration)

### Primary genotype confidence

**The primary genotype confidence should describe target C, the reportable-block genotype agreement.** That is the mass,
within the evaluated pool, of hypotheses carrying the reported unordered allele pair at every reportable (SV) block.
- **Per SV record,** the same family restricted to that block: D at reportable blocks, which is the existing
  `conditional_genotype_posterior` / `conditional_GQ` of those blocks.

Why C:
- It is what SV-focused use consumes.
- It is insensitive to phase and to non-reportable context blocks, which dominate A and B uncertainty here.
- It already exists per block.

How it must be presented:
- **as conditional:** "within the evaluated candidate pool", carrying `candidate_scope`;
- **with its evidence basis beside it:** the per-block `genotype_support` READS / PRIOR_ONLY;
- **named as label agreement:** "reportable-block genotype agreement", not structural equivalence. Similar SV alleles are
  not clustered.

### Phase confidence

**Phase confidence is reported separately,** as E (phase given the reported genotype) and the adjacent het-pair phase
masses, always with `phase_support`.
- In all eight runs it is prior-only, and it should be labelled as such.
- It is never folded into genotype confidence.

### Detailed / debug output

These are diagnostics, not the primary confidence statement:
- **Exact sequence-pair quantities (A):** `winner_candidate_set_probability`, `numeric_separation_gap`, and the
  equivalence set. They belong with the FASTA sequence output.
- **Locus unphased-genotype quantities (B):** `separation` (`caller_genotype_separation_gap`), `carries_called_genotype`
  mass, and the best competitor's fragment, exposure and prior gaps.
- **Candidate-set diagnostics:** `fragment_contenders`, `fragment_term_gap_to_runner_up`, `representation_*` columns,
  `projection_ambiguous`, MIXED counts, pool size and search flags.
- **Naming recommendation:**
  - `numeric_separation_gap` (A) and `separation` (B) should not appear under the unqualified word "separation" in
    headline output;
  - M8-D's misreading came from exactly this.

### Not interchangeable

Conditional posterior mass (prior-weighted, pool-conditional) and likelihood support (read-only margins) are different
quantities. A single confidence number must say which one it is. This audit recommends mass for the primary statement,
with likelihood support shown as its basis label.

## Limitations

- Eight runs from two donors, one locus (ACOT) and one coverage.
- No MIXED or multi-label hypotheses occurred, so interval semantics are validated by fixtures only.
- Masses are conditional on non-exhaustive, run-specific pools.
- No calibration is claimed.

## Qualification (appended by M9-B, 2026-10-20; the text above is unchanged)

**1. A and C describe different propositions.** They are not two measurements of one certainty.
- **A:** the reported exact unordered whole-haplotype sequence pair.
- **C:** joint agreement of the reported unordered block-allele pairs across all reportable blocks.
- **What replaces "understates SV certainty by about 4–5×":** the conditional masses themselves. In HG01975 LOO,
  within each run's evaluated pool:

| run | A | C |
|---|---|---|
| r0.01 | 0.2025 | 0.8378 |
| r0.1 | 0.1551 | 0.8085 |

  Neither value corrects the other.

**2. C is joint, D is marginal.**
- **C** is the mass of hypotheses agreeing at *every* reportable block simultaneously.
- **D** is agreement at *one* block.
- Joint agreement is not the product of block marginals. M9-B fixture J1 has two block marginals of 0.9 with joint
  agreement 0.8; J2 has the same marginals with joint 0.9. Any value in the Fréchet interval [0.8, 0.9] is possible.

**3. Neither C nor D is a calibrated probability of biological correctness.**
- Both are conditional posterior masses within a non-exhaustive, run-specific candidate pool.
- They are labels over block alleles, not structures.
- Example (see M9-B §4): in HG01975 LOO, sv:5 has D = 1.0 while carrying a 10-bp exact-sequence error.

**4. The recommendation is refined in M9-B.** The declared output unit is the reportable block, so the primary value is
per-block D of exact block-allele agreement. C is a locus-level joint summary.
