# M3-H: bounded evidence decomposition of the ctx:0 comparison (diagnostic only)

**Ancestry.** Branch `m3-h-evidence-decomposition` from `54cbad0` (M3-G), plus the M3-G qualifications (`7f02018`).

**No model or caller change.** The binary is the M3-F/M3-G one, `a728eb71…`.

## 1. The frozen matched comparison

**A = R_t0** (ctx:0 = (14,56), the truth pair) against **B = R_c0** (ctx:0 = the called (4,52)).
- **Context:** the repair's final incumbent R, everywhere except ctx:0.
- **Assignment 0:** hap1 carries 14 or 4, and hap2 carries 56 or 52.
- **Same sequences at both coverages:** the same two complete sequence pairs are scored against the 2× reads and the 5×
  twin's reads. No flank, phase or candidate changes.

**Why this context and assignment** (from the M3-G results):
1. **R_t0 is the truth path exactly.** Its bytes are checked against the projected truth by `spell`, so A's flanks and
   phase are the true ones. The only difference from B is the ctx:0 alleles.
2. **R_c0 carries the called likelihood exactly.** Its fragment and exposure terms are bit-identical to the reported
   call's (W1/W2: F −29691.653054385122, E −2129.0207068360523), so the likelihood comparison is the reported one.
3. **The deltas are the reported ones.** The truth-minus-called deltas in this context equal those against the call:
   ΔF −1.332, Δprior −10.164.
4. **Assignment 0 is the better-scoring called assignment** (R_c0 −0.613 vs R_c1 −8.985 against R). The two (14,56)
   assignments score identically.

**Identity checks** (score-set refuses unless they reproduce within 1e-6):
- **2×:** R_c0 equals the saved caller neighbour 8652 (−31872.0131623); R_t0 equals the M3-G probe value
  (−31883.508709682916).
- **5×:** R_t0 is in the saved 5× winning equivalence set (rank 2, seq 2753, −74664.792055680344).

**Libraries.** Each coverage uses its own frozen library file. They differ only in `fragment_start_intensity`
(0.005005 at 2×, 0.012513 at 5×).

**Nesting, verified before freezing.** All 2,145 2× pairs are in the 5× pool, with identical bytes in both mates, and
5× adds 3,229. The subsampler records "fnv1a64 … nested across coverages".

## 2. Likelihood decomposition

**Probe calls.**
- `score.2x` and `score.5x`: `score-set --decompose R_t0,R_c0 --decompose-block ctx:0 --contrast R_t0,R_c0`.
- **Reconstruction checks:** the per-pair sum equals the scorer's fragment term; the canonical-origin reconstruction
  minus the scorer is 0; `placement_evidence_interpretable` is reported.

**`read_detail.py` reports:**
- counts favouring A, favouring B and neutral (exactly ≤ 1e-9; effectively < 0.01 nats);
- concentration (top 1/5/10/20/50) and the delta by placement class;
- the 5× delta split into the shared 2× pairs and the added pairs, and the shared pairs' 2×-vs-5× deltas;
- the discriminating 31-mers at ctx:0 per homologue;
- for the top 10 contributors at 2× and the top 10 added at 5×: per-mate gapless placements on the four haplotypes
  (mismatches, ambiguity, ctx:0 and discriminating overlap), the FR insert, base qualities, and the scorer's origin
  counts.

**Read pairs** are counted once, using the scorer's masses.

**Inconsistency flag, frozen.** A contributor with |Δ| ≥ 0.5 whose best proper-geometry placement has strictly fewer
mismatches on the candidate it scores against. Valid alternative placements, or single reads that favour the wrong
candidate, are not by themselves inconsistencies.

## 3. Prior attribution (current model only; no grid)

- **`prior-pool`** (the M3-A machinery) re-derives every saved pool prior (its gate) and gives the sequence-level priors
  of the pair as counterfactual rows. Its fixed M3-A model set and the path-level flank table are emitted by the tool
  but not used for conclusions; only `const_0.01`, the current model, is reported.
- **`prior_attrib.py`** is an independent forward/backward over the curated panel.
  - **Gate:** it reproduces score-set's sequence-level pair prior within 1e-9, which also shows that no alternative
    spelling contributes.
  - **Attribution, per homologue:** to donor compatibility (templates carrying the ctx:0 allele) and to copying
    transitions (continue into sv:1 vs redraw at the boundary; expected donor changes).
- **Parameter wording:** r = 0.01 is described as a copying-model donor-redraw parameter, not a biological
  recombination rate.

## 4. Decision rule (frozen)

| condition | decision |
|---|---|
| A reconstruction check, identity check or the prior gate fails, or 4 or more of the 2× top-10 contributors have no proper placement on either candidate | **EVIDENCE_INSUFFICIENT** |
| Otherwise, any flagged contributor | **LIKELIHOOD_INCONSISTENCY_TARGETED_TEST** (naming the read pairs) |
| Otherwise | **NO_INCONSISTENCY**: retained as an unresolved low-depth error |

The unrefined-pool-winner issue and the distance-dependent prior/flank questions are recorded separately, for later.

## Bounds

- 300 s counting every attempt, 4 threads, unchanged resource settings.
- **Declared timeouts:** score.2x 60, score.5x 90, spell 40, panel 40, prior.attrib 20, prior.pool 40,
  reads.detail 60 (s).
- Inputs are hashed in INPUTS.tsv and verified before every unit.
