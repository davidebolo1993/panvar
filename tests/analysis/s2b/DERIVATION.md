# S2-B: an admissible pre-placement bound for the implemented fragment likelihood

**What this is:** the target and the argument, written before the prototype. The target is the implemented scorer, defined by `enumerate_placements_exhaustive` and `log_placement_mass`, not the prose description. Every inequality below is checked executably in `s2b_preplacement_bound fixtures`.

## 0. The quantity

For a candidate {a, b} the implemented score is:

```
S(a,b) = −c·(E_a + E_b) + Σ_f ℓ_f( M_fa + M_fb ) + π(a,b)

  ℓ_f(m) = log[ c·m + η·P_bg(f) ]   c = (1−η)·λ
```

- **ℓ_f:** `log_fragment_contribution`. It is non-decreasing in m because c > 0. At η = 0 it is `log(c·m)` and refuses m = 0. Every fixture here and every real profile uses η > 0.
- **π(a,b):** the sequence-level prior, `log_haploid_sequence_prior(a) + (b) + 1[a ≠ b]·log 2`. It is a function of the model and the two sequences only.
- **E_h:** `exposure_for_length(|h|, …)`. It is a function of |h|, the library and the observed mate-length floors only.

So π, E and P_bg are exact before any placement. Only M_fh needs a bound. A bound U_fh ≥ M_fh for every cell gives:

```
UB(a,b) = −c·(E_a+E_b) + Σ_f ℓ_f(U_fa + U_fb) + π(a,b) ≥ S(a,b)
```

This needs only monotonicity of ℓ_f and of addition, with no chimera or per-fragment best-partner argument. That is the I5-A failure mode, and it does not arise because U is per (fragment, haplotype).

## 1. The placement mass, as implemented

From `enumerate_placements_exhaustive` and `log_placement_mass`, for haplotype h (length n) and fragment f with mates m1 and m2:

```
M_fh = Σ_{arm∈{0,1}} Σ_{s=0}^{n−1} Σ_{r=0}^{n} ½ · P(r−s) · Σ_{u∈C_F(s)} Σ_{v∈C_R(r)} 1[valid(s,r,u,v)] · e_u · e_v
```

| term | meaning |
|---|---|
| arm 0 | forward read F = m1, reverse read R = revcomp(m2) |
| arm 1 | forward read F = m2, reverse read R = revcomp(m1) |
| ½ | `kLogMateOrientationPrior` |
| P | the LIBRARY insert distribution (`build_library_insert_distribution`); the insert is r − s whatever the consumed spans |
| C_F(s) | the consumed-span components of `gapped_emission(h, F, s)`, left-anchored at s |
| C_R(r) | the components of `gapped_emission_ending_at(h, R, r)`, right-anchored at r |
| valid | r > s, r − s inside the fragment's support, and r − s ≥ max(consumed_u, consumed_v) |

**Relaxation 1 (drop validity, keep r − s ∈ library support).** All terms are non-negative, so:

```
M_fh ≤ Σ_arm Σ_s Σ_{r: r−s ∈ supp(P)} ½·P(r−s)·F(s)·R(r)

  F(s) = Σ_u e_u = exp(fit.log_emission)   (0 when gapped_emission returns no fit)
  R(r) likewise for the reverse read
```

The fragment support is a subset of the library support, and `P` is zero outside the library support.

**Relation to the production scorer.** The accelerated scorer computes a subset of the oracle's placements: those its pigeonhole seeding reaches, under its caps, which refuse rather than add. It never adds a placement outside the oracle. So M_prod ≤ M_oracle. This is checked on every fixture cell (production ≤ oracle + rounding), and on real data the bound is checked against production masses directly.

## 2. Emission bound without the read–window comparison

`gapped_emission` returns raw(x,s) / Z(L, a):
- raw is the banded affine-gap recursion over paths with at most B = `mismatch_band` edits and net displacement at most Δ = `max_gap_bp`;
- Z is `log_gap_shape_mass(L, a, B)`, the same recursion summed over all 4^L reads;
- a = min(n − s, L + Δ), the clamp under which Z is bit-identical, per the code's own argument.

**Decomposition by edits.** Split raw by the recursion's own edit index e (mismatch, insertion base and deletion base each cost one):

```
raw(x,s) = Σ_{e=0}^{B} raw_e(x,s)
```

**The ceiling W_e.** Define W_e(L, a) as the shape recursion with the same transitions, availability, floor and terminal rules, in which an M step takes either branch: match (weight 1−ε, cost 0) or mismatch (weight ε/3, cost 1). Insertions emit 0.25, as in the emitting pass.

**Claim:** raw_e(x,s) ≤ W_e(L, a) for every read x and window.

**Proof.** A path contributing to raw_e(x,s) has a fixed sequence of M/I/D steps, and each of its M steps is a match or a mismatch, determined by x and the window. The same step sequence with those branch choices is one term of W_e, with identical transition weights and identical emission weights (1−ε, ε/3, 0.25) and the same edit count. W_e also contains the other branch assignments, all non-negative. The availability bound is the same in both: paths never consume more than L + Δ bases, so clamping a changes no path. ∎

**Check.** Running the same recursion with mismatch weight ε (= 3·ε/3) reproduces Z, compared with `log_gap_shape_mass` in the fixtures.

**Consequence.** With an edit lower bound d(s) ≤ the minimum edits of any admissible alignment of x at s:

```
F(s) ≤ F̂(s) = Σ_{e=d(s)}^{B} W_e(L, a) / Z(L, a)        (and F(s) = 0 when d(s) > B)
```

## 3. The edit lower bound from exact pieces

Partition the read into p = B + 1 consecutive pieces (piece j at offset o_j, length ℓ_j ≥ 1). Take any admissible alignment A with e ≤ B edits, starting at s.
- **Broken pieces.** A piece is broken if it contains an edit: a mismatching read base, an inserted read base, or a deletion between two of its own consecutive bases. Each edit breaks at most one piece, so broken ≤ e.
- **Unbroken pieces.** An unbroken piece j is aligned by consecutive matching M steps at constant displacement k_j. So it occurs EXACTLY in h at s + o_j + k_j. The displacement is a cell index of the recursion, so |k_j| ≤ Δ. Each unit of displacement needs one gap base, so |k_j| ≤ e.
- **Threshold count.** Let t* = max_{j unbroken} |k_j|. Let c(t) be the number of pieces with an exact occurrence at s + o_j + k for some |k| ≤ t. Then c(t*) ≥ #unbroken ≥ p − e, so e ≥ max(p − c(t*), t*).
- **The bound.** Therefore:

  ```
  d(s) := min_{0 ≤ t ≤ min(Δ,B)} max(p − c(t), t)   ≤   minimum edits of any admissible alignment at s
  ```

**The right-anchored reverse read.** The argument is the same, mirrored: displacement is measured from r, and an unbroken piece occurs exactly at r − (L − o_j) − k.

**Completeness of the start set.** Every start s with an admissible alignment has c(min(Δ,B)) ≥ p − B ≥ 1, so it receives at least one exact piece occurrence within ±min(Δ,B) of its expected position.
- **What produces the start set.** An exact-occurrence index over the WHOLE haplotype (every position; no cap, no sampling, no external mapper) enumerates every such start.
- **Starts with no occurrence.** These have d ≥ p > B, so F = 0 there. This is not an assumption about hits: it follows from the pigeonhole step above.
- **The scorer's own caps.** `max_piece_hits` and `max_proposed_starts` can only remove placements from production. They are irrelevant to this bound, which covers the oracle.

**Inputs needed before placement:**
- the read bases;
- the haplotype bytes, through an exact substring index;
- the library (ε, gap parameters, Δ, B, insert distribution);
- the precomputed tables W_e and Z per (L, a).

**What is never run:** no `gapped_emission`, no Myers prefilter, no DP over a read–window pair, no pairing of accepted alignments.

## 4. The lazy-evidence variant (exact reuse)

The emission is a pure function of (read bytes, window bytes h[s, s + a)). For the right-anchored read the window is h[r − a, r).
- **What is reused.** For a window whose bytes equal a window already evaluated EXACTLY for the same read — S2-A's early incumbent haplotypes, whose evidence production has already computed — F(s) is that exact value, not F̂(s).
- **Why it stays admissible.** An exact value is its own upper bound.
- **What "already computed" means.** It must have been evaluated for the incumbent and kept. In the prototype the incumbent's exact emissions are recomputed and timed apart, as "already paid".

## 5. Numerical allowance and ties

UB is computed in floating point by an expression order different from production's S. So a candidate c is excluded against an incumbent i only if:

```
UB(c) + τ_c < S(i) − τ_i      τ_x = 64·(F + 8)·2^−53·A_x
```

- **A_x:** for the incumbent, the absolute term sum (`absolute_term_sum`); for c, the sum of |terms| of UB(c), plus |prior|.
- **Ties.** An exact tie, S(c) = S(i), can never be excluded, because UB(c) ≥ S(c) up to the same allowance.
- **Underflow.** U_fh = 0 exactly when no start qualifies; then M_fh = 0 exactly too, by the pigeonhole step. No underflow is involved: all terms are kept as doubles, and the smallest, (ε/3)^B·0.25^B·P, is about 1e−40.

## 6. What this bound does not use, and its known loose components

- **W_e sums every placement of e edits over the read.** The exact emission has one specific pattern. For e ≥ 1 this costs roughly log(L) nats per edit-bearing mate, reduced to log(piece length) only if broken pieces were used to localise edits (not done).
- **The validity constraint is dropped.** It matters only near the insert-support edges.
- **The full insert window is summed** over the qualifying ends.

These are measured, not assumed small.
