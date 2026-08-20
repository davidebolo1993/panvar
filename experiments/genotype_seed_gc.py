#!/usr/bin/env python3
"""The pre-registered seed analysis for docs/reports/genotype-seed-preregistration.md.

    log E[y_k] = alpha_b + beta * GC_k + f_b(u_k) + log(m_k)

Poisson, log link, IRLS. Free intercept per block. Block-specific natural cubic spline in normalized
truth position, 4 df, interior knots at the within-block quartiles and boundary knots at 0 and 1.
Truth multiplicity enters as a fixed offset, so dosage is divided out inside the model rather than in
the response. One beta per seed, fitted over the markers the sample actually carries.

Evidence is the distribution of the per-seed betas and nothing below it. No marker-level standard
error is computed or reported: markers covary through shared fragments, measured at roughly 9-fold for
anchors and 20-fold for markers, so any such interval would be invalid. That is also why the first
draft's stability rule -- comparing across-seed spread against a within-seed OLS standard error -- was
withdrawn as incoherent.

Seeds k = 0..4 are EXPLORATORY: they were analysed at length before this plan existed. The confirmatory
set is k = 5..54.

    genotype_seed_gc.py <seed-dir> [--exploratory]
"""
import sys
from pathlib import Path

import numpy as np
from scipy import stats

MIN_PER_BLOCK = 200
ARRAY_BLOCK = 13
BASE_SEED, SEED_STRIDE = 42, 466
CONFIRMATORY = range(5, 55)


def ns_basis(x, knots):
    """Natural cubic spline basis, truncated-power form, intercept dropped.

    K knots gives K-1 columns: a linear term plus K-2 constrained cubics. The intercept is omitted
    because every block already carries a free one.
    """
    k = np.asarray(knots, float)
    K = len(k)

    def d(j):
        return (np.clip(x - k[j], 0, None) ** 3 - np.clip(x - k[-1], 0, None) ** 3) / (k[-1] - k[j])

    dK1 = d(K - 2)
    return np.column_stack([x] + [d(j) - dK1 for j in range(K - 2)])


def poisson_irls(X, y, offset, iters=100, tol=1e-8, case_weights=None):
    """Poisson GLM by IRLS. Returns (coefficients, converged, rank).

    `case_weights` are OBSERVATION weights and enter the IRLS weight as `mu * w`, leaving X, y and the
    offset untouched. Scaling the design and the response by sqrt(w) instead -- which an earlier
    version did -- is not a weighted GLM at all; it fits a different model.

    Solves the weighted least squares directly on `sqrt(W) X` rather than forming the normal equations
    `X'WX`, which squares the condition number. On this design that spurious ill-conditioning was
    enough to stall one seed of fifty at the iteration cap; solving directly it converges in five, and
    the coefficient moves by 3e-11.

    The tolerance is 1e-8. At 1e-10 the criterion sits below the numerical noise floor of the solve, so
    every fit was flagged non-convergent while being stable to eight decimals from iteration 10.
    """
    beta = np.zeros(X.shape[1])
    eta = np.log(np.maximum(y, 0.5)) - offset
    cw = np.ones_like(y) if case_weights is None else np.asarray(case_weights, float)
    converged = False
    rank = 0
    for _ in range(iters):
        mu = np.exp(np.clip(eta + offset, -30, 30))
        z = eta + (y - mu) / np.maximum(mu, 1e-9)
        s = np.sqrt(np.maximum(mu * cw, 1e-12))
        new, _res, rank, _sv = np.linalg.lstsq(X * s.reshape(-1, 1), z * s, rcond=None)
        if np.max(np.abs(new - beta)) < tol:
            beta, converged = new, True
            break
        beta = new
        eta = X @ beta
    return beta, converged, int(rank)


def design(blk, gc, u, blocks):
    cols, ok = [gc.reshape(-1, 1)], []
    for b in blocks:
        m = blk == b
        cols.append(m.astype(float).reshape(-1, 1))
        ub = u[m]
        knots = list(np.maximum.accumulate([0.0, *np.quantile(ub, [0.25, 0.50, 0.75]), 1.0]))
        if len(set(knots)) < len(knots):
            ok.append((b, False))
            continue
        B = np.zeros((len(u), len(knots) - 1))
        B[m] = ns_basis(ub, knots)
        cols.append(B)
        ok.append((b, True))
    return np.hstack(cols), ok


def fit_seed(path):
    rows = [l.rstrip("\n").split("\t") for l in path.open()]
    c = {n: i for i, n in enumerate(rows[0])}
    body = rows[1:]
    blk = np.array([int(r[c["block"]]) for r in body])
    cls = np.array([r[c["class"]] for r in body])
    gc = np.array([float(r[c["gc"]]) for r in body])
    y = np.array([float(r[c["count"]]) for r in body])
    m = np.array([float(r[c["mult"]]) for r in body])
    u = np.array([float(r[c["upos"]]) for r in body])

    blocks = [b for b in np.unique(blk) if (blk == b).sum() >= MIN_PER_BLOCK]
    keep = np.isin(blk, blocks)
    blk, cls, gc, y, m, u = blk[keep], cls[keep], gc[keep], y[keep], m[keep], u[keep]
    off = np.log(np.maximum(m, 1e-9))

    X, _ = design(blk, gc, u, blocks)
    beta, conv, rank = poisson_irls(X, y, off)

    # Secondary: multiplicity-weighted refit. Same design, same response, same offset -- multiplicity
    # enters as an observation weight and nothing else.
    beta_w = poisson_irls(X, y, off, case_weights=m)[0][0]

    # Secondary: anchors alone, where multiplicity is 1 by construction.
    # Anchors alone, where multiplicity is 1 by construction. The row set must be restricted to the
    # blocks that are actually modelled: a row in a block with no intercept and no spline has its whole
    # count level pushed onto the shared GC coefficient. 336 such rows of 19,330 -- 1.7 percent -- moved
    # this estimate from +0.037 to +0.461, a factor of twelve, and made it look like a large effect
    # reproducing in every seed.
    anc = cls == "anchor"
    ablocks = [b for b in np.unique(blk[anc]) if (blk[anc] == b).sum() >= MIN_PER_BLOCK]
    beta_a = np.nan
    if ablocks:
        sel = anc & np.isin(blk, ablocks)
        Xa, _ = design(blk[sel], gc[sel], u[sel], ablocks)
        beta_a = poisson_irls(Xa, y[sel], off[sel])[0][0]

    # Secondary: per-block slopes, fitted separately.
    per_block = {}
    for b in blocks:
        mk = blk == b
        Xb, _ = design(blk[mk], gc[mk], u[mk], [b])
        per_block[int(b)] = poisson_irls(Xb, y[mk], off[mk])[0][0]

    arr = (blk == ARRAY_BLOCK) & (cls == "informative")
    cpc = y / m
    return dict(
        beta=beta[0], beta_w=beta_w, beta_anchor=beta_a, per_block=per_block,
        converged=conv, rank=rank, ncol=X.shape[1], n=len(y), n_blocks=len(blocks),
        gc_delta=(np.average(gc[arr], weights=m[arr]) - gc[anc].mean()) if arr.any() and anc.any() else np.nan,
        arr_mean=cpc[arr].mean() if arr.any() else np.nan,
        arr_median=float(np.median(cpc[arr])) if arr.any() else np.nan,
        anc_mean=cpc[anc].mean() if anc.any() else np.nan)


def main():
    d = Path(sys.argv[1] if len(sys.argv) > 1 else "seeds")
    want_expl = "--exploratory" in sys.argv
    files = []
    for p in sorted(d.glob("m_*.tsv"), key=lambda q: int(q.stem.split("_")[1])):
        k = (int(p.stem.split("_")[1]) - BASE_SEED) // SEED_STRIDE
        if (k in CONFIRMATORY) != want_expl:
            files.append((k, p))
    if not files:
        sys.exit(f"no {'exploratory' if want_expl else 'confirmatory'} seed dumps in {d}")
    if not want_expl:
        missing = sorted(set(CONFIRMATORY) - {k for k, _ in files})
        if missing:
            sys.exit(f"confirmatory set incomplete: missing k={missing}. The decision rule is defined "
                     f"over exactly k={CONFIRMATORY.start}..{CONFIRMATORY.stop - 1}; a partial set "
                     f"would let the answer depend on which seeds happened to finish.")

    res = [(k, fit_seed(p)) for k, p in files]
    label = "EXPLORATORY (excluded from the decision rule)" if want_expl else "CONFIRMATORY"
    print(f"{label}: {len(res)} seeds\n")

    bad = [k for k, r in res if not r["converged"] or r["rank"] < r["ncol"]]
    if bad:
        print(f"  WARNING: {len(bad)} seeds failed to converge or were rank-deficient: {bad}\n")

    print("per-seed slope")
    for k, r in res:
        print(f"  k={k:<3d} beta {r['beta']:+.5f}  weighted {r['beta_w']:+.5f}  "
              f"anchors {r['beta_anchor']:+.5f}  rank {r['rank']}/{r['ncol']}"
              f"{'' if r['converged'] else '  NOT CONVERGED'}")

    b = np.array([r["beta"] for _, r in res])
    n, sd = len(b), np.std([r["beta"] for _, r in res], ddof=1)
    lo, hi = stats.t.interval(0.95, n - 1, loc=b.mean(), scale=sd / np.sqrt(n))
    print("\nPRIMARY, statistical -- the seed-level slopes and nothing below them")
    print(f"  mean beta       {b.mean():+.5f}")
    print(f"  across-seed sd  {sd:.5f}")
    print(f"  95% t interval  [{lo:+.5f}, {hi:+.5f}]")
    print(f"  positive in     {int((b > 0).sum())} of {n} seeds")
    if want_expl:
        print("  -> exploratory, no verdict")
    else:
        rep = (lo > 0 or hi < 0) and (b > 0).sum() >= 45
        print(f"  -> {'REPRODUCIBLE positive effect' if rep else 'no reproducible positive effect demonstrated in this wgsim experiment'}")

    gd = np.nanmean([r["gc_delta"] for _, r in res])
    lam = np.nanmean([r["anc_mean"] for _, r in res])
    dlam = lam * (np.exp(b.mean() * gd) - 1.0)
    print("\nCO-PRIMARY, practical -- lambda_anchor * (exp(beta * delta_GC) - 1), log link")
    print(f"  multiplicity-weighted GC delta at block {ARRAY_BLOCK}  {gd:+.5f}")
    print(f"  anchor lambda                              {lam:.5f}")
    print(f"  predicted lambda bias                      {dlam:+.5f}")
    print(f"  ratio to the 0.025 half-width              {abs(dlam)/0.025:.2f}   (no verdict attached)")

    am = np.nanmean([r["arr_mean"] for _, r in res])
    amed = np.nanmean([r["arr_median"] for _, r in res])
    print("\nSECONDARY -- reported without thresholds")
    print(f"  anchors count_per_copy   {lam:.4f}")
    print(f"  block {ARRAY_BLOCK} mean            {am:.4f}   excess {100*(am/lam-1):+.2f}%")
    print(f"  block {ARRAY_BLOCK} median          {amed:.4f}   excess {100*(amed/lam-1):+.2f}%")
    print(f"  weighted-refit mean beta {np.mean([r['beta_w'] for _, r in res]):+.5f}")
    print(f"  anchor-only  mean beta   {np.nanmean([r['beta_anchor'] for _, r in res]):+.5f}")
    allb = [(kb, v) for _, r in res for kb, v in r["per_block"].items()]
    print(f"  per-block fits positive  {sum(1 for _, v in allb if v > 0)} of {len(allb)}")

    # The 50 slopes are the evidence, so they belong in a file rather than only in a console line.
    out = d / ("seed_slopes_exploratory.tsv" if want_expl else "seed_slopes.tsv")
    with out.open("w") as f:
        f.write("k\tbeta\tbeta_weighted\tbeta_anchor\tgc_delta\tanchor_lambda"
                "\tarray_mean_cpc\tarray_median_cpc\tconverged\trank\tncol\tn\n")
        for k, r in res:
            f.write(f"{k}\t{r['beta']:.8f}\t{r['beta_w']:.8f}\t{r['beta_anchor']:.8f}\t"
                    f"{r['gc_delta']:.8f}\t{r['anc_mean']:.6f}\t{r['arr_mean']:.6f}\t"
                    f"{r['arr_median']:.6f}\t{int(r['converged'])}\t{r['rank']}\t{r['ncol']}\t{r['n']}\n")
    print(f"\nper-seed values written to {out}")


if __name__ == "__main__":
    main()
