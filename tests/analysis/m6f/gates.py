#!/usr/bin/env python3
"""M6-F gates (frozen with DESIGN.md). Fixture gates G1-G4 and G6 use tiny enumerable libraries; G5 and G7 run the
training-only fitting on the real folds (no evaluation copy is scored).

  G1  normalisation: over all 4^M sequences each model's probabilities sum to 1, for every grid point (incl. monomorphic
      positions handled by the compression)
  G2  compressed forward == exhaustive hidden-path enumeration over every scored position (uncompressed); and
      compressed == uncompressed forward on a larger fixture
  G3  rho = 0 reproduces C1
  G4  redraw at every scored site (r = 1) reproduces C0
  G5  held-out isolation: for each evaluation family, the fitted library, mask, variable sites, inner validation totals
      and selected parameters are identical when that family's copies are replaced by scrambled synthetic components
      with arbitrary weights
  G6  representation-only aliases: splitting a donor into two identical donors with half its weight preserves every
      model's probability (fixture and one real library)
  G7  all models score exactly the same observations (every model call within one score_all sees one digest)

  gates.py <prepare dir> <events tsv> <out prefix>
"""
import builtins
import itertools
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
prep, events_path, prefix = sys.argv[1:4]
import m6f_core as core  # noqa: E402

ALLOW = {os.path.abspath(p) for p in [core.M1D_FA, events_path, f"{prep}/copies.tsv", f"{prep}/prepare.tsv", f"{prep}/families.tsv"]}
_open = builtins.open
OUTDIR = os.path.dirname(os.path.abspath(prefix))


def guarded(path, mode="r", *a, **k):
    p = os.path.abspath(path)
    if "w" in mode or "a" in mode:
        assert p.startswith(OUTDIR), f"write outside the unit directory: {p}"
    else:
        assert p in ALLOW or p.startswith(OUTDIR), f"not an allowlisted input: {p}"
    return _open(path, mode, *a, **k)


builtins.open = guarded
rows, failed = [], 0


def row(g, item, value, expect, ok):
    global failed
    rows.append(f"{g}\t{item}\t{value:.3g}\t{expect}\t{'PASS' if ok else 'FAIL'}")
    failed += 0 if ok else 1


rng = np.random.default_rng(20261020)
# fixture: K = 3 donors, M = 5 positions (two training-monomorphic), irregular distances
pos = np.array([3, 40, 41, 700, 2500])
B = np.array([[0, 1, 2, 3, 0], [0, 1, 0, 3, 1], [2, 1, 2, 3, 2]], dtype=np.uint8)
w = np.array([0.5, 0.3, 0.2])
rho_g = np.array([0.0, 1e-4, 1e-2, 5.0])

# G1 normalisation
tot = {"C0": np.zeros(len(core.EPS)), "C1": np.zeros(len(core.EPS)), "C2": np.zeros((len(core.EPS), len(rho_g)))}
for y in itertools.product(range(4), repeat=len(pos)):
    o = core.Obs(pos, B, np.array(y, dtype=np.uint8), w)
    tot["C0"] += np.exp(core.score_c0(o, core.EPS))
    tot["C1"] += np.exp(core.score_c1(o, core.EPS))
    tot["C2"] += np.exp(core.score_c2(o, core.EPS, rho_g))
for m, t in tot.items():
    row("G1_normalisation", f"{m}: max |sum_y P(y) - 1| over 4^5 sequences and the grid", float(np.max(np.abs(t - 1))), "<= 1e-12", float(np.max(np.abs(t - 1))) <= 1e-12)


def c2_full(pos, B, y, w, eps, rho):   # uncompressed forward over every scored position
    w = w / w.sum()
    e = lambda i: np.where(B[:, i] == y[i], 1 - eps, eps / 3)
    a = w * e(0)
    for i in range(1, len(pos)):
        r = -np.expm1(-rho * (pos[i] - pos[i - 1]))
        a = ((1 - r) * a + r * w * a.sum()) * e(i)
    return np.log(a.sum())


def c2_paths(pos, B, y, w, eps, rho):   # exhaustive hidden-path enumeration
    w = w / w.sum()
    K, total = len(w), 0.0
    for path in itertools.product(range(K), repeat=len(pos)):
        p = w[path[0]] * (1 - eps if B[path[0], 0] == y[0] else eps / 3)
        for i in range(1, len(pos)):
            r = -np.expm1(-rho * (pos[i] - pos[i - 1]))
            p *= ((1 - r) if path[i] == path[i - 1] else 0.0) + r * w[path[i]]
            p *= 1 - eps if B[path[i], i] == y[i] else eps / 3
        total += p
    return np.log(total)


worst = 0.0
for _ in range(20):
    y = rng.integers(0, 4, len(pos)).astype(np.uint8)
    o = core.Obs(pos, B, y, w)
    for ie, eps in enumerate(core.EPS):
        for ir, rho in enumerate(rho_g):
            worst = max(worst, abs(core.score_c2(o, [eps], [rho])[0, 0] - c2_paths(pos, B, y, w, eps, rho)))
row("G2_forward", "compressed forward - exhaustive hidden paths (3^5 paths, 20 sequences x grid)", worst, "<= 1e-10", worst <= 1e-10)
pos2 = np.sort(rng.choice(30000, 60, replace=False))
B2 = rng.integers(0, 4, (4, 60)).astype(np.uint8)
B2[:, ::3] = B2[0, ::3]   # every third position training-monomorphic
w2 = rng.random(4)
worst = 0.0
for _ in range(10):
    y = B2[rng.integers(0, 4)].copy()
    flip = rng.random(60) < 0.1
    y[flip] = rng.integers(0, 4, flip.sum())
    o = core.Obs(pos2, B2, y, w2)
    for eps in core.EPS:
        for rho in core.RHO:
            worst = max(worst, abs(core.score_c2(o, [eps], [rho])[0, 0] - c2_full(pos2, B2, y, w2, eps, rho)))
row("G2_forward", "compressed - uncompressed forward (60 positions, 20 monomorphic, 10 sequences x grid)", worst, "<= 1e-9", worst <= 1e-9)

# G3 / G4
worst3 = worst4 = 0.0
for _ in range(20):
    y = rng.integers(0, 4, len(pos2)).astype(np.uint8)
    o = core.Obs(pos2, B2, y, w2)
    worst3 = max(worst3, float(np.max(np.abs(core.score_c2(o, core.EPS, [0.0])[:, 0] - core.score_c1(o, core.EPS)))))
    worst4 = max(worst4, float(np.max(np.abs(core.score_c2(o, core.EPS, [1e-4], redraw_all=True)[:, 0] - core.score_c0(o, core.EPS)))))
row("G3_zero_switching", "C2(rho = 0) - C1", worst3, "<= 1e-12", worst3 <= 1e-12)
row("G4_redraw_every_site", "C2(r = 1 at every scored site) - C0", worst4, "<= 1e-12", worst4 <= 1e-12)

# G6 aliases (fixture)
worst = 0.0
for _ in range(20):
    y = rng.integers(0, 4, len(pos2)).astype(np.uint8)
    a, b = core.Obs(pos2, B2, y, w2), core.Obs(pos2, np.vstack([B2, B2[:1]]), y, np.r_[w2[0] / 2, w2[1:], w2[0] / 2])
    for f in (lambda o: core.score_c0(o, core.EPS), lambda o: core.score_c1(o, core.EPS), lambda o: core.score_c2(o, core.EPS, core.RHO).ravel()):
        worst = max(worst, float(np.max(np.abs(f(a) - f(b)))))
row("G6_aliases", "fixture: donor split into two identical halves, all models and grid", worst, "<= 1e-12", worst <= 1e-12)

# real data: G5, G6 (real library), G7
copies, refs, comps, families, evaluation = core.load(prep, events_path)
strata = sorted(refs)


def summary(fit):
    if fit is None:
        return None
    return (tuple(fit["library"].items()), core.hashlib.sha256(fit["mask"].tobytes()).hexdigest(), fit["inner_targets"],
            tuple(np.round(np.concatenate([fit["validation"][m].ravel() for m in ("C0", "C1", "C2")]), 9)),
            tuple((m, float(v[0]), None if v[1] is None else float(v[1])) for m, v in sorted(fit["selected"].items())))


iso_ok, n_folds = True, 0
digests_ok, calls = True, 0
seen = []


def recorder(obs, s):
    global digests_ok, calls
    calls += 1
    digests_ok &= s["digest"] == obs.digest


orig = {"c0": core.score_c0, "c1": core.score_c1, "c2": core.score_c2}
log_digest = []
core.score_c0 = lambda o, *a, **k: (log_digest.append(("C0", o.digest)), orig["c0"](o, *a, **k))[1]
core.score_c1 = lambda o, *a, **k: (log_digest.append(("C1", o.digest)), orig["c1"](o, *a, **k))[1]
core.score_c2 = lambda o, *a, **k: (log_digest.append(("C2", o.digest)), orig["c2"](o, *a, **k))[1]
for fam in evaluation:
    for s in strata:
        if not any(c["stratum"] == s and c["family"] == fam for c in copies):
            continue
        a = summary(core.fit_fold(copies, comps[s], s, fam, recorder))
        # perturbed: the held-out family's copies point to scrambled synthetic components with arbitrary weights
        pert_comps = dict(comps[s])
        pert = []
        for i, c in enumerate(copies):
            if c["stratum"] == s and c["family"] == fam:
                key = f"SCRAMBLED_{i}"
                src = comps[s][c["unique_record"]]
                fake = core.Component([], "A" * len(src.base))
                fake.base = rng.integers(0, 4, len(src.base)).astype(np.uint8)
                fake.bad = rng.random(len(src.base)) < 0.3
                pert_comps[key] = fake
                pert.append(dict(c, unique_record=key, weight=str(rng.random() * 5)))
            else:
                pert.append(c)
        b = summary(core.fit_fold(pert, pert_comps, s, fam))
        iso_ok &= a == b
        n_folds += 1
row("G5_heldout_isolation", f"fold fits identical with the held-out family scrambled ({n_folds} family x stratum folds)", 0.0 if iso_ok else 1.0, "identical", iso_ok)
# G7: within each score_all call, the three model calls carry one digest
triples_ok = all(len({d for _, d in log_digest[i:i + 3]}) == 1 and [m for m, _ in log_digest[i:i + 3]] == ["C0", "C1", "C2"] for i in range(0, len(log_digest), 3))
row("G7_same_observations", f"every score_all call: C0, C1, C2 see one observation digest ({calls} calls)", 0.0 if (triples_ok and digests_ok and calls > 0) else 1.0,
    "identical", triples_ok and digests_ok and calls > 0 and len(log_digest) % 3 == 0)
core.score_c0, core.score_c1, core.score_c2 = orig["c0"], orig["c1"], orig["c2"]
# G6 on one real library (training-only fold of the first evaluation family, first stratum it has)
s = next(s for s in strata if any(c["stratum"] == s and c["family"] == evaluation[0] for c in copies))
lib = core.library(copies, s, {evaluation[0]})
mask = core.mask_of(lib, comps[s])
g = sorted({c["family"] for c in copies if c["stratum"] == s and c["family"] != evaluation[0] and float(c["weight"]) > 0})[0]
target = next(c["unique_record"] for c in copies if c["stratum"] == s and c["family"] == g)
lib_g = core.library(copies, s, {evaluation[0], g})
mask_g = core.mask_of(lib_g, comps[s])
k0 = next(iter(lib_g))
alias_comps = dict(comps[s])
alias_comps["ALIAS"] = comps[s][k0]
lib_alias = dict(lib_g)
lib_alias[k0] /= 2
lib_alias["ALIAS"] = lib_g[k0] / 2
o1 = core.observation(lib_g, comps[s], mask_g, target)
o2 = core.observation(lib_alias, alias_comps, mask_g, target)
d = max(float(np.max(np.abs(core.score_c0(o1, core.EPS) - core.score_c0(o2, core.EPS)))), float(np.max(np.abs(core.score_c1(o1, core.EPS) - core.score_c1(o2, core.EPS)))),
        float(np.max(np.abs(core.score_c2(o1, core.EPS, core.RHO) - core.score_c2(o2, core.EPS, core.RHO)))))
row("G6_aliases", f"real library ({s}, training-only target): component split into two aliases", d, "<= 1e-9", d <= 1e-9)

with open(f"{prefix}.gates.tsv", "w") as f:
    f.write("gate\titem\tvalue\texpectation\tresult\n" + "\n".join(rows) + f"\nm6f_gates\t-\t{failed}\t0 failures\t{'PASS' if failed == 0 else 'FAIL'}\n")
print(open(f"{prefix}.gates.tsv").read())
sys.exit(1 if failed else 0)
