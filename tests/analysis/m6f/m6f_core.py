"""M6-F core: three normalised predictive models of a held-out component's bases at declared reference positions.

Observation: the held-out copy's bases y_x at scored positions x (reference coordinates of its stratum), where
scored = training mask (simple correspondence in every positive-weight training component) AND the held-out copy's own
simple correspondence. Library: unique training components with summed positive family-aware weights, normalised.
Residual emission (shared by all models): e(y | b) = 1 - eps if y == b else eps / 3.

  C0  independent sites:      prod_x sum_k w_k e(y_x | b_k(x))
  C1  whole-component mixture: sum_k w_k prod_x e(y_x | b_k(x))
  C2  Li-Stephens copying:     hidden donor path; between consecutive scored positions at distance d the donor is
                               redrawn from w with probability r = 1 - exp(-rho d), else kept
Exactness of the compression: at a training-monomorphic scored position every donor carries the same base, so the
emission is the same for every hidden donor and factors out of every model; the redraw kernels compose
(K(d1) K(d2) = K(d1 + d2)), so only training-variable positions enter the forward pass, with their reference distances.
rho = 0 gives C1; r = 1 at every scored site gives C0.
"""
import hashlib
from collections import defaultdict

import numpy as np

EPS = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2])
RHO = np.array([1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3])
W_BP, LARGE_BP = 10, 50
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


class Component:
    """One unique component aligned to its stratum reference: base code per reference position (4 = absent) and the
    positions without simple correspondence (absent, or within W of an indel / large-deletion edge)."""

    def __init__(self, events, ref):
        L = len(ref)
        self.base = np.array([CODE[c] for c in ref], dtype=np.uint8)
        self.bad = np.zeros(L, dtype=bool)
        for e in events:
            p, n = int(e["ref_pos"]), int(e["ref_len"])
            if e["event"] == "SNV":
                self.base[p] = CODE[e["query_bases"]]
            elif e["event"] == "DEL":
                self.base[p:p + n] = 4
                self.bad[p:p + n] = True
                self.bad[max(0, p - W_BP):min(L, p + W_BP)] = True
                self.bad[max(0, p + n - W_BP):min(L, p + n + W_BP)] = True
            else:   # INS before p
                self.bad[max(0, p - W_BP):min(L, p + W_BP)] = True


class Obs:
    """The scored observation of one held-out copy; every model receives this same object."""

    def __init__(self, pos, B, y, w):
        # pos: scored reference positions (sorted); B: (K, n) library bases at pos; y: (n,) held-out bases; w: (K,) weights
        w = np.asarray(w, dtype=float)
        self.w = w / w.sum()
        var = np.any(B != B[:1], axis=0) if B.shape[1] else np.zeros(0, dtype=bool)
        self.var_pos = pos[var]
        self.Bv = B[:, var]
        self.yv = y[var]
        mono = ~var
        self.mono_match = int(np.sum(y[mono] == B[0, mono])) if B.shape[1] else 0
        self.mono_mismatch = int(np.sum(mono)) - self.mono_match
        self.n_scored = int(len(pos))
        self.digest = hashlib.sha256(np.ascontiguousarray(pos).tobytes() + np.ascontiguousarray(y).tobytes()).hexdigest()


def mono_term(obs, eps):
    return obs.mono_match * np.log1p(-eps) + obs.mono_mismatch * np.log(eps / 3.0)


def emissions(obs, i, eps):   # (P, K)
    match = obs.Bv[:, i] == obs.yv[i]
    return np.where(match[None, :], (1.0 - eps)[:, None], (eps / 3.0)[:, None])


def score_c0(obs, eps):
    eps = np.asarray(eps, dtype=float)
    ll = mono_term(obs, eps)
    for i in range(len(obs.var_pos)):
        ll = ll + np.log(emissions(obs, i, eps) @ obs.w)
    return ll


def score_c1(obs, eps):
    eps = np.asarray(eps, dtype=float)
    acc = np.zeros((len(eps), len(obs.w)))
    for i in range(len(obs.var_pos)):
        acc += np.log(emissions(obs, i, eps))
    a = acc + np.log(obs.w)[None, :]
    m = a.max(axis=1)
    return mono_term(obs, eps) + m + np.log(np.exp(a - m[:, None]).sum(axis=1))


def score_c2(obs, eps, rho, redraw_all=False):
    """Returns (len(eps), len(rho)). redraw_all: r = 1 between every pair of scored positions (gate: reproduces C0)."""
    eps = np.asarray(eps, dtype=float)
    rho = np.asarray(rho, dtype=float)
    E, R = np.meshgrid(eps, rho, indexing="ij")
    er, rr = E.ravel(), R.ravel()
    ll = mono_term(obs, er)
    n = len(obs.var_pos)
    if n == 0:
        return ll.reshape(len(eps), len(rho))
    alpha = emissions(obs, 0, er) * obs.w[None, :]
    c = alpha.sum(axis=1)
    ll = ll + np.log(c)
    alpha /= c[:, None]
    for i in range(1, n):
        if redraw_all:
            r = np.ones_like(rr)
        else:
            r = -np.expm1(-rr * float(obs.var_pos[i] - obs.var_pos[i - 1]))
        alpha = alpha * (1.0 - r)[:, None] + r[:, None] * obs.w[None, :]
        alpha *= emissions(obs, i, er)
        c = alpha.sum(axis=1)
        ll = ll + np.log(c)
        alpha /= c[:, None]
    return ll.reshape(len(eps), len(rho))


def score_all(obs):
    """The single entry point: all three models score the same observation object over their grids."""
    return {"C0": score_c0(obs, EPS), "C1": score_c1(obs, EPS), "C2": score_c2(obs, EPS, RHO), "digest": obs.digest}


# ---- training-only fitting -------------------------------------------------------------------------------------------
def library(copies, stratum, exclude):
    """Unique positive-weight training components of a stratum (families in `exclude` removed) with summed weights."""
    lib = defaultdict(float)
    for c in copies:
        if c["stratum"] == stratum and c["family"] not in exclude and float(c["weight"]) > 0:
            lib[c["unique_record"]] += float(c["weight"])
    return dict(sorted(lib.items()))


def mask_of(lib, comps):
    bad = np.zeros_like(next(iter(comps.values())).bad)
    for k in lib:
        bad |= comps[k].bad
    return ~bad


def observation(lib, comps, mask, held):
    h = comps[held]
    scored = mask & ~h.bad
    pos = np.flatnonzero(scored)
    keys = list(lib)
    B = np.stack([comps[k].base[pos] for k in keys]) if keys else np.zeros((0, len(pos)), dtype=np.uint8)
    obs = Obs(pos, B, h.base[pos], [lib[k] for k in keys])
    L = len(mask)
    obs.excluded_training_mask = int(L - mask.sum())
    obs.excluded_heldout_alignment = int(np.sum(mask & h.bad))
    obs.novel_at_monomorphic = obs.mono_mismatch
    return obs


def pick(scores, shape):
    """Argmax over the grid; ties (within 1e-9) go to smaller eps, then smaller rho (grids are ascending)."""
    flat = scores.ravel()
    best = flat.max()
    i = int(np.flatnonzero(flat >= best - 1e-9)[0])
    return np.unravel_index(i, shape)


def fit_fold(copies, comps, stratum, held_family, recorder=None):
    """Select each model's parameters by inner leave-one-family-out validation over the training families only.
    Returns (library, mask, selected parameters, validation totals). The held-out family is excluded throughout."""
    lib = library(copies, stratum, {held_family})
    if not lib:
        return None
    mask = mask_of(lib, comps)
    inner = sorted({c["family"] for c in copies if c["stratum"] == stratum and c["family"] != held_family and float(c["weight"]) > 0})
    tot = {"C0": np.zeros(len(EPS)), "C1": np.zeros(len(EPS)), "C2": np.zeros((len(EPS), len(RHO)))}
    n_targets = 0
    for g in inner:
        lib_g = library(copies, stratum, {held_family, g})
        if not lib_g:
            continue
        mask_g = mask_of(lib_g, comps)
        for c in copies:
            if c["stratum"] == stratum and c["family"] == g and float(c["weight"]) > 0:
                obs = observation(lib_g, comps, mask_g, c["unique_record"])
                s = score_all(obs)
                if recorder is not None:
                    recorder(obs, s)
                for m in tot:
                    tot[m] += float(c["weight"]) * s[m]
                n_targets += 1
    sel = {"C0": (EPS[pick(tot["C0"], tot["C0"].shape)[0]], None), "C1": (EPS[pick(tot["C1"], tot["C1"].shape)[0]], None)}
    i, j = pick(tot["C2"], tot["C2"].shape)
    sel["C2"] = (EPS[i], RHO[j])
    return {"library": lib, "mask": mask, "selected": sel, "validation": tot, "inner_targets": n_targets}


def score_selected(obs, sel):
    return {"C0": float(score_c0(obs, [sel["C0"][0]])[0]), "C1": float(score_c1(obs, [sel["C1"][0]])[0]),
            "C2": float(score_c2(obs, [sel["C2"][0]], [sel["C2"][1]])[0, 0]), "digest": obs.digest}


# ---- loading (callers install their own input guards) ----------------------------------------------------------------
M1D_FA = "/Users/davide.bolognini/panvar/tests/results/m1d_unit_reachability_20261019/runs/inventory.a1/inventory.fa"


def fasta(path):
    recs, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs[name] = []
        elif name:
            recs[name].append(line)
    return {k: "".join(v) for k, v in recs.items()}


def load(prep, events_path):
    import csv
    copies = list(csv.DictReader(open(f"{prep}/copies.tsv"), delimiter="\t"))
    refs = {l.split("\t")[1]: l.rstrip("\n").split("\t")[2] for l in open(f"{prep}/prepare.tsv") if l.startswith("reference\t")}
    seqs = fasta(M1D_FA)
    ev = defaultdict(list)
    for r in csv.DictReader(open(events_path), delimiter="\t"):
        if r["event"] != "SUMMARY":
            ev[r["label"]].append(r)
    comps = {}
    for s, ref in refs.items():
        units = sorted({c["unique_record"] for c in copies if c["stratum"] == s})
        comps[s] = {u: Component(ev[f"{s}:{u}"], seqs[ref]) for u in units}
    families = list(csv.DictReader(open(f"{prep}/families.tsv"), delimiter="\t"))
    evaluation = [l.rstrip("\n").split("\t")[1].split(",") for l in open(f"{prep}/prepare.tsv") if l.startswith("evaluation_families\t")][0]
    return copies, refs, comps, families, evaluation
