#!/usr/bin/env python3
"""Independent allele-pair oracle for panvar's multiplicity emission.

This is deliberately outside the production genotyper.  It answers four ordered
questions for one block:

  1. Is a certified sequence-nearest pair distinguishable in the marker map?
  2. Does the likelihood recover that pair from its own noiseless multiplicities?
  3. Which panel pair best projects the held-out truth multiplicities?
  4. What changes when the same likelihood sees the observed read counts?

It also reproduces the default 50/50 detected-marker/length candidate pruning.
The HMM, linkage, GQ, ESS discount, and all non-default scoring flags are absent
on purpose.  The ESS factor and full-universe baseline are constant within a
block and therefore do not change pair ranks under the default model.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


MASK64 = (1 << 64) - 1


@dataclass(frozen=True)
class Scenario:
    name: str
    observed: np.ndarray
    detected: np.ndarray


def pair_arg(text: str) -> tuple[int, int]:
    fields = text.split(",")
    if len(fields) != 2:
        raise argparse.ArgumentTypeError("pair must be ALLELE1,ALLELE2")
    a, b = (int(x) for x in fields)
    return (a, b) if a <= b else (b, a)


def read_conf(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    slots: set[int] = set()
    max_allele = -1
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"allele", "slot", "mult", "obs"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: expected columns {sorted(required)}")
        for row in reader:
            allele = int(row["allele"])
            slot = int(row["slot"])
            mult = int(row["mult"])
            obs = float(row["obs"])
            if mult <= 0:
                raise ValueError(f"{path}: non-positive multiplicity for allele {allele}, slot {slot}")
            del obs
            slots.add(slot)
            max_allele = max(max_allele, allele)
    if max_allele < 0 or not slots:
        raise ValueError(f"{path}: no allele markers")

    ordered_slots = np.asarray(sorted(slots), dtype=np.int64)
    slot_to_col = {int(slot): col for col, slot in enumerate(ordered_slots)}
    multiplicity = np.zeros((max_allele + 1, len(ordered_slots)), dtype=np.int16)
    observed = np.full(len(ordered_slots), np.nan, dtype=np.float64)
    # Stream a second time instead of retaining one Python tuple per sparse entry.  The all-kmer KIV
    # diagnostic has 5.3 million rows but only a 457 x 52,492 dense matrix; retaining the rows costs
    # an order of magnitude more memory than the data being analysed.
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            allele = int(row["allele"])
            slot = int(row["slot"])
            mult = int(row["mult"])
            obs = float(row["obs"])
            col = slot_to_col[slot]
            if multiplicity[allele, col] != 0:
                raise ValueError(f"{path}: duplicate allele/slot row {allele}/{slot}")
            multiplicity[allele, col] = mult
            if not math.isnan(observed[col]) and observed[col] != obs:
                raise ValueError(f"{path}: inconsistent observed count for slot {slot}")
            observed[col] = obs
    if np.isnan(observed).any():
        raise AssertionError("internal error: an indexed slot has no observation")
    return multiplicity, observed, ordered_slots


def read_marker_dump(
    path: Path, block: int, ordered_slots: np.ndarray
) -> tuple[np.ndarray, np.ndarray, list[float]]:
    slot_to_col = {int(slot): col for col, slot in enumerate(ordered_slots)}
    truth = np.full(len(ordered_slots), np.nan, dtype=np.float64)
    clumps = np.full(len(ordered_slots), -1, dtype=np.int64)
    anchors: list[float] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"block_index", "marker_class", "slot", "count", "truth_mult", "clump"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: expected columns {sorted(required)}")
        for row in reader:
            if row["marker_class"] == "anchor":
                anchors.append(float(row["count"]))
            if int(row["block_index"]) != block or row["marker_class"] != "informative":
                continue
            slot = int(row["slot"])
            col = slot_to_col.get(slot)
            if col is None:
                continue
            value = row["truth_mult"]
            if value in {"", ".", "NA", "nan"}:
                continue
            truth[col] = float(value)
            clumps[col] = int(row["clump"])
    missing = int(np.isnan(truth).sum())
    if missing:
        raise ValueError(f"{path}: truth multiplicity is missing for {missing} block-{block} markers")
    if np.any(clumps < 0):
        raise ValueError(f"{path}: clump is missing for at least one block-{block} marker")
    return truth, clumps, anchors


def read_lengths(path: Path, n_alleles: int) -> np.ndarray:
    lengths = np.zeros(n_alleles, dtype=np.int64)
    seen = np.zeros(n_alleles, dtype=np.bool_)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"allele", "bp"}.issubset(reader.fieldnames):
            raise ValueError(f"{path}: expected allele and bp columns")
        for row in reader:
            allele = int(row["allele"])
            if 0 <= allele < n_alleles:
                lengths[allele] = int(row["bp"])
                seen[allele] = True
    if not seen.all():
        raise ValueError(f"{path}: missing lengths for {int((~seen).sum())} alleles")
    return lengths


def read_depth(path: Path, block: int) -> tuple[float, float]:
    usable: list[float] = []
    block_lambda: float | None = None
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {"block_index", "lambda_hap", "usable"}.issubset(reader.fieldnames):
            raise ValueError(f"{path}: expected block_index, lambda_hap, usable columns")
        for row in reader:
            value = float(row["lambda_hap"])
            if int(float(row["usable"])) != 0:
                usable.append(value)
            if int(row["block_index"]) == block:
                block_lambda = value
    if block_lambda is None:
        raise ValueError(f"{path}: block {block} is absent")
    if not usable:
        raise ValueError(f"{path}: no usable depths")
    return block_lambda, sum(usable) / len(usable)


def estimate_phi(anchor_counts: Iterable[float]) -> float:
    values = np.asarray(list(anchor_counts), dtype=np.float64)
    if values.size < 20:
        return 0.0
    mean = float(values.mean())
    variance = float(values.var(ddof=1))
    if mean <= 0.0 or variance <= mean:
        return 0.0
    return mean * mean / (variance - mean)


def log_nb(x: float, mean: float, phi: float) -> float:
    mean = max(mean, 1e-9)
    if phi <= 0.0 or not math.isfinite(phi):
        return -mean + x * math.log(mean) - math.lgamma(x + 1.0)
    return (
        math.lgamma(x + phi)
        - math.lgamma(phi)
        - math.lgamma(x + 1.0)
        + phi * math.log(phi / (phi + mean))
        + x * math.log(mean / (phi + mean))
    )


def build_corrections(
    observed: np.ndarray,
    max_multiplicity: int,
    lam: float,
    mu: float,
    phi: float,
    outlier: float = 0.0,
    max_single_multiplicity: int = 1,
) -> np.ndarray:
    """Return correction[multiplicity, marker] used by the production default."""
    if not 0.0 <= outlier <= 0.5:
        raise ValueError("outlier probability must be in [0,0.5]")
    broad_ll = -math.log(max(2.0, lam * max_single_multiplicity + 1.0))

    def mixture(x: float, mean: float) -> float:
        nb = log_nb(x, mean, phi)
        if outlier <= 0.0:
            return nb
        a = math.log1p(-outlier) + nb
        b = math.log(outlier) + broad_ll
        hi = max(a, b)
        return hi + math.log(math.exp(a - hi) + math.exp(b - hi))

    corrections = np.zeros((max_multiplicity + 1, observed.size), dtype=np.float64)
    baseline = np.asarray([mixture(float(x), mu) for x in observed], dtype=np.float64)
    for mult in range(1, max_multiplicity + 1):
        mean = lam * mult + mu
        corrections[mult] = np.asarray(
            [mixture(float(x), mean) for x in observed], dtype=np.float64
        ) - baseline
    return corrections


def all_pair_scores(multiplicity: np.ndarray, corrections: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_alleles, n_markers = multiplicity.shape
    n_pairs = n_alleles * (n_alleles + 1) // 2
    allele1 = np.empty(n_pairs, dtype=np.int32)
    allele2 = np.empty(n_pairs, dtype=np.int32)
    scores = np.empty(n_pairs, dtype=np.float64)
    marker_index = np.arange(n_markers, dtype=np.int64)[None, :]
    offset = 0
    for a in range(n_alleles):
        b = np.arange(a, n_alleles, dtype=np.int32)
        totals = multiplicity[a][None, :] + multiplicity[a:]
        values = corrections[totals, marker_index].sum(axis=1)
        end = offset + b.size
        allele1[offset:end] = a
        allele2[offset:end] = b
        scores[offset:end] = values
        offset = end
    return allele1, allele2, scores


def all_pair_feature_metrics(
    multiplicity: np.ndarray, truth: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Raw count-space distances, independent of the NB likelihood."""
    n_alleles = multiplicity.shape[0]
    n_pairs = n_alleles * (n_alleles + 1) // 2
    l1 = np.empty(n_pairs, dtype=np.float64)
    l2 = np.empty(n_pairs, dtype=np.float64)
    cosine = np.empty(n_pairs, dtype=np.float64)
    truth_norm = max(float(np.linalg.norm(truth)), 1e-300)
    offset = 0
    for a in range(n_alleles):
        totals = multiplicity[a][None, :].astype(np.float64) + multiplicity[a:]
        diff = totals - truth[None, :]
        size = totals.shape[0]
        end = offset + size
        l1[offset:end] = np.abs(diff).sum(axis=1)
        l2[offset:end] = np.square(diff).sum(axis=1)
        denom = np.linalg.norm(totals, axis=1) * truth_norm
        cosine[offset:end] = np.divide(
            totals @ truth,
            denom,
            out=np.zeros(size, dtype=np.float64),
            where=denom > 0.0,
        )
        offset = end
    return l1, l2, cosine


def ascending_rank(values: np.ndarray, index: int) -> tuple[int, float]:
    best = float(values.min())
    value = float(values[index])
    return 1 + int(np.count_nonzero(values < value - 1e-9)), value - best


def descending_rank(values: np.ndarray, index: int) -> tuple[int, float]:
    best = float(values.max())
    value = float(values[index])
    return 1 + int(np.count_nonzero(values > value + 1e-12)), value - best


def residual_summary(pair_vector: np.ndarray, truth: np.ndarray) -> str:
    delta = pair_vector.astype(np.float64) - truth
    return (
        f"exact={int(np.count_nonzero(delta == 0))}"
        f" under={int(np.count_nonzero(delta < 0))}"
        f" over={int(np.count_nonzero(delta > 0))}"
        f" missing={int(np.count_nonzero((pair_vector == 0) & (truth > 0)))}"
        f" extra={int(np.count_nonzero((pair_vector > 0) & (truth == 0)))}"
        f" l1={float(np.abs(delta).sum()):.9g}"
        f" l2={float(np.square(delta).sum()):.9g}"
    )


def likelihood_contrast(
    target_vector: np.ndarray,
    called_vector: np.ndarray,
    truth: np.ndarray,
    clumps: np.ndarray,
    corrections: np.ndarray,
) -> list[tuple[str, int, int, float, float, float]]:
    """Target-minus-called ideal log likelihood by support/multiplicity class."""
    classes = {
        "target_only": (target_vector > 0) & (called_vector == 0),
        "called_only": (called_vector > 0) & (target_vector == 0),
        "both_different": (target_vector > 0) & (called_vector > 0) & (target_vector != called_vector),
        "both_same": target_vector == called_vector,
    }
    columns = np.arange(target_vector.size, dtype=np.int64)
    target_ll = corrections[target_vector.astype(np.int64), columns]
    called_ll = corrections[called_vector.astype(np.int64), columns]
    rows: list[tuple[str, int, int, float, float, float]] = []
    for name, mask in classes.items():
        per_clump: dict[int, float] = {}
        for clump, value in zip(clumps[mask], (target_ll[mask] - called_ll[mask]), strict=True):
            per_clump[int(clump)] = per_clump.get(int(clump), 0.0) + float(value)
        magnitudes = sorted((abs(value) for value in per_clump.values()), reverse=True)
        total_magnitude = sum(magnitudes)
        top3_fraction = sum(magnitudes[:3]) / total_magnitude if total_magnitude > 0.0 else 0.0
        rows.append(
            (
                name,
                int(mask.sum()),
                len(per_clump),
                float(truth[mask].sum()),
                float((target_ll[mask] - called_ll[mask]).sum()),
                top3_fraction,
            )
        )
    return rows


def candidate_set(
    multiplicity: np.ndarray, detected: np.ndarray, lengths: np.ndarray, budget: int
) -> np.ndarray:
    n_alleles = multiplicity.shape[0]
    keep_n = min(n_alleles, max(2, budget))
    marker_count = np.count_nonzero(multiplicity, axis=1)
    seen = np.count_nonzero((multiplicity > 0) & detected[None, :], axis=1)
    fraction = np.divide(seen, marker_count, out=np.zeros(n_alleles), where=marker_count > 0)
    score_order = np.lexsort((np.arange(n_alleles), -fraction))

    kept: list[int] = []
    taken = np.zeros(n_alleles, dtype=np.bool_)
    by_score = keep_n - keep_n // 2
    for allele in score_order[:by_score]:
        ai = int(allele)
        kept.append(ai)
        taken[ai] = True

    if len(kept) < keep_n:
        length_order = np.lexsort((np.arange(n_alleles), lengths))
        want = keep_n - len(kept)
        for k in range(want):
            start = (k * (n_alleles - 1)) // (want - 1) if want > 1 else n_alleles // 2
            for step in range(n_alleles):
                index = (start + step) % n_alleles
                ai = int(length_order[index])
                if not taken[ai]:
                    kept.append(ai)
                    taken[ai] = True
                    break
    return np.asarray(sorted(kept), dtype=np.int32)


def splitmix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def exact_equivalence_class(multiplicity: np.ndarray, target: tuple[int, int]) -> list[tuple[int, int]]:
    """Hash all pairs additively, then verify every hash match byte-for-byte."""
    n_alleles, n_markers = multiplicity.shape
    weights1 = [splitmix64(j + 0x13579BDF) for j in range(n_markers)]
    weights2 = [splitmix64(j + 0x2468ACE0) for j in range(n_markers)]
    hashes1: list[int] = []
    hashes2: list[int] = []
    for row in multiplicity:
        h1 = 0
        h2 = 0
        for col in np.flatnonzero(row):
            mult = int(row[col])
            h1 = (h1 + weights1[int(col)] * mult) & MASK64
            h2 = (h2 + weights2[int(col)] * mult) & MASK64
        hashes1.append(h1)
        hashes2.append(h2)

    ta, tb = target
    target_h1 = (hashes1[ta] + hashes1[tb]) & MASK64
    target_h2 = (hashes2[ta] + hashes2[tb]) & MASK64
    target_vector = multiplicity[ta].astype(np.int32) + multiplicity[tb]
    matches: list[tuple[int, int]] = []
    for a in range(n_alleles):
        needed1 = (target_h1 - hashes1[a]) & MASK64
        needed2 = (target_h2 - hashes2[a]) & MASK64
        for b in range(a, n_alleles):
            if hashes1[b] != needed1 or hashes2[b] != needed2:
                continue
            if np.array_equal(multiplicity[a].astype(np.int32) + multiplicity[b], target_vector):
                matches.append((a, b))
    return matches


def pair_mask(allele1: np.ndarray, allele2: np.ndarray, kept: np.ndarray) -> np.ndarray:
    allowed = np.zeros(max(int(allele2.max()) + 1, 1), dtype=np.bool_)
    allowed[kept] = True
    return allowed[allele1] & allowed[allele2]


def pair_index(allele1: np.ndarray, allele2: np.ndarray, pair: tuple[int, int]) -> int:
    hits = np.flatnonzero((allele1 == pair[0]) & (allele2 == pair[1]))
    if hits.size != 1:
        raise ValueError(f"pair {pair} is outside the allele space")
    return int(hits[0])


def optional_pair_index(
    allele1: np.ndarray, allele2: np.ndarray, pair: tuple[int, int]
) -> int | None:
    hits = np.flatnonzero((allele1 == pair[0]) & (allele2 == pair[1]))
    if hits.size > 1:
        raise AssertionError(f"pair {pair} was enumerated more than once")
    return int(hits[0]) if hits.size == 1 else None


def rank_of(scores: np.ndarray, index: int, mask: np.ndarray | None = None) -> tuple[int, float]:
    active = np.ones(scores.size, dtype=np.bool_) if mask is None else mask
    if not active[index]:
        return -1, math.nan
    value = scores[index]
    best = float(scores[active].max())
    rank = 1 + int(np.count_nonzero(active & (scores > value + 1e-9)))
    return rank, value - best


def top_pairs(
    scores: np.ndarray, allele1: np.ndarray, allele2: np.ndarray, mask: np.ndarray | None, n: int
) -> list[tuple[int, int, int, float, float]]:
    active = np.ones(scores.size, dtype=np.bool_) if mask is None else mask
    indices = np.flatnonzero(active)
    order = indices[np.argsort(-scores[indices], kind="stable")[:n]]
    best = float(scores[order[0]])
    return [
        (rank, int(allele1[idx]), int(allele2[idx]), float(scores[idx]), float(scores[idx] - best))
        for rank, idx in enumerate(order, start=1)
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conf", required=True, type=Path, help="BLOCK.conf.tsv from --dump-block")
    parser.add_argument("--markers", required=True, type=Path, help="all-marker dump carrying truth_mult")
    parser.add_argument("--alleles", required=True, type=Path, help="BLOCK.tsv allele summary")
    parser.add_argument("--depth", required=True, type=Path, help="reads.depth.tsv from the same run")
    parser.add_argument("--block", type=int, required=True)
    parser.add_argument("--target", type=pair_arg, required=True, help="certified pair A,B")
    parser.add_argument("--called", type=pair_arg, help="final called pair A,B")
    parser.add_argument("--max-alleles", type=int, default=64)
    parser.add_argument(
        "--candidate-only",
        action="store_true",
        help="rank only the reproduced candidate set (useful for very large all-kmer dumps)",
    )
    parser.add_argument(
        "--outliers",
        default="0",
        help="comma-separated --marker-outlier values to sweep (default: 0)",
    )
    parser.add_argument(
        "--clump-powers",
        default="0",
        help="comma-separated powers p for marker weight 1/clump_size^p (default: 0)",
    )
    parser.add_argument("--lambda-hap", type=float, help="override block lambda")
    parser.add_argument("--mu", type=float, help="override error background")
    parser.add_argument("--phi", type=float, help="override NB shape")
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--output", type=Path, help="write ranked summaries as TSV")
    args = parser.parse_args()

    multiplicity, observed, slots = read_conf(args.conf)
    truth, clumps, anchors = read_marker_dump(args.markers, args.block, slots)
    lengths = read_lengths(args.alleles, multiplicity.shape[0])
    depth_lambda, mean_lambda = read_depth(args.depth, args.block)
    lam = args.lambda_hap if args.lambda_hap is not None else depth_lambda
    mu = args.mu if args.mu is not None else max(0.01, 0.02 * mean_lambda)
    phi = args.phi if args.phi is not None else estimate_phi(anchors)
    outliers = [float(value) for value in args.outliers.split(",")]
    if not outliers or any(value < 0.0 or value > 0.5 for value in outliers):
        raise ValueError("--outliers must contain values in [0,0.5]")
    clump_powers = [float(value) for value in args.clump_powers.split(",")]
    if not clump_powers or any(value < 0.0 for value in clump_powers):
        raise ValueError("--clump-powers must contain non-negative values")
    unique_clumps, clump_inverse, clump_sizes = np.unique(
        clumps, return_inverse=True, return_counts=True
    )

    for name, pair in (("target", args.target), ("called", args.called)):
        if pair is not None and (pair[0] < 0 or pair[1] >= multiplicity.shape[0]):
            raise ValueError(f"{name} pair {pair} is outside 0..{multiplicity.shape[0] - 1}")

    target_vector = multiplicity[args.target[0]].astype(np.float64) + multiplicity[args.target[1]]
    scenarios = [
        Scenario("panel_self", lam * target_vector + mu, target_vector > 0),
        Scenario("truth_projection", lam * truth + mu, truth > 0),
        Scenario("observed_reads", observed, observed > 0),
    ]

    equivalence = exact_equivalence_class(multiplicity, args.target)
    print(f"alleles={multiplicity.shape[0]} markers={multiplicity.shape[1]} pairs={multiplicity.shape[0] * (multiplicity.shape[0] + 1) // 2}")
    print(f"lambda={lam:.9g} mu={mu:.9g} phi={phi:.9g} anchors={len(anchors)}")
    print(f"target={args.target[0]},{args.target[1]} exact_marker_equivalents={len(equivalence)} pairs={equivalence}")

    output_rows: list[dict[str, object]] = []
    max_pair_mult = int(multiplicity.max()) * 2
    allele1: np.ndarray | None = None
    allele2: np.ndarray | None = None
    for outlier in outliers:
      for clump_power in clump_powers:
       marker_weights = np.power(clump_sizes[clump_inverse].astype(np.float64), -clump_power)
       marker_weights /= marker_weights.mean()
       for scenario in scenarios:
        kept = candidate_set(multiplicity, scenario.detected, lengths, args.max_alleles)
        corrections = build_corrections(
            scenario.observed,
            max_pair_mult,
            lam,
            mu,
            phi,
            outlier,
            int(multiplicity.max()),
        )
        corrections *= marker_weights[None, :]
        score_matrix = multiplicity
        if args.candidate_only:
            score_matrix = multiplicity[kept]
        pa, pb, scores = all_pair_scores(score_matrix, corrections)
        if args.candidate_only:
            pa = kept[pa]
            pb = kept[pb]
        if allele1 is None:
            allele1, allele2 = pa, pb
        elif not args.candidate_only and not (
            np.array_equal(allele1, pa) and np.array_equal(allele2, pb)
        ):
            raise AssertionError("pair enumeration changed across scenarios")
        kept_mask = pair_mask(pa, pb, kept)
        target_idx = optional_pair_index(pa, pb, args.target)
        target_survives = target_idx is not None and bool(kept_mask[target_idx])
        if target_idx is None:
            full_rank, full_delta = -1, math.nan
            kept_rank, kept_delta = -1, math.nan
        else:
            full_rank, full_delta = rank_of(scores, target_idx)
            kept_rank, kept_delta = rank_of(scores, target_idx, kept_mask)
        called_summary = ""
        called_idx: int | None = None
        if args.called is not None:
            called_idx = optional_pair_index(pa, pb, args.called)
            called_survives = called_idx is not None and bool(kept_mask[called_idx])
            if called_idx is None:
                called_rank, called_delta = -1, math.nan
            else:
                called_rank, called_delta = rank_of(scores, called_idx)
            called_summary = (
                f" called_rank={called_rank} called_delta={called_delta:.6g}"
                f" called_survives={int(called_survives)}"
            )
        else:
            called_rank, called_delta, called_survives = -1, math.nan, False

        best = top_pairs(scores, pa, pb, None, args.top)
        best_kept = top_pairs(scores, pa, pb, kept_mask, min(args.top, int(kept_mask.sum())))
        print(
            f"outlier={outlier:.9g} clump_power={clump_power:.9g} {scenario.name}:"
            f" best={best[0][1]},{best[0][2]}"
            f" target_rank={full_rank} target_delta={full_delta:.6g}"
            f" target_survives={int(target_survives)}"
            f" kept_best={best_kept[0][1]},{best_kept[0][2]} kept_target_rank={kept_rank}"
            f" kept_target_delta={kept_delta:.6g}{called_summary}"
        )
        if (
            scenario.name == "truth_projection"
            and outlier == outliers[0]
            and clump_power == clump_powers[0]
            and target_idx is not None
        ):
            l1, l2, cosine = all_pair_feature_metrics(score_matrix, truth)
            t_l1 = ascending_rank(l1, target_idx)
            t_l2 = ascending_rank(l2, target_idx)
            t_cos = descending_rank(cosine, target_idx)
            print(
                "truth_feature_rank:"
                f" target_l1={t_l1[0]} delta={t_l1[1]:.9g}"
                f" target_l2={t_l2[0]} delta={t_l2[1]:.9g}"
                f" target_cosine={t_cos[0]} delta={t_cos[1]:.9g}"
            )
            target_pair_vector = (
                multiplicity[args.target[0]].astype(np.int32) + multiplicity[args.target[1]]
            )
            print(f"truth_residual target: {residual_summary(target_pair_vector, truth)}")
            if args.called is not None and called_idx is not None:
                called_pair_vector = (
                    multiplicity[args.called[0]].astype(np.int32) + multiplicity[args.called[1]]
                )
                c_l1 = ascending_rank(l1, called_idx)
                c_l2 = ascending_rank(l2, called_idx)
                c_cos = descending_rank(cosine, called_idx)
                print(
                    "truth_feature_rank:"
                    f" called_l1={c_l1[0]} delta={c_l1[1]:.9g}"
                    f" called_l2={c_l2[0]} delta={c_l2[1]:.9g}"
                    f" called_cosine={c_cos[0]} delta={c_cos[1]:.9g}"
                )
                print(f"truth_residual called: {residual_summary(called_pair_vector, truth)}")
                for name, count, n_clumps, truth_mass, ll_delta, top3_fraction in likelihood_contrast(
                    target_pair_vector, called_pair_vector, truth, clumps, corrections
                ):
                    print(
                        f"truth_contrast {name}: markers={count} clumps={n_clumps}"
                        f" truth_mult={truth_mass:.9g} target_minus_called_ll={ll_delta:.9g}"
                        f" top3_clump_abs_fraction={top3_fraction:.6g}"
                    )
        for scope, rows in (("full", best), ("kept", best_kept)):
            for rank, a, b, score, delta in rows:
                output_rows.append(
                    {
                        "scenario": scenario.name,
                        "outlier": outlier,
                        "clump_power": clump_power,
                        "scope": scope,
                        "rank": rank,
                        "allele1": a,
                        "allele2": b,
                        "score": score,
                        "delta_from_best": delta,
                        "target": int((a, b) == args.target),
                        "called": int(args.called is not None and (a, b) == args.called),
                        "total_bp": int(lengths[a] + lengths[b]),
                    }
                )
        output_rows.extend(
            [
                {
                    "scenario": scenario.name,
                    "outlier": outlier,
                    "clump_power": clump_power,
                    "scope": "audit_target",
                    "rank": full_rank,
                    "allele1": args.target[0],
                    "allele2": args.target[1],
                    "score": scores[target_idx] if target_idx is not None else math.nan,
                    "delta_from_best": full_delta,
                    "target": 1,
                    "called": int(args.called == args.target),
                    "total_bp": int(lengths[args.target[0]] + lengths[args.target[1]]),
                },
                {
                    "scenario": scenario.name,
                    "outlier": outlier,
                    "clump_power": clump_power,
                    "scope": "audit_target_kept",
                    "rank": kept_rank,
                    "allele1": args.target[0],
                    "allele2": args.target[1],
                    "score": scores[target_idx] if target_idx is not None else math.nan,
                    "delta_from_best": kept_delta,
                    "target": 1,
                    "called": int(args.called == args.target),
                    "total_bp": int(lengths[args.target[0]] + lengths[args.target[1]]),
                },
            ]
        )
        if args.called is not None and called_idx is not None:
            output_rows.append(
                {
                    "scenario": scenario.name,
                    "outlier": outlier,
                    "clump_power": clump_power,
                    "scope": "audit_called",
                    "rank": called_rank,
                    "allele1": args.called[0],
                    "allele2": args.called[1],
                    "score": scores[called_idx],
                    "delta_from_best": called_delta,
                    "target": int(args.called == args.target),
                    "called": 1,
                    "total_bp": int(lengths[args.called[0]] + lengths[args.called[1]]),
                }
            )

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fields = [
            "scenario", "outlier", "clump_power", "scope", "rank", "allele1", "allele2", "score",
            "delta_from_best", "target", "called", "total_bp",
        ]
        with args.output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(output_rows)
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
