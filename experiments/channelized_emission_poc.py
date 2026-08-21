#!/usr/bin/env python3
"""Proof of concept: separate support mismatch from multiplicity fit.

The production negative-binomial emission lets a marker absent from a candidate
contribute an arbitrarily large penalty.  In an off-panel tandem array, a small
number of novel unit-variant markers can therefore veto the pair whose shared
markers have the best multiplicities.

This diagnostic replaces that with two explicit channels:

  * support: a bounded penalty for observed/predicted presence disagreement,
    divided across markers in the same fragment-scale clump;
  * multiplicity: negative-binomial deviance from the best integer copy state,
    evaluated only where both observation and candidate support the marker.

It is not a production model.  The sweep establishes whether this separation
can rescue a sequence-certified pair before investing in C++ implementation.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from ideal_multiplicity_oracle import (
    Scenario,
    all_pair_scores,
    candidate_set,
    estimate_phi,
    pair_arg,
    pair_index,
    pair_mask,
    rank_of,
    read_conf,
    read_depth,
    read_lengths,
    read_marker_dump,
    top_pairs,
)
from ideal_multiplicity_oracle import log_nb


def channel_corrections(
    observed: np.ndarray,
    active: np.ndarray,
    max_multiplicity: int,
    lam: float,
    mu: float,
    phi: float,
    support_penalty: float,
    support_weights: np.ndarray,
) -> np.ndarray:
    """Per-marker channel score indexed by candidate pair multiplicity."""
    log_likelihood = np.empty((max_multiplicity + 1, observed.size), dtype=np.float64)
    for mult in range(1, max_multiplicity + 1):
        mean = lam * mult + mu
        log_likelihood[mult] = np.asarray(
            [log_nb(float(value), mean, phi) for value in observed], dtype=np.float64
        )

    # Zero is not part of the integer copy fit.  It is handled by the support channel.
    best_supported = log_likelihood[1:].max(axis=0)
    correction = np.zeros_like(log_likelihood)
    correction[0, active] = -support_penalty * support_weights[active]
    for mult in range(1, max_multiplicity + 1):
        correction[mult, active] = log_likelihood[mult, active] - best_supported[active]
        correction[mult, ~active] = -support_penalty * support_weights[~active]
    return correction


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conf", required=True, type=Path)
    parser.add_argument("--markers", required=True, type=Path)
    parser.add_argument("--alleles", required=True, type=Path)
    parser.add_argument("--depth", required=True, type=Path)
    parser.add_argument("--block", required=True, type=int)
    parser.add_argument("--target", required=True, type=pair_arg)
    parser.add_argument("--called", type=pair_arg)
    parser.add_argument("--max-alleles", type=int, default=64)
    parser.add_argument("--penalties", default="0,0.1,0.3,1,3,10,30")
    parser.add_argument(
        "--observed-thresholds",
        default="0.1,0.25,0.5",
        help="observed support thresholds as fractions of haploid lambda",
    )
    parser.add_argument(
        "--support-clump-power",
        type=float,
        default=1.0,
        help="presence penalty weight is 1/clump_size^p (default: 1)",
    )
    parser.add_argument("--lambda-hap", type=float)
    parser.add_argument("--mu", type=float)
    parser.add_argument("--phi", type=float)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    multiplicity, observed, slots = read_conf(args.conf)
    truth, clumps, anchors = read_marker_dump(args.markers, args.block, slots)
    lengths = read_lengths(args.alleles, multiplicity.shape[0])
    depth_lambda, mean_lambda = read_depth(args.depth, args.block)
    lam = args.lambda_hap if args.lambda_hap is not None else depth_lambda
    mu = args.mu if args.mu is not None else max(0.01, 0.02 * mean_lambda)
    phi = args.phi if args.phi is not None else estimate_phi(anchors)
    penalties = [float(value) for value in args.penalties.split(",")]
    thresholds = [float(value) for value in args.observed_thresholds.split(",")]
    if any(value < 0.0 for value in penalties + thresholds):
        raise ValueError("penalties and thresholds must be non-negative")

    _, inverse, sizes = np.unique(clumps, return_inverse=True, return_counts=True)
    support_weights = np.power(sizes[inverse].astype(np.float64), -args.support_clump_power)
    # Only relative support weight versus the multiplicity deviance matters.  A mean-one
    # normalization makes penalty values comparable between clump powers.
    support_weights /= support_weights.mean()

    target_vector = multiplicity[args.target[0]].astype(np.float64) + multiplicity[args.target[1]]
    base_scenarios = [
        Scenario("panel_self", lam * target_vector + mu, target_vector > 0),
        Scenario("truth_projection", lam * truth + mu, truth > 0),
    ]
    scenarios = list(base_scenarios)
    for threshold in thresholds:
        scenarios.append(
            Scenario(
                f"observed_reads_t{threshold:g}",
                observed,
                observed >= threshold * lam,
            )
        )

    max_pair_mult = 2 * int(multiplicity.max())
    rows: list[dict[str, object]] = []
    for penalty in penalties:
        for scenario in scenarios:
            corrections = channel_corrections(
                scenario.observed,
                scenario.detected,
                max_pair_mult,
                lam,
                mu,
                phi,
                penalty,
                support_weights,
            )
            allele1, allele2, scores = all_pair_scores(multiplicity, corrections)
            kept = candidate_set(multiplicity, scenario.detected, lengths, args.max_alleles)
            kept_pairs = pair_mask(allele1, allele2, kept)
            target_index = pair_index(allele1, allele2, args.target)
            target_rank, target_delta = rank_of(scores, target_index)
            kept_rank, kept_delta = rank_of(scores, target_index, kept_pairs)
            best = top_pairs(scores, allele1, allele2, None, 1)[0]
            kept_best = top_pairs(scores, allele1, allele2, kept_pairs, 1)[0]

            called_rank = -1
            called_delta = float("nan")
            if args.called is not None:
                called_index = pair_index(allele1, allele2, args.called)
                called_rank, called_delta = rank_of(scores, called_index)
            print(
                f"penalty={penalty:.9g} {scenario.name}: best={best[1]},{best[2]}"
                f" target_rank={target_rank} target_delta={target_delta:.6g}"
                f" kept_best={kept_best[1]},{kept_best[2]} kept_target_rank={kept_rank}"
                f" kept_target_delta={kept_delta:.6g} called_rank={called_rank}"
                f" called_delta={called_delta:.6g}"
            )
            rows.append(
                {
                    "penalty": penalty,
                    "scenario": scenario.name,
                    "best_a": best[1],
                    "best_b": best[2],
                    "target_rank": target_rank,
                    "target_delta": target_delta,
                    "target_survives": int(kept_pairs[target_index]),
                    "kept_best_a": kept_best[1],
                    "kept_best_b": kept_best[2],
                    "kept_target_rank": kept_rank,
                    "kept_target_delta": kept_delta,
                    "called_rank": called_rank,
                    "called_delta": called_delta,
                }
            )

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fields = list(rows[0])
        with args.output.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
