#!/usr/bin/env python3
"""Standalone proof-of-concept for dosage-first graph genotyping.

This file deliberately does not import or modify panvar's genotyper.  It consumes
the diagnostic files produced by ``panvar genotype --dump-block`` and asks a
falsifiable question: does estimating marker dosage before choosing an allele
pair rescue a real held-out tandem-array failure?

The model has two stages.

1. Turn read counts into a latent diploid marker-dosage vector

       z_j = max(0, (count_j - background) / haploid_depth)

2. At an array whose graph paths demonstrate a reliable multiplicity/length
   calibration, choose the supported diploid dosage stratum nearest z, then use
   a bounded variance-stabilised residual to choose sequence composition inside
   that stratum.

The dosage stratum is intentionally conditional.  It is not applied unless
marker multiplicity predicts path length with high R^2 and the candidate totals
form separated length bands.  Simple bubbles instead use internal and junction
components directly.  Across bubbles, the same local emissions can be joined by
a switching HMM, allowing mosaic haplotypes instead of forcing one panel pair
through the whole locus.

Only NumPy is required.  If GFA truth paths are supplied, a small shared edlib
library is built in the system temporary directory from panvar's vendored edlib
source; production build files and sources remain untouched.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import os
from pathlib import Path
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass
class Dump:
    dosage_by_allele: np.ndarray
    observed: np.ndarray
    allele_bp: np.ndarray
    allele_seq: list[str]


@dataclass(order=True)
class PairScore:
    loss: float
    allele1: int
    allele2: int
    total_bp: int


def load_dump(conf_path: Path, fasta_path: Path) -> Dump:
    rows: list[tuple[int, int, int, int]] = []
    slots: set[int] = set()
    max_allele = -1
    observation: dict[int, int] = {}
    with conf_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            allele = int(row["allele"])
            slot = int(row["slot"])
            mult = int(row["mult"])
            obs = int(row["obs"])
            old = observation.setdefault(slot, obs)
            if old != obs:
                raise ValueError(f"slot {slot} has inconsistent observations {old} and {obs}")
            rows.append((allele, slot, mult, obs))
            slots.add(slot)
            max_allele = max(max_allele, allele)
    if max_allele < 0 or not slots:
        raise ValueError("empty genotype diagnostic dump")

    ordered_slots = sorted(slots)
    slot_index = {slot: i for i, slot in enumerate(ordered_slots)}
    observed = np.asarray([observation[s] for s in ordered_slots], dtype=np.float64)
    dosage = np.zeros((max_allele + 1, len(ordered_slots)), dtype=np.float32)
    for allele, slot, mult, _ in rows:
        dosage[allele, slot_index[slot]] = mult

    seqs: list[str] = []
    current: int | None = None
    with fasta_path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = int(line[1:].split()[0])
                while len(seqs) <= current:
                    seqs.append("")
            elif current is None:
                raise ValueError(f"sequence before FASTA header in {fasta_path}")
            else:
                seqs[current] += line.upper()
    if len(seqs) != max_allele + 1 or any(not s for s in seqs):
        raise ValueError("the diagnostic FASTA does not contain every dumped allele")
    bp = np.asarray([len(s) for s in seqs], dtype=np.int64)
    return Dump(dosage, observed, bp, seqs)


def huber(x: np.ndarray, c: float) -> np.ndarray:
    absolute = np.abs(x)
    return np.where(absolute <= c, 0.5 * x * x, c * (absolute - 0.5 * c))


def robust_count_loss(
    pair_dosage: np.ndarray,
    observed: np.ndarray,
    haploid_depth: float,
    background: float,
    huber_c: float,
) -> np.ndarray:
    """Anscombe-like square-root residual followed by bounded-influence loss."""
    predicted = haploid_depth * pair_dosage + background
    residual = 2.0 * np.sqrt(predicted + 0.375) - 2.0 * np.sqrt(observed + 0.375)
    return np.sum(huber(residual, huber_c), axis=1)


def all_pair_scores(
    dump: Dump,
    haploid_depth: float,
    background: float,
    huber_c: float,
) -> list[PairScore]:
    out: list[PairScore] = []
    matrix = dump.dosage_by_allele
    for allele1 in range(matrix.shape[0]):
        pair = matrix[allele1:] + matrix[allele1]
        losses = robust_count_loss(pair, dump.observed, haploid_depth, background, huber_c)
        for offset, loss in enumerate(losses):
            allele2 = allele1 + offset
            out.append(PairScore(float(loss), allele1, allele2,
                                 int(dump.allele_bp[allele1] + dump.allele_bp[allele2])))
    out.sort()
    return out


def fit_graph_length_calibration(dump: Dump, latent: np.ndarray) -> tuple[float, float, float, float]:
    marker_mass = np.sum(dump.dosage_by_allele, axis=1, dtype=np.float64)
    slope, intercept = np.polyfit(marker_mass, dump.allele_bp.astype(float), 1)
    fitted = slope * marker_mass + intercept
    residual = dump.allele_bp - fitted
    total = dump.allele_bp - np.mean(dump.allele_bp)
    r2 = 1.0 - float(np.dot(residual, residual) / max(1e-12, np.dot(total, total)))
    diploid_bp = float(slope * np.sum(latent) + 2.0 * intercept)
    return float(slope), float(intercept), r2, diploid_bp


def supported_length_bands(allele_bp: np.ndarray, minimum_gap: int) -> list[tuple[int, int]]:
    totals = allele_bp[:, None] + allele_bp[None, :]
    values = np.unique(totals[np.triu_indices(len(allele_bp))])
    if len(values) == 0:
        return []
    bands: list[tuple[int, int]] = []
    start = int(values[0])
    previous = int(values[0])
    for value0 in values[1:]:
        value = int(value0)
        if value - previous >= minimum_gap:
            bands.append((start, previous))
            start = value
        previous = value
    bands.append((start, previous))
    return bands


def interval_distance(value: float, interval: tuple[int, int]) -> float:
    lo, hi = interval
    if value < lo:
        return lo - value
    if value > hi:
        return value - hi
    return 0.0


def select_dosage_band(
    scores: Sequence[PairScore],
    bands: Sequence[tuple[int, int]],
    target_bp: float,
) -> tuple[tuple[int, int], list[PairScore]]:
    if not bands:
        raise ValueError("no supported length bands")
    band = min(bands, key=lambda x: (interval_distance(target_bp, x), x[0]))
    selected = [score for score in scores if band[0] <= score.total_bp <= band[1]]
    if not selected:
        raise AssertionError("selected a length band containing no candidate pairs")
    return band, selected


def rank_of(scores: Sequence[PairScore], pair: tuple[int, int]) -> int | None:
    wanted = tuple(sorted(pair))
    for rank, score in enumerate(scores, 1):
        if (score.allele1, score.allele2) == wanted:
            return rank
    return None


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def truth_walks_from_gfa(
    gfa_path: Path,
    bubbles_path: Path,
    bubble_id: int,
    path_names: Sequence[str],
) -> list[str]:
    bubble: dict[str, str] | None = None
    with bubbles_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if int(row["bubble_id"]) == bubble_id:
                bubble = row
                break
    if bubble is None:
        raise ValueError(f"bubble {bubble_id} is absent from {bubbles_path}")
    source = bubble["source"]
    sink = bubble["sink"]
    inside = set(filter(None, bubble.get("inside_nodes", "").split(";")))

    nodes: dict[str, str] = {}
    paths: dict[str, list[tuple[str, bool]]] = {}
    wanted = set(path_names)
    with gfa_path.open() as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if fields[0] == "S":
                nodes[fields[1]] = fields[2]
            elif fields[0] == "P" and fields[1] in wanted:
                paths[fields[1]] = [(token[:-1], token[-1] == "-")
                                    for token in fields[2].split(",")]
    missing = wanted.difference(paths)
    if missing:
        raise ValueError(f"truth paths missing from GFA: {sorted(missing)}")

    def best_interval(steps: Sequence[tuple[str, bool]]) -> tuple[int, int, bool]:
        src = [i for i, step in enumerate(steps) if step[0] == source]
        snk = [i for i, step in enumerate(steps) if step[0] == sink]
        candidates: list[tuple[int, int, int, bool]] = []
        for a in src:
            for b in snk:
                if a == b:
                    continue
                lo, hi = sorted((a, b))
                count = sum(steps[i][0] in inside for i in range(lo + 1, hi))
                if count:
                    candidates.append((-count, hi - lo, lo, a < b))
        if not candidates:
            raise ValueError("truth path does not traverse the requested bubble")
        _, span, lo, forward = min(candidates)
        return lo, lo + span, forward

    result: list[str] = []
    for name in path_names:
        steps = paths[name]
        lo, hi, forward = best_interval(steps)
        chosen = list(steps[lo:hi + 1])
        if not forward:
            chosen = [(node, not reverse) for node, reverse in reversed(chosen)]
        if chosen[0][0] != source or chosen[-1][0] != sink:
            raise ValueError(f"failed to orient bubble walk for {name}")
        result.append("".join(reverse_complement(nodes[node]) if reverse else nodes[node]
                              for node, reverse in chosen))
    return result


class EdlibAlignConfig(ctypes.Structure):
    _fields_ = [("k", ctypes.c_int), ("mode", ctypes.c_int), ("task", ctypes.c_int),
                ("additionalEqualities", ctypes.c_void_p),
                ("additionalEqualitiesLength", ctypes.c_int)]


class EdlibAlignResult(ctypes.Structure):
    _fields_ = [("status", ctypes.c_int), ("editDistance", ctypes.c_int),
                ("endLocations", ctypes.POINTER(ctypes.c_int)),
                ("startLocations", ctypes.POINTER(ctypes.c_int)),
                ("numLocations", ctypes.c_int),
                ("alignment", ctypes.POINTER(ctypes.c_ubyte)),
                ("alignmentLength", ctypes.c_int), ("alphabetLength", ctypes.c_int)]


def load_edlib(repo_root: Path) -> ctypes.CDLL:
    source = repo_root / "external/edlib/edlib/src/edlib.cpp"
    include = repo_root / "external/edlib/edlib/include"
    if not source.exists():
        raise FileNotFoundError("vendored edlib is absent; initialise the external/edlib submodule")
    suffix = ".dylib" if os.uname().sysname == "Darwin" else ".so"
    library = Path(tempfile.gettempdir()) / f"panvar_multiplicity_poc_edlib{suffix}"
    if not library.exists() or library.stat().st_mtime < source.stat().st_mtime:
        subprocess.run(["c++", "-O3", "-std=c++17", "-shared", "-fPIC", str(source),
                        "-I", str(include), "-o", str(library)], check=True)
    edlib = ctypes.CDLL(str(library))
    edlib.edlibNewAlignConfig.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int,
                                          ctypes.c_void_p, ctypes.c_int]
    edlib.edlibNewAlignConfig.restype = EdlibAlignConfig
    edlib.edlibAlign.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
                                 ctypes.c_int, EdlibAlignConfig]
    edlib.edlibAlign.restype = EdlibAlignResult
    edlib.edlibFreeAlignResult.argtypes = [EdlibAlignResult]
    return edlib


def nw_stats(edlib: ctypes.CDLL, first: str, second: str) -> tuple[int, int]:
    # EDLIB_MODE_NW=0, EDLIB_TASK_PATH=2.  PATH supplies alignmentLength, the
    # denominator used by panvar's graded genotype accuracy calculation.
    config = edlib.edlibNewAlignConfig(-1, 0, 2, None, 0)
    a, b = first.encode(), second.encode()
    result = edlib.edlibAlign(a, len(a), b, len(b), config)
    try:
        if result.status != 0 or result.editDistance < 0:
            raise RuntimeError("edlib global alignment failed")
        return int(result.editDistance), int(result.alignmentLength)
    finally:
        edlib.edlibFreeAlignResult(result)


def pair_identity(edlib: ctypes.CDLL, candidate: Sequence[str], truth: Sequence[str]) -> tuple[float, float]:
    stats = [[nw_stats(edlib, candidate[c], truth[t]) for t in range(2)] for c in range(2)]
    direct = stats[0][0][0] + stats[1][1][0]
    crossed = stats[1][0][0] + stats[0][1][0]
    order = (0, 1) if direct <= crossed else (1, 0)
    chosen = [stats[order[t]][t] for t in range(2)]
    identities = [1.0 - edits / max(1, aligned) for edits, aligned in chosen]
    weighted = 1.0 - sum(x[0] for x in chosen) / max(1, sum(x[1] for x in chosen))
    return float(sum(identities) / 2.0), float(weighted)


def synthetic_cn_case() -> None:
    # Every allele contains the same unit markers; only traversal multiplicity differs.
    alleles = np.asarray([[copy] * 12 for copy in range(1, 7)], dtype=float)
    truth = alleles[1] + alleles[4]  # CN 2 + 5 = total CN 7
    observed = 15.0 * truth + 0.3
    presence = (alleles > 0).astype(float)
    presence_losses = []
    dosage_losses = []
    for i in range(len(alleles)):
        for j in range(i, len(alleles)):
            presence_losses.append((float(np.sum(((presence[i] + presence[j]) > 0) != (observed > 1))), i, j))
            dosage_losses.append((float(np.sum((alleles[i] + alleles[j] - truth) ** 2)), i, j))
    presence_losses.sort()
    dosage_losses.sort()
    presence_totals = {i + j + 2 for loss, i, j in presence_losses if loss == presence_losses[0][0]}
    dosage_totals = {i + j + 2 for loss, i, j in dosage_losses if loss == dosage_losses[0][0]}
    print("synthetic/copy-number")
    print(f"  truth diploid CN: 7")
    print(f"  presence-only optimum permits total CN: {sorted(presence_totals)}")
    print(f"  dosage optimum permits total CN: {sorted(dosage_totals)}")
    assert dosage_totals == {7}


def synthetic_flank_case() -> None:
    # Columns: repeat dosage, left-A, right-A, left-B, right-B.  Internally A and
    # B are identical; their source/repeat and repeat/sink junctions identify them.
    allele_a = np.asarray([3, 1, 1, 0, 0], dtype=float)
    allele_b = np.asarray([3, 0, 0, 1, 1], dtype=float)
    candidates = [allele_a, allele_b]
    observed = 20.0 * (allele_a + allele_b) + 0.2

    def score(pair: tuple[int, int], columns: slice) -> float:
        predicted = 20.0 * (candidates[pair[0]] + candidates[pair[1]]) + 0.2
        return float(np.sum((np.sqrt(predicted[columns]) - np.sqrt(observed[columns])) ** 2))

    pairs = [(0, 0), (0, 1), (1, 1)]
    internal = [score(pair, slice(0, 1)) for pair in pairs]
    with_flanks = [score(pair, slice(None)) for pair in pairs]
    print("synthetic/flank-junctions")
    print(f"  internal-only losses AA/AB/BB: {internal} (unidentifiable)")
    print(f"  with junction losses AA/AB/BB: {with_flanks} (AB recovered)")
    assert len(set(internal)) == 1 and int(np.argmin(with_flanks)) == 1


def viterbi(emission: np.ndarray, switch_penalty: float) -> list[int]:
    blocks, states = emission.shape
    cost = emission[0].copy()
    back = np.zeros((blocks, states), dtype=int)
    for block in range(1, blocks):
        nxt = np.empty(states)
        for state in range(states):
            choices = cost + switch_penalty * (np.arange(states) != state)
            previous = int(np.argmin(choices))
            nxt[state] = choices[previous] + emission[block, state]
            back[block, state] = previous
        cost = nxt
    path = [int(np.argmin(cost))]
    for block in range(blocks - 1, 0, -1):
        path.append(int(back[block, path[-1]]))
    return list(reversed(path))


def synthetic_mosaic_case() -> None:
    # The panel contains only all-A and all-B donor haplotypes.  The sample is
    # homozygous A in the first half and homozygous B in the second: neither
    # whole donor pair is correct, but all local graph components exist.
    states = np.asarray([2.0, 1.0, 0.0])  # AA, AB, BB dosage of A marker
    truth = np.asarray([2.0, 2.0, 0.0, 0.0])
    observed = 18.0 * truth + 0.2
    predicted = 18.0 * states[None, :] + 0.2
    emission = (np.sqrt(observed[:, None]) - np.sqrt(predicted)) ** 2
    fixed = int(np.argmin(np.sum(emission, axis=0)))
    path = viterbi(emission, switch_penalty=0.5)
    print("synthetic/mosaic")
    print(f"  best forced locus-wide state: {['AA', 'AB', 'BB'][fixed]}")
    print(f"  block HMM states: {[['AA', 'AB', 'BB'][state] for state in path]}")
    assert path == [0, 0, 2, 2]


def run_synthetic(_: argparse.Namespace) -> int:
    synthetic_cn_case()
    synthetic_flank_case()
    synthetic_mosaic_case()
    print("synthetic verdict: all dosage/flank/mosaic assertions passed")
    return 0


def print_scores(label: str, scores: Sequence[PairScore], truth_bp: int | None, top: int) -> None:
    print(label)
    print("  rank\talleles\ttotal_bp\tcount_loss" + ("\t|bp-truth|" if truth_bp is not None else ""))
    for rank, score in enumerate(scores[:top], 1):
        suffix = f"\t{abs(score.total_bp - truth_bp)}" if truth_bp is not None else ""
        print(f"  {rank}\t{score.allele1},{score.allele2}\t{score.total_bp}\t{score.loss:.6f}{suffix}")


def run_real(args: argparse.Namespace) -> int:
    dump = load_dump(args.conf, args.alleles)
    if args.haploid_depth <= 0:
        raise ValueError("--haploid-depth must be positive")
    latent = np.maximum(0.0, (dump.observed - args.background) / args.haploid_depth)
    scores = all_pair_scores(dump, args.haploid_depth, args.background, args.huber_c)
    slope, intercept, r2, target_bp = fit_graph_length_calibration(dump, latent)
    bands = supported_length_bands(dump.allele_bp, args.minimum_band_gap)
    eligible = r2 >= args.minimum_r2 and len(bands) >= 2
    if eligible:
        band, constrained = select_dosage_band(scores, bands, target_bp)
    else:
        band, constrained = (0, 0), scores

    print("real-data dosage projection")
    print(f"  alleles={dump.dosage_by_allele.shape[0]} markers={dump.dosage_by_allele.shape[1]}")
    print(f"  latent diploid marker mass={np.sum(latent):.3f}")
    print(f"  graph calibration: bp = {slope:.6f} * marker_mass + {intercept:.3f}; R^2={r2:.8f}")
    print(f"  graph-calibrated diploid span={target_bp:.3f} bp")
    print(f"  separated supported bands={len(bands)} (gap >= {args.minimum_band_gap} bp)")
    if eligible:
        print(f"  selected dosage band={band[0]}..{band[1]} bp; pairs in band={len(constrained)}")
        if args.depth_sensitivity > 0.0:
            sensitivity: list[tuple[float, tuple[int, int]]] = []
            for scale in (1.0 - args.depth_sensitivity, 1.0 + args.depth_sensitivity):
                depth = args.haploid_depth * scale
                z = np.maximum(0.0, (dump.observed - args.background) / depth)
                _, _, _, span = fit_graph_length_calibration(dump, z)
                alternate, _ = select_dosage_band(scores, bands, span)
                sensitivity.append((depth, alternate))
            stable = all(alternate == band for _, alternate in sensitivity)
            print("  depth sensitivity: " + ", ".join(
                f"lambda={depth:.4f} -> {alternate[0]}..{alternate[1]}"
                for depth, alternate in sensitivity))
            print("  dosage-band stability=" + ("PASS" if stable else "FAIL (posterior depth must be integrated)"))
    else:
        print("  dosage-band gate DECLINED: calibration or band separation is insufficient")
    print_scores("unconstrained robust count projection", scores, args.truth_bp, args.top)
    print_scores("dosage-stratified robust count projection", constrained, args.truth_bp, args.top)

    if args.current_pair is not None:
        current = tuple(int(x) for x in args.current_pair.split(","))
        if len(current) != 2:
            raise ValueError("--current-pair must be A,B")
        print(f"current pair {current[0]},{current[1]}: unrestricted rank={rank_of(scores, current)}, "
              f"dosage-band rank={rank_of(constrained, current)}")

    if args.gfa is not None:
        if args.bubbles is None or args.bubble_id is None or len(args.truth_path) != 2:
            raise ValueError("identity check needs --gfa, --bubbles, --bubble-id and two --truth-path values")
        truth = truth_walks_from_gfa(args.gfa, args.bubbles, args.bubble_id, args.truth_path)
        repo_root = Path(__file__).resolve().parents[1]
        edlib = load_edlib(repo_root)
        for rank, selected in enumerate(constrained[:args.identity_top], 1):
            selected_identity = pair_identity(
                edlib, [dump.allele_seq[selected.allele1], dump.allele_seq[selected.allele2]], truth)
            print(f"selected rank {rank} ({selected.allele1},{selected.allele2}) sequence identity: "
                  f"mean={selected_identity[0]:.6f} length-weighted={selected_identity[1]:.6f}")
        if args.current_pair is not None:
            current = tuple(int(x) for x in args.current_pair.split(","))
            current_identity = pair_identity(
                edlib, [dump.allele_seq[current[0]], dump.allele_seq[current[1]]], truth)
            print(f"current pair sequence identity:  mean={current_identity[0]:.6f} "
                  f"length-weighted={current_identity[1]:.6f}")

    if args.truth_bp is not None:
        selected_error = abs(constrained[0].total_bp - args.truth_bp)
        print(f"primary held-out check: selected |bp error|={selected_error} bp")
    return 0


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)
    syn = sub.add_parser("synthetic", help="run copy-number, flank and mosaic assertions")
    syn.set_defaults(function=run_synthetic)

    real = sub.add_parser("real", help="score a panvar --dump-block diagnostic")
    real.add_argument("--conf", type=Path, required=True, help="BLOCK.conf.tsv")
    real.add_argument("--alleles", type=Path, required=True, help="BLOCK.fa")
    real.add_argument("--haploid-depth", type=float, required=True)
    real.add_argument("--background", type=float, default=1.0)
    real.add_argument("--huber-c", type=float, default=1.0)
    real.add_argument("--minimum-r2", type=float, default=0.995)
    real.add_argument("--minimum-band-gap", type=int, default=1000)
    real.add_argument("--depth-sensitivity", type=float, default=0.01,
                      help="relative +/- depth perturbation used to audit dosage-band stability")
    real.add_argument("--truth-bp", type=int, help="evaluation only; never enters scoring")
    real.add_argument("--current-pair", help="production alleles A,B, for comparison only")
    real.add_argument("--top", type=int, default=5)
    real.add_argument("--identity-top", type=int, default=1,
                      help="exactly align this many dosage-stratified candidates")
    real.add_argument("--gfa", type=Path, help="optional GFA for exact sequence-identity comparison")
    real.add_argument("--bubbles", type=Path, help="bubble CSV paired with --gfa")
    real.add_argument("--bubble-id", type=int)
    real.add_argument("--truth-path", action="append", default=[],
                      help="truth path name; provide exactly twice")
    real.set_defaults(function=run_real)
    return ap


def main() -> int:
    args = parser().parse_args()
    return int(args.function(args))


if __name__ == "__main__":
    raise SystemExit(main())
