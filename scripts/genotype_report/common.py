"""Shared vocabulary, table I/O and small statistics for the panvar genotype reporting tool.

Everything here is deterministic: TSVs are read as strings (no NA inference), written with a fixed float format and
line terminator, and every ordering is explicit.
"""
import hashlib
import math
import os
from decimal import ROUND_HALF_UP, Decimal

import pandas as pd

# The headline filters' values, in the fail-closed precedence order of genotype_headline_filters.md section 2
# (worst first), and the call classes of the locus table.
STATUSES = ["NO_CALL", "NOT_EVALUATED", "LOW_QUALITY", "PASS", "NOT_APPLICABLE"]
CALL_CLASSES = ["strict_call", "genotype_equivalent_call", "no_call", "NA"]
CALLED = ("strict_call", "genotype_equivalent_call")

# The documented usage policies (genotype_headline_filters.md section 7). They are downstream selections, not filters.
EXACT_GENOTYPE_USABLE_RULE = "GENOTYPE_FILTER == PASS and DOSAGE_FILTER in {PASS, NOT_APPLICABLE}"
DOSAGE_CONTENT_USABLE_RULE = "DOSAGE_FILTER == PASS"
PHASED_GENOTYPE_USABLE_RULE = "EXACT_GENOTYPE_USABLE and PHASE_FILTER == PASS"

# The normalised-excess denominator, as the sealed C4 evaluation package defines it (evaluate_integration.py): a
# deterministic diploid length scale from the evaluator's recorded allele lengths, the larger of the two haplotypes'
# summed truth lengths and summed called lengths. It is not an observed alignment length.
NORMALIZED_EXCESS_FORMULA = "excess / max(truth1_bp + truth2_bp, called1_bp + called2_bp)"
NORMALIZATION_COLUMNS = ["truth1_bp", "truth2_bp", "called1_bp", "called2_bp"]

FLOAT_FORMAT = None   # the shortest representation that round-trips (repr): exact and deterministic

# The caption of every model-order block figure whose block labels carry the observed truth-allele length range.
BLOCK_LENGTH_CAPTION = ("Length below each block ID is the observed truth-allele length range across truth-defined "
                        "evaluation observations. It is not a reference span and may not include alleles absent from "
                        "the evaluated dataset. Exact minimum, median and maximum values are in data/block_lengths.tsv.")


class Refusal(Exception):
    """An input the tool will not report on: a malformed schema, mixed incompatible runs, or a request that needs an
    input that was not supplied. Raised before any output is written where possible."""


def read_tsv(path, what=None):
    if not os.path.isfile(path):
        raise Refusal(f"{what or 'table'} not found: {path}")
    try:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_filter=False)
    except pd.errors.EmptyDataError:
        raise Refusal(f"{what or 'table'} is empty: {path}")


def require_columns(frame, columns, what):
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise Refusal(f"{what}: missing required column(s) {', '.join(missing)}")


def write_tsv(frame, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, float_format=FLOAT_FORMAT, lineterminator="\n", na_rep="NA")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def attach_run_id(frame, what):
    """Every joined table names its run by `run_id`, or by `dataset` and `job` (`cell` is accepted for `job`), which
    become run_id = dataset/job: the relative prefix path a directory of dataset subdirectories yields."""
    frame = frame.copy()
    if "run_id" in frame.columns:
        return frame
    job = "job" if "job" in frame.columns else "cell" if "cell" in frame.columns else None
    if job is None:
        raise Refusal(f"{what}: no run key (run_id, or dataset and job)")
    frame["run_id"] = (frame["dataset"] + "/" + frame[job]) if "dataset" in frame.columns else frame[job]
    return frame


def attach_block_id(frame, what):
    frame = frame.copy()
    if "block_id" not in frame.columns:
        if "block" not in frame.columns:
            raise Refusal(f"{what}: no block_id column")
        frame["block_id"] = frame["block"]
    return frame


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def to_number(series):
    """Exact decimal-to-float parsing (Python's float): pandas' fast string parser can round the last digits, which
    would make a table read back from disk differ from the one written. NA, empty and 'none' read as NaN."""
    return series.map(_float).astype(float)


def to_bool(series):
    """True/False/NA from the spellings the evaluators write; anything else is a schema error."""
    mapping = {"True": True, "true": True, "1": True, "False": False, "false": False, "0": False,
               "NA": None, "": None, "nan": None}
    bad = sorted(set(series) - set(mapping))
    if bad:
        raise Refusal(f"not a boolean: {', '.join(bad[:5])}")
    return series.map(mapping)


def wilson(k, n, z=1.959963984540054):
    """Wilson score interval for k successes of n (95% by default); (nan, nan) when n is 0."""
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def ordered_unique(values):
    seen, out = set(), []
    for v in values:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def coverage_order(values):
    """Coverages sort numerically when every value is a number, otherwise in order of first appearance."""
    vals = ordered_unique(values)
    try:
        return sorted(vals, key=lambda v: float(v))
    except ValueError:
        return vals


def stable_jitter(key, width):
    """A jitter offset in [-width/2, width/2) from a hash of the point's identity, so a point's position never depends
    on row order or a random state."""
    h = int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
    return (h / 0xFFFFFFFF - 0.5) * width


def _plain(d):
    """A Decimal in plain notation without a trailing exponent or superfluous zeros ('117', '117.5')."""
    return format(d.normalize(), "f")


def format_length_range(lo, hi):
    """A compact length range for a block label; the exact values stay in the table.
    Equal ends: one value ('117 bp', '48.0 kb'). Both below 1,000 bp: bp ('1–69 bp'). Otherwise kb with one decimal,
    rounded half up, with an exact 0 written '0' ('12.3–12.4 kb', '0–7.7 kb'); when two different ends would
    print identically at one decimal, decimals are added (up to three, exact for whole bp) until they differ."""
    a, b = Decimal(repr(float(lo))), Decimal(repr(float(hi)))

    def kb(v, places):
        return "0" if v == 0 else str((v / 1000).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP))

    if a == b:
        return f"{_plain(a)} bp" if a < 1000 else f"{kb(a, 1)} kb"
    if a < 1000 and b < 1000:
        return f"{_plain(a)}–{_plain(b)} bp"
    places = 1
    while places < 3 and kb(a, places) == kb(b, places):
        places += 1
    return f"{kb(a, places)}–{kb(b, places)} kb"
