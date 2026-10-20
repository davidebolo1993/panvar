#!/usr/bin/env python3
"""M9-F junction-aware proposal score (diagnostic prototype; frozen rule JP1, see DESIGN.md). Not a likelihood, not a
confidence value, not a final-score term: a search heuristic for blocks whose allele-pair proposal emissions all tie.

Rule JP1 (one declared version):
  window   W(h, a) = L_h + allele_a + R_h, with L_h / R_h the last / first K-1 bases of the incumbent homologue h around the
           block (spelled through neighbouring blocks; shorter only at the locus ends). Every K-mer of W overlaps the
           allele or spans its boundary; for an empty allele they span the L|R junction; for a short allele some span both.
  markers  ALL canonical K-mers (K = 31; canonical = lexicographic min of the K-mer and its reverse complement, equivalent
           to the caller's 2-bit min), not closed syncmers -- a DECLARED deviation from the caller's marker scheme.
  inventory  the distinct canonical K-mers over W(h, a) for every allele a and both homologues h (fixed per comparison).
  copies   copies(m | x on h1, y on h2) = n(m, W(1, x)) + n(m, W(2, y)) + e(m), where e(m) counts m in both incumbent
           haplotypes outside the two window spans (occurrences elsewhere in the locus are kept, never assumed unique).
  reads    y(m) = occurrences of canonical m in every read of every pair (both mates), skipping K-mers with non-ACGT bases:
           the caller's count_marker_observations convention applied to all K-mers.
  expected lambda(m) = rate * copies(m) + background, with the caller's own depth fit (mosaic_status
           marker_depth_rate_per_copy, marker_background); applicability to all K-mers is gated empirically per cell.
  score    S(x on h1, y on h2) = sum over the inventory of [y(m) log lambda(m) - lambda(m)] (Poisson kernel; constants
           dropped), computed exactly as a delta over the K-mers the candidate's windows touch.
  pair     the unordered pair {x, y} is ranked by max(S(x on h1, y on h2), S(y on h1, x on h2)); both assignments are kept.
  ties     |S - S'| <= 1e-9 is a tie; ties are reported; the deterministic order inside a tie is (x, y) ascending and is
           NOT evidence that the reads distinguish the tied pairs.
"""
from collections import Counter
import math

K = 31
TIE = 1e-9
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def canon(w):
    r = w.translate(COMP)[::-1]
    return w if w <= r else r


def kmers(s):
    for i in range(len(s) - K + 1):
        w = s[i:i + K]
        if set(w) <= set("ACGT"):
            yield canon(w)


def windows(hap_prefix, hap_suffix, alleles):
    L, R = hap_prefix[-(K - 1):], hap_suffix[:K - 1]
    return [Counter(kmers(L + a + R)) for a in alleles]


def elsewhere(haps, inventory):
    """haps: list of (prefix, incumbent_allele, suffix); counts inventory K-mers outside the window span of each."""
    e = Counter()
    for pre, a, suf in haps:
        full = pre + a + suf
        lo = max(0, len(pre) - (K - 1))           # first K-mer start touching the block region
        hi = len(pre) + len(a)                    # K-mers starting at >= hi lie entirely in the suffix... unless they start
        for i in range(len(full) - K + 1):        # inside L; the window covers starts in [lo, hi) exactly
            if lo <= i < hi:
                continue
            w = full[i:i + K]
            if set(w) <= set("ACGT"):
                c = canon(w)
                if c in inventory:
                    e[c] += 1
    return e


def read_counts(reads, inventory):
    y = Counter()
    for s in reads:
        for c in kmers(s):
            if c in inventory:
                y[c] += 1
    return y


def make_scorer(e, y, rate, bg):
    f = lambda m, c: y.get(m, 0) * math.log(rate * c + bg) - (rate * c + bg)

    def score(w1, w2):
        touched = set(w1) | set(w2)
        return sum(f(m, e.get(m, 0) + w1.get(m, 0) + w2.get(m, 0)) - f(m, e.get(m, 0)) for m in touched)
    return score


def rank_pairs(W1, W2, score):
    """W1[a], W2[a]: window counters per allele on homologue 1 / 2. Returns rows sorted by best score, then (x, y)."""
    rows = []
    n = len(W1)
    for x in range(n):
        for y in range(x, n):
            ab = score(W1[x], W2[y])
            ba = score(W1[y], W2[x]) if x != y else ab
            rows.append(dict(x=x, y=y, score_ab=ab, score_ba=ba, best=max(ab, ba), best_assignment="ab" if ab >= ba else "ba"))
    rows.sort(key=lambda r: (-r["best"], r["x"], r["y"]))
    # explicit tie classes (chained within TIE of the class head)
    cls, head = 0, None
    for i, r in enumerate(rows):
        if head is None or head - r["best"] > TIE:
            cls += 1
            head = r["best"]
        r["rank"], r["tie_class"] = i + 1, cls
    sizes = Counter(r["tie_class"] for r in rows)
    for r in rows:
        r["tie_class_size"] = sizes[r["tie_class"]]
    return rows


def oracle_score(full_haps, inventory, reads, rate, bg):
    """Independent straightforward implementation: spell whole haplotypes, count every inventory K-mer in them and in the
    reads by direct scanning, and sum the Poisson kernel over the whole inventory."""
    copies = {m: 0 for m in inventory}
    for h in full_haps:
        for i in range(len(h) - K + 1):
            w = h[i:i + K]
            r = w.translate(COMP)[::-1]
            c = min(w, r)
            if c in copies and set(w) <= set("ACGT"):
                copies[c] += 1
    obs = {m: 0 for m in inventory}
    for s in reads:
        for i in range(len(s) - K + 1):
            w = s[i:i + K]
            if not set(w) <= set("ACGT"):
                continue
            r = w.translate(COMP)[::-1]
            c = min(w, r)
            if c in obs:
                obs[c] += 1
    return sum(obs[m] * math.log(rate * copies[m] + bg) - (rate * copies[m] + bg) for m in inventory)
