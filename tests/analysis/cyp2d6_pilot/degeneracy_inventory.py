#!/usr/bin/env python3
"""CYP2D6/D7 pilot diagnostic: block-representation degeneracies of a genotype model (truth-free, model-only).

A model haplotype is one allele per block; its sequence is the concatenation of the block alleles. Spelling is
non-injective when two different allele tuples over the same consecutive blocks i..j spell the same sequence. Then a
caller can identify the sequence without a unique per-block genotype.

Search (exhaustive over every allele combination, not only panel rows): at the first differing block i the two tuples
hold different alleles, one a proper prefix of the other (equal sequences within one block are reported separately as
within-block duplicates). The longer side's surplus (the overhang) is carried block by block: both sides append one
allele; the shorter side must agree with the overhang, which leaves a new overhang on one side or the other. A tuple
pair is an equivalence over i..j when the overhang is empty after block j. Closing is minimal: the search does not extend
a closed pair. Every state is explored; the state count is reported and nothing is capped silently (--max-states
aborts with an error rather than truncating).

Usage: degeneracy_inventory.py <out tsv> <summary tsv> <model name>=<dump prefix> ...
(<dump prefix> from cyp2d6_block_alleles: <prefix>.block_alleles.tsv, .panel.tsv, .summary.tsv)
"""
import csv
import hashlib
import sys

MAX_STATES = 50_000_000


def load(prefix):
    blocks, ids, kinds = [], [], []
    with open(prefix + ".block_alleles.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            b, a = int(r["block_index"]), int(r["allele_index"])
            while len(blocks) <= b:
                blocks.append([])
                ids.append(r["block_id"])
                kinds.append(r["kind"])
            if a != len(blocks[b]):
                sys.exit(f"{prefix}: allele indices out of order at block {b}")
            if len(r["sequence"]) != int(r["length"]):
                sys.exit(f"{prefix}: length mismatch at block {b} allele {a}")
            blocks[b].append(r["sequence"])
    panel = {}
    with open(prefix + ".panel.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            panel[r["template"]] = [int(x) for x in r["alleles"].split(",")]
    fp = dict(l.rstrip("\n").split("\t") for l in open(prefix + ".summary.tsv") if not l.startswith("key"))["model_fingerprint"]
    return blocks, ids, kinds, panel, fp


def search(blocks):
    """Yields (i, j, tuple_a, tuple_b) for every minimal cross-block equivalence; returns the explored state count."""
    found, states = [], 0
    nb = len(blocks)
    for i in range(nb):
        seqs = blocks[i]
        stack = []
        for x, sx in enumerate(seqs):
            for y, sy in enumerate(seqs):
                if x != y and len(sx) < len(sy) and sy.startswith(sx):
                    # side "b" ahead by the overhang: tuple_a is behind
                    stack.append((i + 1, sy[len(sx):], "b", (x,), (y,)))
        while stack:
            b, w, ahead, ta, tb = stack.pop()
            states += 1
            if states > MAX_STATES:
                raise RuntimeError(f"state budget {MAX_STATES} exhausted at start block {i}; inventory incomplete")
            if b >= nb:
                continue          # the chain ends with an unresolved overhang: no equivalence
            for xa, s in enumerate(blocks[b]):          # s: allele appended by the side that is behind
                if len(s) <= len(w):
                    if not w.startswith(s):
                        continue
                    rest_behind = None
                else:
                    if not s.startswith(w):
                        continue
                    rest_behind = s[len(w):]
                for yb, t in enumerate(blocks[b]):      # t: allele appended by the side that is ahead
                    if rest_behind is None:
                        w2, ahead2 = w[len(s):] + t, ahead
                    elif t.startswith(rest_behind):
                        w2, ahead2 = t[len(rest_behind):], ahead
                    elif rest_behind.startswith(t):
                        w2, ahead2 = rest_behind[len(t):], ("a" if ahead == "b" else "b")
                    else:
                        continue
                    na, nb_ = (ta + (xa,), tb + (yb,)) if ahead == "b" else (ta + (yb,), tb + (xa,))
                    if w2 == "":
                        found.append((i, b, na, nb_))
                    else:
                        stack.append((b + 1, w2, ahead2, na, nb_))
    return found, states


def main():
    out_path, summary_path = sys.argv[1:3]
    rows, summary = [], []
    for arg in sys.argv[3:]:
        name, prefix = arg.split("=", 1)
        blocks, ids, kinds, panel, fp = load(prefix)
        within = []
        for b, seqs in enumerate(blocks):
            seen = {}
            for a, s in enumerate(seqs):
                if s in seen:
                    within.append((b, b, (seen[s],), (a,)))
                seen.setdefault(s, a)
        cross, states = search(blocks)
        eq = 0
        for kind, items in (("within_block", within), ("cross_block_compensating", cross)):
            for i, j, ta, tb in items:
                eq += 1
                spelled = "".join(blocks[k][a] for k, a in zip(range(i, j + 1), ta))
                assert spelled == "".join(blocks[k][a] for k, a in zip(range(i, j + 1), tb))
                carriers = {}
                for label, tup in (("a", ta), ("b", tb)):
                    carriers[label] = sorted(t for t, v in panel.items() if tuple(v[i:j + 1]) == tup)
                rows.append({
                    "model": name, "model_fingerprint": fp, "equivalence_id": f"{name}:{eq}", "type": kind,
                    "block_interval": f"{ids[i]}..{ids[j]}", "first_block_index": i, "last_block_index": j,
                    "differing_blocks": ",".join(ids[k] for k, (p, q) in zip(range(i, j + 1), zip(ta, tb)) if p != q),
                    "tuple_a": ",".join(map(str, ta)), "tuple_b": ",".join(map(str, tb)),
                    "tuple_a_lengths": ",".join(str(len(blocks[k][a])) for k, a in zip(range(i, j + 1), ta)),
                    "tuple_b_lengths": ",".join(str(len(blocks[k][a])) for k, a in zip(range(i, j + 1), tb)),
                    "spelled_length": len(spelled), "spelled_sha256": hashlib.sha256(spelled.encode()).hexdigest(),
                    "templates_tuple_a": len(carriers["a"]), "templates_tuple_b": len(carriers["b"]),
                    "template_names_tuple_a": ";".join(carriers["a"]) or "none",
                    "template_names_tuple_b": ";".join(carriers["b"]) or "none",
                    "minimal_compound_interval": f"{ids[i]}..{ids[j]}"})
        intervals = sorted({(r["first_block_index"], r["last_block_index"]) for r in rows if r["model"] == name})
        merged = []
        for lo, hi in intervals:
            if merged and lo <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], hi)
            else:
                merged.append([lo, hi])
        summary.append({"model": name, "model_fingerprint": fp, "blocks": len(blocks),
                        "alleles": sum(len(s) for s in blocks), "within_block_duplicates": len(within),
                        "cross_block_equivalences": len(cross), "search_states": states,
                        "coalescence_intervals": ";".join(f"{ids[lo]}..{ids[hi]}" for lo, hi in merged) or "none"})
    with open(out_path, "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else ["model"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(summary_path, "w") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(summary)
    for s in summary:
        print("\t".join(f"{k}={v}" for k, v in s.items()))


if __name__ == "__main__":
    main()
