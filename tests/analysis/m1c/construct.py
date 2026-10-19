#!/usr/bin/env python3
"""M1-C: construct the frozen C4 allocation alternatives (diagnostic evaluator only; truth allowed here, never in a caller).

Input: one `genotype_probe spell` output (FASTA + spelled.tsv) holding `reported.{1,2}`, `saved_2+2.{1,2}` (the M1-B
selected 2+2 comparator) and `truth.{1,2}` (the donor's truth haplotypes projected onto the block chain), all in the model
frame with the target-block interval. Output: one FASTA and one `score-set --explicit` table with at most 12 declared
alternatives, plus a construction record.

Alternatives (fixed before scoring):
  reported, saved_2+2         spelled panel paths, with their saved totals as expectations (must reproduce)
  realloc_<phase>_p<k>        the reported pair's own units reallocated 3+1 -> 2+2 (below), k = 1, 2; phases A, B
  truth                       the actual whole-locus truth pair (diagnostic reference)
  truthblock_in_reported_a/b  the reported flanks with the truth target-block bytes (truth1 into the 3-unit side; swapped)
  reportedblock_in_truth_a/b  the truth flanks with the reported target-block bytes (the 3-unit allele into truth1; swapped)

Reallocation (an exact partition; nothing duplicated, discarded or reversed):
  * an ANCHOR is a 31-mer (forward strand of the model frame) occurring exactly 3 times in the 3-unit side's target block,
    exactly once in the 1-unit side's block, exactly twice in each truth block, and nowhere else on either reported
    haplotype; its consecutive spacings on the 3-unit side must agree within 5 %. Phase A = the first anchor by position in
    the 1-unit block; phase B = the first anchor whose position is at least half that block's length beyond phase A's.
  * with anchor starts a1 < a2 < a3 on the 3-unit haplotype H and q on the 1-unit haplotype G, the full periods are
    P1 = H[a1:a2] and P2 = H[a2:a3]. realloc moves P_k: H' = H[:a_k] + H[a_(k+1):], G' = G[:q] + P_k + G[q:].
  * before any change, H[:a_k] + P_k + H[a_(k+1):] == H and G[:q] + G[q:] == G are verified; afterwards the base
    composition of H' + G' equals that of H + G, the total length is preserved, and the anchor counts become 2 + 2.
  * new junctions are recorded (31 bp either side of each cut); their content differences are measured, not assumed away.
If no anchor qualifies, no reallocation is constructed and the limitation is recorded.

  construct.py --spell PREFIX --out OUT_PREFIX --expected-reported X --expected-comparator Y
"""
import argparse
import collections
import hashlib

K = 31


def fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line)
    return {k: "".join(v) for k, v in out.items()}


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def kcount(s):
    c = collections.Counter()
    for i in range(len(s) - K + 1):
        c[s[i:i + K]] += 1
    return c


def positions(s, k):
    out, i = [], s.find(k)
    while i >= 0:
        out.append(i)
        i = s.find(k, i + 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spell", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--expected-reported", required=True)
    ap.add_argument("--expected-comparator", required=True)
    a = ap.parse_args()
    seq = fasta(a.spell + ".fa")
    iv = {}
    for line in open(a.spell + ".spelled.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "record":
            iv[f[0]] = (int(f[3]), int(f[4]))
    blk = {r: seq[r][iv[r][0]:iv[r][1]] for r in seq}
    log = [f"source\t{a.spell}"]
    # Which reported side carries more units is decided by anchors, not lengths: try both orientations.
    found = None
    for h, g in (("reported.1", "reported.2"), ("reported.2", "reported.1")):
        H, G = seq[h], seq[g]
        hb, gb = blk[h], blk[g]
        whole_h, whole_g = kcount(H), kcount(G)
        t1, t2 = kcount(blk["truth.1"]), kcount(blk["truth.2"])
        bh = kcount(hb)
        anchors = []
        for i in range(len(gb) - K + 1):
            k = gb[i:i + K]
            if whole_g[k] != 1 or bh[k] != 3 or whole_h[k] != 3 or t1[k] != 2 or t2[k] != 2:
                continue
            p = positions(hb, k)
            gaps = [p[1] - p[0], p[2] - p[1]]
            if abs(gaps[0] - gaps[1]) > 0.05 * max(gaps):
                continue
            anchors.append((i, k, p))
        if anchors:
            found = (h, g, anchors)
            break
    alts = []   # (label, rec1, seq1, b1, e1, rec2, seq2, b2, e2, expected)

    def add(label, s1, i1, s2, i2, expected="NA"):
        alts.append((label, s1, i1, s2, i2, expected))

    add("reported", seq["reported.1"], iv["reported.1"], seq["reported.2"], iv["reported.2"], a.expected_reported)
    add("saved_2+2", seq["saved_2+2.1"], iv["saved_2+2.1"], seq["saved_2+2.2"], iv["saved_2+2.2"], a.expected_comparator)
    records = []
    if found is None:
        log.append("reallocation\tNOT CONSTRUCTED: no anchor satisfies the declared unit-boundary rule")
    else:
        h, g, anchors = found
        H, G = seq[h], seq[g]
        (hb0, he0), (gb0, ge0) = iv[h], iv[g]
        first = anchors[0]
        second = next((x for x in anchors if x[0] >= first[0] + len(blk[g]) // 2), None)
        log.append(f"three_unit_side\t{h}\none_unit_side\t{g}\nqualifying_anchors\t{len(anchors)}")
        for phase, anchor in (("A", first), ("B", second)):
            if anchor is None:
                log.append(f"phase_{phase}\tNOT CONSTRUCTED: no anchor at least half a unit beyond phase A")
                continue
            gi, k, p = anchor
            cuts = [hb0 + x for x in p]
            q = gb0 + gi
            log.append(f"phase_{phase}_anchor\t{k}\nphase_{phase}_anchor_offset_in_one_unit_block\t{gi}\n"
                       f"phase_{phase}_cuts_three_unit_haplotype\t{','.join(map(str, cuts))}\nphase_{phase}_cut_one_unit_haplotype\t{q}\n"
                       f"phase_{phase}_period_bp\t{cuts[1] - cuts[0]},{cuts[2] - cuts[1]}")
            for kk in (1, 2):
                lo, hi = cuts[kk - 1], cuts[kk]
                moved = H[lo:hi]
                assert H[:lo] + moved + H[hi:] == H, "partition of the 3-unit haplotype does not reconstruct it"
                assert G[:q] + G[q:] == G
                H2, G2 = H[:lo] + H[hi:], G[:q] + moved + G[q:]
                assert len(H2) + len(G2) == len(H) + len(G)
                assert collections.Counter(H2 + G2) == collections.Counter(H + G), "base composition changed"
                assert H2.count(k) == 2 and G2.count(k) == 2, "anchor counts are not 2 + 2 after the move"
                ih, ig = (hb0, he0 - len(moved)), (gb0, ge0 + len(moved))
                label = f"realloc_{phase}_p{kk}"
                junctions = [("three_unit_side", H2[lo - K:lo + K]), ("one_unit_side_left", G2[q - K:q + K]),
                             ("one_unit_side_right", G2[q + len(moved) - K:q + len(moved) + K])]
                records.append(f"{label}\tmoved\t{h}[{lo}:{hi}) -> {g} at {q}\n{label}\tmoved_bp\t{len(moved)}\n"
                               f"{label}\tmoved_sha256\t{sha(moved)}\n{label}\tlengths\t{len(H2)},{len(G2)} (from {len(H)},{len(G)})\n"
                               + "".join(f"{label}\tnew_junction_{n}\t{j}\n" for n, j in junctions)
                               + f"{label}\tsha256\t{sha(H2)},{sha(G2)}\n")
                # Keep homologue order as reported: side 1 stays side 1.
                if h == "reported.1":
                    add(label, H2, ih, G2, ig)
                else:
                    add(label, G2, ig, H2, ih)
    add("truth", seq["truth.1"], iv["truth.1"], seq["truth.2"], iv["truth.2"])
    r1, r2, tr1, tr2 = "reported.1", "reported.2", "truth.1", "truth.2"
    three = found[0] if found else r1
    one = r2 if three == r1 else r1

    def splice(base, block_bytes):
        b, e = iv[base]
        s = seq[base][:b] + block_bytes + seq[base][e:]
        return s, (b, b + len(block_bytes))

    for suffix, (x, y) in (("a", (tr1, tr2)), ("b", (tr2, tr1))):
        s3, i3 = splice(three, blk[x])
        s1, i1 = splice(one, blk[y])
        if three == r1:
            add(f"truthblock_in_reported_{suffix}", s3, i3, s1, i1)
        else:
            add(f"truthblock_in_reported_{suffix}", s1, i1, s3, i3)
    for suffix, (x, y) in (("a", (tr1, tr2)), ("b", (tr2, tr1))):
        sx, ix = splice(x, blk[three])
        sy, iy = splice(y, blk[one])
        if x == tr1:
            add(f"reportedblock_in_truth_{suffix}", sx, ix, sy, iy)
        else:
            add(f"reportedblock_in_truth_{suffix}", sy, iy, sx, ix)
    assert len(alts) <= 12, "more than 12 alternatives"
    with open(a.out + ".fa", "w") as fa, open(a.out + ".explicit.tsv", "w") as tsv:
        tsv.write("# label\tfasta\trecord1\trecord2\tbegin1\tend1\tbegin2\tend2\texpected_total\n")
        for label, s1, i1, s2, i2, expected in alts:
            fa.write(f">{label}.1\n{s1}\n>{label}.2\n{s2}\n")
            tsv.write(f"{label}\t{a.out}.fa\t{label}.1\t{label}.2\t{i1[0]}\t{i1[1]}\t{i2[0]}\t{i2[1]}\t{expected}\n")
    with open(a.out + ".construction.tsv", "w") as f:
        f.write("\n".join(log) + "\n" + "".join(records))
        for label, s1, i1, s2, i2, _ in alts:
            f.write(f"{label}\tbp\t{len(s1)},{len(s2)}\n{label}\tblock\t{i1[0]}-{i1[1]},{i2[0]}-{i2[1]}\n"
                    f"{label}\tsha256\t{sha(s1)},{sha(s2)}\n")
    print(open(a.out + ".construction.tsv").read(), end="")
    print(f"declared {len(alts)} alternatives")


if __name__ == "__main__":
    main()
