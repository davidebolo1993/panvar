#!/usr/bin/env python3
"""M1-D: the frozen diagnostic candidate set for one coverage (C4 NA19650 sv:7), reported flanks fixed.

The ORACLE label applies to every composition here: components were selected as nearest to truth's units. This is a
representation-feasibility test, not a truth-free candidate generator, and not evidence that a caller would find them.

Candidates (at most 12; side 1 = the reported 3-unit homologue, side 2 = the reported 1-unit homologue):
  reported                 the reported 3+1 pair (M1-B rank-1), saved total as expectation
  saved_2+2                the M1-B selected panel 2+2 comparator, saved total as expectation
  truthblock_a             truth sv:7 bytes in the reported flanks (truth1 -> side 1), the truth-block reference
  whole_nearest_a          the nearest existing whole sv:7 alleles W(truth1) -> side 1, W(truth2) -> side 2
  periodonly_a             W(truthH) with only its period replaced by the nearest panel period (if constructed)
  oracle_compose_a         nearest prefix + period + tail for truth1 -> side 1, for truth2 -> side 2 (2+2)
  oracle_compose_b         the same two composed blocks exchanged between homologues (2+2)
  oracle_compose_c         the two composed blocks with their periods exchanged (2+2)
  oracle_compose_d31       truth2's composed period moved to side 1 after side 1's period (3+1)
  oracle_compose_e13       truth1's composed period moved to side 2 before side 2's period (1+3)
c/d31/e13 use exactly the composed a/b unit inventory (asserted). Every candidate side = reported left flank + block +
reported right flank (asserted); every component equals its frozen inventory record (SHA-256 asserted); every join's
62-bp context is recorded.

  construct.py --spell PREFIX --inventory PREFIX --reach PREFIX --truth-blocks FA --out OUT_PREFIX --expected-reported X --expected-comparator Y
"""
import argparse
import collections
import csv
import hashlib


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


def main():
    ap = argparse.ArgumentParser()
    for x in ("spell", "inventory", "reach", "truth_blocks", "out", "expected_reported", "expected_comparator"):
        ap.add_argument("--" + x.replace("_", "-"), required=True)
    a = ap.parse_args()
    seq = fasta(a.spell + ".fa")
    iv = {}
    for line in open(a.spell + ".spelled.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "record":
            iv[f[0]] = (int(f[3]), int(f[4]))
    inv = fasta(a.inventory + ".fa")
    inv_sha = {}
    for r in csv.DictReader(open(a.inventory + ".tsv"), delimiter="\t"):
        inv_sha.setdefault(r["unique_record"], set()).add(r["sha256"])
    for rec, shas in inv_sha.items():
        assert shas == {sha(inv[rec])}, f"inventory record {rec} does not match its source components"
    anchor = next(r["anchor"] for r in csv.DictReader(open(a.inventory + ".anchors.tsv"), delimiter="\t") if r["phase"] == "A")
    composed = fasta(a.reach + ".composed.fa")
    structure = {r["record"]: r["components"].split("+") for r in csv.DictReader(open(a.reach + ".composed_structure.tsv"), delimiter="\t")}
    blk = {r: seq[r][iv[r][0]:iv[r][1]] for r in seq}
    # Amendment 2: the reachability step read truth from the M1-B block file (records truth1/truth2); the spelled truth
    # (truth.1/truth.2) must carry the identical sv:7 bytes.
    reach_truth = fasta(a.truth_blocks)
    assert blk["truth.1"] == reach_truth["truth1"] and blk["truth.2"] == reach_truth["truth2"], "the two truth sources differ"
    assert blk["reported.1"].count(anchor) == 3 and blk["reported.2"].count(anchor) == 1, "side 1 must be the 3-unit homologue"
    flank = {s: (seq[f"reported.{s}"][:iv[f"reported.{s}"][0]], seq[f"reported.{s}"][iv[f"reported.{s}"][1]:]) for s in (1, 2)}
    log, alts = [], []

    def put(label, b1, b2, parts1="-", parts2="-", expected="NA"):
        sides = []
        for s, b in ((1, b1), (2, b2)):
            left, right = flank[s]
            full = left + b + right
            assert full[:len(left)] == left and full[len(left) + len(b):] == right and full[len(left):len(left) + len(b)] == b
            sides.append((full, (len(left), len(left) + len(b))))
        alts.append((label, sides, expected, b1, b2))
        log.append(f"{label}\tblocks_bp\t{len(b1)},{len(b2)}\n{label}\tcomponents\t{parts1} | {parts2}\n"
                   f"{label}\tsha256\t{sha(sides[0][0])},{sha(sides[1][0])}")

    def comp(record_names):
        for r in record_names:
            assert r in inv, f"{r} is not a frozen inventory record"
        return "".join(inv[r] for r in record_names)

    # reported and the comparator as spelled (the comparator keeps its own flanks: it is a whole panel diplotype)
    alts.append(("reported", [(seq["reported.1"], iv["reported.1"]), (seq["reported.2"], iv["reported.2"])],
                 a.expected_reported, blk["reported.1"], blk["reported.2"]))
    alts.append(("saved_2+2", [(seq["saved_2+2.1"], iv["saved_2+2.1"]), (seq["saved_2+2.2"], iv["saved_2+2.2"])],
                 a.expected_comparator, blk["saved_2+2.1"], blk["saved_2+2.2"]))
    put("truthblock_a", blk["truth.1"], blk["truth.2"], "truth.1 block", "truth.2 block")
    put("whole_nearest_a", composed["whole.truth1"], composed["whole.truth2"], structure["whole.truth1"][0],
        structure["whole.truth2"][0])
    if "periodonly.truth1" in composed and "periodonly.truth2" in composed:
        put("periodonly_a", composed["periodonly.truth1"], composed["periodonly.truth2"],
            "+".join(structure["periodonly.truth1"]), "+".join(structure["periodonly.truth2"]))
    else:
        log.append("periodonly_a\tNOT CONSTRUCTED: a nearest whole allele has a different unit count")
    c1, c2 = structure["compose.truth1"], structure["compose.truth2"]
    if len(c1) == 3 and len(c2) == 3:
        p1, u1, t1 = c1
        p2, u2, t2 = c2
        assert comp(c1) == composed["compose.truth1"] and comp(c2) == composed["compose.truth2"]
        variants = {"oracle_compose_a": ([p1, u1, t1], [p2, u2, t2]), "oracle_compose_b": ([p2, u2, t2], [p1, u1, t1]),
                    "oracle_compose_c": ([p1, u2, t1], [p2, u1, t2]), "oracle_compose_d31": ([p1, u1, u2, t1], [p2, t2]),
                    "oracle_compose_e13": ([p1, t1], [p2, u1, u2, t2])}
        inventory = collections.Counter(c1 + c2)
        for label, (x, y) in variants.items():
            assert collections.Counter(x + y) == inventory, f"{label} changes the unit inventory"
            b1, b2 = comp(x), comp(y)
            for parts, b in ((x, b1), (y, b2)):
                joins, at = [], 0
                for r in parts[:-1]:
                    at += len(inv[r])
                    joins.append(b[at - 31:at + 31])
                log.append(f"{label}\tjoins\t{';'.join(joins)}")
            put(label, b1, b2, "+".join(x), "+".join(y))
    else:
        log.append("oracle_compose\tNOT CONSTRUCTED: a truth haplotype does not have exactly two units under the anchor")
    assert len(alts) <= 12
    with open(a.out + ".fa", "w") as fa, open(a.out + ".explicit.tsv", "w") as tsv, open(a.out + ".blocks.fa", "w") as bf:
        tsv.write("# label\tfasta\trecord1\trecord2\tbegin1\tend1\tbegin2\tend2\texpected_total\n")
        for label, sides, expected, b1, b2 in alts:
            fa.write(f">{label}.1\n{sides[0][0]}\n>{label}.2\n{sides[1][0]}\n")
            tsv.write(f"{label}\t{a.out}.fa\t{label}.1\t{label}.2\t{sides[0][1][0]}\t{sides[0][1][1]}\t{sides[1][1][0]}\t"
                      f"{sides[1][1][1]}\t{expected}\n")
            bf.write(f">{label}.1\n{b1}\n>{label}.2\n{b2}\n")
        bf.write(f">truth.1\n{blk['truth.1']}\n>truth.2\n{blk['truth.2']}\n")
    with open(a.out + ".pairs.tsv", "w") as p:
        for label, *_ in alts:
            for s in (1, 2):
                for t in (1, 2):
                    p.write(f"{label}.{s}:truth.{t}\t{a.out}.blocks.fa\t{label}.{s}\t{a.out}.blocks.fa\ttruth.{t}\n")
    open(a.out + ".construction.tsv", "w").write("\n".join(log) + "\n")
    print("\n".join(log))
    print(f"declared {len(alts)} candidates")


if __name__ == "__main__":
    main()
