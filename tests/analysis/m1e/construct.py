#!/usr/bin/env python3
"""M1-E: the frozen residual-sequence restoration set for one coverage (C4 NA19650 sv:7; ORACLE sensitivity interventions,
not inferred corrections). Reported flanks and the shared prefix are fixed.

Starting point: the M1-D oracle composition, components P (prefix_u6, = both truth prefixes exactly), U1 = period_u13,
T1 = tail_u44 (truth1) and U2 = period_u35, T2 = tail_u33 (truth2). A patch REPLACES a whole component with the truth
component's bytes (M1-D truth partition with the panel anchor); nothing else changes.
Patch sets: none, U1, T1, U2, T2, all. For each, the identical patched inventory is arranged two ways:
  a22_<patch>   side 1 = P + U1 + T1, side 2 = P + U2 + T2   (2+2)
  e13_<patch>   side 1 = P + T1,      side 2 = P + U1 + U2 + T2   (1+3, the M1-D matched-inventory arrangement)
Controls: reported (saved total as expectation) and truthblock (truth sv:7 bytes; must equal a22_all byte for byte, so the
scorer interns it once). 14 declared diplotypes. Asserted: component identities (inventory SHA-256, truth components),
a22_none == the M1-D composition, matched inventories per patch set, every side = reported left flank + block + reported
right flank, a22_all == truthblock. Recorded: lengths, SHA-256 and the 62-bp context of every join.

  construct.py --inventory P --structure F --truth-components F --truth-blocks F --composed F --spell P --expected-reported X --out P
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
    for x in ("inventory", "structure", "truth_components", "truth_blocks", "composed", "spell", "expected_reported", "out"):
        ap.add_argument("--" + x.replace("_", "-"), required=True)
    a = ap.parse_args()
    inv = fasta(a.inventory + ".fa")
    for r in csv.DictReader(open(a.inventory + ".tsv"), delimiter="\t"):
        assert sha(inv[r["unique_record"]]) == r["sha256"], "inventory record does not match its source component"
    struct = {r["record"]: r["components"].split("+") for r in csv.DictReader(open(a.structure), delimiter="\t")}
    p1, u1, t1 = struct["compose.truth1"]
    p2, u2, t2 = struct["compose.truth2"]
    assert p1 == p2, "the two compositions use different prefixes"
    tc = fasta(a.truth_components)
    tb = fasta(a.truth_blocks)
    composed = fasta(a.composed)
    P = inv[p1]
    assert tc["truth1.prefix"] == P and tc["truth2.prefix"] == P, "truth prefixes are not the shared panel prefix"
    panel = {"U1": inv[u1], "T1": inv[t1], "U2": inv[u2], "T2": inv[t2]}
    truth = {"U1": tc["truth1.period1"], "T1": tc["truth1.tail"], "U2": tc["truth2.period1"], "T2": tc["truth2.tail"]}
    names = {"U1": u1, "T1": t1, "U2": u2, "T2": t2}
    seq = fasta(a.spell + ".fa")
    iv = {}
    for line in open(a.spell + ".spelled.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "record":
            iv[f[0]] = (int(f[3]), int(f[4]))
    flank = {s: (seq[f"reported.{s}"][:iv[f"reported.{s}"][0]], seq[f"reported.{s}"][iv[f"reported.{s}"][1]:]) for s in (1, 2)}
    log, alts = [], []

    def put(label, b1, b2, parts, expected="NA"):
        sides = []
        for s, b in ((1, b1), (2, b2)):
            left, right = flank[s]
            full = left + b + right
            assert full.startswith(left) and full.endswith(right) and full[len(left):len(full) - len(right)] == b
            sides.append((full, (len(left), len(left) + len(b))))
        alts.append((label, sides, b1, b2, expected))
        log.append(f"{label}\tblocks_bp\t{len(b1)},{len(b2)}\n{label}\tcomponents\t{parts}\n{label}\tsha256\t{sha(sides[0][0])},{sha(sides[1][0])}")

    alts.append(("reported", [(seq["reported.1"], iv["reported.1"]), (seq["reported.2"], iv["reported.2"])],
                 seq["reported.1"][iv["reported.1"][0]:iv["reported.1"][1]], seq["reported.2"][iv["reported.2"][0]:iv["reported.2"][1]],
                 a.expected_reported))
    for comp in ("U1", "T1", "U2", "T2"):
        log.append(f"component_{comp}\tpanel {names[comp]} {len(panel[comp])} bp sha {sha(panel[comp])[:16]}; truth {len(truth[comp])} bp "
                   f"sha {sha(truth[comp])[:16]}")
    for patch in ("none", "U1", "T1", "U2", "T2", "all"):
        c = {k: (truth[k] if patch in (k, "all") else panel[k]) for k in panel}
        src = {k: ("truth" if patch in (k, "all") else names[k]) for k in panel}
        b22 = (P + c["U1"] + c["T1"], P + c["U2"] + c["T2"])
        b13 = (P + c["T1"], P + c["U1"] + c["U2"] + c["T2"])
        assert collections.Counter([P, c["U1"], c["T1"], P, c["U2"], c["T2"]]) == collections.Counter([P, c["T1"], P, c["U1"], c["U2"], c["T2"]])
        assert len(b22[0]) + len(b22[1]) == len(b13[0]) + len(b13[1])
        if patch == "none":
            assert b22 == (composed["compose.truth1"], composed["compose.truth2"]), "a22_none is not the M1-D composition"
        if patch == "all":
            assert b22 == (tb["truth1"], tb["truth2"]), "fully corrected 2+2 blocks differ from the truth-block reference"
        for arr, (x, y), parts in (("a22", b22, f"P+U1[{src['U1']}]+T1[{src['T1']}] | P+U2[{src['U2']}]+T2[{src['T2']}]"),
                                   ("e13", b13, f"P+T1[{src['T1']}] | P+U1[{src['U1']}]+U2[{src['U2']}]+T2[{src['T2']}]")):
            seams = []
            for block, pieces in ((x, ([P, c["U1"], c["T1"]] if arr == "a22" else [P, c["T1"]])),
                                  (y, ([P, c["U2"], c["T2"]] if arr == "a22" else [P, c["U1"], c["U2"], c["T2"]]))):
                at = 0
                for piece in pieces[:-1]:
                    at += len(piece)
                    seams.append(block[at - 31:at + 31])
            log.append(f"{arr}_{patch}\tjoins\t{';'.join(seams)}")
            put(f"{arr}_{patch}", x, y, parts)
    put("truthblock", tb["truth1"], tb["truth2"], "truth1 | truth2 (reference; = a22_all)")
    assert len(alts) <= 14
    with open(a.out + ".fa", "w") as fa, open(a.out + ".explicit.tsv", "w") as tsv, open(a.out + ".blocks.fa", "w") as bf:
        tsv.write("# label\tfasta\trecord1\trecord2\tbegin1\tend1\tbegin2\tend2\texpected_total\n")
        for label, sides, b1, b2, expected in alts:
            fa.write(f">{label}.1\n{sides[0][0]}\n>{label}.2\n{sides[1][0]}\n")
            tsv.write(f"{label}\t{a.out}.fa\t{label}.1\t{label}.2\t{sides[0][1][0]}\t{sides[0][1][1]}\t{sides[1][1][0]}\t{sides[1][1][1]}\t{expected}\n")
            bf.write(f">{label}.1\n{b1}\n>{label}.2\n{b2}\n")
        bf.write(f">truth.1\n{tb['truth1']}\n>truth.2\n{tb['truth2']}\n")
    open(a.out + ".construction.tsv", "w").write("\n".join(log) + "\n")
    print("\n".join(l for l in log if "joins" not in l and "sha256" not in l))
    print(f"declared {len(alts)} diplotypes")


if __name__ == "__main__":
    main()
