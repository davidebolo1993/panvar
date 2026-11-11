#!/usr/bin/env python3
"""SAM files holding the pairs of a FASTQ pair, for checking that recruit reads the same fragments from
both: one name-sorted, one coordinate-sorted. Positions are made up (the reads are not aligned), placed to
exercise the pairing: most mates 250 bp apart, some far apart on the same contig, some on another contig,
some with an unmapped mate placed beside it, some pairs unmapped at the end; reverse-strand records hold
the reverse complement; some records are duplicates (kept); extra secondary, supplementary and QC-failed
records (skipped) are added.

    reads_sam.py <reads_1.fq> <reads_2.fq> <name_sorted.sam> <coordinate_sorted.sam> [seed]
"""
import random
import sys

COMP = str.maketrans("ACGTN", "TGCAN")


def fastq(path):
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline().strip()
            f.readline()
            f.readline()
            name = h[1:].split()[0]
            if name[-2:] in ("/1", "/2"):
                name = name[:-2]
            yield name, s


def record(name, flag, rname, pos, mapq, rnext, pnext, seq):
    if flag & 16:
        seq = seq.translate(COMP)[::-1]
    cigar = "*" if flag & 4 else f"{len(seq)}M"
    return [name, flag, rname, pos, mapq, cigar, rnext, pnext, 0, seq, "I" * len(seq)]


def main():
    r1, r2, out_name, out_coord = sys.argv[1:5]
    rng = random.Random(int(sys.argv[5]) if len(sys.argv) > 5 else 7)
    contigs = [("chrA", 20_000_000), ("chrB", 20_000_000)]
    recs = []
    for (n1, s1), (n2, s2) in zip(fastq(r1), fastq(r2)):
        assert n1 == n2
        kind = rng.random()
        c1 = rng.choice(contigs)
        p1 = rng.randrange(1, c1[1] - 2_000_000)
        rev = rng.random() < 0.5
        f1, f2 = 1 | 64 | (16 if rev else 32), 1 | 128 | (32 if rev else 16)
        if rng.random() < 0.05:
            f1 |= 1024
            f2 |= 1024
        if kind < 0.70:        # mates near each other
            c2, p2 = c1, p1 + 250
        elif kind < 0.80:      # far on the same contig
            c2, p2 = c1, p1 + rng.randrange(150_000, 1_500_000)
        elif kind < 0.90:      # another contig
            c2 = contigs[1] if c1 == contigs[0] else contigs[0]
            p2 = rng.randrange(1, c2[1] - 1000)
        elif kind < 0.95:      # mate 2 unmapped, placed with mate 1
            c2, p2 = c1, p1
            f2 = (f2 | 4) & ~16 & ~32
            f1 |= 8
        else:                  # both unmapped: the end of the file
            f1, f2 = 1 | 4 | 8 | 64, 1 | 4 | 8 | 128
            recs.append(record(n1, f1, "*", 0, 0, "*", 0, s1))
            recs.append(record(n2, f2, "*", 0, 0, "*", 0, s2))
            continue
        same = "=" if c1 == c2 else c2[0]
        same_back = "=" if c1 == c2 else c1[0]
        recs.append(record(n1, f1, c1[0], p1, 60, same, p2, s1))
        recs.append(record(n2, f2, c2[0], p2, 0 if f2 & 4 else 60, same_back, p1, s2))
        if rng.random() < 0.03:   # a secondary and a supplementary copy, with altered sequence: skipped
            recs.append(record(n1, (f1 | 256) & ~1024, c1[0], p1 + 5000, 0, same, p2, s1[::-1]))
            recs.append(record(n2, (f2 | 2048) & ~1024, c2[0], p2 + 7000, 0, same_back, p1, s2[::-1]))
    for i in range(50):   # QC-failed reads: skipped, they must not appear
        c = rng.choice(contigs)
        p = rng.randrange(1, c[1] - 1000)
        recs.append(record(f"qcfail{i}", 512 | 4 | 8 if False else 512, c[0], p, 60, "*", 0, "ACGT" * 40))
    header = "".join(f"@SQ\tSN:{c}\tLN:{n}\n" for c, n in contigs) + "@RG\tID:x\tSM:synthetic\n"
    order = {c: i for i, (c, _) in enumerate(contigs)}

    def line(r):
        return "\t".join(str(x) for x in r) + "\n"

    with open(out_name, "w") as f:
        f.write("@HD\tVN:1.6\tSO:queryname\n" + header)
        for r in sorted(recs, key=lambda r: (r[0], 0 if r[1] & 64 else 1, r[1] & 2304)):
            f.write(line(r))
    with open(out_coord, "w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n" + header)
        for r in sorted(recs, key=lambda r: (order.get(r[2], 99), r[3] if r[2] != "*" else 0, r[0], r[1])):
            f.write(line(r))


if __name__ == "__main__":
    main()
