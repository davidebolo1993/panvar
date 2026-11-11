#!/usr/bin/env python3
"""Truth for dev/recruit_assembly_truth.sh, from alignments to the sample's own assembly.

    recruit_assembly_score.py intervals <loci.fa> <loci.paf>      each locus on each haplotype (BED)
    recruit_assembly_score.py score <loci.bed> <aligned.tsv> <recruit.tsv>

intervals: per locus and haplotype (contigs prefixed h1_, h2_), from primary alignments only (a
secondary one may be a nearby paralog): the contig holding most of the locus's aligned bases, the span of
the alignments on it near the best one, and the share of the locus they cover (query bases counted once);
a locus covered under half is flagged (a deletion, or the assembly broken there).

score: aligned.tsv holds primary and supplementary alignments (name, mate, contig, start, end, MAPQ). Per locus, each
fragment aligned into it, or screened by recruit, is classed by its mates:
    locus_q        a mate (or part of one) overlaps the locus on either haplotype at MAPQ >= 20
    locus_low      a mate overlaps it, all such mates below MAPQ 20 (repeats the aligner cannot place)
    elsewhere_q    no mate there, every aligned mate MAPQ >= 20 (expect never certain)
    elsewhere_low  no mate there, an aligned mate below MAPQ 20
    unaligned      screened, no mate aligned
and counted by recruit's fate (unscreened: no row for the locus). False certain fragments are written to
false_certain.tsv beside the output for inspection.
"""
import collections
import csv
import sys

MAPQ = 20


def intervals(fa, paf):
    qlen = {}
    with open(fa) as f:
        name = None
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                qlen[name] = 0
            else:
                qlen[name] += len(line.strip())
    alns = collections.defaultdict(list)
    with open(paf) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            if "tp:A:S" in x[12:]:
                continue
            q, t, ts, te, matches = x[0], x[5], int(x[7]), int(x[8]), int(x[9])
            alns[(q, t[:3])].append((t, ts, te, matches, int(x[2]), int(x[3])))
    for q in qlen:
        for hap in ("h1_", "h2_"):
            a = alns.get((q, hap), [])
            if not a:
                print(f"# {q} {hap[:2]}: not found", file=sys.stderr)
                continue
            by = collections.Counter()
            for t, _, _, m, _, _ in a:
                by[t] += m
            contig = by.most_common(1)[0][0]
            best = max((x for x in a if x[0] == contig), key=lambda x: x[3])
            near = [x for x in a if x[0] == contig and x[1] < best[2] + qlen[q] and x[2] > best[1] - qlen[q]]
            s, e = min(x[1] for x in near), max(x[2] for x in near)
            covered, end = 0, 0   # query bases in the union of the alignments' query ranges
            for qs, qe in sorted((x[4], x[5]) for x in near):
                if qe > end:
                    covered += qe - max(qs, end)
                    end = qe
            cover = covered / qlen[q]
            flag = "" if cover >= 0.5 else "\tlow_cover"
            print(f"{contig}\t{s}\t{e}\t{q}\t{cover:.3f}{flag}")


def score(bed, aligned, table):
    loci = collections.defaultdict(list)
    with open(bed) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            loci[x[3]].append((x[0], int(x[1]), int(x[2])))
    mates = collections.defaultdict(list)   # name -> [(contig, start, end, mapq)]
    with open(aligned) as f:
        for line in f:
            n, _, c, s, e, q = line.rstrip("\n").split("\t")
            mates[n].append((c, int(s), int(e), int(q)))
    fates = collections.defaultdict(dict)   # locus -> name -> fate
    with open(table) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fates[r["locus"]][r["fragment"]] = r["fate"]

    def inside(m, ivs):
        return any(m[0] == c and m[1] < e and m[2] > s for c, s, e in ivs)

    fates_cols = ["certain", "ambiguous", "unrecruited", "unscreened"]
    print("\t".join(["locus", "class", "fragments"] + fates_cols + ["certain_share", "recruited_share"]))
    with open("false_certain.tsv", "w") as fc:
        fc.write("locus\tfragment\talignments\n")
        for locus in sorted(set(loci) | set(fates)):
            ivs = loci.get(locus, [])
            names = set(fates.get(locus, {}))
            names |= {n for n, ms in mates.items() if any(inside(m, ivs) for m in ms)}
            counts = collections.defaultdict(collections.Counter)
            for n in names:
                ms = mates.get(n, [])
                if any(inside(m, ivs) and m[3] >= MAPQ for m in ms):
                    cls = "locus_q"
                elif any(inside(m, ivs) for m in ms):
                    cls = "locus_low"
                elif not ms:
                    cls = "unaligned"
                elif all(m[3] >= MAPQ for m in ms):
                    cls = "elsewhere_q"
                else:
                    cls = "elsewhere_low"
                fate = fates.get(locus, {}).get(n, "unscreened")
                counts[cls][fate] += 1
                if fate == "certain" and not cls.startswith("locus"):
                    fc.write(f"{locus}\t{n}\t{cls}:" + ";".join(f"{c}:{s}-{e}:q{q}" for c, s, e, q in ms) + "\n")
            for cls in ("locus_q", "locus_low", "elsewhere_q", "elsewhere_low", "unaligned"):
                c = counts[cls]
                total = sum(c.values())
                if total == 0:
                    continue
                share = c["certain"] / total
                rec = (c["certain"] + c["ambiguous"]) / total
                print("\t".join([locus, cls, str(total)] + [str(c[f]) for f in fates_cols] + [f"{share:.4f}", f"{rec:.4f}"]))


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "intervals":
        intervals(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 5 and sys.argv[1] == "score":
        score(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
