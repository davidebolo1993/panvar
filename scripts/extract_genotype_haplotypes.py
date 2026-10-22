#!/usr/bin/env python3
"""Extract a panel haplotype from the graph, independently of the caller, from the extraction reference that
sequence_first_v2 writes in <prefix>.sequence_haplotypes.tsv (column `extraction`):

    gfa_path=<path name>;interval=0-<length>;orientation=<forward|reverse_complement>;whole_locus_path

The path's steps are concatenated from the GFA's S-lines (each step reverse-complemented when its orientation is '-'),
the interval is taken (the whole locus path), and the result is reverse-complemented when the reference says so. The
SHA-256 and length must equal the row's `sha256` and `length_bp`.

  extract_genotype_haplotypes.py <graph.gfa[.gz]> <extraction reference> [--fasta out.fa]
      -> prints "<sha256>\t<length>" for the extracted sequence
  extract_genotype_haplotypes.py --check <graph.gfa[.gz]> <prefix>.sequence_haplotypes.tsv
      -> extracts EVERY panel path of EVERY PANEL row and compares; exit 0 iff all match
"""
import gzip
import hashlib
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def load(gfa):
    opener = gzip.open if gfa.endswith(".gz") else open
    segments, paths = {}, {}
    with opener(gfa, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] == "S":
                segments[p[1]] = p[2]
            elif p[0] == "P":
                paths[p[1]] = p[2]
    return segments, paths


def path_sequence(segments, paths, name):
    out = []
    for step in paths[name].split(","):
        node, orient = step[:-1], step[-1]
        out.append(segments[node] if orient == "+" else revcomp(segments[node]))
    return "".join(out)


def extract(segments, paths, reference):
    fields = dict(item.split("=", 1) for item in reference.split(";") if "=" in item)
    seq = path_sequence(segments, paths, fields["gfa_path"])
    start, end = (int(x) for x in fields["interval"].split("-"))
    seq = seq[start:end]
    return revcomp(seq) if fields["orientation"] == "reverse_complement" else seq


def main(argv):
    if argv[0] == "--check":
        segments, paths = load(argv[1])
        rows = [l.rstrip("\n").split("\t") for l in open(argv[2])]
        head, bad, checked = rows[0], 0, 0
        for row in rows[1:]:
            r = dict(zip(head, row))
            if r["origin"] != "PANEL":
                continue
            interval = r["extraction"].split("interval=")[1].split(";")[0]
            for name, orient in zip(r["panel_paths"].split(","), r["path_orientations"].split(",")):
                seq = extract(segments, paths, f"gfa_path={name};interval={interval};orientation={orient}")
                ok = hashlib.sha256(seq.encode()).hexdigest() == r["sha256"] and str(len(seq)) == r["length_bp"]
                checked += 1
                bad += not ok
                print(f"{r['haplotype']}\t{name}\t{orient}\t{'MATCH' if ok else 'MISMATCH'}")
        print(f"{checked} panel path extractions, {bad} mismatches")
        return 1 if bad else 0
    segments, paths = load(argv[0])
    seq = extract(segments, paths, argv[1])
    if "--fasta" in argv:
        open(argv[argv.index("--fasta") + 1], "w").write(">extracted\n" + seq + "\n")
    print(f"{hashlib.sha256(seq.encode()).hexdigest()}\t{len(seq)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
