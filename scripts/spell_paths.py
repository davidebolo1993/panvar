#!/usr/bin/env python3
"""spell_paths.py - write each GFA path out as FASTA.

The graph's P lines already spell their haplotypes exactly at the `bubble` stage, so this is how you
get per-haplotype sequence for read simulation. Simulate from these, not from a post-`panphorte`
graph: approximate folding rewrites tandem arrays, so reads drawn from a folded graph would carry the
collapsed sequence the genotyper is trying to be tested against.

Usage:
  spell_paths.py -i <graph.gfa[.gz]> -o <out_dir> [--paths NAME[,NAME...]] [--list]
  spell_paths.py -i <graph.gfa> --list                 # names only, one per line
  spell_paths.py -i <graph.gfa> -o fa --paths HG002#1#chr1  # just these

Options:
  -i, --gfa <path>      Input GFA, plain or gzipped (required)
  -o, --out-dir <dir>   Directory for the per-path FASTA files
      --paths <names>   Comma-separated path names to emit (default: all)
      --list            Print path names and lengths, write nothing
      --width <N>       FASTA line width (default 60)
  -h, --help            Show this help
"""

import argparse
import gzip
import os
import sys

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load(gfa):
    segments, paths = {}, {}
    for line in _open(gfa):
        f = line.rstrip("\n").split("\t")
        if f[0] == "S":
            segments[f[1]] = f[2]
        elif f[0] == "P":
            paths[f[1]] = [(s[:-1], s[-1]) for s in f[2].split(",") if s]
    return segments, paths


def spell(segments, steps):
    out = []
    for node, orient in steps:
        seq = segments.get(node, "")
        out.append(seq if orient == "+" else seq.translate(_COMP)[::-1])
    return "".join(out)


def safe_name(name):
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)


def main():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("-i", "--gfa", required=True)
    ap.add_argument("-o", "--out-dir")
    ap.add_argument("--paths")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--width", type=int, default=60)
    ap.add_argument("-h", "--help", action="store_true")
    args = ap.parse_args()
    if args.help:
        print(__doc__)
        return 0

    segments, paths = load(args.gfa)
    wanted = list(paths)
    if args.paths:
        requested = [p for p in args.paths.split(",") if p]
        missing = [p for p in requested if p not in paths]
        if missing:
            sys.exit(f"spell_paths: not in graph: {', '.join(missing)}")
        wanted = requested

    if args.list:
        for name in wanted:
            print(f"{name}\t{len(spell(segments, paths[name]))}")
        return 0

    if not args.out_dir:
        sys.exit("spell_paths: -o/--out-dir is required unless --list")
    os.makedirs(args.out_dir, exist_ok=True)
    for name in wanted:
        seq = spell(segments, paths[name])
        out = os.path.join(args.out_dir, safe_name(name) + ".fa")
        with open(out, "w") as fh:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), args.width):
                fh.write(seq[i : i + args.width] + "\n")
    print(f"wrote {len(wanted)} FASTA files to {args.out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
