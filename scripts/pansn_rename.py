#!/usr/bin/env python3
"""Rewrite GFA path names to PanSN (sample#hap#contig[:start-end]) IN PLACE.

panvar's gene annotation (`--gtf`) only fires when the reference path name is PanSN, so it can read
the chromosome + absolute start. Some graphs name paths with underscores instead
(`GRCh38_0_chr6:160509252-160734894`, `HG00097_1_JBIRDD010000043.1:...`); this converts the first two
underscores of each P-line name to '#' (`GRCh38#0#chr6:160509252-160734894`, etc.).

Idempotent: names that already contain '#' are left unchanged, so re-running is safe. W-lines already
carry sample/hap/contig in separate columns (PanSN by construction) and are left untouched.

  pansn_rename.py <graph.gfa[.gz]> [<graph2.gfa[.gz]> ...]
"""
import gzip
import os
import shutil
import sys
import tempfile


def open_text(path, mode):
    return gzip.open(path, mode + "t") if path.endswith(".gz") else open(path, mode)


def to_pansn(name):
    """SAMPLE_HAP_CONTIG[:coords] -> SAMPLE#HAP#CONTIG[:coords] (first two '_' -> '#')."""
    if "#" in name:
        return name  # already PanSN
    parts = name.split("_", 2)
    if len(parts) < 3:
        return name  # not the expected sample_hap_contig shape; leave as-is
    return parts[0] + "#" + parts[1] + "#" + parts[2]


def convert(path):
    changed = 0
    fd, tmp = tempfile.mkstemp(suffix=".gfa.gz" if path.endswith(".gz") else ".gfa",
                               dir=os.path.dirname(os.path.abspath(path)))
    os.close(fd)
    with open_text(path, "r") as fin, open_text(tmp, "w") as fout:
        for line in fin:
            if line.startswith("P\t"):
                f = line.rstrip("\n").split("\t")
                new = to_pansn(f[1])
                if new != f[1]:
                    f[1] = new
                    changed += 1
                    fout.write("\t".join(f) + "\n")
                    continue
            fout.write(line)
    if changed:
        shutil.move(tmp, path)
    else:
        os.remove(tmp)  # leave the file byte-identical when nothing changed
    return changed


def main(argv):
    if len(argv) < 2 or argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    for path in argv[1:]:
        n = convert(path)
        print(f"{path}: renamed {n} P-line path name(s) to PanSN" if n
              else f"{path}: already PanSN (no change)")


if __name__ == "__main__":
    main(sys.argv)
