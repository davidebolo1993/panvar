"""Helpers shared by the transcripts checks: GFA walks, reverse complement, edit distances."""
import gzip
import os
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def walks(gfa, names=None):
    """Each path's sequence, as written (only `names` when given)."""
    seg, paths = {}, {}
    with (gzip.open if gfa.endswith(".gz") else open)(gfa, "rt") as fh:
        for line in fh:
            if line[0] == "S":
                f = line.rstrip("\n").split("\t")
                seg[f[1]] = f[2].upper()
            elif line[0] == "P":
                f = line.split("\t", 3)
                if names is None or f[1] in names:
                    paths[f[1]] = f[2].split(",")
    return {n: "".join(seg[s[:-1]] if s[-1] == "+" else rc(seg[s[:-1]]) for s in p) for n, p in paths.items()}


def edit_distances(pairs):
    """{key: distance} for (key, a, b) pairs, with dev/edit_distance.cpp (built into build/ when missing)."""
    if not pairs:
        return {}
    exe = os.path.join(REPO, "build", "edit_distance")
    if not os.path.exists(exe):
        subprocess.run(["c++", "-std=c++17", "-O2", "-I" + os.path.join(REPO, "external/edlib/edlib/include"),
                        os.path.join(REPO, "dev/edit_distance.cpp"), os.path.join(REPO, "build/libedlib.a"), "-o", exe],
                       check=True)
    out = subprocess.run([exe], input="".join(f"{k}\t{a}\t{b}\n" for k, a, b in pairs), capture_output=True,
                         text=True, check=True).stdout
    return {k: int(d) for k, d in (line.split("\t") for line in out.splitlines())}
