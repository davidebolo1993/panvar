#!/usr/bin/env python3
"""mosaic-truth on minus-strand paths (GFA step orientation in the two-route self-check).

A synthetic locus whose every path, the reference included, traverses its nodes on the '-' strand (as the CYP2D6
graph does), with non-palindromic node sequences. Checks:
  1. a switched mosaic is accepted by mosaic-truth's two-route self-check, and the emitted MOSAIC#1 path spells
     (reverse-complementing each '-' step) exactly the sequence expected independently here: the reverse complements of
     the intended nodes, in walk order; its steps are donor A#1's own up to an interior cut and B#1's after it, with
     '-' steps among them; the unswitched homologue spells donor A's second template exactly;
  2. an unswitched control spells both donor rows exactly;
  3. the comparison is sensitive to orientation: the same graph with one node read on the '+' strand by every path
     (its sequence unchanged) is still a valid locus, and its emitted mosaic no longer matches the expected sequence;
  4. two runs are byte-identical.
Usage: genotype_mosaic_truth_orientation_test.py <genotype_probe> [<pre-fix genotype_probe>]
With a pre-fix binary, also requires that it REFUSES case 1 (the self-check concatenated '-' steps unreversed).
"""
import os
import random
import shutil
import subprocess
import sys
import tempfile

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def node_sequences():
    """Ten non-palindromic nodes: contexts 1, 4, 7, 10 and bubble alleles {2, 3}, {5, 6}, {8, 9}."""
    rng = random.Random(20261014)
    lengths = {1: 80, 2: 70, 3: 12, 4: 80, 5: 64, 6: 9, 7: 80, 8: 66, 9: 11, 10: 80}
    seqs = {}
    for n, length in lengths.items():
        while True:
            s = "".join(rng.choice("ACGT") for _ in range(length))
            if s != rc(s):
                break
        seqs[n] = s
    return seqs


# panel rows: allele at bubbles (2|3, 5|6, 8|9); every path is written on the minus strand, from node 10 down to 1
ROWS = {"ref#0#chain": (2, 5, 8), "A#1#h1": (2, 6, 8), "A#2#h2": (3, 5, 9), "B#1#h1": (3, 6, 9), "B#2#h2": (2, 5, 9)}


def write_gfa(path, seqs, plus_node=None):
    """All paths traverse 10 -> 1 on the '-' strand; plus_node (if any) is read on the '+' strand by every path."""
    def step(n):
        return f"{n}{'+' if n == plus_node else '-'}"
    lines = ["H\tVN:Z:1.0"] + [f"S\t{n}\t{s}" for n, s in sorted(seqs.items())]
    edges = set()
    for alleles in ROWS.values():
        walk = [10, alleles[2], 7, alleles[1], 4, alleles[0], 1]
        for a, b in zip(walk, walk[1:]):
            edges.add((a, b))
    for a, b in sorted(edges):
        oa, ob = ("+" if a == plus_node else "-"), ("+" if b == plus_node else "-")
        lines.append(f"L\t{a}\t{oa}\t{b}\t{ob}\t0M")
    for name, alleles in ROWS.items():
        walk = [10, alleles[2], 7, alleles[1], 4, alleles[0], 1]
        lines.append(f"P\t{name}\t{','.join(step(n) for n in walk)}\t*")
    open(path, "w").write("\n".join(lines) + "\n")


def spell(gfa, name):
    seq, steps = {}, None
    for line in open(gfa):
        f = line.rstrip("\n").split("\t")
        if f[0] == "S":
            seq[f[1]] = f[2]
        elif f[0] == "P" and f[1] == name:
            steps = f[2].split(",")
    assert steps is not None, f"no path {name} in {gfa}"
    return "".join(rc(seq[s[:-1]]) if s.endswith("-") else seq[s[:-1]] for s in steps), steps


def run(probe, gfa, out, switches_hap1, require):
    os.makedirs(out, exist_ok=True)
    cmd = [probe, "mosaic-truth", "--gfa", gfa, "--reference", "ref#0#chain", "--donor-a", "A", "--donor-b", "B",
           "--switches-hap1", switches_hap1, "--switches-hap2", "", "--require-off-panel", require, "--min-variant-bp", "5",
           "--out-gfa", os.path.join(out, "staged.gfa"), "--out", os.path.join(out, "truth.tsv")]
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr).strip()


failures = []


def check(ok, what):
    print(("PASS " if ok else "FAIL ") + what)
    if not ok:
        failures.append(what)


def main():
    probe = sys.argv[1]
    prefix_probe = sys.argv[2] if len(sys.argv) > 2 else None
    tmp = tempfile.mkdtemp(prefix="mosaic_truth_orientation_")
    try:
        seqs = node_sequences()
        gfa = os.path.join(tmp, "minus.gfa")
        write_gfa(gfa, seqs)
        # the chain follows the reference, node 10 to node 1, as five blocks: ctx:0 (node 10), sv:1 (8|9), sv:2 (5|6),
        # sv:3 (2|3), ctx:4 (node 1), the shared contexts 7 and 4 absorbed into the bubble blocks. Switching hap1
        # before block 3 takes sv:3 and ctx:4 from B#1: alleles (3, 6, 8), in no panel row.
        rcode, log = run(probe, gfa, os.path.join(tmp, "switched"), "3", "one")
        check(rcode == 0, f"switched mosaic accepted by the two-route self-check ({log.splitlines()[-1] if log else ''})")
        expected1 = "".join(rc(seqs[n]) for n in (10, 8, 7, 6, 4, 3, 1))
        expected2 = "".join(rc(seqs[n]) for n in (10, 9, 7, 5, 4, 3, 1))
        if rcode == 0:
            staged = os.path.join(tmp, "switched", "staged.gfa")
            got1, steps1 = spell(staged, "MOSAIC#1#synthetic:0-0")
            got2, _ = spell(staged, "MOSAIC#2#synthetic:0-0")
            check(got1 == expected1, "MOSAIC#1 spells exactly the independently expected sequence (reverse complements, walk order)")
            # the staged graph is the reference-sorted one (the reference reads '+' there, alternative alleles keep
            # '-'): MOSAIC#1 must be donor A#1's own steps up to an interior cut and B#1's after it, orientation kept
            a1, b1 = spell(staged, "A#1#h1")[1], spell(staged, "B#1#h1")[1]
            cuts = [k for k in range(1, len(steps1)) if steps1[:k] == a1[:k] and steps1[k:] == b1[k:] and a1 != steps1 != b1]
            check(bool(cuts), f"MOSAIC#1 joins A#1's own steps to B#1's at an interior cut, orientations kept (cuts {cuts})")
            check(any(s.endswith("-") for s in steps1), "the mosaic walk contains '-' steps, so reverse complementing is exercised")
            check(got2 == expected2, "MOSAIC#2 (unswitched) spells donor A's second row exactly")
            check(got1 != "".join(seqs[n] for n in (10, 8, 7, 6, 4, 3, 1)), "the expected sequence is not the unreversed concatenation")
        # unswitched control: both homologues are donor A's rows
        rcode, log = run(probe, gfa, os.path.join(tmp, "control"), "", "none")
        check(rcode == 0, "unswitched control accepted")
        if rcode == 0:
            staged = os.path.join(tmp, "control", "staged.gfa")
            check(spell(staged, "MOSAIC#1#synthetic:0-0")[0] == spell(staged, "A#1#h1")[0], "control MOSAIC#1 == A#1 spelling")
            check(spell(staged, "MOSAIC#2#synthetic:0-0")[0] == spell(staged, "A#2#h2")[0], "control MOSAIC#2 == A#2 spelling")
        # negative control: node 6 read on the '+' strand by every path, sequences unchanged
        flipped = os.path.join(tmp, "flipped.gfa")
        write_gfa(flipped, seqs, plus_node=6)
        rcode, log = run(probe, flipped, os.path.join(tmp, "flipped"), "3", "one")
        check(rcode == 0, f"flipped-orientation locus is still a valid locus ({log.splitlines()[-1] if log else ''})")
        if rcode == 0:
            got1, _ = spell(os.path.join(tmp, "flipped", "staged.gfa"), "MOSAIC#1#synthetic:0-0")
            check(got1 != expected1, "flipping one step's orientation (sequence unchanged) makes the expected-sequence comparison fail")
            check(got1 == "".join(rc(seqs[n]) if n != 6 else seqs[n] for n in (10, 8, 7, 6, 4, 3, 1)),
                  "and the flipped mosaic spells node 6 forward, as its path now reads it")
        # determinism
        r1 = run(probe, gfa, os.path.join(tmp, "again"), "3", "one")[0]
        same = r1 == 0 and all(open(os.path.join(tmp, "again", f), "rb").read() == open(os.path.join(tmp, "switched", f), "rb").read()
                               for f in ("staged.gfa", "truth.tsv"))
        check(same, "two runs are byte-identical (staged.gfa, truth.tsv)")
        if prefix_probe:
            rcode, log = run(prefix_probe, gfa, os.path.join(tmp, "prefix"), "3", "one")
            check(rcode != 0 and "the splice does not reproduce the truth" in log,
                  "the pre-fix binary refuses the same minus-strand mosaic at its self-check")
    finally:
        shutil.rmtree(tmp)
    if failures:
        print(f"{len(failures)} failure(s)")
        return 1
    print("mosaic-truth orientation: all checks pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
