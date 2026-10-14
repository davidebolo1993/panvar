#!/usr/bin/env python3
"""Re-express staged MOSAIC truth paths in the ORIGINAL graph's node frame (truth construction glue; no science).

mosaic-truth writes its staged graph in the reference-SORTED node frame of the graph it was given, which renumbers
nodes and may store a node reverse-complemented. Calls use the original graph, and the evaluator places a truth graph
in the model's frame only when both share node identities, so the mosaic paths must be written in original node ids.

Method: every template path present in both graphs gives a step-by-step correspondence (sorted node, strand) ->
(original node, strand). The mapping must be one-to-one and consistent across all paths, and each mapped node's
sequence must equal the original's (or its reverse complement when the strand flips). The MOSAIC walks are then
translated step by step. Refuses unless: every MOSAIC node is mapped; every junction of the translated walk is an
original-graph edge in the walked orientation; and each translated path spells exactly the staged MOSAIC spelling.

Usage: translate_mosaic_frame.py <original.gfa[.gz]> <staged.gfa> <out truth.gfa> <out report.tsv>
Writes the original graph's lines unchanged plus the two MOSAIC P lines (original ids).
"""
import gzip
import sys

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


def read(path):
    op = gzip.open if path.endswith(".gz") else open
    seq, paths, edges, lines = {}, {}, set(), []
    with op(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            lines.append(line)
            p = line.split("\t")
            if p[0] == "S":
                seq[p[1]] = p[2]
            elif p[0] == "P":
                paths[p[1]] = [(s[:-1], s[-1] == "-") for s in p[2].split(",")]
            elif p[0] == "L":
                a, ra, b, rb = p[1], p[2] == "-", p[3], p[4] == "-"
                edges.add((a, ra, b, rb))
                edges.add((b, not rb, a, not ra))   # the same edge read backwards
    return seq, paths, edges, lines


def spell(seq, steps):
    return "".join(rc(seq[n]) if r else seq[n] for n, r in steps)


def main():
    orig_path, staged_path, out_gfa, out_report = sys.argv[1:5]
    oseq, opaths, oedges, olines = read(orig_path)
    sseq, spaths, _, _ = read(staged_path)
    mapping, conflicts, shared = {}, 0, 0
    for name, ssteps in spaths.items():
        if name.startswith("MOSAIC#") or name not in opaths:
            continue
        osteps = opaths[name]
        if len(osteps) != len(ssteps):
            # a sorted path may be stored in the opposite direction
            raise SystemExit(f"REFUSED: path {name} has {len(ssteps)} steps staged but {len(osteps)} original")
        if spell(sseq, ssteps) == spell(oseq, osteps):
            pairs = zip(ssteps, osteps)
        elif spell(sseq, ssteps) == rc(spell(oseq, osteps)):
            pairs = zip(ssteps, [(n, not r) for n, r in reversed(osteps)])
        else:
            raise SystemExit(f"REFUSED: path {name} spells differently in the two graphs")
        shared += 1
        for (sn, sr), (on, orr) in pairs:
            flip = sr != orr
            if sn in mapping and mapping[sn] != (on, flip):
                conflicts += 1
            mapping[sn] = (on, flip)
    if conflicts:
        raise SystemExit(f"REFUSED: {conflicts} inconsistent node correspondences")
    if len({v[0] for v in mapping.values()}) != len(mapping):
        raise SystemExit("REFUSED: the node correspondence is not one-to-one")
    for sn, (on, flip) in mapping.items():
        if (rc(oseq[on]) if flip else oseq[on]) != sseq[sn]:
            raise SystemExit(f"REFUSED: staged node {sn} does not carry original node {on}'s sequence")
    out_lines, report = list(olines), []
    for name in ("MOSAIC#1#synthetic:0-0", "MOSAIC#2#synthetic:0-0"):
        ssteps = spaths[name]
        missing = [sn for sn, _ in ssteps if sn not in mapping]
        if missing:
            raise SystemExit(f"REFUSED: {name} uses {len(missing)} node(s) no shared path maps, e.g. {missing[0]}")
        osteps = [(mapping[sn][0], sr != mapping[sn][1]) for sn, sr in ssteps]
        bad_edges = sum((a, ra, b, rb) not in oedges for (a, ra), (b, rb) in zip(osteps, osteps[1:]))
        if bad_edges:
            raise SystemExit(f"REFUSED: {name} translated walk leaves the original graph at {bad_edges} junction(s)")
        staged_seq, orig_seq = spell(sseq, ssteps), spell(oseq, osteps)
        if staged_seq != orig_seq:
            raise SystemExit(f"REFUSED: {name} translated walk does not spell the staged sequence")
        out_lines.append(f"P\t{name}\t{','.join(n + ('-' if r else '+') for n, r in osteps)}\t*")
        report.append((name, len(osteps), sum(r for _, r in osteps), len(orig_seq), "identical"))
    with open(out_gfa, "w") as f:
        f.write("\n".join(out_lines) + "\n")
    with open(out_report, "w") as f:
        f.write(f"shared_paths_used\t{shared}\nmapped_nodes\t{len(mapping)}\noriginal_nodes\t{len(oseq)}\n"
                f"flipped_nodes\t{sum(v[1] for v in mapping.values())}\n")
        f.write("path\tsteps\tminus_steps\tbp\tspelling_vs_staged\n")
        for r in report:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"translated 2 MOSAIC paths through {len(mapping)} mapped nodes ({shared} shared paths); spellings identical")


if __name__ == "__main__":
    main()
