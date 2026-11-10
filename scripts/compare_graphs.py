#!/usr/bin/env python3
"""compare_graphs.py - compare the structure of graphs of one or several loci.

Usage:
  compare_graphs.py -i <first.gfa[.gz]> -i <other.gfa[.gz]> [-i ...] -o <prefix> [options]
  compare_graphs.py --list <graphs.tsv> -o <prefix> [options]

Options:
  -i, --gfa <path>       A graph, plain or gzipped; at least two. The others are compared with the first
  -l, --labels <a,b,..>  One label per -i graph (default: the file names)
      --list <tsv>       Several loci: locus, stage and path per line, tab-separated (optional header);
                         each locus needs two or more graphs, compared with its first
  -o, --out <prefix>     Output prefix (required): <prefix>.graphs.tsv and <prefix>.paths.tsv
      --hub-degree <N>   Distinct neighbours on one side that make a node a hub, as in rebuild (default 50)
  -h, --help             Show this help
"""

import argparse
import gzip
import hashlib
import os
import re
import statistics
import sys
from collections import defaultdict

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_WALK_STEP = re.compile(r"([<>])([^<>]+)")


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def walk_name(f):
    # W: sample, haplotype, sequence, 0-based start, end.
    name = f"{f[1]}#{f[2]}#{f[3]}"
    if f[4].isdigit():
        name += f":{int(f[4]) + 1}-{f[5] if f[5].isdigit() else '*'}"
    return name


_INTERVAL = re.compile(r":\d+-(\d+|\*)$")


def match_paths(rows, first):
    """Each path's counterpart in the first graph: same name, else the only one with the same name less
    its :start-end interval."""
    by_base = defaultdict(list)
    for name in first:
        by_base[_INTERVAL.sub("", name)].append(name)
    match = {}
    for name in rows:
        if name in first:
            match[name] = name
        else:
            cands = by_base.get(_INTERVAL.sub("", name), [])
            if len(cands) == 1:
                match[name] = cands[0]
    return match


def load(path):
    segments, links, paths = {}, set(), {}
    for line in _open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == "S":
            segments[f[1]] = f[2]
        elif f[0] == "L":
            fwd, rev = (f[1], f[2], f[3], f[4]), (f[3], flip(f[4]), f[1], flip(f[2]))
            links.add(min(fwd, rev))
        elif f[0] == "P":
            paths[f[1]] = [(s[:-1], s[-1]) for s in f[2].split(",") if s]
        elif f[0] == "W":
            paths[walk_name(f)] = [(n, "+" if o == ">" else "-") for o, n in _WALK_STEP.findall(f[6])]
    return segments, links, paths


def flip(o):
    return "-" if o == "+" else "+"


def degrees(links):
    # Distinct neighbours on each side of a segment; the right side is the end a '+' step leaves by.
    sides = defaultdict(set)
    loops = set()
    for a, ao, b, bo in links:
        if a == b:
            loops.add(a)
            continue
        sides[(a, ao == "+")].add(b)
        sides[(b, bo == "-")].add(a)
    deg = defaultdict(int)
    for (node, _), nbrs in sides.items():
        deg[node] = max(deg[node], len(nbrs))
    return deg, loops


def spell_digest(segments, steps, rc_cache):
    parts = []
    for node, orient in steps:
        seq = segments.get(node, "")
        if orient == "-":
            if node not in rc_cache:
                rc_cache[node] = seq.translate(_COMP)[::-1]
            seq = rc_cache[node]
        parts.append(seq)
    seq = "".join(parts).upper()
    return len(seq), hashlib.blake2b(seq.encode(), digest_size=16).digest()


def summarise(path, label, hub_degree):
    segments, links, paths = load(path)
    if not segments:
        sys.exit(f"compare_graphs: {path} has no nodes")
    deg, loops = degrees(links)
    lengths = [len(s) for s in segments.values()]
    rc_cache = {}
    rows = {}
    for name, steps in paths.items():
        distinct = len({n for n, _ in steps})
        length, digest = spell_digest(segments, steps, rc_cache)
        rows[name] = {"steps": len(steps), "distinct": distinct, "length": length, "digest": digest}
    graph = {
        "graph": label, "file": path, "nodes": len(segments), "edges": len(links),
        "sequence_bp": sum(lengths), "median_node_bp": statistics.median(lengths),
        "single_base_nodes": sum(1 for x in lengths if x == 1) / len(lengths),
        "self_loops": len(loops), "hubs": sum(1 for d in deg.values() if d >= hub_degree),
        "max_degree": max(deg.values(), default=0), "paths": len(paths),
        "revisiting_paths": sum(1 for r in rows.values() if r["steps"] > r["distinct"]),
        "median_steps_per_kb": statistics.median(
            [r["steps"] / (r["length"] / 1000) for r in rows.values() if r["length"] > 0] or [0]),
    }
    # As rebuild logs them: nodes per kb of the longest haplotype, and the graph's sequence over the median one.
    lens = sorted(r["length"] for r in rows.values()) or [0]
    graph["nodes_per_kb"] = len(segments) / (lens[-1] / 1000) if lens[-1] else 0.0
    graph["graph_over_median_haplotype"] = graph["sequence_bp"] / lens[len(lens) // 2] if lens[len(lens) // 2] else 0.0
    return graph, rows


def fmt(x):
    if x is None:
        return "NA"
    if isinstance(x, float):
        return f"{x:.4g}"
    return str(x)


def read_list(path):
    entries = []
    for n, line in enumerate(open(path), 1):
        f = line.rstrip("\n").split("\t")
        if not line.strip() or line.startswith("#") or (n == 1 and f[:3] == ["locus", "stage", "path"]):
            continue
        if len(f) < 3:
            sys.exit(f"compare_graphs: {path}:{n}: expected locus, stage and path separated by tabs")
        entries.append((f[0], f[1], f[2]))
    return entries


def main():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("-i", "--gfa", action="append", default=[])
    ap.add_argument("-l", "--labels")
    ap.add_argument("--list")
    ap.add_argument("-o", "--out")
    ap.add_argument("--hub-degree", type=int, default=50)
    ap.add_argument("-h", "--help", action="store_true")
    args = ap.parse_args()
    if args.help or not (args.gfa or args.list):
        print(__doc__)
        return 0
    if args.gfa and args.list:
        sys.exit("compare_graphs: give -i or --list, not both")
    if not args.out:
        sys.exit("compare_graphs: -o/--out is required")
    if args.list:
        entries = read_list(args.list)
    else:
        if len(args.gfa) < 2:
            sys.exit("compare_graphs: give at least two graphs (-i twice)")
        labels = args.labels.split(",") if args.labels else [os.path.basename(p) for p in args.gfa]
        if len(labels) != len(args.gfa):
            sys.exit(f"compare_graphs: {len(labels)} labels for {len(args.gfa)} graphs")
        entries = [("locus", lab, p) for lab, p in zip(labels, args.gfa)]
    seen = set()
    for locus, stage, _ in entries:
        if (locus, stage) in seen:
            sys.exit(f"compare_graphs: {locus} {stage} is listed twice (give distinct labels with -l)")
        seen.add((locus, stage))
    loci = list(dict.fromkeys(e[0] for e in entries))
    for locus in loci:
        if sum(1 for e in entries if e[0] == locus) < 2:
            sys.exit(f"compare_graphs: locus {locus} has one graph; each needs two or more")

    graphs, per_path = [], []
    for locus, stage, path in entries:
        print(f"compare_graphs: reading {path}", file=sys.stderr)
        g, rows = summarise(path, stage, args.hub_degree)
        g["locus"] = locus
        graphs.append(g)
        per_path.append(rows)

    # Each graph against the first one listed for its locus, path by path.
    first_of = {}
    for g, rows in zip(graphs, per_path):
        first_of.setdefault(g["locus"], rows)
    matches = []
    for g, rows in zip(graphs, per_path):
        first = first_of[g["locus"]]
        match = match_paths(rows, first)
        matches.append(match)
        g["paths_in_first"] = len(match)
        g["paths_same_sequence"] = sum(1 for n, m in match.items() if rows[n]["digest"] == first[m]["digest"])
        renamed = sum(1 for n, m in match.items() if n != m)
        unmatched = len(rows) - len(match)
        if renamed or unmatched:
            print(f"compare_graphs: {g['locus']} {g['graph']}: {renamed} paths matched to the first graph by "
                  f"name less their interval, {unmatched} with no counterpart there", file=sys.stderr)

    gcols = ["locus", "graph", "file", "nodes", "edges", "sequence_bp", "median_node_bp", "single_base_nodes",
             "self_loops", "hubs", "max_degree", "nodes_per_kb", "graph_over_median_haplotype", "paths",
             "revisiting_paths", "median_steps_per_kb", "paths_in_first", "paths_same_sequence"]
    with open(args.out + ".graphs.tsv", "w") as out:
        out.write("\t".join(gcols) + "\n")
        for g in graphs:
            out.write("\t".join(fmt(g[c]) for c in gcols) + "\n")
    with open(args.out + ".paths.tsv", "w") as out:
        out.write("locus\tgraph\tpath\tsteps\tdistinct_nodes\trevisits\tlength_bp\tlength_diff_bp\tsame_sequence\n")
        for g, rows, match in zip(graphs, per_path, matches):
            first = first_of[g["locus"]]
            for name, r in rows.items():
                f = first[match[name]] if name in match else None
                diff = None if f is None else r["length"] - f["length"]
                same = None if f is None else int(r["digest"] == f["digest"])
                out.write("\t".join(fmt(x) for x in (g["locus"], g["graph"], name, r["steps"], r["distinct"],
                          r["steps"] - r["distinct"], r["length"], diff, same)) + "\n")

    show = ["nodes", "edges", "sequence_bp", "hubs", "max_degree", "self_loops", "revisiting_paths",
            "paths_same_sequence"]
    width = max(len(c) for c in show)
    heads = [g["graph"] if len(loci) == 1 else f"{g['locus']}:{g['graph']}" for g in graphs]
    print(f"{'':{width}}  " + "  ".join(f"{h:>14}" for h in heads))
    for c in show:
        print(f"{c:{width}}  " + "  ".join(f"{fmt(g[c]):>14}" for g in graphs))
    print(f"wrote {args.out}.graphs.tsv and {args.out}.paths.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
