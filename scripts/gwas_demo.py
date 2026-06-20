#!/usr/bin/env python3
"""Pangenome association harness: PRESENCE/ABSENCE vs COUNT (multiplicity) on one feature substrate.

Reads a panvar describe sample-level fsm file (`<feature> | sample:count ...`) and a phenotype
table, and tests each discriminative feature two ways:

  presence/absence : binarize count>0  (what a presence/absence-based association method sees)
  count            : use the per-sample feature count = local copy number (a count-based method)

The substrate is chosen with --substrate:

  kmer  : k-mer multiplicity (describe fsm_kmers.samples.txt.gz) -- the primary, tested layer
  graph : node/edge dosage   (describe fsm_graph.samples.txt.gz) -- the complementary layer

For a copy-number locus (e.g. a tandem repeat) the unit feature is present in every sample, so the
presence/absence test has no contrast and cannot find it, while the count test recovers it. Multiple
testing is Benjamini-Hochberg (q) + Bonferroni. Features are annotated with their graph node(s), a
genomic position (call node_track.tsv) and the called variant (call variant_nodes.tsv) for a
Manhattan plot and traceback.

Outputs <out>.assoc.tsv and prints a summary. numpy + scipy only.
"""
import argparse
import glob
import gzip
import os
import sys

import numpy as np
from scipy import stats


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_fsm(path):
    """fsm-lite: '<feature> | strain:count strain:count ...' -> {feature: {strain: count}}."""
    out = {}
    with open_maybe_gz(path) as fh:
        for line in fh:
            if "|" not in line:
                continue
            left, right = line.split("|", 1)
            feature = left.split()[0]
            carriers = {}
            for tok in right.split():
                s, _, c = tok.rpartition(":")
                if s:
                    carriers[s] = float(c)
            out[feature] = carriers
    return out


def feature_nodes(label, substrate, kn):
    """Graph node id(s) backing a feature, for position/variant annotation.

    kmer  : look up the decoded k-mer in the describe feature map (kmer -> nodes).
    graph : the feature label IS a graph coordinate -- a node id, or an 'a+>b+' edge key whose two
            endpoints are the nodes (orientation stripped)."""
    if substrate == "graph":
        if ">" in label:
            return [p.rstrip("+-") for p in label.split(">") if p]
        return [label.rstrip("+-")]
    return sorted(kn.get(label, set()), key=lambda x: (len(x), x))


def parse_phenotypes(path, id_col, col):
    strains, y = [], []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if id_col not in header or col not in header:
            sys.exit(f"phenotype file needs columns '{id_col}' and '{col}' (have: {header})")
        pn_i, y_i = header.index(id_col), header.index(col)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(pn_i, y_i) or f[y_i] == "" or f[y_i] == "NA":
                continue
            strains.append(f[pn_i])
            y.append(float(f[y_i]))
    return strains, np.asarray(y, float)


def load_kmer_nodes(feature_glob):
    """kmer (decoded seq) -> set(node_ids) from describe kmer_features.tsv.gz files."""
    kn = {}
    for fp in glob.glob(feature_glob):
        with open_maybe_gz(fp) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            if "kmer" not in hdr or "nodes" not in hdr:
                continue
            ki, ni = hdr.index("kmer"), hdr.index("nodes")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= max(ki, ni):
                    continue
                # describe writes node ids as a ';'-separated list
                kn.setdefault(f[ki], set()).update(n for n in f[ni].replace(",", ";").split(";") if n)
    return kn


def load_node_pos(node_track):
    pos, bub = {}, {}
    if not node_track or not os.path.exists(node_track):
        return pos, bub
    with open(node_track) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            nid = f[idx["node_id"]]
            gp = f[idx["genomic_pos"]]
            if gp not in ("", "NA", "."):
                pos[nid] = float(gp)
            bub[nid] = f[idx["bubble_id"]]
    return pos, bub


def load_variant_nodes(vn_path):
    """node_id -> 'variant_id(svtype)' from call variant_nodes.tsv."""
    nv = {}
    if not vn_path or not os.path.exists(vn_path):
        return nv
    with open(vn_path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            vid, svt, nodes = f[0], f[2], f[3]
            for n in nodes.split(","):
                if n:
                    nv.setdefault(n, f"{vid}({svt})")
    return nv


def load_node_genes(ng_path):
    """node_id -> gene-name string (';'-joined) from call node_genes.tsv (from --gtf)."""
    ng = {}
    if not ng_path or not os.path.exists(ng_path):
        return ng
    with open(ng_path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0]:
                ng[f[0]] = f[1]
    return ng


def bh_q(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fsm", required=True)
    ap.add_argument("--phenotypes", required=True)
    ap.add_argument("--phenotype-col", required=True)
    ap.add_argument("--id-col", default="sample",
                    help="phenotype id column matched to the fsm strain ids (default: sample)")
    ap.add_argument("--mode", choices=["continuous", "binary"], required=True)
    ap.add_argument("--substrate", choices=["kmer", "graph"], default="kmer",
                    help="feature type in --fsm: kmer (fsm_kmers.samples) or graph node/edge "
                         "(fsm_graph.samples) (default: kmer)")
    ap.add_argument("--feature-map", help="glob for describe kmer_features.tsv.gz (traceback/position)")
    ap.add_argument("--node-track", help="call node_track.tsv (genomic position)")
    ap.add_argument("--variant-nodes", help="call variant_nodes.tsv (variant annotation)")
    ap.add_argument("--node-genes", help="call node_genes.tsv (gene traceback, from call --gtf)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fsm = parse_fsm(args.fsm)
    strains, y = parse_phenotypes(args.phenotypes, args.id_col, args.phenotype_col)
    sidx = {s: i for i, s in enumerate(strains)}
    n = len(strains)
    if n == 0:
        sys.exit("no phenotyped samples")

    kn = load_kmer_nodes(args.feature_map) if args.feature_map else {}
    npos, nbub = load_node_pos(args.node_track)
    nv = load_variant_nodes(args.variant_nodes)
    ng = load_node_genes(args.node_genes)
    # bubble -> median reference position, so off-reference k-mers still land at their locus
    bub_pos = {}
    for nd, p in npos.items():
        bub_pos.setdefault(nbub.get(nd, "."), []).append(p)
    bub_pos = {b: float(np.median(v)) for b, v in bub_pos.items()}

    recs = []
    gene_anns = []   # parallel to recs (pre-sort), appended as the trailing 'genes' column
    for kmer, carriers in fsm.items():
        c = np.zeros(n)
        for s, cnt in carriers.items():
            if s in sidx:
                c[sidx[s]] = cnt
        if np.all(c == c[0]):  # invariant across the cohort
            continue
        b = (c > 0).astype(float)

        # presence/absence
        if np.all(b == b[0]):
            pa_p = 1.0  # present (or absent) in everyone -> no contrast, untestable
        elif args.mode == "continuous":
            pa_p = stats.linregress(b, y).pvalue
        else:
            tab = [[int(((b == 1) & (y == 1)).sum()), int(((b == 1) & (y == 0)).sum())],
                   [int(((b == 0) & (y == 1)).sum()), int(((b == 0) & (y == 0)).sum())]]
            pa_p = stats.fisher_exact(tab)[1]

        # count (multiplicity)
        if args.mode == "continuous":
            count_p = stats.linregress(c, y).pvalue
        else:
            g1, g0 = c[y == 1], c[y == 0]
            count_p = stats.mannwhitneyu(g1, g0, alternative="two-sided").pvalue if len(g1) and len(g0) else 1.0

        nodes = feature_nodes(kmer, args.substrate, kn)
        ps = [npos[nd] for nd in nodes if nd in npos]
        bub = next((nbub[nd] for nd in nodes if nd in nbub), ".")
        pos = float(np.median(ps)) if ps else bub_pos.get(bub, float("nan"))
        var = next((nv[nd] for nd in nodes if nd in nv), ".")
        # graph-order x: smallest numeric node id (odgi-sorted graph -> ~reference order),
        # placing EVERY k-mer including those with no reference coordinate.
        num_nodes = [int(nd) for nd in nodes if nd.isdigit()]
        node_min = float(min(num_nodes)) if num_nodes else float("nan")
        gset = []
        for nd in nodes:
            for g in ng.get(nd, "").split(";"):
                if g and g not in gset:
                    gset.append(g)
        gene_anns.append(";".join(sorted(gset)) or ".")
        recs.append([kmer, bub, pos, node_min, ",".join(nodes) or ".", var,
                     int(b.sum()), float(c.max()), pa_p, count_p])

    if not recs:
        sys.exit("no testable (variable) k-mers")

    # record layout: 0 kmer 1 bubble 2 pos 3 node_min 4 nodes 5 variant
    #                6 n_carriers 7 max_count 8 pa_p 9 count_p
    pa_q = bh_q([r[8] for r in recs])
    count_q = bh_q([r[9] for r in recs])
    m = len(recs)
    for i, r in enumerate(recs):
        r += [min(1.0, r[8] * m), min(1.0, r[9] * m), pa_q[i], count_q[i], gene_anns[i]]
    #                10 pa_bonf 11 count_bonf 12 pa_q 13 count_q 14 genes

    recs.sort(key=lambda r: r[13])  # by count_q
    with open(args.out + ".assoc.tsv", "w") as fh:
        fh.write("kmer\tbubble_id\tpos\tnode_min\tnodes\tvariant\tn_carriers\tmax_count\t"
                 "pa_p\tcount_p\tpa_bonf\tcount_bonf\tpa_q\tcount_q\tgenes\n")
        for r in recs:
            fh.write("{}\t{}\t{:.1f}\t{:.0f}\t{}\t{}\t{}\t{:.0f}\t{:.3e}\t{:.3e}\t{:.3e}\t{:.3e}\t{:.3e}\t{:.3e}\t{}\n".format(*r))

    sig = lambda q: q < 0.05
    n_count = sum(sig(r[13]) for r in recs)
    n_pa = sum(sig(r[12]) for r in recs)
    print(f"[{args.mode}/{args.substrate}] {m} testable features | count q<0.05: {n_count} | P/A q<0.05: {n_pa}")
    print("  top 5 by count:")
    for r in recs[:5]:
        print(f"    {r[0][:24]}.. nodes={r[4]} var={r[5]} count_q={r[13]:.2e} pa_q={r[12]:.2e} maxCN={r[7]:.0f} carriers={r[6]}/{n}")
    print(f"  wrote {args.out}.assoc.tsv")


if __name__ == "__main__":
    main()
