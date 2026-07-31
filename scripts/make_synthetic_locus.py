#!/usr/bin/env python3
"""Generate a synthetic locus -- haplotypes AND the graph -- with every variant known by construction.

The point is a dataset where the expected answer is known exactly, so a failure is a bug rather than a
measurement. Real loci cannot do this: their truth is whatever the graph happens to encode, and a wrong
call is indistinguishable from an unrepresentable one.

What it emits:

  * a reference backbone, and haplotype DESIGNS differing from it at four well-separated SV sites
    (DEL, INS, DUP with variable copy number, INV) plus shared-position SNPs in between, which are what
    give the inter-bubble backbone blocks their identifying signal;
  * every design TWICE under different names -- an exact twin. That is what makes leave-one-out
    interpretable: hold out a sample's two haplotypes and its twins remain, so the panel can still
    represent the sample exactly. Expected accuracy is then 100%, and any shortfall is ours;
  * the GFA itself, with a W line per haplotype, built from the same edits. Emitting the graph rather
    than inferring it removes any doubt about what the true bubbles are -- there are exactly four SV
    bubbles plus one per SNP position.

The inversion is represented the way a real graph represents one: the same node traversed in reverse
orientation, which is the case that exercises orientation handling downstream.

Outputs (into --out):
  ref.fa, haplotypes.fa, hap_<name>.fa   sequences
  graph.gfa                              the graph, W line per haplotype, reference path named "ref"
  truth.variants.tsv                     per design: state at each SV site
  truth.haplotypes.tsv                   haplotype -> design, twin
  truth.bubbles.tsv                      the SV bubbles that `bubble` is expected to find
"""

import argparse
import os
import random

COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s):
    return s.translate(COMP)[::-1]


def random_seq(rng, n, gc=0.41):
    # Not uniform: pure-random sequence makes every k-mer unique and the problem unrealistically easy.
    return "".join(rng.choices("ACGT", weights=[(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2], k=n))


def write_fasta(path, name, seq, width=60):
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i + width] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--backbone-bp", type=int, default=50000)
    ap.add_argument("--designs", type=int, default=8)
    ap.add_argument("--snps", type=int, default=60,
                    help="SNP sites, shared positions with per-design alleles")
    ap.add_argument("--paralog-del", action="store_true",
                    help="plant a SECOND deletion site whose sequence (and flanking context) is an "
                         "exact copy of the first, in a different block, polymorphic independently. "
                         "Every syncmer of either site then occurs in both AND varies in both, so "
                         "confinement must strip them all and both bubbles are left with no local "
                         "evidence. That forces the call to come from linkage -- the C4A/C4B case, and "
                         "the one mechanism the clean synthetic locus never exercises.")
    ap.add_argument("--segdups", type=int, default=0,
                    help="segmental duplication PAIRS to plant. Each copies a stretch of backbone to "
                         "another position in a DIFFERENT block, so its syncmers occur in two blocks "
                         "at once. That is what marker confinement strips, and what makes a real "
                         "locus hard -- without it every marker is trivially block-unique.")
    ap.add_argument("--segdup-len", type=int, default=4000)
    ap.add_argument("--segdup-divergence", type=float, default=0.01,
                    help="per-base divergence between the two copies. At 1%% a 31-mer survives in "
                         "both with probability ~0.73, so most markers are shared.")
    ap.add_argument("--dup-divergence", type=float, default=0.0,
                    help="per-base divergence between tandem copies. 0 keeps copies identical, so "
                         "copy number shows up purely as marker multiplicity -- the clean case. Raise "
                         "it to make the array progressively harder.")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.out, exist_ok=True)
    backbone = random_seq(rng, args.backbone_bp)

    # SV sites at fixed FRACTIONS of the backbone, so the block structure scales with it. At the
    # default 50 kb these are the familiar 8000/18000/28000/38000.
    DEL_POS, DEL_LEN = int(args.backbone_bp * 0.16), 400
    INS_POS, INS_LEN = int(args.backbone_bp * 0.36), 600
    DUP_POS, DUP_LEN = int(args.backbone_bp * 0.56), 500
    INV_POS, INV_LEN = int(args.backbone_bp * 0.76), 700
    ins_payload = random_seq(rng, INS_LEN)
    # The paralogue sits inside the backbone stretch between INS and DUP, i.e. a different block.
    DEL2_POS = int(args.backbone_bp * 0.46)
    PARA_FLANK = 400
    sv_spans = [(DEL_POS, DEL_POS + DEL_LEN), (INS_POS, INS_POS + INS_LEN),
                (DUP_POS, DUP_POS + DUP_LEN), (INV_POS, INV_POS + INV_LEN)]
    if args.paralog_del:
        sv_spans.append((DEL2_POS, DEL2_POS + DEL_LEN))

    def in_sv(p):
        return any(a - 50 <= p < b + 50 for a, b in sv_spans)


    # Segmental duplications, planted before anything else so every haplotype inherits them. Each
    # pair puts near-identical sequence in two different blocks of the chain.
    segdup_pairs = []
    if args.segdups > 0:
        L = args.segdup_len
        # Deliberately spanning different blocks: flank<->flank, then backbone<->backbone.
        # One copy per block, chosen so the two copies land in DIFFERENT blocks of the chain -- that
        # is the whole point. Block midpoints, in order: leading flank, the three backbones, trailing
        # flank.
        mids = [DEL_POS // 2,
                (DEL_POS + DEL_LEN + INS_POS) // 2,
                (INS_POS + INS_LEN + DUP_POS) // 2,
                (DUP_POS + DUP_LEN + INV_POS) // 2,
                (INV_POS + INV_LEN + args.backbone_bp) // 2]
        starts = [max(200, m - L // 2) for m in mids]
        candidates = [(starts[0], starts[4]), (starts[1], starts[3]), (starts[2], starts[4] + L)]
        for k in range(min(args.segdups, len(candidates))):
            a, b = candidates[k]
            if any(a - 200 < e and s_ - 200 < a + L + 200 for s_, e in sv_spans) or \
               any(b - 200 < e and s_ - 200 < b + L + 200 for s_, e in sv_spans):
                continue
            src = backbone[a:a + L]
            cp = list(src)
            for i in range(L):
                if rng.random() < args.segdup_divergence:
                    cp[i] = rng.choice([x for x in "ACGT" if x != cp[i]])
            backbone = backbone[:b] + "".join(cp) + backbone[b + L:]
            segdup_pairs.append((a, b, L))

    # Copy the first deletion site, with flanking context, onto the paralogue position. Both sites then
    # spell the same sequence, so a syncmer from either occurs in both.
    if args.paralog_del:
        a0, a1 = DEL_POS - PARA_FLANK, DEL_POS + DEL_LEN + PARA_FLANK
        b0 = DEL2_POS - PARA_FLANK
        backbone = backbone[:b0] + backbone[a0:a1] + backbone[b0 + (a1 - a0):]

    # Read the tandem unit only now: the segdups above may have rewritten this stretch.
    dup_unit = backbone[DUP_POS:DUP_POS + DUP_LEN]

    # Draw from a pre-built grid of legal slots rather than rejection-sampling. With a minimum spacing
    # there are only backbone_bp/spacing slots, so rejection sampling stalls once the request approaches
    # that bound -- the acceptance probability goes to zero and the loop never finishes.
    SNP_SPACING = 100
    def in_para(q):
        if not args.paralog_del:
            return False
        return (DEL_POS - PARA_FLANK - 50 <= q < DEL_POS + DEL_LEN + PARA_FLANK + 50 or
                DEL2_POS - PARA_FLANK - 50 <= q < DEL2_POS + DEL_LEN + PARA_FLANK + 50)
    slots = [q for q in range(200, args.backbone_bp - 200, SNP_SPACING)
             if not in_sv(q) and not in_para(q)]
    if args.snps > len(slots):
        raise SystemExit(f"--snps {args.snps} exceeds the {len(slots)} slots available at "
                         f"{SNP_SPACING} bp spacing in a {args.backbone_bp} bp backbone")
    snp_pos = sorted(rng.sample(slots, args.snps))

    # Per design: state at each SV site, and an allele at each SNP position.
    designs = []
    for d in range(args.designs):
        drng = random.Random(args.seed * 1000 + d)
        alt = {}
        for p in snp_pos:
            ref_b = backbone[p]
            alt[p] = ref_b if drng.random() < 0.5 else drng.choice([b for b in "ACGT" if b != ref_b])
        designs.append({
            "DEL": (d % 2 == 1),
            "INS": (d // 2) % 2 == 1,
            "DUP": [1, 2, 3, 2][d % 4],
            "INV": (d // 4) % 2 == 1,
            "DEL2": (d // 8) % 2 == 1,      # independent of DEL, so the two sites are not correlated
            "snp": alt,
        })

    # Mirror each SNP inside a segdup copy onto its paralogue, with the SAME per-design allele. Without
    # this the duplicated stretch is identical across every haplotype, so its syncmers are constant and
    # become depth anchors rather than markers -- and confinement, which only sees markers, never fires.
    # A real paralogue pair carries the same variation twice, and that is what makes a marker ambiguous.
    snp_set = set(snp_pos)
    for a, b, L in segdup_pairs:
        for p_ in [q for q in snp_pos if a <= q < a + L]:
            q_ = b + (p_ - a)
            if in_sv(q_) or q_ in snp_set or q_ >= args.backbone_bp - 200:
                continue
            snp_pos.append(q_)
            snp_set.add(q_)
            for d in range(args.designs):
                designs[d]["snp"][q_] = designs[d]["snp"][p_]
    snp_pos.sort()

    def dup_copies(n):
        if args.dup_divergence <= 0.0:
            return [dup_unit] * n
        out = []
        crng = random.Random(args.seed * 77 + n)
        for _ in range(n):
            u = list(dup_unit)
            for i in range(len(u)):
                if crng.random() < args.dup_divergence:
                    u[i] = crng.choice([b for b in "ACGT" if b != u[i]])
            out.append("".join(u))
        return out

    # ---- graph -----------------------------------------------------------------------------
    # Walk the backbone left to right emitting shared nodes, one bubble per SNP position, and the four
    # SV structures. Node ids are assigned in reference order, which is what the downstream sort assumes.
    segs = []          # ("shared", seq) | ("snp", pos) | ("sv", which)
    cur = 0
    events = [(p, "snp", p) for p in snp_pos]
    events += [(DEL_POS, "sv", "DEL"), (INS_POS, "sv", "INS"),
               (DUP_POS, "sv", "DUP"), (INV_POS, "sv", "INV")]
    if args.paralog_del:
        events.append((DEL2_POS, "sv", "DEL2"))
    events.sort()
    for pos, kind, payload in events:
        if pos > cur:
            segs.append(("shared", backbone[cur:pos]))
        if kind == "snp":
            segs.append(("snp", pos))
            cur = pos + 1
        else:
            segs.append(("sv", payload))
            cur = pos + (DEL_LEN if payload in ("DEL", "DEL2") else INS_LEN if payload == "INS"
                         else DUP_LEN if payload == "DUP" else INV_LEN)
    if cur < len(backbone):
        segs.append(("shared", backbone[cur:]))

    nodes = []                 # (id, seq)
    walks = {f"design{d}#{t}": [] for d in range(args.designs) for t in ("a", "b")}
    ref_walk = []

    def add_node(seq):
        nodes.append((len(nodes) + 1, seq))
        return len(nodes)

    def push(hap_names, nid, orient="+"):
        for h in hap_names:
            walks[h].append((nid, orient))

    all_haps = list(walks.keys())
    for kind, payload in segs:
        if kind == "shared":
            nid = add_node(payload)
            ref_walk.append((nid, "+"))
            push(all_haps, nid)
        elif kind == "snp":
            p = payload
            bases = sorted({designs[d]["snp"][p] for d in range(args.designs)} | {backbone[p]})
            ids = {b: add_node(b) for b in bases}
            ref_walk.append((ids[backbone[p]], "+"))
            for d in range(args.designs):
                push([f"design{d}#a", f"design{d}#b"], ids[designs[d]["snp"][p]])
        else:
            if payload in ("DEL", "DEL2"):
                p0 = DEL_POS if payload == "DEL" else DEL2_POS
                nid = add_node(backbone[p0:p0 + DEL_LEN])
                ref_walk.append((nid, "+"))
                for d in range(args.designs):
                    if not designs[d][payload]:
                        push([f"design{d}#a", f"design{d}#b"], nid)
            elif payload == "INS":
                base = add_node(backbone[INS_POS:INS_POS + INS_LEN])
                alt = add_node(ins_payload)
                ref_walk.append((base, "+"))
                for d in range(args.designs):
                    hs = [f"design{d}#a", f"design{d}#b"]
                    push(hs, base)
                    if designs[d]["INS"]:
                        push(hs, alt)
            elif payload == "DUP":
                # One node per copy index; a haplotype with N copies traverses the first N. Identical
                # copies (default) make copy number show up purely as marker multiplicity.
                maxcn = max(designs[d]["DUP"] for d in range(args.designs))
                copies = dup_copies(maxcn)
                ids = [add_node(c) for c in copies]
                ref_walk.append((ids[0], "+"))
                for d in range(args.designs):
                    hs = [f"design{d}#a", f"design{d}#b"]
                    for k in range(designs[d]["DUP"]):
                        push(hs, ids[k])
            else:  # INV -- the same node traversed in reverse, as a real graph represents it
                nid = add_node(backbone[INV_POS:INV_POS + INV_LEN])
                ref_walk.append((nid, "+"))
                for d in range(args.designs):
                    push([f"design{d}#a", f"design{d}#b"], nid,
                         "-" if designs[d]["INV"] else "+")

    # ---- sequences, straight from the walks so graph and FASTA cannot disagree ---------------
    seq_of = {nid: s for nid, s in nodes}

    def spell(walk):
        return "".join(seq_of[n] if o == "+" else revcomp(seq_of[n]) for n, o in walk)

    ref_seq = spell(ref_walk)
    hap_seq = {h: spell(w) for h, w in walks.items()}

    write_fasta(os.path.join(args.out, "ref.fa"), "ref", ref_seq)
    with open(os.path.join(args.out, "haplotypes.fa"), "w") as allf:
        for h in all_haps:
            allf.write(f">{h}\n")
            for i in range(0, len(hap_seq[h]), 60):
                allf.write(hap_seq[h][i:i + 60] + "\n")
            write_fasta(os.path.join(args.out, f"hap_{h.replace('#', '_')}.fa"), h, hap_seq[h])

    # ---- GFA -------------------------------------------------------------------------------
    edges = set()
    for w in [ref_walk] + list(walks.values()):
        for i in range(1, len(w)):
            edges.add((w[i - 1][0], w[i - 1][1], w[i][0], w[i][1]))
    with open(os.path.join(args.out, "graph.gfa"), "w") as f:
        f.write("H\tVN:Z:1.1\n")
        for nid, s in nodes:
            f.write(f"S\t{nid}\t{s}\n")
        for a, ao, b, bo in sorted(edges):
            f.write(f"L\t{a}\t{ao}\t{b}\t{bo}\t0M\n")
        f.write("P\tref\t" + ",".join(f"{n}{o}" for n, o in ref_walk) + "\t*\n")
        for h in all_haps:
            f.write(f"P\t{h}\t" + ",".join(f"{n}{o}" for n, o in walks[h]) + "\t*\n")

    # ---- truth -----------------------------------------------------------------------------
    with open(os.path.join(args.out, "truth.variants.tsv"), "w") as f:
        f.write("design\tDEL\tINS\tDUP_copies\tINV\n")
        for d in range(args.designs):
            s = designs[d]
            f.write(f"design{d}\t{int(s['DEL'])}\t{int(s['INS'])}\t{s['DUP']}\t{int(s['INV'])}\n")
    with open(os.path.join(args.out, "truth.haplotypes.tsv"), "w") as f:
        f.write("haplotype\tdesign\ttwin\tbp\n")
        for h in all_haps:
            d = h.split("#")[0]
            twin = d + "#" + ("b" if h.endswith("a") else "a")
            f.write(f"{h}\t{d}\t{twin}\t{len(hap_seq[h])}\n")
    with open(os.path.join(args.out, "truth.bubbles.tsv"), "w") as f:
        f.write("site\ttype\tref_pos\tsize\tn_alleles_expected\n")
        f.write(f"DEL\tDEL\t{DEL_POS}\t{DEL_LEN}\t2\n")
        f.write(f"INS\tINS\t{INS_POS}\t{INS_LEN}\t2\n")
        f.write(f"DUP\tDUP\t{DUP_POS}\t{DUP_LEN}\t3\n")
        f.write(f"INV\tINV\t{INV_POS}\t{INV_LEN}\t2\n")

    if segdup_pairs:
        print("  segdups: " + ", ".join(f"{a}<->{b} ({L}bp, {100*args.segdup_divergence:.1f}% div)"
                                        for a, b, L in segdup_pairs))
    print(f"{len(all_haps)} haplotypes ({args.designs} designs x 2 twins), "
          f"{len(nodes)} nodes, {len(edges)} edges, {len(snp_pos)} SNP sites")
    print(f"  ref {len(ref_seq)} bp; DEL@{DEL_POS}/{DEL_LEN} INS@{INS_POS}/{INS_LEN} "
          f"DUP@{DUP_POS}/{DUP_LEN} INV@{INV_POS}/{INV_LEN}")
    for d in range(args.designs):
        s = designs[d]
        print(f"  design{d}: DEL={int(s['DEL'])} INS={int(s['INS'])} DUP={s['DUP']} "
              f"INV={int(s['INV'])}  {len(hap_seq[f'design{d}#a'])} bp")


if __name__ == "__main__":
    main()
