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
    ap.add_argument("--abut", action="store_true",
                    help="place the DEL and INS sites back to back with no backbone between them, so "
                         "two bubbles share a boundary node. build_block_chain drops the backbone when "
                         "source == sink, a branch nothing else exercises.")
    ap.add_argument("--mosaics", type=int, default=0,
                    help="emit N extra haplotypes that are NOT in the panel, each switching donor "
                         "design at every SV site and every inter-SV stretch. They test the chain the "
                         "hard way: the backbones that carry a marker-poor bubble are themselves from "
                         "different donors, so nothing along the chain is stable and the model must "
                         "switch almost everywhere. Truth is written per region.")
    ap.add_argument("--nested-del", action="store_true",
                    help="plant a large polymorphic deletion spanning ALL the inner SV sites. It "
                         "creates a bubble that CONTAINS the others, which is ankrd36c's structure and "
                         "the case the block chain assumes cannot happen.")
    ap.add_argument("--twin-divergence", type=int, default=0,
                    help="SNPs by which a design's second copy differs from its first. 0 (default) "
                         "makes them exact twins, so leave-one-out has a perfectly representable "
                         "answer. Raise it to test the realistic case: the sample's own haplotypes are "
                         "gone and only a NEAR-identical one remains, where the right answer is no "
                         "longer an exact allele match but the most similar allele available.")
    ap.add_argument("--dup-unit-bp", type=int, default=500,
                    help="tandem repeat unit length. Raise it with --dup-max-cn for a VNTR-scale array "
                         "like LPA's KIV-2 or ANKRD36C's, where copy number IS the variant and marker "
                         "multiplicity is the only signal that can measure it.")
    ap.add_argument("--dup-fold-divergence", type=float, default=0.0,
                    help="graph keeps ONE consensus repeat unit traversed N times while the emitted "
                         "haplotypes carry copies diverging from it at this per-base rate, which is what "
                         "panphorte normalization produces. Reads then come from real copies and the "
                         "graph from a consensus, so k-mer and alignment evidence diverge sharply")
    ap.add_argument("--dup-node-per-copy", action="store_true",
                    help="emit one node per array copy instead of revisiting a single node. That is "
                         "the opposite of what a real pangenome graph does (lpa: 75%% of path steps "
                         "revisit a node, no >100bp sequence under two ids) and it makes the array "
                         "easier, so it exists only to reproduce results measured before the fix")
    ap.add_argument("--dup-per-design-divergence", type=float, default=0.0,
                    help="give EVERY design its own array copies, diverged from the base unit by this "
                         "per-base rate. Without it all designs share one set of copy nodes and the "
                         "array bubble has only as many alleles as there are distinct copy numbers -- "
                         "heavily shared, so leave-one-out always leaves a carrier. With it each "
                         "design's array is unique, reproducing what a real VNTR looks like: lpa's "
                         "KIV-2 block carries 457 distinct alleles among 466 haplotypes, so removing a "
                         "haplotype removes the only carrier of its allele.")
    ap.add_argument("--dup-min-cn", type=int, default=1)
    ap.add_argument("--dup-max-cn", type=int, default=3)
    ap.add_argument("--paralog-ins", action="store_true",
                    help="a second insertion site carrying the SAME payload as the first, in a "
                         "different block, polymorphic independently. The payload's syncmers then "
                         "occur in both sites and vary in both, so confinement strips them; junction "
                         "syncmers differ, so the bubble is left partly rather than wholly stripped.")
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
    # --abut puts INS immediately after DEL, so the two bubbles share a boundary and no backbone lies
    # between them.
    INS_POS, INS_LEN = (int(args.backbone_bp * 0.16) + 400 if args.abut
                        else int(args.backbone_bp * 0.36)), 600
    DUP_POS, DUP_LEN = int(args.backbone_bp * 0.56), args.dup_unit_bp
    INV_POS, INV_LEN = int(args.backbone_bp * 0.76), 700
    ins_payload = random_seq(rng, INS_LEN)
    # The paralogue sits inside the backbone stretch between INS and DUP, i.e. a different block.
    DEL2_POS = int(args.backbone_bp * 0.46)
    INS2_POS = int(args.backbone_bp * 0.66)
    PARA_FLANK = 400
    sv_spans = [(DEL_POS, DEL_POS + DEL_LEN), (INS_POS, INS_POS + INS_LEN),
                (DUP_POS, DUP_POS + DUP_LEN), (INV_POS, INV_POS + INV_LEN)]
    if args.paralog_del:
        sv_spans.append((DEL2_POS, DEL2_POS + DEL_LEN))
    if args.paralog_ins:
        sv_spans.append((INS2_POS, INS2_POS + INS_LEN))

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
            "DUP": args.dup_min_cn + (d % max(1, args.dup_max_cn - args.dup_min_cn + 1)),
            "INV": (d // 4) % 2 == 1,
            "DEL2": (d // 8) % 2 == 1,      # independent of DEL, so the two sites are not correlated
            "INS2": (d // 16) % 2 == 1,
            "NEST": (d % 5 == 4),          # a minority carry the container deletion
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

    # Near-twin mode: give the "b" copy a different allele at a few SNP sites, so holding out a
    # sample leaves a haplotype that is CLOSE but not identical. The right answer is then no longer an
    # exact allele match -- it is the most similar allele available -- which is the realistic case.
    twin_snp = {}
    if args.twin_divergence > 0:
        trng = random.Random(args.seed * 31337)
        for d in range(args.designs):
            for p_ in trng.sample(snp_pos, min(args.twin_divergence, len(snp_pos))):
                cur_b = designs[d]["snp"][p_]
                twin_snp[(f"design{d}#b", p_)] = trng.choice([x for x in "ACGT" if x != cur_b])

    # ---- graph -----------------------------------------------------------------------------
    # Walk the backbone left to right emitting shared nodes, one bubble per SNP position, and the four
    # SV structures. Node ids are assigned in reference order, which is what the downstream sort assumes.
    NEST_A = DEL_POS - 1000
    NEST_B = INV_POS + INV_LEN + 1000
    segs = []          # (start_pos, kind, payload)
    cur = 0
    events = [(p, "snp", p) for p in snp_pos]
    events += [(DEL_POS, "sv", "DEL"), (INS_POS, "sv", "INS"),
               (DUP_POS, "sv", "DUP"), (INV_POS, "sv", "INV")]
    if args.paralog_del:
        events.append((DEL2_POS, "sv", "DEL2"))
    if args.paralog_ins:
        events.append((INS2_POS, "sv", "INS2"))
    events.sort()
    for pos, kind, payload in events:
        if pos > cur:
            segs.append((cur, "shared", backbone[cur:pos]))
        if kind == "snp":
            segs.append((pos, "snp", pos))
            cur = pos + 1
        else:
            segs.append((pos, "sv", payload))
            cur = pos + (DEL_LEN if payload in ("DEL", "DEL2")
                         else INS_LEN if payload in ("INS", "INS2")
                         else DUP_LEN if payload == "DUP" else INV_LEN)
    if cur < len(backbone):
        segs.append((cur, "shared", backbone[cur:]))

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
    # Haplotypes carrying the nested deletion skip every segment inside its span, which is what makes
    # the enclosing bubble contain the inner ones.
    def carriers(exclude_nested):
        if not (args.nested_del and exclude_nested):
            return all_haps
        return [h for h in all_haps if not designs[int(h.split("#")[0][6:])]["NEST"]]

    # Where each seg's steps land in each haplotype's walk, so a mosaic can be spliced out of walks
    # that are already correct rather than rebuilt from scratch. A donor that bypasses a seg simply
    # contributed no steps there, and the slice is empty -- which is the right answer.
    seg_marks = {h: [] for h in walks}
    for seg_pos, kind, payload in segs:
        for h in seg_marks:
            seg_marks[h].append(len(walks[h]))
        inside_nest = args.nested_del and NEST_A <= seg_pos < NEST_B
        here = carriers(inside_nest)
        if kind == "shared":
            nid = add_node(payload)
            ref_walk.append((nid, "+"))
            push(here, nid)
        elif kind == "snp":
            p = payload
            bases = sorted({designs[d]["snp"][p] for d in range(args.designs)} | {backbone[p]}
                           | {v for (h, q), v in twin_snp.items() if q == p})
            ids = {b: add_node(b) for b in bases}
            ref_walk.append((ids[backbone[p]], "+"))
            for d in range(args.designs):
                for h in (f"design{d}#a", f"design{d}#b"):
                    if h in here:
                        push([h], ids[twin_snp.get((h, p), designs[d]["snp"][p])])
        else:
            if payload in ("DEL", "DEL2"):
                p0 = DEL_POS if payload == "DEL" else DEL2_POS
                nid = add_node(backbone[p0:p0 + DEL_LEN])
                ref_walk.append((nid, "+"))
                for d in range(args.designs):
                    if not designs[d][payload]:
                        push([h for h in (f"design{d}#a", f"design{d}#b") if h in here], nid)
            elif payload in ("INS", "INS2"):
                p0 = INS_POS if payload == "INS" else INS2_POS
                base = add_node(backbone[p0:p0 + INS_LEN])
                # A distinct node carrying IDENTICAL sequence: two independent insertions of the same
                # payload, which is what makes the payload's syncmers ambiguous between the two sites.
                alt = add_node(ins_payload)
                ref_walk.append((base, "+"))
                for d in range(args.designs):
                    hs = [h for h in (f"design{d}#a", f"design{d}#b") if h in here]
                    push(hs, base)
                    if designs[d][payload]:
                        push(hs, alt)
            elif payload == "DUP":
                # A real pangenome graph collapses identical sequence into ONE node and expresses
                # repetition as revisiting it, so a haplotype with N copies traverses the same node N
                # times and copy number is that node's traversal multiplicity. Measured on the lpa
                # graph: 75% of path steps revisit a node already used, and no sequence over 100 bp is
                # stored under two node ids.
                #
                # This generator used to emit one node per copy index instead, which is the opposite
                # structure -- N distinct nodes each traversed once, 0% revisits. That made the array
                # easier than the real thing in a way that matters: read placement is unambiguous when
                # the copies are separate nodes, and per-node coverage then needs the opposite
                # weighting. Any result measured on the old shape should be re-measured on this one.
                # `--dup-node-per-copy` restores it for exactly that comparison.
                maxcn = max(designs[d]["DUP"] for d in range(args.designs))
                copies = dup_copies(maxcn)
                ids = [add_node(c) for c in copies] if args.dup_node_per_copy else [add_node(copies[0])]
                ref_walk.append((ids[0], "+"))
                if args.dup_per_design_divergence > 0.0:
                    # Each design walks its OWN copies, so every array allele is distinct. The unit is
                    # diverged ONCE per design and then repeated, so copies within a haplotype stay
                    # identical and copy number still shows up as marker MULTIPLICITY -- which is what a
                    # real tandem array looks like, the copies being a recent duplication of one unit.
                    # Diverging each copy independently instead would make every copy's syncmers unique,
                    # so multiplicity would never accumulate and copy number would be carried by marker
                    # count alone. That is a different and much harder problem than the real one.
                    for d in range(args.designs):
                        drng2 = random.Random(args.seed * 991 + d)
                        u = list(dup_unit)
                        for t in range(len(u)):
                            if drng2.random() < args.dup_per_design_divergence:
                                u[t] = drng2.choice([b for b in "ACGT" if b != u[t]])
                        unit_d = "".join(u)
                        # One node per DESIGN (its own diverged unit), traversed as many times as
                        # that design has copies -- distinct unit variants are distinct nodes, and
                        # repetition of a variant is revisiting its node, which is what a real graph does.
                        own = ([add_node(unit_d) for _ in range(designs[d]["DUP"])]
                               if args.dup_node_per_copy
                               else [add_node(unit_d)] * designs[d]["DUP"])
                        hs = [h for h in (f"design{d}#a", f"design{d}#b") if h in here]
                        for nid2 in own:
                            push(hs, nid2)
                else:
                    for d in range(args.designs):
                        hs = [h for h in (f"design{d}#a", f"design{d}#b") if h in here]
                        for k in range(designs[d]["DUP"]):
                            push(hs, ids[k] if args.dup_node_per_copy else ids[0])
            else:  # INV -- the same node traversed in reverse, as a real graph represents it
                nid = add_node(backbone[INV_POS:INV_POS + INV_LEN])
                ref_walk.append((nid, "+"))
                for d in range(args.designs):
                    push([h for h in (f"design{d}#a", f"design{d}#b") if h in here], nid,
                         "-" if designs[d]["INV"] else "+")

    for h in seg_marks:
        seg_marks[h].append(len(walks[h]))

    # ---- mosaic haplotypes: a different donor design in every region of the chain ------------
    # These are NOT in the panel. Each switches donor at every SV site and every stretch between two,
    # so the backbones that carry a marker-poor bubble are themselves from different donors and nothing
    # along the chain is stable. This is the hard case for a linkage-based model.
    mosaic_walks = {}
    mosaic_truth = []      # (name, region, kind, donor)
    if args.mosaics > 0:
        mrng = random.Random(args.seed * 4242)
        for m in range(args.mosaics):
            name = f"mosaic{m}"
            walk = []
            region = 0
            donors = {}
            for i, (seg_pos, kind, payload) in enumerate(segs):
                if kind == "sv":
                    region += 1
                if region not in donors:
                    donors[region] = mrng.randrange(args.designs)
                src = f"design{donors[region]}#a"
                walk.extend(walks[src][seg_marks[src][i]:seg_marks[src][i + 1]])
                if kind == "sv":
                    region += 1
                    donors.setdefault(region, mrng.randrange(args.designs))
            mosaic_walks[name] = walk
            for r, d in sorted(donors.items()):
                mosaic_truth.append((name, r, "sv" if r % 2 == 1 else "between", f"design{d}"))

    # ---- sequences, straight from the walks so graph and FASTA cannot disagree ---------------
    seq_of = {nid: s for nid, s in nodes}

    def spell(walk):
        return "".join(seq_of[n] if o == "+" else revcomp(seq_of[n]) for n, o in walk)

    # --dup-fold-divergence reproduces what panphorte's normalization does to a tandem array, which is
    # the ONE structure a real folded graph has that no other fixture here reproduces: the GRAPH holds a
    # single consensus repeat unit traversed N times, while the SAMPLE's copies each differ from that
    # consensus by a percent or two. Graph and haplotypes deliberately disagree, which nothing else in
    # this generator does -- everywhere else they are emitted from the same edits precisely so they
    # cannot.
    #
    # It matters because it separates two things that lpa conflates. A k-mer spanning a copy that
    # diverges 2% from the consensus fails to match with probability 1 - 0.98^31, about 47%, so the
    # marker path should collapse here. Coverage should not: an aligner is untroubled by 2%, and copy
    # number is the consensus node's traversal count. If coverage works here and fails on lpa, the lpa
    # failure is something else and this fixture localises it.
    fold_rng = random.Random(args.seed * 7717)
    def spell_sample(walk):
        out = []
        seen = {}
        for n, o in walk:
            base = seq_of[n]
            if args.dup_fold_divergence > 0.0 and len(base) >= args.dup_unit_bp:
                # Each traversal of the folded node is an independent copy in the sample, so successive
                # visits diverge independently -- which is what makes them real copies rather than a
                # repeat of one consensus.
                seen[n] = seen.get(n, 0) + 1
                u = list(base)
                for t in range(len(u)):
                    if fold_rng.random() < args.dup_fold_divergence:
                        u[t] = fold_rng.choice([b for b in "ACGT" if b != u[t]])
                base = "".join(u)
            out.append(base if o == "+" else revcomp(base))
        return "".join(out)

    ref_seq = spell(ref_walk)
    hap_seq = {h: spell_sample(w) for h, w in walks.items()}
    mosaic_seq = {n: spell(w) for n, w in mosaic_walks.items()}

    write_fasta(os.path.join(args.out, "ref.fa"), "ref", ref_seq)
    with open(os.path.join(args.out, "haplotypes.fa"), "w") as allf:
        for h in all_haps:
            allf.write(f">{h}\n")
            for i in range(0, len(hap_seq[h]), 60):
                allf.write(hap_seq[h][i:i + 60] + "\n")
            write_fasta(os.path.join(args.out, f"hap_{h.replace('#', '_')}.fa"), h, hap_seq[h])

    # ---- GFA -------------------------------------------------------------------------------
    edges = set()
    for w in [ref_walk] + list(walks.values()) + list(mosaic_walks.values()):
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
        # Mosaics go in as paths so their per-block truth can be resolved by spelling, exactly as for a
        # panel haplotype. They must then be excluded from the panel when genotyping -- which is what
        # leave-one-out already does -- or the test is trivial.
        for n, w in mosaic_walks.items():
            f.write(f"P\t{n}\t" + ",".join(f"{a}{o}" for a, o in w) + "\t*\n")

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
    for n, seq in mosaic_seq.items():
        write_fasta(os.path.join(args.out, f"hap_{n}.fa"), n, seq)
    if mosaic_truth:
        with open(os.path.join(args.out, "truth.mosaics.tsv"), "w") as f:
            f.write("mosaic\tregion\tregion_kind\tdonor\n")
            for name, r, kind, d in mosaic_truth:
                f.write(f"{name}\t{r}\t{kind}\t{d}\n")

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
              f"INV={int(s['INV'])}"
              + (f" DEL2={int(s['DEL2'])}" if args.paralog_del else "")
              + (f" INS2={int(s['INS2'])}" if args.paralog_ins else "")
              + f"  {len(hap_seq[f'design{d}#a'])} bp")


if __name__ == "__main__":
    main()
