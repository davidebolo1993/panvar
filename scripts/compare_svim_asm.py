#!/usr/bin/env python3
"""Compare panvar's per-haplotype CN/CNBP against svim-asm, an assembly-based SV caller.

panvar reads structure off the graph; svim-asm reads it off a minimap2 alignment of one assembly to
another. They share no code and no assumptions, so agreement is evidence and disagreement is a lead.

The comparison is per HAPLOTYPE, against that haplotype's own CNBP -- not against the region VCF's
representative SVLEN, which is one merged event's size and is not what any particular carrier holds.

Coordinates stay in the reference PATH's own frame throughout: the FASTA handed to minimap2 is the
spelled reference path, so svim-asm's positions and the bubble spans computed here are already the same
frame and nothing has to be projected.

Its own uncertainty is reported rather than absorbed, because this script is as capable of manufacturing
a disagreement as of finding one: calls no bubble could claim, calls equidistant between two spans, how
many calls were summed into each comparison, and comparisons whose net cancels a gain against a loss.
A number here is only worth acting on once those read zero.

The ceiling, stated because it is easy to overclaim: both sides start from the same graph haplotypes.
Agreement shows CNBP corresponds to an assembly-level SV length, not that the graph is biologically
right. That would need the original assemblies, which are not in the repository.

  compare_svim_asm.py --gfa bubble.sorted.gfa --bubbles bubble.bubbles.csv --vcf call.region.vcf \
                      --reference <path name> --out-dir out [--samples N] [--threads N]
"""
import argparse
import csv
import gzip
import os
import random
import re
import shutil
import subprocess
import statistics
import sys


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_gfa(path):
    """segment sequences, and every path as a list of (node, is_reverse)."""
    seq, paths = {}, {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("S\t"):
                f = line.rstrip("\n").split("\t")
                seq[f[1]] = f[2]
            elif line.startswith("P\t"):
                f = line.rstrip("\n").split("\t")
                paths[f[1]] = [(s[:-1], s[-1] == "-") for s in f[2].split(",")]
    return seq, paths


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def spell(seq, steps):
    out = []
    for node, rev in steps:
        s = seq.get(node, "")
        out.append(s.translate(_COMP)[::-1] if rev else s)
    return "".join(out)


def ref_offsets(seq, steps):
    """cumulative start offset of every step, and a node -> [starts] index."""
    starts, pos, by_node = [], 0, {}
    for node, _rev in steps:
        starts.append(pos)
        by_node.setdefault(node, []).append(pos)
        pos += len(seq.get(node, ""))
    return starts, by_node, pos


def read_vcf(path):
    hdr, rows = None, []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                hdr = line.rstrip("\n").split("\t")
            elif not line.startswith("#"):
                rows.append(line.rstrip("\n").split("\t"))
    return hdr, rows


def info_of(field):
    d = {}
    for kv in field.split(";"):
        k, _, v = kv.partition("=")
        d[k] = v if _ else True
    return d


def safe(name):
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, **kw)


def svim_events(vcf_path):
    """(pos, signed_len) for every length-changing call svim-asm made."""
    out = []
    if not os.path.exists(vcf_path):
        return out
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            inf = info_of(f[7])
            t = inf.get("SVTYPE")
            if t not in ("DEL", "INS", "DUP", "DUP:TANDEM"):
                continue
            try:
                out.append((int(f[1]), int(inf.get("SVLEN", 0))))
            except ValueError:
                continue
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gfa", required=True, help="the graph call was run on")
    ap.add_argument("--bubbles", required=True, help="the bubbles CSV call was run with")
    ap.add_argument("--vcf", required=True, help="call's region VCF")
    ap.add_argument("--reference", required=True, help="reference path name")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--samples", type=int, default=40,
                    help="how many carrier haplotypes to compare (0 = all)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--bin-dir", default="",
                    help="directory holding minimap2, samtools and svim-asm")
    ap.add_argument("--margin", type=int, default=2000,
                    help="bp of slack when attributing an svim-asm call to a bubble span; breakpoints "
                         "from an alignment do not land where graph boundaries do")
    ap.add_argument("--min-size-bp", type=int, default=1000,
                    help="relative agreement (within 5%%/10%%) is reported only for haplotypes whose "
                         "CNBP reaches this; a 40 bp disagreement on a 50 bp delta is not a 5%% question")
    ap.add_argument("--keep", action="store_true", help="keep per-haplotype BAM/VCF working files")
    args = ap.parse_args()

    # --bin-dir pins a directory holding the three external tools; empty means look on PATH, which is
    # what a machine other than the author's will have.
    env_bin = args.bin_dir
    for tool in ("minimap2", "samtools", "svim-asm"):
        if env_bin:
            if not os.path.exists(os.path.join(env_bin, tool)):
                sys.exit(f"compare_svim_asm: {tool} not found in {env_bin}")
        elif shutil.which(tool) is None:
            sys.exit(f"compare_svim_asm: {tool} not found on PATH (or pass --bin-dir)")

    seq, paths = load_gfa(args.gfa)
    if args.reference not in paths:
        sys.exit(f"compare_svim_asm: reference path not in the graph: {args.reference}")

    _starts, ref_by_node, ref_len = ref_offsets(seq, paths[args.reference])

    # bubble -> [lo, hi) in the reference path's own coordinates
    span = {}
    with open(args.bubbles) as fh:
        for row in csv.DictReader(fh):
            src, snk = row["source"], row["sink"]
            if src not in ref_by_node or snk not in ref_by_node:
                continue
            lo = min(ref_by_node[src] + ref_by_node[snk])
            hi = max(ref_by_node[src] + ref_by_node[snk]) + len(seq.get(snk, ""))
            span[row["bubble_id"]] = (lo, hi)

    hdr, rows = read_vcf(args.vcf)
    dup_rows = [r for r in rows if info_of(r[7]).get("SVTYPE") == "DUP"]
    if not dup_rows:
        sys.exit("compare_svim_asm: the region VCF has no DUP records")

    # a bubble is "pure" when the DUP is the only record there; anything else is a compound site whose
    # CNBP legitimately carries content beyond the copy-number event
    per_bubble = {}
    for r in rows:
        per_bubble.setdefault(info_of(r[7]).get("BUBBLE_ID"), []).append(r)

    # every haplotype that has a CNBP on some DUP record, with the CN it was called at
    cn_of = {}
    for r in dup_rows:
        keys = r[8].split(":")
        if "CNBP" not in keys or "CN" not in keys:
            continue
        bi, ci = keys.index("CNBP"), keys.index("CN")
        for i in range(9, len(hdr)):
            v = r[i].split(":")
            if len(v) > bi and v[bi] not in (".", "") and v[ci].lstrip("-").isdigit():
                cn_of.setdefault(hdr[i], {})[info_of(r[7]).get("BUBBLE_ID")] = int(v[ci])
    cn_of.pop(args.reference, None)
    wanted = sorted(cn_of)

    # Stratified, not uniform. A random draw from a cohort where two thirds sit at one CN can easily
    # contain no gain, no zero and no reference-like haplotype at all -- and those are the classes worth
    # checking. The extremes and one representative of each behaviour are taken first, random fill after.
    if args.samples and len(wanted) > args.samples:
        ref_cn = {}
        for r in dup_rows:
            inf = info_of(r[7])
            try:
                ref_cn[inf.get("BUBBLE_ID")] = int(inf.get("REF_CN", "0"))
            except ValueError:
                ref_cn[inf.get("BUBBLE_ID")] = 0
        must = []
        for bid in sorted({b for m in cn_of.values() for b in m}):
            here = {h: m[bid] for h, m in cn_of.items() if bid in m}
            if not here:
                continue
            rc = ref_cn.get(bid, 0)
            lo = min(here, key=lambda h: here[h])
            hi = max(here, key=lambda h: here[h])
            zero = next((h for h in sorted(here) if here[h] == 0), None)
            like = next((h for h in sorted(here) if here[h] == rc), None)
            loss = next((h for h in sorted(here) if here[h] < rc), None)
            gain = next((h for h in sorted(here) if here[h] > rc), None)
            must += [h for h in (lo, hi, zero, like, loss, gain) if h]
        must = sorted(dict.fromkeys(must))[:args.samples]
        rest = [h for h in wanted if h not in set(must)]
        random.Random(args.seed).shuffle(rest)
        wanted = sorted(must + rest[:max(0, args.samples - len(must))])
        print(f"stratified selection: {len(must)} required (per-bubble min/max CN, zero, "
              f"reference-like, a loss and a gain) + {len(wanted) - len(must)} random")

    os.makedirs(args.out_dir, exist_ok=True)
    ref_fa = os.path.join(args.out_dir, "reference.fa")
    with open(ref_fa, "w") as fh:
        fh.write(">reference\n" + spell(seq, paths[args.reference]) + "\n")
    run([os.path.join(env_bin, "samtools"), "faidx", ref_fa])

    records = []
    unmatched = []     # svim-asm calls no bubble span could claim
    ambiguous = []     # calls equidistant from two spans, attributed by an arbitrary tie-break
    for n, hap in enumerate(wanted, 1):
        if hap not in paths:
            continue
        tag = safe(hap)
        work = os.path.join(args.out_dir, tag)
        os.makedirs(work, exist_ok=True)
        hap_fa = os.path.join(work, "hap.fa")
        with open(hap_fa, "w") as fh:
            fh.write(">" + tag + "\n" + spell(seq, paths[hap]) + "\n")
        bam = os.path.join(work, "aln.bam")
        with open(os.path.join(work, "aln.sam"), "wb") as sam:
            sam.write(run([os.path.join(env_bin, "minimap2"), "-a", "-x", "asm5", "--cs", "-r2k",
                           "-t", str(args.threads), ref_fa, hap_fa]).stdout)
        run([os.path.join(env_bin, "samtools"), "sort", "-o", bam,
             os.path.join(work, "aln.sam")])
        run([os.path.join(env_bin, "samtools"), "index", bam])
        run([os.path.join(env_bin, "svim-asm"), "haploid", work, bam, ref_fa])
        events = svim_events(os.path.join(work, "variants.vcf"))
        sys.stderr.write(f"\r  [{n}/{len(wanted)}] {tag[:44]:44s} {len(events)} svim-asm calls")
        sys.stderr.flush()

        # Each svim-asm call belongs to ONE bubble: the span containing it, else the nearest span
        # within the margin. Summing every call within a margin of each span instead double-counted a
        # neighbouring site -- at GSTM1 that pulled bubble 4's separate 4236 bp deletion into bubble 8
        # and looked like a 4 kb disagreement.
        attributed = {}
        for p, l in events:
            ranked = sorted(((0 if lo <= p <= hi else min(abs(p - lo), abs(p - hi)), bid)
                             for bid, (lo, hi) in span.items()))
            if not ranked or ranked[0][0] > args.margin:
                unmatched.append((hap, p, l, ranked[0][0] if ranked else -1))
                continue
            best_d, best = ranked[0]
            # A call equidistant from two spans, or inside neither and near both, is attributed by an
            # arbitrary tie-break. Counting those rather than hiding them is the difference between a
            # disagreement that means something and one that is an artefact of this script.
            if len(ranked) > 1 and ranked[1][0] == best_d:
                ambiguous.append((hap, p, l, best, ranked[1][1], best_d))
            attributed.setdefault(best, []).append((p, l))

        for r in dup_rows:
            bid = info_of(r[7]).get("BUBBLE_ID")
            if bid not in span:
                continue
            keys = r[8].split(":")
            if "CNBP" not in keys or "CN" not in keys:
                continue
            bi, ci = keys.index("CNBP"), keys.index("CN")
            col = hdr.index(hap)
            v = r[col].split(":")
            if len(v) <= bi or v[bi] in (".", ""):
                continue
            calls = attributed.get(bid, [])
            svim_bp = sum(l for _p, l in calls)
            # Signed summation can hide a deletion and an insertion cancelling, so the events are kept
            # alongside the total: a net of zero from one 0 bp call and a net of zero from +5k/-5k are
            # very different findings.
            records.append({
                "bubble": bid, "variant": r[2], "hap": hap,
                "cn": v[ci], "cnbp": int(v[bi]), "svim_bp": svim_bp,
                "n_svim_calls": len(calls),
                "svim_gain": sum(l for _p, l in calls if l > 0),
                "svim_loss": sum(l for _p, l in calls if l < 0),
                "svim_detail": ";".join(f"{p}:{l:+d}" for p, l in sorted(calls)) or ".",
                "n_records_in_bubble": len(per_bubble.get(bid, [])),
                "class": "pure_module" if len(per_bubble.get(bid, [])) == 1 else "compound_site",
            })
        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)
    sys.stderr.write("\n")

    out_tsv = os.path.join(args.out_dir, "svim_asm_comparison.tsv")
    with open(out_tsv, "w") as fh:
        fh.write("bubble_id\tvariant_id\thaplotype\tcn\tcnbp\tsvim_asm_bp\tdiff_bp\tn_svim_calls\t"
                 "svim_gain_bp\tsvim_loss_bp\trecords_in_bubble\tclass\tsvim_calls\n")
        for r in records:
            fh.write(f"{r['bubble']}\t{r['variant']}\t{r['hap']}\t{r['cn']}\t{r['cnbp']}\t"
                     f"{r['svim_bp']}\t{r['svim_bp'] - r['cnbp']}\t{r['n_svim_calls']}\t"
                     f"{r['svim_gain']}\t{r['svim_loss']}\t{r['n_records_in_bubble']}\t"
                     f"{r['class']}\t{r['svim_detail']}\n")

    def report(label, subset):
        subset = [r for r in subset if r["cnbp"] != 0 or r["svim_bp"] != 0]
        if not subset:
            print(f"{label:22s} (no non-zero comparisons)")
            return
        diffs = [r["svim_bp"] - r["cnbp"] for r in subset]
        # Relative agreement is only a meaningful question at SV scale: a 40 bp difference on a 50 bp
        # delta is 80% and says nothing, while the same 40 bp on an 18 kb deletion is 0.2%.
        big = [(r, d) for r, d in zip(subset, diffs) if abs(r["cnbp"]) >= args.min_size_bp]
        line = (f"{label:22s} n={len(subset):5d}  "
                f"median|err|={statistics.median(abs(d) for d in diffs):7.0f} bp  "
                f"signed bias={statistics.mean(diffs):+9.1f} bp")
        if big:
            rel = [abs(d) / abs(r["cnbp"]) for r, d in big]
            line += (f"  | >={args.min_size_bp}bp: n={len(big):4d}"
                     f"  within 5%={100.0 * sum(1 for x in rel if x <= 0.05) / len(rel):5.1f}%"
                     f"  within 10%={100.0 * sum(1 for x in rel if x <= 0.10) / len(rel):5.1f}%")
        print(line)

    print(f"\nwrote {out_tsv}  ({len(records)} comparisons over {len(wanted)} haplotypes)")
    n_calls = sum(r["n_svim_calls"] for r in records)
    print(f"svim-asm calls summed into a comparison: {n_calls}"
          f"   attributed to no bubble: {len(unmatched)}"
          f"   equidistant between two: {len(ambiguous)}")
    if unmatched:
        u = sorted(unmatched, key=lambda x: -abs(x[2]))[:4]
        print("  largest unattributed: " +
              ", ".join(f"{l:+d} bp at {p} ({d} bp from the nearest span)" for _h, p, l, d in u))
    if ambiguous:
        a = ambiguous[:3]
        print("  ambiguous: " +
              ", ".join(f"{l:+d} bp at {p} -> bubble {b1} or {b2}, both {d} bp away"
                        for _h, p, l, b1, b2, d in a))
    cancel = [r for r in records if r["svim_gain"] > 0 and r["svim_loss"] < 0]
    if cancel:
        print(f"  {len(cancel)} comparison(s) sum a gain and a loss together; the net can hide both "
              f"(see svim_calls in the TSV)")
    report("all", records)
    for cls in ("pure_module", "compound_site"):
        report(cls, [r for r in records if r["class"] == cls])
    for bid in sorted({r["bubble"] for r in records}, key=lambda x: int(x)):
        sub = [r for r in records if r["bubble"] == bid]
        report(f"bubble {bid}", sub)


main()
