#!/usr/bin/env python3
"""Compare panvar's per-haplotype CN/CNBP against svim-asm, an assembly-based SV caller.

panvar reads structure off the graph; svim-asm reads it off a minimap2 alignment of one assembly to
another. They share no code and no assumptions, so agreement is evidence and disagreement is a lead.

The comparison is per HAPLOTYPE, against that haplotype's own CNBP -- not against the region VCF's
representative SVLEN, which is one merged event's size and is not what any particular carrier holds.

Coordinates stay in the reference PATH's own frame throughout: the FASTA handed to minimap2 is the
spelled reference path, so svim-asm's positions and the bubble spans computed here are already the same
frame and nothing has to be projected.

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
    ap.add_argument("--bin-dir", default=os.path.expanduser("~/miniconda3/envs/svimasm/bin"),
                    help="directory holding minimap2, samtools and svim-asm")
    ap.add_argument("--margin", type=int, default=2000,
                    help="bp of slack when attributing an svim-asm call to a bubble span; breakpoints "
                         "from an alignment do not land where graph boundaries do")
    ap.add_argument("--min-size-bp", type=int, default=1000,
                    help="relative agreement (within 5%%/10%%) is reported only for haplotypes whose "
                         "CNBP reaches this; a 40 bp disagreement on a 50 bp delta is not a 5%% question")
    ap.add_argument("--keep", action="store_true", help="keep per-haplotype BAM/VCF working files")
    args = ap.parse_args()

    env_bin = args.bin_dir
    for tool in ("minimap2", "samtools", "svim-asm"):
        if not os.path.exists(os.path.join(env_bin, tool)):
            sys.exit(f"compare_svim_asm: {tool} not found in {env_bin}")

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

    # every haplotype that has a CNBP on some DUP record
    wanted = set()
    for r in dup_rows:
        keys = r[8].split(":")
        if "CNBP" not in keys:
            continue
        bi = keys.index("CNBP")
        for i in range(9, len(hdr)):
            v = r[i].split(":")
            if len(v) > bi and v[bi] not in (".", ""):
                wanted.add(hdr[i])
    wanted.discard(args.reference)
    wanted = sorted(wanted)
    if args.samples and len(wanted) > args.samples:
        random.Random(args.seed).shuffle(wanted)
        wanted = sorted(wanted[:args.samples])

    os.makedirs(args.out_dir, exist_ok=True)
    ref_fa = os.path.join(args.out_dir, "reference.fa")
    with open(ref_fa, "w") as fh:
        fh.write(">reference\n" + spell(seq, paths[args.reference]) + "\n")
    run([os.path.join(env_bin, "samtools"), "faidx", ref_fa])

    records = []
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
            best, best_d = None, None
            for bid, (lo, hi) in span.items():
                d = 0 if lo <= p <= hi else min(abs(p - lo), abs(p - hi))
                if best_d is None or d < best_d:
                    best, best_d = bid, d
            if best is not None and best_d <= args.margin:
                attributed.setdefault(best, []).append(l)

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
            svim_bp = sum(attributed.get(bid, []))
            records.append({
                "bubble": bid, "variant": r[2], "hap": hap,
                "cn": v[ci], "cnbp": int(v[bi]), "svim_bp": svim_bp,
                "n_records_in_bubble": len(per_bubble.get(bid, [])),
                "class": "pure_module" if len(per_bubble.get(bid, [])) == 1 else "compound_site",
            })
        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)
    sys.stderr.write("\n")

    out_tsv = os.path.join(args.out_dir, "svim_asm_comparison.tsv")
    with open(out_tsv, "w") as fh:
        fh.write("bubble_id\tvariant_id\thaplotype\tcn\tcnbp\tsvim_asm_bp\tdiff_bp\t"
                 "records_in_bubble\tclass\n")
        for r in records:
            fh.write(f"{r['bubble']}\t{r['variant']}\t{r['hap']}\t{r['cn']}\t{r['cnbp']}\t"
                     f"{r['svim_bp']}\t{r['svim_bp'] - r['cnbp']}\t{r['n_records_in_bubble']}\t"
                     f"{r['class']}\n")

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
    report("all", records)
    for cls in ("pure_module", "compound_site"):
        report(cls, [r for r in records if r["class"] == cls])
    for bid in sorted({r["bubble"] for r in records}, key=lambda x: int(x)):
        sub = [r for r in records if r["bubble"] == bid]
        report(f"bubble {bid}", sub)


main()
