#!/usr/bin/env python3
"""Call SNVs and SVs from a graph's haplotypes, to test beside panvar's records.

Each haplotype path is aligned to the reference path (minimap2 asm5); SNVs come from paftools.js call and
SVs from svim-asm (haploid), merged across haplotypes by bcftools merge and truvari collapse (default
thresholds, no size limit).

A haplotype reads 1 where it carries the allele; NA where it is not aligned over the site, or aligned there
more than once; else 0 when it has a base aligned at a SNV, or bases aligned to at least half of an SV's
reference span (both bases around an insertion). A haplotype missing most of a deletion's bases carries a
deletion called elsewhere, so it reads NA.

Writes in --out-dir, for snv and svim, per haplotype and, with --samples, per sample (the sum of its
haplotypes, NA if any is NA):
  haplotype/<layer>/bimbam_<layer>.bimbam.gz, samples.txt.gz, feature_annot.<layer>.tsv.gz
  sample/<layer>/...
AF and AN are over each level's own haplotypes, as describe writes them. SNV ids are
snv_<chrom>_<pos>_<ref>_<alt>, on the chromosome when the reference path's name carries its range.

  assembly_variants.py --gfa graph.gfa --reference GRCh38 --tools <dir> --merge-tools <dir> \\
                       --out-dir out [--samples samples.tsv] [--threads N] [--min-sv-bp 50]
"""
import argparse
import bisect
import collections
import concurrent.futures
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def w_path_name(sample, hap, seqid, start, end):
    """A W line's path name as panvar gives it: sample#hap#seqid:start-end, the start made 1-based; no
    interval without a numeric start, and an end that is not a number written '*'."""
    numeric = lambda t: re.fullmatch(r"[0-9]{1,18}", t) is not None
    if not numeric(start):
        return "%s#%s#%s" % (sample, hap, seqid)
    return "%s#%s#%s:%d-%s" % (sample, hap, seqid, int(start) + 1, end if numeric(end) else "*")


def load_gfa(path):
    """Segment sequences, every P and W path as a list of (node, is_reverse), and the names of the W paths."""
    seq, paths, walks = {}, {}, set()
    with _open(path) as fh:
        for line in fh:
            if line.startswith("S\t"):
                f = line.rstrip("\n").split("\t")
                seq[f[1]] = f[2]
            elif line.startswith("P\t"):
                f = line.rstrip("\n").split("\t")
                paths[f[1]] = [(s[:-1], s[-1] == "-") for s in f[2].split(",")]
            elif line.startswith("W\t"):
                f = line.rstrip("\n").split("\t")
                name = w_path_name(*f[1:6])
                paths[name] = [(n, o == "<") for o, n in re.findall(r"([<>])([^<>]+)", f[6])]
                walks.add(name)
    return seq, paths, walks


def spell(seq, steps):
    return "".join(seq[n].translate(_COMP)[::-1] if rev else seq[n] for n, rev in steps)


def info_of(field):
    d = {}
    for kv in field.split(";"):
        k, _, v = kv.partition("=")
        d[k] = v
    return d


def aligned_bases(tstart, tend, cs):
    """The reference intervals, half-open, whose bases a cs string aligns to a base of the query: matches
    and substitutions, not deletions or introns. Stops if the string does not end at the alignment's end."""
    out, t = [], tstart
    for op in re.findall(r"(:\d+|=[A-Za-z]+|\*[a-z]{2}|[+-][A-Za-z]+|~[a-z]{2}\d+[a-z]{2})", cs):
        if op[0] == ":":
            n = int(op[1:])
        elif op[0] == "=":
            n = len(op) - 1
        elif op[0] == "*":
            n = 1
        elif op[0] == "-":
            t += len(op) - 1
            continue
        elif op[0] == "~":
            t += int(op[3:-2])
            continue
        else:
            continue
        if out and out[-1][1] == t:
            out[-1][1] = t + n
        else:
            out.append([t, t + n])
        t += n
    if t != tend:
        sys.exit("cs string ends at %d, its alignment at %d" % (t, tend))
    return [tuple(x) for x in out]


class Intervals:
    """Reference intervals, half-open, with the number of alignments over each base."""

    def __init__(self, spans):
        events = collections.Counter()
        for s, e in spans:
            events[s] += 1
            events[e] -= 1
        self.bounds, self.depth = [], []
        d = 0
        for p in sorted(events):
            d += events[p]
            self.bounds.append(p)
            self.depth.append(d)

    def depth_at(self, pos):
        i = bisect.bisect_right(self.bounds, pos) - 1
        return self.depth[i] if i >= 0 else 0

    def covered(self, start, end):
        """Bases of [start, end) under at least one interval."""
        total = 0
        i = max(0, bisect.bisect_right(self.bounds, start) - 1)
        while i < len(self.bounds) and self.bounds[i] < end:
            lo = max(start, self.bounds[i])
            hi = min(end, self.bounds[i + 1]) if i + 1 < len(self.bounds) else end
            if self.depth[i] > 0 and hi > lo:
                total += hi - lo
            i += 1
        return total

    def single_over(self, start, end):
        """True when [start, end) lies under exactly one alignment throughout."""
        i = bisect.bisect_right(self.bounds, start) - 1
        if i < 0 or self.depth[i] != 1:
            return False
        j = i + 1
        while j < len(self.bounds) and self.bounds[j] < end:
            if self.depth[j] != 1:
                return False
            j += 1
        return True


def call_haplotype(name, label, fasta_text, ref_fa, tools, work_root, min_sv_bp):
    """SNVs as {(pos0, ref, alt)}, svim-asm's VCF lines (sample `label`), the alignments' reference spans
    and the reference bases they align."""
    env = dict(os.environ, PATH=tools + os.pathsep + os.environ.get("PATH", ""))
    work = tempfile.mkdtemp(prefix="hap_", dir=work_root)
    try:
        hap_fa = os.path.join(work, "hap.fa")
        with open(hap_fa, "w") as fh:
            fh.write(fasta_text)
        mm = [os.path.join(tools, "minimap2"), "-x", "asm5", "--cs", "-r2k", "-t", "1"]
        paf = subprocess.run(mm + ["-c", ref_fa, hap_fa], check=True, capture_output=True, text=True,
                             env=env).stdout
        # primary alignments of 1 kb or more: a secondary places the same sequence again, a shorter one is unreliable
        spans, bases, rows = [], [], []
        for line in paf.splitlines():
            f = line.split("\t")
            if int(f[10]) < 1000 or "tp:A:P" not in f[12:]:
                continue
            cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
            spans.append((int(f[7]), int(f[8])))
            bases.extend(aligned_bases(int(f[7]), int(f[8]), cs) if cs is not None else [(int(f[7]), int(f[8]))])
            rows.append(line)
        cover = (Intervals(spans), Intervals(bases))
        snvs = set()
        if rows:
            rows.sort(key=lambda r: (r.split("\t")[5], int(r.split("\t")[7])))
            sorted_paf = os.path.join(work, "aln.paf")
            with open(sorted_paf, "w") as fh:
                fh.write("\n".join(rows) + "\n")
            out = subprocess.run([os.path.join(tools, "k8"), os.path.join(tools, "paftools.js"), "call",
                                  "-l", "1000", "-L", "1000", sorted_paf],
                                 check=True, capture_output=True, text=True, env=env).stdout
            for line in out.splitlines():
                f = line.split("\t")
                # V chrom start end depth mapq ref alt ...: a single-base change under one alignment.
                if f[0] != "V" or f[4] != "1" or len(f[6]) != 1 or len(f[7]) != 1 or f[6] == "-" or f[7] == "-":
                    continue
                snvs.add((int(f[2]), f[6].upper(), f[7].upper()))
        sam = subprocess.run(mm + ["-a", ref_fa, hap_fa], check=True, capture_output=True, env=env).stdout
        with open(os.path.join(work, "aln.sam"), "wb") as fh:
            fh.write(sam)
        bam = os.path.join(work, "aln.bam")
        subprocess.run([os.path.join(tools, "samtools"), "sort", "-o", bam, os.path.join(work, "aln.sam")],
                       check=True, capture_output=True, env=env)
        subprocess.run([os.path.join(tools, "samtools"), "index", bam], check=True, capture_output=True, env=env)
        subprocess.run([os.path.join(tools, "svim-asm"), "haploid", "--min_sv_size", str(min_sv_bp),
                        "--sample", label, work, bam, ref_fa], check=True, capture_output=True, env=env)
        vcf_lines = []
        vcf = os.path.join(work, "variants.vcf")
        if os.path.exists(vcf):
            for line in open(vcf):
                if not line.startswith("#"):
                    t = info_of(line.split("\t")[7]).get("SVTYPE", "").split(":")[0]
                    if t not in ("DEL", "INS", "DUP", "INV"):
                        continue
                vcf_lines.append(line)
        return name, snvs, vcf_lines, cover
    finally:
        shutil.rmtree(work, ignore_errors=True)


def merge_svs(vcfs, labels, ref_fa, merge_tools, work, min_sv_bp):
    """One record per SV across haplotypes: the haplotypes' svim-asm VCFs merged with bcftools merge,
    then truvari collapse with its default thresholds (same type, within 500 bp, 95% size and sequence
    similarity, repeats unrolled) but no upper size limit: the default skips SVs over 50 kb, which would
    leave each long call alone. Calls of one event placed at different breakpoints, or of different
    sizes, stay separate records. Returns [(type, pos0, length, {label})]."""
    tool = lambda name: os.path.join(merge_tools, name)
    files = []
    for label, lines in zip(labels, vcfs):
        if not any(not l.startswith("#") for l in lines):
            continue
        path = os.path.join(work, label + ".vcf")
        with open(path, "w") as fh:
            fh.writelines(lines)
        subprocess.run([tool("bcftools"), "sort", "-Oz", "-o", path + ".gz", path], check=True, capture_output=True)
        subprocess.run([tool("tabix"), "-p", "vcf", path + ".gz"], check=True, capture_output=True)
        files.append(path + ".gz")
    if not files:
        return []
    merged = os.path.join(work, "merged.vcf.gz")
    if len(files) == 1:
        shutil.copy(files[0], merged)
    else:
        listing = os.path.join(work, "files.txt")
        with open(listing, "w") as fh:
            fh.write("\n".join(files) + "\n")
        subprocess.run([tool("bcftools"), "merge", "-m", "none", "-l", listing, "-Oz", "-o", merged],
                       check=True, capture_output=True)
    subprocess.run([tool("tabix"), "-f", "-p", "vcf", merged], check=True, capture_output=True)
    collapsed = os.path.join(work, "collapsed.vcf")
    subprocess.run([tool("truvari"), "collapse", "-i", merged, "-o", collapsed, "-c", os.path.join(work, "removed.vcf"),
                    "-f", ref_fa, "-s", str(min_sv_bp), "-S", "-1"], check=True, capture_output=True)
    records, samples = [], []
    for line in open(collapsed):
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            samples = f[9:]
            continue
        inf = info_of(f[7])
        t = inf.get("SVTYPE", "").split(":")[0]
        pos0 = int(f[1]) - 1
        end0 = int(inf["END"]) if "END" in inf else pos0 + 1
        length = abs(int(inf["SVLEN"])) if "SVLEN" in inf else max(end0 - pos0, abs(len(f[4]) - len(f[3])))
        if length < min_sv_bp:
            continue
        gi = f[8].split(":").index("GT")
        carriers = {smp for smp, val in zip(samples, f[9:]) if "1" in re.split(r"[/|]", val.split(":")[gi])}
        if carriers:
            records.append((t, pos0, length, carriers))
    return records


def write_layer(root, layer, ids, annot_rows, columns, dosage):
    """dosage(row_index, column) -> value or None (NA)."""
    d = os.path.join(root, layer)
    os.makedirs(d, exist_ok=True)
    with gzip.open(os.path.join(d, "samples.txt.gz"), "wt") as fh:
        fh.write("".join(c + "\n" for c in columns))
    with gzip.open(os.path.join(d, "bimbam_%s.bimbam.gz" % layer), "wt") as fh:
        for i, rid in enumerate(ids):
            vals = (dosage(i, c) for c in columns)
            fh.write(rid + ", A, B, " + ", ".join("NA" if v is None else str(v) for v in vals) + "\n")
    with gzip.open(os.path.join(d, "feature_annot.%s.tsv.gz" % layer), "wt") as fh:
        fh.write("feature_id\tlayer\tencoding\tbubbles\tnodes\tsvtype\tgene\tAF\tAN\tchrom\tpos\tend\n")
        for row in annot_rows:
            fh.write("\t".join(row) + "\n")


def main():
    ap = argparse.ArgumentParser(description="SNVs and SVs called from a graph's haplotypes, as BIMBAM matrices.",
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=32, width=100))
    ap.add_argument("--gfa", required=True, metavar="GFA", help="the graph")
    ap.add_argument("--reference", required=True, metavar="NAME", help="reference path, or part of its name")
    ap.add_argument("--tools", required=True, metavar="DIR", help="minimap2, k8, paftools.js, samtools, svim-asm")
    ap.add_argument("--merge-tools", required=True, metavar="DIR", help="bcftools, tabix, truvari")
    ap.add_argument("--out-dir", required=True, metavar="DIR", help="output folder")
    ap.add_argument("--samples", metavar="TSV", help="sample, hap1, hap2: also per-sample matrices")
    ap.add_argument("--sv-tables", metavar="PREFIX", help="also svim-asm's SVs per haplotype and merged, as tables")
    ap.add_argument("--threads", type=int, default=2, metavar="N", help="threads (default 2)")
    ap.add_argument("--min-sv-bp", type=int, default=50, metavar="N", help="smallest SV (default 50)")
    a = ap.parse_args()
    for tool in ("minimap2", "k8", "paftools.js", "samtools", "svim-asm"):
        if not os.path.exists(os.path.join(a.tools, tool)):
            sys.exit("no %s in --tools %s" % (tool, a.tools))
    for tool in ("bcftools", "tabix", "truvari"):
        if not os.path.exists(os.path.join(a.merge_tools, tool)):
            sys.exit("no %s in --merge-tools %s" % (tool, a.merge_tools))

    seq, paths, walks = load_gfa(a.gfa)
    # as panvar resolves -r: the exact name, else the paths containing it, else the same ignoring case
    refs = ([p for p in paths if p == a.reference] or [p for p in paths if a.reference in p]
            or [p for p in paths if a.reference.lower() in p.lower()])
    if len(refs) != 1:
        sys.exit("--reference %s matches %d paths" % (a.reference, len(refs)))
    ref = refs[0]
    # Positions from the reference path's chrom:start-end, as panvar reads it: 1-based inclusive when the
    # interval is the path's length, else half-open; a W path's name is always 1-based.
    m = re.search(r"(?:^|#)([^#:]+):(\d+)-(\d+|\*)$", ref)
    if m and m.group(3) == "*" and ref not in walks:
        m = None
    chrom, offset = (m.group(1), int(m.group(2))) if m else (ref, 0)
    if m and offset > 0:
        if ref in walks:
            offset -= 1
        elif int(m.group(3)) - offset + 1 == sum(len(seq[node]) for node, _ in paths[ref]):
            offset -= 1

    os.makedirs(a.out_dir, exist_ok=True)
    work_root = tempfile.mkdtemp(prefix="assembly_variants_", dir=a.out_dir)
    try:
        ref_seq = spell(seq, paths[ref])
        ref_fa = os.path.join(work_root, "reference.fa")
        with open(ref_fa, "w") as fh:
            fh.write(">reference\n" + ref_seq + "\n")
        subprocess.run([os.path.join(a.tools, "samtools"), "faidx", ref_fa], check=True, capture_output=True)

        haps = [p for p in paths if p != ref]
        # Sample names in the VCFs: the path name made safe, one per haplotype.
        label = {h: "h%d_%s" % (i, re.sub(r"[^A-Za-z0-9._-]", "_", h)) for i, h in enumerate(haps)}
        results = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, a.threads)) as pool:
            futs = [pool.submit(call_haplotype, h, label[h], ">%s\n%s\n" % (label[h], spell(seq, paths[h])),
                                ref_fa, a.tools, work_root, a.min_sv_bp) for h in haps]
            for n, fut in enumerate(concurrent.futures.as_completed(futs), 1):
                name, snvs, vcf_lines, cover = fut.result()
                results[name] = (snvs, vcf_lines, cover)
                sys.stderr.write("\r  aligned %d/%d haplotypes" % (n, len(haps)))
                sys.stderr.flush()
        sys.stderr.write("\n")
        hap_of = {label[h]: h for h in haps}
        if a.sv_tables:
            # each haplotype's calls before merging: VCF POS, type without subtype, |SVLEN| or END - POS
            with open(a.sv_tables + ".haplotype_svs.tsv", "w") as fh:
                fh.write("haplotype\ttype\tpos\tend\tlength\n")
                for h in haps:
                    for line in results[h][1]:
                        if line.startswith("#"):
                            continue
                        f = line.rstrip("\n").split("\t")
                        inf = info_of(f[7])
                        t = inf.get("SVTYPE", "").split(":")[0]
                        pos = int(f[1])
                        end = int(inf["END"]) if "END" in inf else pos
                        length = abs(int(inf["SVLEN"])) if "SVLEN" in inf else max(end - pos, abs(len(f[4]) - len(f[3])))
                        if length >= a.min_sv_bp:
                            fh.write("%s\t%s\t%d\t%d\t%d\n" % (h, t, offset + pos, offset + end, length))
        svrec = [(t, p, ln, {hap_of[c] for c in carriers})
                 for t, p, ln, carriers in merge_svs([results[h][1] for h in haps], [label[h] for h in haps],
                                                     ref_fa, a.merge_tools, work_root, a.min_sv_bp)]
    finally:
        shutil.rmtree(work_root, ignore_errors=True)

    columns = list(paths)   # every path, the reference included, as describe's haplotype matrices
    # SNVs: one row per (position, alt).
    sites = sorted({s for snvs, _, _ in results.values() for s in snvs})
    snv_ids = ["snv_%s_%d_%s_%s" % (chrom, offset + p + 1, r, alt) for p, r, alt in sites]

    def snv_dose(i, hap):
        if hap == ref:
            return 0
        snvs, _, (spans, bases) = results[hap]
        p, r, alt = sites[i]
        if (p, r, alt) in snvs:
            return 1
        return 0 if spans.depth_at(p) == 1 and bases.depth_at(p) == 1 else None

    # SVs: merged across haplotypes; a haplotype without the SV must align at least half of its span.
    sv_ids = ["svim_%s_%s_%d_%d" % (t, chrom, offset + p + 1, ln) for t, p, ln, _ in svrec]
    seen = collections.Counter()
    for i, rid in enumerate(sv_ids):   # two records at one position and size: kept apart by a suffix
        seen[rid] += 1
        if seen[rid] > 1:
            sv_ids[i] = "%s_%d" % (rid, seen[rid])
    if a.sv_tables:
        with open(a.sv_tables + ".merged_svs.tsv", "w") as fh:
            fh.write("id\ttype\tpos\tend\tlength\tcarriers\n")
            for rid, (t, p, ln, carriers) in zip(sv_ids, svrec):
                end = offset + p + 1 + (0 if t == "INS" else ln if t == "DEL" else ln - 1)
                fh.write("%s\t%s\t%d\t%d\t%d\t%s\n" % (rid, t, offset + p + 1, end, ln, ",".join(sorted(carriers))))

    def sv_dose(i, hap):
        if hap == ref:
            return 0
        t, p, ln, carriers = svrec[i]
        if hap in carriers:
            return 1
        spans, bases = results[hap][2]
        if t == "INS":   # between the anchor base p and the next one, absent at the path's end
            flanks = [b for b in (p, p + 1) if b < len(ref_seq)]
            if not spans.single_over(p, flanks[-1] + 1):
                return None
            return 0 if all(bases.depth_at(b) for b in flanks) else None
        # A deletion's VCF position is the base before it; other types start at their position.
        lo = p + 1 if t == "DEL" else p
        hi = lo + ln
        if not spans.single_over(max(0, lo - 1), min(len(ref_seq), hi + 1)):
            return None
        return 0 if 2 * bases.covered(lo, hi) >= hi - lo else None

    # Reference span: a SNV its base, an insertion its anchor base, a deletion POS to its last base, an
    # inversion or duplication its own bases.
    snv_span = ["%s\t%d\t%d" % (chrom, offset + p + 1, offset + p + 1) for p, _, _ in sites]
    sv_span = ["%s\t%d\t%d" % (chrom, offset + p + 1, offset + p + 1 + (0 if t == "INS" else ln if t == "DEL" else ln - 1))
               for t, p, ln, _ in svrec]

    def annot(layer, ids, types, cols, dose, alleles=lambda c: 1):
        """AF and AN over the columns written: AN counts `alleles(c)` haplotypes per observed column, AF is the
        dosage summed over them divided by AN, as describe computes them at each level."""
        rows = []
        for i, rid in enumerate(ids):
            an = dosage = 0
            for c in cols:
                v = dose(i, c)
                if v is not None:
                    an += alleles(c)
                    dosage += v
            # no observed column leaves AF undefined: '.', as describe writes it
            af = "%.6g" % (dosage / an) if an else "."
            span = (snv_span if layer == "snv" else sv_span)[i]
            rows.append([rid, layer, "dosage", ".", ".", types[i], ".", af, str(an), span])
        return rows

    sv_types = [t for t, _, _, _ in svrec]
    snv_annot = annot("snv", snv_ids, ["SNV"] * len(snv_ids), columns, snv_dose)
    sv_annot = annot("svim", sv_ids, sv_types, columns, sv_dose)

    hap_root = os.path.join(a.out_dir, "haplotype")
    write_layer(hap_root, "snv", snv_ids, snv_annot, columns, snv_dose)
    write_layer(hap_root, "svim", sv_ids, sv_annot, columns, sv_dose)
    print("haplotypes: %d; SNVs: %d; svim-asm SVs: %d" % (len(columns), len(snv_ids), len(sv_ids)))

    if a.samples:
        # Read as panvar's read_cosigt_table does: optional header, comma-separated contigs per haplotype,
        # empty parts skipped, a repeated contig counted once, every sample with at least one haplotype.
        table, seen, first = [], set(), True
        # lines split on "\n" only, one trailing "\r" removed, as panvar does
        with open(a.samples, "rb") as fh:
            raw = fh.read()
        if raw[:2] == b"\x1f\x8b":
            sys.exit("--samples %s is gzip-compressed; give it uncompressed" % a.samples)
        for line in raw.decode("utf-8", "surrogateescape").split("\n"):
            if line.endswith("\r"):
                line = line[:-1]
            if not line:
                continue
            f = line.split("\t")
            if first:
                first = False
                if f[0].lower() in ("sample", "sample_id", "sample.id") or f[0].startswith("#"):
                    continue
            if not f[0]:
                sys.exit("--samples has a row with an empty sample id")
            if f[0] in seen:
                sys.exit("--samples lists sample '%s' twice" % f[0])
            seen.add(f[0])
            haps = []
            for cell in f[1:]:
                parts = []
                for part in cell.split(","):
                    if part and part not in parts:
                        parts.append(part)
                if parts:
                    haps.append(parts)
            if not haps:
                sys.exit("--samples: sample '%s' has no haplotype" % f[0])
            table.append((f[0], haps))
        if not table:
            sys.exit("--samples lists no sample")
        names = [s for s, _ in table]
        by_name = dict(table)
        for s, haps in table:
            for hap in haps:
                for part in hap:
                    if part not in paths:
                        sys.exit("--samples names '%s' (sample %s), which is not a path of the graph" % (part, s))

        def summed(dose):
            def d(i, sample):
                total = 0
                for hap in by_name[sample]:
                    # a haplotype in several contigs is known where exactly one is observed (two would
                    # count the site twice), as describe and associate read it
                    obs = [v for v in (dose(i, part) for part in hap) if v is not None]
                    if len(obs) != 1:
                        return None
                    total += obs[0]
                return total
            return d

        # AF and AN over the samples' haplotypes, not the panel's: the sample level's own denominator.
        smp_root = os.path.join(a.out_dir, "sample")
        slots = lambda sample: len(by_name[sample])
        for layer, ids, types, dose in (("snv", snv_ids, ["SNV"] * len(snv_ids), summed(snv_dose)),
                                        ("svim", sv_ids, sv_types, summed(sv_dose))):
            write_layer(smp_root, layer, ids, annot(layer, ids, types, names, dose, slots), names, dose)
        print("samples: %d" % len(names))


if __name__ == "__main__":
    main()
