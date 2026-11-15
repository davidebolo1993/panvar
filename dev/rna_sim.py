#!/usr/bin/env python3
"""Sources and reads for the RNA recruitment benchmark, one locus at a time.

    rna_sim.py only --gfa <locus.gfa[.gz]> --sample <a,b,...> --keep <ref,...> --out <only.gfa>
        the graph with only the held-out samples' paths and the reference paths (the truth is projected
        onto it, which is faster than onto the whole graph and gives the same transcripts)
    rna_sim.py transcriptome --out <tx.fa.gz> <cdna.fa.gz> [<ncrna.fa.gz> ...]
        Ensembl transcript FASTAs cut to the primary chromosomes: transcripts on patches and alternate
        haplotypes carry other gene ids, so a locus gene's copy there would count as a paralog
    rna_sim.py simulate --gfa <only.gfa> --hst <truth prefix> --samples <a,b> --tx <tx.fa.gz>
                        --interval chr:start-end --out <dir> [options]
        expressed transcripts of the held-out haplotypes, with perturbations the annotation lacks, and
        negatives from elsewhere in the transcriptome; reads of three classes from both

Sources (the truth table's `source`):
    locus       the held-out haplotypes' transcripts (truth HSTs), and perturbed molecules made from their
                walks: `skip` (an internal exon skipped, at a junction no transcript of the haplotype has),
                `retained_intron`, `pre_mrna` (the whole unspliced span), `alt_end` (the first or last exon
                extended into its flank)
    paralog     reference transcripts outside the locus whose gene symbol matches --paralog-regex, or that
                share 31-mers with the locus transcripts
    pseudogene  processed pseudogenes sharing 31-mers with the locus transcripts
    transcriptome  a sample of all other reference transcripts, with a log-normal expression
Classes: short (2x150 pairs, 0.2% substitutions), kinnex (whole molecules, 0.1% errors), ont (pbsim3's
ONT model on molecules cut at their 5' end).
"""
import argparse
import gzip
import math
import os
import random
import re
import subprocess
import sys

K = 31
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


def opener(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def fasta(path):
    name, seq = None, []
    with opener(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def sample_of(path):
    return path.split("#")[0]


# ---------------------------------------------------------------- only / transcriptome

def cmd_only(a):
    keep = set(a.sample.split(",")) | set(a.keep.split(","))
    n = 0
    with opener(a.gfa) as f, open(a.out, "w") as out:
        for line in f:
            if line[0] in "PW":
                fld = line.split("\t", 3)
                s = fld[1] if line[0] == "W" else sample_of(fld[1])
                if s not in keep:
                    continue
                n += 1
            out.write(line)
    print(f"kept {n} paths", file=sys.stderr)


def tx_header(h):
    """(id, chrom, start, end, gene id, symbol, biotype) of an Ensembl cDNA/ncRNA header."""
    f = h.split()
    loc = f[2].split(":")
    kv = dict(x.split(":", 1) for x in f[3:] if ":" in x and not x.startswith("description"))
    return (f[0], loc[0], loc[2], int(loc[3]), int(loc[4]), kv.get("gene", "").split(".")[0],
            kv.get("gene_symbol", ""), kv.get("transcript_biotype", ""))


def cmd_transcriptome(a):
    n = kept = 0
    with gzip.open(a.out, "wt", compresslevel=1) as out:
        for path in a.inputs:
            for h, s in fasta(path):
                n += 1
                if tx_header(h)[1] != "chromosome":
                    continue
                kept += 1
                out.write(f">{h}\n{s}\n")
    print(f"{kept} of {n} transcripts on primary chromosomes", file=sys.stderr)


# ---------------------------------------------------------------- truth transcripts and walks

def spell_paths(gfa, samples):
    seq, walks = {}, {}
    with opener(gfa) as f:
        for line in f:
            if line[0] == "S":
                fld = line.rstrip("\n").split("\t")
                seq[fld[1]] = fld[2]
            elif line[0] == "P":
                fld = line.rstrip("\n").split("\t")
                if sample_of(fld[1]) in samples:
                    walks[fld[1]] = [(s[:-1], s[-1] == "-") for s in fld[2].split(",")]
    return {name: "".join(rc(seq[n]) if r else seq[n] for n, r in steps) for name, steps in walks.items()}


def read_truth(prefix, samples, paths):
    """Expressible HSTs of the held-out haplotypes: rows with a copy, their exons on the walk as
    (start, end) in transcript order, and the direction the transcript runs along the walk."""
    seqs = {}
    for h, s in fasta(prefix + ".hst.fa"):
        seqs[h.split()[0]] = s
    rows = []
    bad = 0
    with open(prefix + ".hst.tsv") as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["copy"] == "0" or sample_of(r["haplotype"]) not in samples or r["source"] == "annotated":
                continue
            name = f"{r['transcript_id']}|{r['gene_name']}|{r['haplotype']}|{r['copy']}"
            if name not in seqs:
                continue
            starts = [int(x) for x in r["exon_starts"].split(",")]
            lens = [int(x) for x in r["exon_lengths"].split(",")]
            exons = [(s, s + n) for s, n in zip(starts, lens)]
            walk = paths[r["haplotype"]]
            fwd = "".join(walk[s:e] for s, e in exons)
            rev = "".join(rc(walk[s:e]) for s, e in exons)
            if fwd == seqs[name]:
                d = 1
            elif rev == seqs[name]:
                d = -1
            else:
                bad += 1
                continue
            rows.append(dict(name=name, tid=r["transcript_id"], gid=r["gene_id"], gene=r["gene_name"],
                             hap=r["haplotype"], exons=exons, dir=d, seq=seqs[name], status=r["status"]))
    print(f"{len(rows)} truth HSTs, {bad} not rebuilt from their walk", file=sys.stderr)
    return rows


def molecule(walk, exons, d):
    """A molecule's sequence from walk intervals in transcript order, and its junctions (molecule
    coordinates where one interval ends and the next begins)."""
    parts = [walk[s:e] if d > 0 else rc(walk[s:e]) for s, e in exons]
    junctions, pos = [], 0
    for p in parts[:-1]:
        pos += len(p)
        junctions.append(pos)
    return "".join(parts), junctions


def perturb(kind, row, walk, rng, known):
    """A molecule the annotation lacks, made from the HST's walk, or None where it does not apply."""
    ex, d = list(row["exons"]), row["dir"]
    n = len(ex)
    if kind == "skip":
        if n < 3:
            return None
        for i in rng.sample(range(1, n - 1), n - 2):
            seq, junc = molecule(walk, ex[:i] + ex[i + 1:], d)
            j = junc[i - 1]
            if seq[max(0, j - 20):j + 20] not in known:   # a junction no expressed transcript has
                return seq, junc, f"e{i + 1}"
        return None
    if kind == "retained_intron":
        if n < 2:
            return None
        i = rng.randrange(n - 1)
        a, b = ex[i], ex[i + 1]
        merged = (min(a[0], b[0]), max(a[1], b[1]))
        if merged[1] - merged[0] > 20000:
            return None
        seq, junc = molecule(walk, ex[:i] + [merged] + ex[i + 2:], d)
        return seq, junc, f"i{i + 1}"
    if kind == "pre_mrna":
        if n < 2:
            return None
        lo, hi = min(s for s, _ in ex), max(e for _, e in ex)
        if hi - lo > 60000:
            return None
        return molecule(walk, [(lo, hi)], d)[0], [], "span"
    if kind == "alt_end":
        five = rng.random() < 0.5
        x = rng.randint(100, 500) if five else rng.randint(200, 1000)
        i = 0 if five else n - 1
        s, e = ex[i]
        # the 5' end lies upstream along the transcript: before the exon on a forward walk
        grow_left = (five and d > 0) or (not five and d < 0)
        s, e = (max(0, s - x), e) if grow_left else (s, min(len(walk), e + x))
        ex[i] = (s, e)
        seq, junc = molecule(walk, ex, d)
        return seq, junc, "5p" if five else "3p"
    raise ValueError(kind)


def kmers(s, step=1):
    s = s.upper()
    out = set()
    for i in range(0, len(s) - K + 1, step):
        k = s[i:i + K]
        if "N" not in k:
            out.add(min(k, rc(k)))
    return out


# ---------------------------------------------------------------- errors and reads

def add_errors(s, rate, rng, sub=1.0):
    """Substitutions (share `sub`) and 1-base indels at `rate` per base."""
    out = []
    for c in s:
        if rng.random() < rate:
            r = rng.random()
            if r < sub:
                out.append(rng.choice([b for b in "ACGT" if b != c]))
            elif r < sub + (1 - sub) / 2:
                out.append(c + rng.choice("ACGT"))
            # else a deletion
        else:
            out.append(c)
    return "".join(out)


def crosses(junctions, a, b, flank):
    return any(a + flank <= j <= b - flank for j in junctions)


def path_contig(path):
    """(contig, 0-based start) of a PanSN path name `sample#hap#contig:start-end`."""
    c = path.split("#", 2)[2]
    contig, span = c.rsplit(":", 1)
    return contig, int(span.split("-")[0])


def read_cat(cat_dir, paths):
    """Transcripts of the haplotypes' own (CAT) annotations, cut from their walks by the annotation's
    coordinates: one record per transcript wholly on the walk."""
    out = []
    if not cat_dir or not os.path.isdir(cat_dir):
        return out
    files = os.listdir(cat_dir)
    for hap, walk in paths.items():
        s, h = hap.split("#")[:2]
        contig, off = path_contig(hap)
        pick = [f for f in files if f.startswith(f"{s}_{h}_{contig}_")]
        if not pick:
            continue
        tx = {}
        with open(os.path.join(cat_dir, pick[0])) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fld = line.rstrip("\n").split("\t")
                if len(fld) < 9 or fld[0] != contig or fld[2] not in ("transcript", "exon"):
                    continue
                att = dict(x.split("=", 1) for x in fld[8].split(";") if "=" in x)
                key = att.get("transcript_id") or att.get("Parent") or att.get("ID")
                if fld[2] == "transcript":
                    key = att.get("transcript_id") or att.get("ID")
                    t = tx.setdefault(key, {"exons": []})
                    t.update(strand=fld[6], gene=att.get("source_gene_common_name") or att.get("gene_name", "."),
                             tid=att.get("source_transcript", key).split(".")[0],
                             cls=att.get("transcript_class", "."))
                else:
                    tx.setdefault(key, {"exons": []})["exons"].append((int(fld[3]) - 1 - off, int(fld[4]) - off))
        for key, t in tx.items():
            ex = sorted(t["exons"])
            if not ex or "strand" not in t or ex[0][0] < 0 or ex[-1][1] > len(walk):
                continue
            d = 1 if t["strand"] == "+" else -1
            if d < 0:
                ex = ex[::-1]
            out.append(dict(name=key, tid=t["tid"], gene=t["gene"], hap=hap, exons=ex, dir=d, cls=t["cls"]))
    return out


def write_fq(f, name, seq, q):
    f.write(f"@{name}\n{seq}\n+\n{q * len(seq)}\n")


def cmd_simulate(a):
    samples = set(a.samples.split(","))
    paths = spell_paths(a.gfa, samples)
    truth = read_truth(a.hst, samples, paths)
    if not truth:
        sys.exit("rna_sim: no truth transcripts")
    known = set()
    for r in truth:   # every junction of the haplotypes' transcripts, 20 bases on each side
        _, junc = molecule(paths[r["hap"]], r["exons"], r["dir"])
        for j in junc:
            known.add(r["seq"][max(0, j - 20):j + 20])
    cat = read_cat(a.cat_dir, paths)
    if cat:   # the annotation's transcripts against the projection, for the log
        proj = {(r["tid"].split(".")[0], r["hap"]): r["seq"] for r in truth}
        same = sum(1 for c in cat if proj.get((c["tid"], c["hap"])) == molecule(paths[c["hap"]], c["exons"], c["dir"])[0])
        both = sum(1 for c in cat if (c["tid"], c["hap"]) in proj)
        print(f"annotated: {len(cat)} transcripts on {len({c['hap'] for c in cat})} haplotypes; "
              f"{both} of a projected transcript, {same} with its sequence", file=sys.stderr)
    locus_genes = {r["gid"] for r in truth}
    locus_names = {r["gene"] for r in truth}

    # negatives from the reference transcriptome, outside the locus
    rng0 = random.Random(a.seed)
    chrom, span = a.interval.rsplit(":", 1)
    lo, hi = (int(x) for x in span.split("-"))
    chrom = chrom.removeprefix("chr")
    keep_chroms = set(c.removeprefix("chr") for c in a.chroms.split(",")) if a.chroms else None
    para_re = re.compile(a.paralog_regex) if a.paralog_regex else None
    locus_k = set()
    for r in truth:
        locus_k |= kmers(r["seq"])
    pool, para, pseudo = [], [], []
    for h, s in fasta(a.tx):
        tid, _, c, s0, e0, gid, sym, bt = tx_header(h)
        if gid in locus_genes or (c == chrom and s0 < hi + a.margin and e0 > lo - a.margin):
            continue
        if keep_chroms is not None and c not in keep_chroms:
            continue
        if len(s) < 200:
            continue
        rec = (tid, sym, s)
        is_pseudo = bt in ("processed_pseudogene", "transcribed_processed_pseudogene")
        if is_pseudo or (para_re and para_re.search(sym)):
            ks = kmers(s, 3)
            share = len(ks & locus_k) / max(1, len(ks))
            if is_pseudo and share >= a.pseudo_share:
                pseudo.append(rec)
                continue
            if not is_pseudo and (para_re and para_re.search(sym) or share >= a.pseudo_share):
                para.append(rec)
                continue
        pool.append(rec)
    sample = rng0.sample(pool, min(a.tx_transcripts, len(pool)))
    # other loci's haplotype-specific transcripts: a paralog's own alleles, which no reference holds
    hap_para = []
    for path in a.hap_negatives:
        for h, s in fasta(path):
            f = h.split()[0].split("|")
            if len(f) >= 3 and f[1] not in locus_names and len(s) >= 200:
                hap_para.append((f[0], f[1], s))
    print(f"negatives: {len(para)} paralog, {len(hap_para)} haplotype paralog, {len(pseudo)} pseudogene, "
          f"{len(sample)} of {len(pool)} transcriptome transcripts", file=sys.stderr)
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "negatives.tsv"), "w") as f:
        f.write("source\ttranscript\tgene\tlength\n")
        for src, recs in (("paralog", para), ("paralog_hap", hap_para), ("pseudogene", pseudo),
                          ("transcriptome", sample)):
            for tid, sym, s in recs:
                f.write(f"{src}\t{tid}\t{sym}\t{len(s)}\n")
    for k in range(1, a.seeds + 1):
        simulate_seed(a, k, paths, truth, known, cat, para, hap_para, pseudo, sample)


def simulate_seed(a, k, paths, truth, known, cat, para, hap_para, pseudo, sample):
    """One expression draw and its reads, in <out>/s<k>."""
    rng = random.Random(a.seed * 1000 + k)
    out_dir = os.path.join(a.out, f"s{k}")
    os.makedirs(out_dir, exist_ok=True)
    prefix = f"{a.prefix}s{k}_"

    # expression: genes on or off, transcripts log-normal, the same on every haplotype and copy
    genes = sorted({r["gene"] for r in truth} | {c["gene"] for c in cat})
    on = {g: rng.random() < a.gene_on for g in genes}
    if not any(on.values()):
        on[rng.choice(genes)] = True
    gexp = {g: math.exp(rng.gauss(0, 1.0)) for g in genes}
    texp = {}
    for r in truth:
        tid = r["tid"].split(".")[0]
        if tid not in texp:
            texp[tid] = math.exp(rng.gauss(0, a.sigma)) * gexp[r["gene"]] * on[r["gene"]]
    cexp = {}   # an annotated transcript takes its source transcript's expression, else its own draw
    for c in cat:
        if c["tid"] not in texp and c["tid"] not in cexp:
            cexp[c["tid"]] = math.exp(rng.gauss(0, a.sigma)) * gexp[c["gene"]] * on[c["gene"]]
    texp_all = {**cexp, **texp}
    expressed = [r for r in truth if texp[r["tid"].split(".")[0]] > 0]
    cat_on = [c for c in cat if texp_all[c["tid"]] > 0]
    sexp = [math.exp(rng.gauss(0, 2.0)) for _ in sample]
    with open(os.path.join(out_dir, "expression.tsv"), "w") as f:
        f.write("transcript\tgene\ton\texpression\n")
        first = {}
        for r in truth:
            first.setdefault(r["tid"].split(".")[0], r["gene"])
        for c in cat:
            first.setdefault(c["tid"], c["gene"])
        for tid in sorted(texp_all):
            f.write(f"{tid}\t{first[tid]}\t{int(on[first[tid]])}\t{texp_all[tid]:.4g}\n")

    kinds = ["annotated", "skip", "retained_intron", "pre_mrna", "alt_end"]
    share = [1 - 4 * a.perturbed] + [a.perturbed] * 4

    def locus_molecules(n, by_length):
        """n locus molecules: (seq, junctions, truth fields)."""
        out = []
        if not expressed:
            return out
        w = [texp[x["tid"].split(".")[0]] * (len(x["seq"]) if by_length else 1) for x in expressed]
        for kind in rng.choices(kinds, weights=share, k=n):
            for _ in range(200):
                r = rng.choices(expressed, weights=w)[0]
                walk = paths[r["hap"]]
                if kind == "annotated":
                    seq, junc = molecule(walk, r["exons"], r["dir"])
                    what = "."
                else:
                    m = perturb(kind, r, walk, rng, known)
                    if m is None:
                        continue
                    seq, junc, what = m
                out.append((seq, junc, ["locus", kind, what, r["gene"], r["tid"], r["hap"], r["status"]]))
                break
        return out

    def cat_molecules(n, by_length):
        if not cat_on or n == 0:
            return []
        mols = [(c, molecule(paths[c["hap"]], c["exons"], c["dir"])) for c in cat_on]
        w = [texp_all[c["tid"]] * (len(m[0]) if by_length else 1) for c, m in mols]
        return [(m[0], m[1], ["cat", "annotated", c["cls"], c["gene"], c["tid"], c["hap"], "."])
                for c, m in rng.choices(mols, weights=w, k=n)]

    def negative_molecules(src, recs, n, weights=None):
        if not recs or n == 0:
            return []
        picks = rng.choices(recs, weights=weights, k=n)
        return [(s, [], [src, ".", ".", sym, tid, ".", "."]) for tid, sym, s in picks]

    def all_molecules(n_locus, n_cat, by_length):
        mols = locus_molecules(n_locus, by_length)
        mols += cat_molecules(n_cat, by_length)
        mols += negative_molecules("paralog", para, a.n_paralog)
        mols += negative_molecules("paralog_hap", hap_para, a.n_paralog)
        mols += negative_molecules("pseudogene", pseudo, a.n_pseudo)
        w = [e * len(r[2]) if by_length else e for e, r in zip(sexp, sample)]
        mols += negative_molecules("transcriptome", sample, a.n_tx, w)
        return mols

    truth_f = open(os.path.join(out_dir, "truth.tsv"), "w")
    truth_f.write("read\tclass\tsource\tperturbation\tdetail\tgene\ttranscript\thaplotype\tstatus\t"
                  "molecule_length\tjunctions\tread_start\tread_end\tjunction_read\tfull_length\n")

    def record(name, cls, fields, mlen, junc, a0, b0, jr, full):
        truth_f.write("\t".join([name, cls] + fields + [str(mlen), str(len(junc)), str(a0), str(b0),
                                                        str(int(jr)), str(int(full))]) + "\n")

    # short pairs: one fragment per molecule, molecules drawn by expression times length
    if a.n_short:
        with gzip.open(os.path.join(out_dir, "short_1.fq.gz"), "wt", compresslevel=1) as f1, \
                gzip.open(os.path.join(out_dir, "short_2.fq.gz"), "wt", compresslevel=1) as f2:
            for i, (seq, junc, fields) in enumerate(all_molecules(a.n_short, a.n_cat_short, True)):
                L = len(seq)
                fl = min(L, max(120, int(rng.gauss(a.fragment, a.fragment_sd))))
                s0 = rng.randint(0, L - fl)
                frag = seq[s0:s0 + fl]
                if rng.random() < 0.5:
                    frag = rc(frag)
                rl = min(a.read_length, fl)
                m1 = add_errors(frag[:rl], a.short_error, rng)
                m2 = add_errors(rc(frag)[:rl], a.short_error, rng)
                name = f"{prefix}s{i:07d}"
                write_fq(f1, name, m1, "?")
                write_fq(f2, name, m2, "?")
                # a junction counts when either mate crosses it with `flank` bases on each side
                jr = (crosses(junc, s0, s0 + rl, a.flank) or crosses(junc, s0 + fl - rl, s0 + fl, a.flank))
                record(name, "short", fields, L, junc, s0, s0 + fl, jr, False)

    # Kinnex-like: whole molecules
    if a.n_long:
        with gzip.open(os.path.join(out_dir, "kinnex.fq.gz"), "wt", compresslevel=1) as f:
            for i, (seq, junc, fields) in enumerate(all_molecules(a.n_long, a.n_cat_long, False)):
                s = seq if rng.random() < 0.5 else rc(seq)
                name = f"{prefix}k{i:07d}"
                write_fq(f, name, add_errors(s, a.kinnex_error, rng, sub=0.5), "?")
                record(name, "kinnex", fields, len(seq), junc, 0, len(seq), crosses(junc, 0, len(seq), a.flank), True)

    # ONT cDNA-like: molecules cut at their 5' end, then pbsim3's ONT error model
    if a.n_ont:
        mols = all_molecules(a.n_ont, a.n_cat_long, False)
        tdir = os.path.join(out_dir, "ont_tmp")
        os.makedirs(tdir, exist_ok=True)
        cut = {}
        with open(os.path.join(tdir, "mol.txt"), "w") as f:
            for i, (seq, junc, fields) in enumerate(mols):
                L = len(seq)
                s0 = 0
                if rng.random() > a.ont_full:
                    s0 = L - max(min(L, 200), int(L * rng.uniform(0.3, 0.95)))
                cut[f"m{i}"] = s0
                sense = rng.random() < 0.5
                f.write(f"m{i}\t{int(sense)}\t{int(not sense)}\t{seq[s0:]}\n")
        subprocess.run([a.pbsim, "--strategy", "trans", "--method", "errhmm", "--errhmm", a.pbsim_model,
                        "--transcript", os.path.join(tdir, "mol.txt"), "--accuracy-mean", str(a.ont_accuracy),
                        "--prefix", os.path.join(tdir, "out"), "--id-prefix", "o", "--seed", str(a.seed * 1000 + k)],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        spans = {}
        with gzip.open(os.path.join(tdir, "out.maf.gz"), "rt") as f:
            src = None
            for line in f:
                if line.startswith("s "):
                    fld = line.split()
                    if src is None:
                        src = (fld[1], int(fld[2]), int(fld[3]), fld[4], int(fld[5]))
                    else:
                        spans[fld[1]] = src
                        src = None
        n = 0
        with gzip.open(os.path.join(tdir, "out.fq.gz"), "rt") as fin, \
                gzip.open(os.path.join(out_dir, "ont.fq.gz"), "wt", compresslevel=1) as f:
            while True:
                h = fin.readline()
                if not h:
                    break
                s = fin.readline().rstrip("\n")
                fin.readline()
                fin.readline()
                rid = h[1:].split()[0]
                mid, start, size, strand, srclen = spans[rid]
                if strand == "-":   # MAF places minus-strand spans on the reverse complement
                    start = srclen - start - size
                seq, junc, fields = mols[int(mid[1:])]
                a0 = cut[mid] + start
                b0 = a0 + size
                name = f"{prefix}o{n:07d}"
                n += 1
                write_fq(f, name, s, "+")
                record(name, "ont", fields, len(seq), junc, a0, b0, crosses(junc, a0, b0, a.flank),
                       a0 <= 50 and b0 >= len(seq) - 50)
        subprocess.run(["rm", "-rf", tdir])
    truth_f.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("only")
    p.add_argument("--gfa", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--keep", default="GRCh38")
    p.add_argument("--out", required=True)
    p = sub.add_parser("transcriptome")
    p.add_argument("--out", required=True)
    p.add_argument("inputs", nargs="+")
    p = sub.add_parser("simulate")
    p.add_argument("--gfa", required=True, help="graph holding the held-out samples' paths")
    p.add_argument("--hst", required=True, help="prefix of the truth projection (.hst.fa, .hst.tsv)")
    p.add_argument("--samples", required=True)
    p.add_argument("--tx", required=True, help="reference transcripts (rna_sim.py transcriptome)")
    p.add_argument("--interval", required=True, help="the locus on the reference, chr:start-end")
    p.add_argument("--margin", type=int, default=10000, help="transcripts this close to the locus are not negatives")
    p.add_argument("--chroms", default="", help="negatives only from these chromosomes (a partial background)")
    p.add_argument("--paralog-regex", default="")
    p.add_argument("--pseudo-share", type=float, default=0.05, help="share of 31-mers with the locus transcripts")
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--prefix", default="", help="read names start with it (several runs scored together)")
    p.add_argument("--seeds", type=int, default=1, help="expression draws, each with its reads in <out>/s<k>")
    p.add_argument("--cat-dir", default="", help="haplotypes' own annotations (GFF3 per haplotype, named "
                   "<sample>_<hap>_<contig>_...): their transcripts are a second truth, source `cat`")
    p.add_argument("--n-cat-short", type=int, default=1000)
    p.add_argument("--n-cat-long", type=int, default=300)
    p.add_argument("--hap-negatives", action="append", default=[],
                   help="haplotype-specific transcripts of other loci (panvar transcripts FASTA), source paralog_hap")
    p.add_argument("--gene-on", type=float, default=0.75)
    p.add_argument("--sigma", type=float, default=1.5)
    p.add_argument("--perturbed", type=float, default=0.1, help="share of locus molecules of each perturbation")
    p.add_argument("--n-short", type=int, default=3000, help="locus fragments")
    p.add_argument("--n-long", type=int, default=1000, help="locus Kinnex reads")
    p.add_argument("--n-ont", type=int, default=1000, help="locus ONT reads")
    p.add_argument("--n-paralog", type=int, default=1000)
    p.add_argument("--n-pseudo", type=int, default=500)
    p.add_argument("--n-tx", type=int, default=3000)
    p.add_argument("--tx-transcripts", type=int, default=20000, help="transcripts in the transcriptome sample")
    p.add_argument("--fragment", type=int, default=250)
    p.add_argument("--fragment-sd", type=int, default=60)
    p.add_argument("--read-length", type=int, default=150)
    p.add_argument("--short-error", type=float, default=0.002)
    p.add_argument("--kinnex-error", type=float, default=0.001)
    p.add_argument("--ont-accuracy", type=float, default=0.96)
    p.add_argument("--ont-full", type=float, default=0.4, help="share of ONT molecules not cut at the 5' end")
    p.add_argument("--flank", type=int, default=10, help="bases on each side for a read to cross a junction")
    p.add_argument("--pbsim", default="pbsim")
    p.add_argument("--pbsim-model", default="")
    a = ap.parse_args()
    {"only": cmd_only, "transcriptome": cmd_transcriptome, "simulate": cmd_simulate}[a.cmd](a)


if __name__ == "__main__":
    main()
