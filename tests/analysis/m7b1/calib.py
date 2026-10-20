#!/usr/bin/env python3
"""M7-B background and calibration.

  calib.py background <unit dir>             spell + verify the frozen background; simulate one per cell (verified lineage)
  calib.py reference <locus> <unit dir>      reference FASTA (m7bg + locus reference path), faidx, jellyfish
  calib.py preproc <cell> <unit dir>         locityper preproc on the cell's background reads; calibration summary
"""
import gzip
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m7b_common import *  # noqa: E402,F401,F403

mode = sys.argv[1]


def lineage_bin(lineage):
    for r in table(f"{OLD}/runs/provenance/lineage.tsv"):
        if r["lineage"] == lineage:
            assert r["verdict"] == "PASS", f"lineage {lineage} is {r['verdict']}: the affected comparison stops"
            return r["simulator"]


if mode == "background":
    d = sys.argv[2]
    log = f"{d}/commands.log"
    parts, rows = [], []
    for gfa, name in BG_PARTS:
        s = spell_gfa_path(gfa, name)
        n_bad = sum(1 for ch in s if ch not in "ACGT")
        rows.append(dict(check=f"spell {name}", value=f"{len(s)} bp (region {region_length(name)}), non-ACGT {n_bad}",
                         ok=len(s) == region_length(name) and n_bad == 0))
        parts.append(s)
    bg = parts[0] + "N" * SPACER + parts[1]
    for (gfa, name), s_ in zip(BG_PARTS, parts):
        rows.append(dict(check=f"source spelling sha256 {name}", value=hashlib.sha256(s_.encode()).hexdigest(), ok=True))
    rows.append(dict(check="complete background contains only A/C/G/T", value=f"non-ACGT {sum(1 for ch in bg if ch not in 'ACGT')}",
                     ok=all(ch in "ACGT" for ch in bg)))
    rows.append(dict(check="artificial junction (1-based): last ACOT base | first GSTM1 base", value=f"{len(parts[0])} | {len(parts[0]) + 1}", ok=True))
    rows.append(dict(check="background sha256", value=hashlib.sha256(bg.encode()).hexdigest(), ok=True))
    open(f"{d}/m7bg.fa", "w").write(f">m7bg\n{bg}\n")
    with open(f"{d}/m7bg.gfa", "w") as f:
        f.write("H\tVN:Z:1.0\n")
        f.write(f"S\t1\t{bg}\n")
        for p in ("M7REF#0#m7bg", "M7BG#1#m7bg", "M7BG#2#m7bg"):
            f.write(f"P\t{p}\t1+\t*\n")
    spelled = {p: spell_gfa_path(f"{d}/m7bg.gfa", p) for p in ("M7REF#0#m7bg", "M7BG#1#m7bg", "M7BG#2#m7bg")}
    donors = [p for p in spelled if p.startswith("M7BG#")]
    rows.append(dict(check="ploidy: exactly two donor paths, identical spelling, reference present", value=f"{len(donors)} donor paths",
                     ok=len(donors) == 2 and len(set(spelled.values())) == 1 and spelled["M7REF#0#m7bg"] == bg))
    rows.append(dict(check="background length", value=f"{len(bg)} bp (no spacer)", ok=True))
    for cell, c in CELLS.items():
        binary = lineage_bin(c["lineage"])
        o = f"{d}/bg_{cell}"
        os.makedirs(o)
        code = run([binary, "simulate", "--gfa", f"{d}/m7bg.gfa", "--reference", "M7REF#0#m7bg", "--donor", "M7BG", "--depth", c["depth"],
                    "--read-length", 100, "--insert-mean", 350, "--insert-sd", 40, "--sigmas", 4, "--error", 0.002, "--outlier", 0.001,
                    "--seed", c["seed"], "--seed-salt", f"m7b-bg-{cell}", "--out", o], log)
        kv = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(f"{o}/simulation.tsv") if "\t" in l} if code == 0 else {}
        lam_expected = (c["depth"] / 2) / (2 * 100) / (1 - 0.001)
        lam_ok = code == 0 and abs(float(kv["lambda"]) - lam_expected) <= 1e-5 * lam_expected   # Amendment 1: simulation.tsv prints 6 significant digits
        hap_ok = code == 0 and kv.get("haplotype_0_bp") == kv.get("haplotype_1_bp") == str(len(bg))
        z = [(float(kv[f"haplotype_{i}_templates"]) - float(kv[f"haplotype_{i}_expected_templates"])) / float(kv[f"haplotype_{i}_expected_templates"]) ** 0.5
             for i in (0, 1)] if code == 0 else [99]
        rows.append(dict(check=f"simulate background {cell} ({binary.split('/')[-3]}; depth {c['depth']}; salt m7b-bg-{cell})",
                         value=f"exit {code}; lambda {kv.get('lambda')} vs {lam_expected:.7g}; templates z {','.join(f'{x:.2f}' for x in z)}; "
                               f"fragments {kv.get('fragments_written')}, skipped non-ACGT {kv.get('fragments_skipped_non_acgt')}",
                         ok=lam_ok and hap_ok and max(abs(x) for x in z) < 5))
    write_tsv(f"{d}/background.tsv", rows)
    print(open(f"{d}/background.tsv").read())
    sys.exit(0 if all(r["ok"] for r in rows) else 1)

if mode == "reference":
    locus, d = sys.argv[2], sys.argv[3]
    log = f"{d}/commands.log"
    L = LOCI[locus]
    ref = spell_gfa_path(L["gfa"], L["ref_path"])
    bg = fasta(f"{OUT}/runs/background/m7bg.fa")["m7bg"]
    with open(f"{d}/ref.fa", "w") as f:
        f.write(f">m7bg\n{bg}\n>{L['contig']}\n{ref}\n")
    open(f"{d}/locus.tsv", "w").write(f"contig\t{L['contig']}\nlength\t{len(ref)}\nregion_length\t{region_length(L['ref_path'])}\nbackground_length\t{len(bg)}\n")
    assert len(ref) == region_length(L["ref_path"]), "reference path spelling length differs from its region"
    c1 = run([SAMTOOLS, "faidx", f"{d}/ref.fa"], log)
    c2 = run([JELLYFISH, "count", "-C", "-m", 25, "-s", "5M", "-t", 4, "--out-counter-len", 2, "--lower-count", 2, "-o", f"{d}/ref.jf", f"{d}/ref.fa"], log)
    print(open(f"{d}/locus.tsv").read())
    sys.exit(0 if c1 == 0 and c2 == 0 else 1)

if mode == "preproc":
    cell, d = sys.argv[2], sys.argv[3]
    c = CELLS[cell]
    log = f"{d}/commands.log"
    lineage_bin(c["lineage"])
    refdir = f"{OUT}/runs/reference.{c['locus']}"
    bg = f"{OUT}/runs/background/bg_{cell}"
    bglen = len(fasta(f"{OUT}/runs/background/m7bg.fa")["m7bg"])
    code = run([LOCITYPER, "preproc", "-i", f"{bg}/reads_1.fq", f"{bg}/reads_2.fq", "-r", f"{refdir}/ref.fa", "-j", f"{refdir}/ref.jf",
                "-o", f"{d}/preproc", "-b", f"m7bg:1-{bglen}", "--seed", PREPROC_SEED, "-@", 4, "--debug",
                "--strobealign", STROBEALIGN, "--samtools", SAMTOOLS, "--jellyfish", JELLYFISH], log)
    if code != 0:
        sys.exit(1)
    distr = json.load(gzip.open(f"{d}/preproc/distr.gz", "rt"))
    ins, err, dep = distr["insert_distr"], distr["error_profile"], distr["bg_depth"]
    n, p = ins.get("n"), ins.get("p")
    rows = [dict(key="read_len", value=distr["seq_info"].get("read_len")), dict(key="background_total_reads", value=distr["seq_info"].get("total_reads")),
            dict(key="insert_nbinom_n_p", value=f"{n} {p}")]
    if n is not None and p is not None:
        rows.append(dict(key="insert_mean_var_implied", value=f"{n * (1 - p) / p:.2f} {n * (1 - p) / p ** 2:.2f}"))
    rows += [dict(key=f"error_{k}", value=v) for k, v in err.items()]
    rows += [dict(key="depth_window", value=dep.get("window")), dict(key="depth_neighb", value=dep.get("neighb")), dict(key="depth_ploidy", value=dep.get("ploidy"))]
    nn, pp = dep.get("n", []), dep.get("p", [])
    means = [a * (1 - b) / b for a, b in zip(nn, pp)] if nn and pp else []
    rows.append(dict(key="depth_nbinom_mean_by_gc(0..100)", value=",".join(f"{m:.3f}" for m in means)))
    # per-GC support from the debug summary (brotli)
    dbg = f"{d}/preproc/depth.csv.br"
    support = {}
    if os.path.exists(dbg):
        txt = subprocess.run([BROTLI, "-dc", dbg], capture_output=True, text=True).stdout.splitlines()
        hdr = None
        for line in txt:
            f = line.split("\t")
            if f[0] == "gc":
                hdr = f
                continue
            if hdr and len(f) == len(hdr):
                r = dict(zip(hdr, f))
                support[int(r["gc"])] = int(r["nwindows"])
        open(f"{d}/gc_support.tsv", "w").write("gc\tnwindows\n" + "".join(f"{g}\t{support[g]}\n" for g in sorted(support)))
    supported = [g for g, k in support.items() if k >= 100]
    rows.append(dict(key="background_windows_total", value=sum(support.values()) if support else "n/a (no depth.csv.br)"))
    rows.append(dict(key="gc_bins_with_>=100_windows", value=f"{min(supported)}-{max(supported)} ({len(supported)} bins)" if supported else "none"))
    # explicit calibration validation (M7-B.1): finite, valid fitted values; usable windows; Locityper-supported GC bins
    fin = lambda x: isinstance(x, (int, float)) and math.isfinite(x)
    val = []
    val.append(("insert NB n, p finite with n > 0 and 0 < p < 1", fin(n) and fin(p) and n > 0 and 0 < p < 1))
    val.append(("error-profile values finite", all(fin(v) for v in err.values())))
    val.append(("depth NB arrays: 101 GC bins, all finite with n > 0 and 0 < p < 1", len(nn) == 101 and len(pp) == 101 and
                all(fin(a) and fin(b) and a > 0 and 0 < b < 1 for a, b in zip(nn, pp))))
    val.append(("depth window and neighbourhood sizes positive", fin(dep.get("window")) and fin(dep.get("neighb")) and dep.get("window") > 0 and dep.get("neighb") > 0))
    val.append(("background windows counted (depth.csv.br present, total > 0)", bool(support) and sum(support.values()) > 0))
    val.append(("at least one GC bin with >= 100 background windows (Locityper min_tail_obs)", bool(supported)))
    for k, ok in val:
        rows.append(dict(key=f"VALIDATION: {k}", value="PASS" if ok else "FAIL"))
    write_tsv(f"{d}/calibration.tsv", rows)
    open(f"{d}/calibration_ok", "w").write("PASS\n" if all(ok for _, ok in val) else "FAIL\n")
    print(open(f"{d}/calibration.tsv").read()[:3500])
    sys.exit(0 if all(ok for _, ok in val) else 1)
