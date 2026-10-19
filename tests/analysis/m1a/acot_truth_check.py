#!/usr/bin/env python3
"""M1-A: ACOT HG00097 (LZO) -- compare the saved simulation truth (truth.fa, cov30 directory; cov5 reads are a nested
subsample of the same simulation) with the called haplotypes of the S2-L combined 5x and 30x runs, as whole sequences
(SHA-256 of uppercased bytes), allowing a global homologue exchange. Read-only; prints one row per run.
  acot_truth_check.py <out.tsv>
"""
import hashlib
import sys

TRUTH = "/Users/davide.bolognini/panvar/tests/results/genotype_perf_bounded_20261015/inputs/acot/HG00097/cov30/truth.fa"
RUNS = {"acot_5": "/Users/davide.bolognini/panvar/tests/results/s2l_cross_locus_20261019/runs/combined.acot_5",
        "acot_30": "/Users/davide.bolognini/panvar/tests/results/s2l_cross_locus_20261019/runs/combined.acot_30"}


def records(path):
    name, seq, out = None, [], {}
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                out[name] = "".join(seq).upper()
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line)
    if name:
        out[name] = "".join(seq).upper()
    return out


def main(out):
    truth = {k: hashlib.sha256(v.encode()).hexdigest() for k, v in records(TRUTH).items()}
    rows = ["run\ttruth_records\tcalled_hap1_sha\tcalled_hap2_sha\twhole_sequence_agreement\tcalled_paths"]
    for tag, d in RUNS.items():
        hap = {}
        header = None
        for line in open(f"{d}/call.sequence_haplotypes.tsv"):
            cells = line.rstrip("\n").split("\t")
            if header is None:
                header = cells
                continue
            row = dict(zip(header, cells))
            hap[row.get("haplotype", str(len(hap)))] = row
        shas = [r.get("sha256", "NA") for r in hap.values()]
        paths = [r.get("panel_path", r.get("panel_paths", "NA")) for r in hap.values()]
        agree = sorted(shas) == sorted(truth.values())
        rows.append(f"{tag}\t{len(truth)}\t{shas[0] if shas else 'NA'}\t{shas[1] if len(shas) > 1 else 'NA'}\t"
                    f"{'EXACT (both haplotypes, homologue exchange allowed)' if agree else 'DIFFERENT'}\t{';'.join(paths)}")
    open(out, "w").write("\n".join(rows) + "\n")
    print("\n".join(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
